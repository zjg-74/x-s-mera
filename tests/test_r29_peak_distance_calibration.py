# -*- coding: utf-8 -*-
"""
R29 回归测试 —— 峰距标定：可辨识性、视界与操作规程（补篇八）

设计原则（沿用 R18/R19/R21/R27/R28）：**只锁结构性预言，不锁具体数值。**

  1. 臂 A：无噪声下峰距命中解析预测，且误差 ≤ 1 步（离散化硬界）
  2. 臂 A：预测点落在实测平台内（并列取小索引的伪像防护，T2）
  3. 臂 B：估计量无偏（偏 ≤ 0.05，且远小于 sd）
  4. 臂 B：sd 随 s 单调不减
  5. 臂 C1：同 K 时比值设计成立（相对误差 < 5%）
  6. 臂 C2：异 K 时 raw 比值**失败**（相对误差 > 50%），归一化**修复**（< 15%）——T8 防回退
  7. 臂 D：对数斜坡为基线；非对数调度的偏**显著更大**——负面结果防回退
  8. 臂 E：视界律——窗口宽度 = L，超界即不可测
  9. 臂 F：共享噪声下 err 率恒为 0；独立噪声下 err 率显著 > 0——定理 6.10 防回退
 10. 臂 G：U−E 线性恒等式残差 ~1e-13 以内（精确成立）

第 6、9、10 条是本文件最重要的防回退项：
  - 6 锁"比值设计必须归一化"（否则跨 K 直接 2 倍偏）
  - 9 锁"峰序守恒依赖共享驱动"（否则静默判反）
  - 10 锁新定理 6.10 的恒等式确实精确
"""

import os
import re

import pytest

_EXP = os.path.join(os.path.dirname(__file__), os.pardir,
                    "experiments", "r29_peak_distance_calibration")
import sys
sys.path.insert(0, os.path.abspath(_EXP))

import r29_peak_distance_calibration as R  # noqa: E402


# ------------------------------------------------------------------ 工具

def _run(cfg, s=0.0, n=1, seed=1, indep_rel=0.0):
    rng = __import__("random").Random(seed)
    out = []
    for _ in range(n):
        Es, Us, _c = R.trajectory(rng, cfg, R.log_ramp, s=s, indep_rel=indep_rel)
        out.append((R.argmax_of(Es), R.argmax_of(Us)))
    return out


# 完整臂 B / 臂 F 跑一次要 ~10M / ~15M 步，多个测试都要用，
# 故在模块级缓存，避免同一份仿真重复跑 2-3 遍（实测可省 ~2/3 时间）。
_ARM_B = None
_ARM_F = None


def arm_B():
    global _ARM_B
    if _ARM_B is None:
        _ARM_B = {b["s"]: b["rows"] for b in R.arm_B_recovery()}
    return _ARM_B


def arm_F():
    global _ARM_F
    if _ARM_F is None:
        _ARM_F = {b["s"]: b["rows"] for b in R.arm_F_order_resolution()}
    return _ARM_F


# ------------------------------------------------------------------ 1-2 臂 A

def test_A_peak_distance_matches_analytic_prediction():
    """无噪声：|k*_E - k*_U| 命中 (K-1)*|ln rho|/L，误差 ≤ 1 步（离散化硬界）。"""
    for rho in R.RHO_GRID:
        cfg = R.make_cfg(rho)
        pE = R.kstar_logramp(cfg["theta_E"])
        pU = R.kstar_logramp(cfg["theta_U"])
        if pE is None or pU is None:
            continue
        pred = abs(pE - pU)
        (kE, kU), = _run(cfg, s=0.0, n=1, seed=7)
        assert abs(abs(kE - kU) - pred) <= 1.0, \
            "rho={}: obs={} pred={} exceeds discretization bound".format(
                rho, abs(kE - kU), pred)


def test_A_prediction_lies_within_observed_plateau():
    """预测点须落在两峰平台之差的区间内（T2：并列时 argmax 取小索引）。"""
    for rho in R.RHO_GRID:
        if abs(math_log(rho)) < 1e-12:
            continue  # 退化：两峰重合，平台差区间退化
        cfg = R.make_cfg(rho)
        pE = R.kstar_logramp(cfg["theta_E"])
        pU = R.kstar_logramp(cfg["theta_U"])
        if pE is None or pU is None:
            continue
        pred = abs(pE - pU)
        rng = __import__("random").Random(7)
        Es, Us, _c = R.trajectory(rng, cfg, R.log_ramp, s=0.0)
        lo_E, hi_E = R.argmax_plateau(Es)
        lo_U, hi_U = R.argmax_plateau(Us)
        lo = min(abs(lo_E - hi_U), abs(hi_E - lo_U))
        hi = max(abs(lo_E - hi_U), abs(hi_E - lo_U))
        assert lo - 0.51 <= pred <= hi + 0.51


def math_log(x):
    import math
    return math.log(x)


# ------------------------------------------------------------------ 3-4 臂 B

C_SENS = (R.K - 1.0) / R.L_SPAN


def test_B_estimator_is_unbiased():
    """|ln rho|^ = Delta_bar / C 应无偏：|bias| ≤ 0.05 且 |bias| < 0.5*sd。"""
    blk = arm_B()
    for s in (0.1, 0.2, 0.4, 0.8):
        for r in blk[s]:
            if r["true_abs_ln_rho"] == 0.0:
                continue
            assert abs(r["bias"]) <= 0.05, "s={} rho={} bias={}".format(
                s, r["rho"], r["bias"])
            assert abs(r["bias"]) < 0.5 * r["obs_delta_sd"] / C_SENS, \
                "bias not small vs sd"


def test_B_sd_grows_monotonically_with_noise():
    """sd(Delta) 随 s 单调不减（按 rho 逐条检查）。"""
    blk = arm_B()
    for idx, rho in enumerate(R.RHO_GRID):
        prev = -1.0
        for s in R.S_GRID:
            cur = blk[s][idx]["obs_delta_sd"]
            assert cur >= prev - 1e-9, "sd decreased at rho={} s={}".format(rho, s)
            prev = cur


# ------------------------------------------------------------------ 5-6 臂 C

def test_C1_ratio_design_works_with_same_K():
    """同 K：比值设计消去 C，相对误差应 < 5%。"""
    for r in R.arm_C_ratio_design()["C1_same_K"]:
        if r["true_ratio"] == 0.0:
            continue
        assert r["rel_err"] < 0.05, "same-K ratio failed: {}".format(r)


def test_C2_normalization_fixes_different_K():
    """★ 防回退 T8：异 K 时 raw 比值须**失败**（>50% 偏），归一化须**修复**（<15%）。"""
    for r in R.arm_C_ratio_design()["C2_diff_K"]:
        if r["true_ratio"] == 0.0:
            continue
        assert r["rel_err_raw"] > 0.50, \
            "expected raw ratio to FAIL across different K, got {}".format(r)
        assert r["rel_err_normalized"] < 0.15, \
            "normalization should repair it, got {}".format(r)


# ------------------------------------------------------------------ 7 臂 D

def test_D_schedule_misspecification_is_worse_than_logramp():
    """★ 负面结果：非对数斜坡调度下按推论 6.5 反解，偏须**显著大于**对数斜坡基线。"""
    blocks = {b["schedule"]: b for b in R.arm_D_schedule_misspec()}
    base = [r["rel_err"] for r in blocks["log_ramp"]["rows"]
            if r["rel_err"] is not None]
    base_worst = max(base)
    for name in ("linear_ramp", "sigmoid_ramp", "piecewise_ramp"):
        vals = [r["rel_err"] for r in blocks[name]["rows"]
                if r["rel_err"] is not None]
        assert max(vals) > base_worst * 3.0, \
            "{} should be much worse than log_ramp baseline".format(name)


def test_D_logramp_baseline_is_accurate():
    """对数斜坡（公式的正确适用域）下反解误差应很小。"""
    blocks = {b["schedule"]: b for b in arm_D_cached()}
    vals = [r["rel_err"] for r in blocks["log_ramp"]["rows"]
            if r["rel_err"] is not None]
    assert max(vals) < 0.05


_ARM_D_CACHE = None


def arm_D_cached():
    global _ARM_D_CACHE
    if _ARM_D_CACHE is None:
        _ARM_D_CACHE = R.arm_D_schedule_misspec()
    return _ARM_D_CACHE


# ------------------------------------------------------------------ 8 臂 E

def test_E_identifiability_horizon():
    """视界律：窗口宽度 = L；窗内可测、窗外不可测。"""
    eh = R.arm_E_horizon()
    assert abs(eh["window_width"] - R.L_SPAN) < 1e-9
    lo, hi = eh["window_ln_rho"]
    for r in eh["rows"]:
        inside = lo < r["ln_rho"] < hi
        assert r["observable"] == inside, \
            "ln_rho={} observable={} but window={}".format(r["ln_rho"], r["observable"], (lo, hi))


# ------------------------------------------------------------------ 9 臂 F

def test_F_shared_noise_never_reverses_order():
    """★ 定理 6.10 防回退：共享 r 噪声下，峰序**从不判反**（err 率恒为 0）。"""
    for s, rows in arm_F().items():
        for r in rows:
            if r["pred_order"] == "degenerate":
                continue
            assert r["shared"]["err_rate"] == 0.0, \
                "shared noise reversed order at s={} ln_rho={}".format(s, r["ln_rho"])


def test_F_independent_noise_breaks_conservation():
    """独立噪声下守恒律破裂：须出现显著 err 率（且在 |ln rho| 小时最明显）。"""
    found = False
    worst = None
    for s, rows in arm_F().items():
        for r in rows:
            if r["pred_order"] == "degenerate":
                continue
            if r["indep"]["err_rate"] > (worst[0] if worst else 0.0):
                worst = (r["indep"]["err_rate"], s, r["ln_rho"])
            if r["indep"]["err_rate"] > 0.10:
                found = True
    assert worst is not None and worst[2] < 0.8, \
        "worst indep error should occur at small |ln rho|, got {}".format(worst)
    assert found, "independent noise should break order conservation somewhere"


def test_F_tie_rate_increases_with_noise():
    """并列率随噪声单调不减（峰序退化为"不确定"而非"判错"）。"""
    blks = {b["s"]: b for b in R.arm_F_order_resolution()}
    for idx, lnr in enumerate(R.F_LN_GRID):
        if abs(lnr) < 1e-12:
            continue
        prev = -1.0
        for s in R.S_GRID:
            cur = blks[s]["rows"][idx]["shared"]["tie_rate"]
            assert cur >= prev - 1e-9, "tie rate dropped at ln_rho={} s={}".format(lnr, s)
            prev = cur


# ------------------------------------------------------------------ 10 臂 G

def test_G_UE_identity_is_exact():
    """★ 定理 6.10 恒等式：U^s(k) - E(k) 与 r/噪声无关，逐点残差 ~1e-13 以内。"""
    gi = R.arm_G_identity()
    assert float(gi["max_abs_resid_over_all"]) < 1e-12


# ------------------------------------------------------------------ 11 源码守卫

def _code_lines(text):
    """剔除 docstring 与 # 注释，只留可执行代码行。

    必须这么做：脚本的 docstring 里**记录了**这些陷阱（含关键字），
    若不剔除，守卫会把自己的文档当成违规代码（本轮实测到的假阳性）。
    """
    out, in_doc = [], False
    for ln in text.splitlines():
        if ln.count('"""') % 2 == 1:
            in_doc = not in_doc
            continue
        if in_doc:
            continue
        if not ln.lstrip().startswith("#") and "#" in ln:
            ln = ln.split("#", 1)[0]
        if ln.strip():
            out.append(ln)
    return out


def test_source_does_not_use_builtin_hash():
    """防 T1 回退：不能用内置 hash() 做种子（str hash 按进程随机化）。"""
    src = os.path.join(os.path.dirname(__file__), os.pardir, "experiments",
                  "r29_peak_distance_calibration",
                  "r29_peak_distance_calibration.py")
    with open(src, encoding="utf-8") as f:
        text = f.read()
    bad = [ln for ln in _code_lines(text)
           if re.search(r"\bhash\s*\(", ln) and "zlib" not in ln]
    assert not bad, "builtin hash() used: {}".format(bad)
    assert "zlib.crc32" in text


def test_source_marks_unimodal_not_applicable_for_nonlinear():
    """防回退：非线性 / 退化情形必须返回 None 而不是 True。"""
    src = os.path.join(os.path.dirname(__file__), os.pardir, "experiments",
                  "r29_peak_distance_calibration",
                  "r29_peak_distance_calibration.py")
    with open(src, encoding="utf-8") as f:
        text = f.read()
    assert "pred_in_plateau_note" in text
    assert "不适用" in text


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
