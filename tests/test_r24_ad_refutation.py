# -*- coding: utf-8 -*-
"""
R24 回归测试 —— A_D 证否签名的判别力分解（补篇九）

设计原则（沿用 R18/R19/R21/R27/R28/R29）：**只锁结构性预言，不锁具体数值。**

  1.  臂 A：签名只由 rho 决定 ⟹ 跨任务（不同 theta_E）一致
  2.  臂 A：归一化峰距命中解析预测，误差 ≤ 离散化硬界
  3.  臂 A：任务内线性失配残差 ≡ 0（定理 6.11）
  4.  臂 A：最优缩放 v* = 1/rho_exp
  5.  臂 A：E 单峰（无内部极小）
  6.  臂 B：eta = 0 时跨任务判据**无误报**
  7.  臂 B ★：eta > 0 时跨任务判据**假阳**，误判率随 eta 上升 —— T12 防回退
  8.  臂 B ★：同一批数据上，任务内残差判据**始终报 A_D 成立**
  9.  臂 B ★：同一批数据上，任务内签名唯一性**始终成立**
 10.  臂 C：残差关于 eps 单调不减
 11.  臂 C ★：残差/eps 关于 eps 为常数（R38）
 12.  臂 C ★：内部极小法在全部网格点漏检（N3 防回退）
 13.  臂 C ★：任务内双象限法在全部网格点漏检（N3 防回退）
 14.  臂 E2 ★：残差/eps 关于 theta_E 为常数（R38 标度不变性）
 15.  臂 E2 ★：翻转点**不是** w = 1（N4 防回退）
 16.  臂 D ★：共享噪声下残差恒为 0，与 s 无关（定理 6.11 防回退）
 17.  臂 D ★：独立噪声下残差 > 0 且随 s 增大（T16 的来源）
 18.  臂 F ★：三分类 SOP 在全部噪声水平下全部判对
 19.  臂 F ★：同一批数据上，R24 原判据误判 D3
 20.  臂 G ★：r 可观测 ⟹ 四系数精确复原；r 不可观测 ⟹ 误差随 s 增长
 21.  臂 G2 ★：规范变换下 (E, U) 逐点不变（定理 6.12 防回退）
 22.  臂 G2 ★：规范变换下四个系数确实改变
 23.  臂 H ★：eps = 0 时 log-log 斜率 ≈ -1/2（R39 防回退）
 24.  臂 H ★：eps > 0 时残差收敛到系统地板、斜率 ≈ 0
 25.  源码守卫：不使用内置 hash()

第 7–9、16–19、21–24 条是本文件最重要的防回退项：
  - 7–9  锁「R24 原签名对 rho 任务依赖假阳，残差判据不假阳」
  - 16–17 锁「守恒依赖共享驱动」及其破裂条件
  - 18–19 锁「三分类 SOP 有效，而 R24 原判据在 D3 上失效」
  - 21–22 锁「A_D 的规范自由度」（B-B23 不闭合的结构原因）
  - 23–24 锁「重复测量平均律」（T16 的修正手段）
"""

import math
import os
import random
import re
import sys

import pytest

_EXP = os.path.join(os.path.dirname(__file__), os.pardir,
                    "experiments", "r24_ad_refutation")
sys.path.insert(0, os.path.abspath(_EXP))

import r24_ad_refutation_signature as R  # noqa: E402


# ------------------------------------------------------------------ 缓存

_CACHE = {}


def _cached(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


def arm_A():
    return _cached("A", R.arm_A_prediction)


def arm_B():
    return _cached("B", lambda: R.arm_B_confound(n_rep=60, s_grid=(0.0, 0.2)))


def arm_C():
    return _cached("C", lambda: R.arm_C_separation(n_rep=40, s_grid=(0.0, 0.2)))


def arm_D():
    return _cached("D", lambda: R.arm_D_noise_regime(n_rep=40))


def arm_E2():
    return _cached("E2", lambda: R.arm_E2_threshold(
        n_rep=40, s_grid=(0.0,), theta_grid=(0.5, 1.0, 2.0)))


def arm_F():
    return _cached("F", lambda: R.arm_F_sop(n_rep=40, s_grid=(0.0, 0.2, 0.8)))


def arm_G():
    return _cached("G", lambda: R.arm_G_identification(
        n_rep=30, s_grid=(0.0, 0.4, 0.8), ))


def arm_G2():
    return _cached("G2", R.arm_G2_gauge)


def arm_H():
    return _cached("H", lambda: R.arm_H_obs_floor(n_rep=40))


# ------------------------------------------------------------------ 1-5 臂 A

def test_A_signature_is_determined_by_rho_alone():
    """签名只由 rho 决定 ⟹ 跨任务（不同 theta_E）必须一致。"""
    for r in arm_A():
        assert r["cross_task_consistent"], \
            "u_plus={} rho={}: signature varies across tasks".format(
                r["u_plus"], r["rho"])


def test_A_peak_distance_matches_analytic_prediction():
    """无噪声：归一化峰距命中 ln(rho)/L，误差 ≤ 离散化硬界 1/(K-1)。"""
    bnd = 1.0 / (R.K - 1.0)
    for r in arm_A():
        assert r["max_abs_d_err"] <= bnd, \
            "rho={}: d_err={} exceeds discretization bound {}".format(
                r["rho"], r["max_abs_d_err"], bnd)


def test_A_linear_mismatch_residual_is_identically_zero():
    """★ 定理 6.11：A_D + 共享驱动 ⟹ 任务内线性失配残差逐样本为 0。"""
    for r in arm_A():
        for t in r["tasks"]:
            assert t["lin_resid_rms"] < R.RESID_TOL, \
                "u_plus={} rho={} theta_E={}: resid={} not zero".format(
                    r["u_plus"], r["rho"], t["theta_E"], t["lin_resid_rms"])


def test_A_optimal_scaling_recovers_rho_exp():
    """★ 定理 6.11 的副产品：最优缩放 v* = 1/rho_exp。"""
    for r in arm_A():
        assert r["max_vstar_err"] < 1e-9, \
            "u_plus={}: vstar error {} too large".format(
                r["u_plus"], r["max_vstar_err"])


def test_A_E_trajectory_is_unimodal():
    """A_D ⟹ E 单峰 ⟹ 无内部极小。"""
    for r in arm_A():
        assert r["max_interior_min"] == 0, \
            "rho={}: found {} interior minima under A_D".format(
                r["rho"], r["max_interior_min"])


# ------------------------------------------------------------------ 6-9 臂 B

def test_B_no_false_positive_when_rho_is_task_invariant():
    """eta = 0（rho 真的任务不变）⟹ 跨任务判据必须**不**报证否。"""
    for blk in arm_B():
        if blk["eta"] != 0.0:
            continue
        for r in blk["rows"]:
            assert r["r24_cross_task_both_rate"] == 0.0, \
                "false positive at eta=0, s={}: {}".format(
                    r["s"], r["r24_cross_task_both_rate"])


def test_B_r24_criterion_false_positives_on_rho_task_dependence():
    """★ 防回退 T12：eta 大时跨任务判据的误判率必须显著 > 0。

    这是 R24 被取代的**直接证据**：A_D 在每个任务内部严格成立，
    只因 rho 跨任务变化，R24 原签名就报「证否」。
    """
    rows_by_eta = {}
    for blk in arm_B():
        rows_by_eta.setdefault(blk["eta"], []).extend(blk["rows"])
    eta_max = max(rows_by_eta)
    worst = max(r["r24_cross_task_both_rate"] for r in rows_by_eta[eta_max])
    assert eta_max > 0.0 and worst > 0.5, \
        "expected R24 criterion to false-positive at eta={}, got {}".format(
            eta_max, worst)


def test_B_within_task_residual_never_fires_under_rho_task_dependence():
    """★ 同一批数据：任务内残差判据必须**始终**报 A_D 成立。"""
    for blk in arm_B():
        for r in blk["rows"]:
            assert r["within_task_max_resid_rms"] < R.RESID_TOL, \
                "eta={} s={}: residual={} should stay ~0".format(
                    blk["eta"], r["s"], r["within_task_max_resid_rms"])


def test_B_within_task_signature_stays_unique_under_rho_task_dependence():
    """★ 同一批数据：任务内签名必须始终唯一（不并存两种象限）。"""
    for blk in arm_B():
        for r in blk["rows"]:
            assert r["within_task_multi_sig_rate"] == 0.0, \
                "eta={} s={}: multi-signature rate {} should be 0".format(
                    blk["eta"], r["s"], r["within_task_multi_sig_rate"])


# ------------------------------------------------------------------ 10-13 臂 C

def _c_at(eps, s):
    for blk in arm_C():
        if abs(blk["eps"] - eps) < 1e-12:
            for r in blk["rows"]:
                if abs(r["s"] - s) < 1e-12:
                    return r
    return None


def test_C_residual_is_monotone_in_violation_strength():
    """残差关于 eps 单调不减（灵敏度可排序）。"""
    for s in (0.0, 0.2):
        prev = -1.0
        for eps in R.EPS_GRID:
            r = _c_at(eps, s)
            assert r is not None
            assert r["lin_resid_rms_median"] >= prev - 1e-12, \
                "residual decreased at s={} eps={}".format(s, eps)
            prev = r["lin_resid_rms_median"]


def test_C_residual_is_linear_in_violation_strength():
    """★ R38：残差/eps 关于 eps 为常数（线性响应）。"""
    ratios = []
    for eps in R.EPS_GRID:
        if eps <= 0.0:
            continue
        r = _c_at(eps, 0.0)
        assert r is not None
        ratios.append(r["lin_resid_rms_median"] / eps)
    assert len(ratios) >= 5
    lo, hi = min(ratios), max(ratios)
    assert hi - lo < 1e-6 * max(1.0, hi), \
        "residual/eps not constant: {} .. {}".format(lo, hi)


def test_C_interior_min_criterion_never_fires():
    """★ N3 防回退：二次注入下内部极小法在本网格上恒漏检。"""
    for blk in arm_C():
        for r in blk["rows"]:
            assert r["within_task_interior_min_rate"] <= 0.05, \
                "eps={} s={}: unexpected interior-min detection {}".format(
                    blk["eps"], r["s"],
                    r["within_task_interior_min_rate"])


def test_C_within_task_both_quadrant_criterion_never_fires():
    """★ N3 防回退：任务内双象限法在本网格上恒漏检。"""
    for blk in arm_C():
        for r in blk["rows"]:
            assert r["within_task_both_sig_rate"] == 0.0, \
                "eps={} s={}: unexpected both-quadrant detection".format(
                    blk["eps"], r["s"])


# ------------------------------------------------------------------ 14-15 臂 E2

def test_E2_residual_sensitivity_is_theta_E_independent():
    """★ R38 标度不变性：残差/eps 关于 theta_E 为常数。

    这是本轮**最强的结构量**：eps 跨 20 倍、theta_E 跨 4 倍，
    灵敏度系数六位有效数字一致 ⟹ 锁的是结构，不是数字。
    """
    per_theta = {}
    for blk in arm_E2():
        vals = [r["resid_over_eps"] for r in blk["rows"] if r["eps"] > 0]
        assert vals
        per_theta[blk["theta_E"]] = sum(vals) / len(vals)
    assert len(per_theta) >= 2
    lo, hi = min(per_theta.values()), max(per_theta.values())
    assert hi - lo < 1e-6 * max(1.0, hi), \
        "residual/eps varies with theta_E: {}".format(per_theta)


def test_E2_flip_point_is_not_at_w_equals_one():
    """★ N4 防回退：翻转点**不是** w = 1（不是 eps_c），且随 theta_E 变化。

    锁这条的意义：防止后续有人把「峰退到左边界」当成
    一个与模型无关的可辨识阈值来用。
    """
    flips = {}
    for blk in arm_E2():
        rows = sorted(blk["rows"], key=lambda r: r["w"])
        flip = None
        for r in rows:
            if r["peak_at_left_edge_rate"] > 0.5:
                flip = r["w"]
                break
        assert flip is not None, \
            "no flip found for theta_E={}".format(blk["theta_E"])
        flips[blk["theta_E"]] = flip
    assert len(set(flips.values())) > 1, \
        "flip point should depend on theta_E, got {}".format(flips)
    assert any(abs(w - 1.0) > 1e-9 for w in flips.values()), \
        "flip point should not sit at w=1, got {}".format(flips)


# ------------------------------------------------------------------ 16-17 臂 D

def test_D_shared_noise_keeps_conservation_exact():
    """★ 定理 6.11 防回退：共享驱动下残差恒为 0，与噪声水平 s 无关。"""
    for r in arm_D():
        if r["regime"] != "shared":
            continue
        assert r["lin_resid_rms_median"] < R.RESID_TOL, \
            "s={} rho={}: shared-noise residual {} not zero".format(
                r["s"], r["rho"], r["lin_resid_rms_median"])


def test_D_independent_noise_breaks_conservation_and_grows_with_s():
    """★ T16：独立通道噪声破坏守恒 ⟹ 残差 > 0 且随 s 单调不减。"""
    for rho in (1.0 / 3.0, 1.0, 3.0):
        prev = -1.0
        for s in R.S_GRID:
            rows = [r for r in arm_D()
                    if r["regime"] == "indep" and r["s"] == s
                    and abs(r["rho"] - rho) < 1e-9]
            assert len(rows) == 1
            cur = rows[0]["lin_resid_rms_median"]
            assert cur >= prev - 1e-12
            if s > 0.0:
                assert cur > 1e-4, \
                    "independent noise at s={} should break conservation".format(s)
            prev = cur


# ------------------------------------------------------------------ 18-19 臂 F

def test_F_sop_classifies_all_three_datasets_correctly():
    """★ 三分类 SOP：D1⟶accept, D2⟶reject_A_D, D3⟶rho_task_dependent。"""
    want = {"D1_A_D_invariant": "accept",
            "D2_A_D_violated": "reject_A_D",
            "D3_rho_task_dependent": "rho_task_dependent"}
    for r in arm_F():
        for key, exp in want.items():
            assert r[key]["verdict"] == exp, \
                "s={}: {} judged {} (expected {})".format(
                    r["s"], key, r[key]["verdict"], exp)


def test_F_r24_alone_misjudges_rho_task_dependence():
    """★ 同一批数据上 R24 原判据必须误判 D3 为「证否」。"""
    for r in arm_F():
        assert r["r24_alone_would_misjudge_D3_as_reject"], \
            "s={}: R24 criterion should misjudge D3".format(r["s"])


# ------------------------------------------------------------------ 20 臂 G

def test_G_regression_recovers_four_coefficients_when_r_observable():
    """★ P8a：r 可观测 ⟹ 回归精确复原四系数，与 s 无关。"""
    for r in arm_G():
        if r["obs"] != "r":
            continue
        assert max(r["rel_err_rho_exp"], r["rel_err_rho_res"]) < 1e-9, \
            "u_plus={} rho={} s={}: error too large".format(
                r["u_plus"], r["rho"], r["s"])


def test_G_r_unobservable_degrades_with_noise_but_stays_consistent():
    """★ P8b：r 不可观测 ⟹ 抖动随 s 增长；但 s=0 时仍精确。"""
    rows = [r for r in arm_G() if r["obs"] == "rbar"]
    assert rows
    for r in rows:
        if r["s"] == 0.0:
            assert max(r["rel_err_rho_exp"], r["rel_err_rho_res"]) < 1e-9
    worst_lo = max(r["rel_err_rho_res"] for r in rows if r["s"] <= 0.4)
    worst_hi = max(r["rel_err_rho_res"] for r in rows if r["s"] >= 0.8)
    assert worst_hi >= worst_lo, \
        "error should not shrink with more noise"


# ------------------------------------------------------------------ 21-22 臂 G2

def test_G2_gauge_transformation_leaves_trajectories_invariant():
    """★ 定理 6.12 防回退：G_lambda 下 (E(k), U(k)) 逐点不变。"""
    for r in arm_G2():
        assert r["max_abs_E_diff"] == 0.0, \
            "lambda={}: E trajectory moved by {}".format(
                r["lambda"], r["max_abs_E_diff"])
        assert r["max_abs_U_diff"] == 0.0, \
            "lambda={}: U trajectory moved by {}".format(
                r["lambda"], r["max_abs_U_diff"])


def test_G2_gauge_transformation_does_change_the_coefficients():
    """★ 定理 6.12 的另一半：四个系数（的绝对值）确实变了 ⟹ 不可辨识。"""
    rows = {r["lambda"]: r for r in arm_G2()}
    assert 1.0 in rows
    base = rows[1.0]
    changed = []
    for lam, r in sorted(rows.items()):
        if lam == 1.0:
            continue
        same = (r["e_plus"], r["u_plus"], r["theta_E"], r["theta_U"]) == \
               (base["e_plus"], base["u_plus"],
                base["theta_E"], base["theta_U"])
        if not same:
            changed.append(lam)
        # 规范变换的代数：e+ -> e+/lambda, theta_E -> lambda*theta_E
        assert abs(r["e_plus"] - base["e_plus"] / lam) < 1e-12
        assert abs(r["theta_E"] - base["theta_E"] * lam) < 1e-12
    assert changed, "coefficients must change under the gauge transformation"


# ------------------------------------------------------------------ 23-24 臂 H

def test_H_averaging_drives_residual_to_zero_when_no_violation():
    """★ R39 防回退：eps=0 时残差按 N^{-1/2} 衰减 ⟹ log-log 斜率 ≈ -1/2。"""
    for sg in (0.10, 0.40):
        sl = R.fit_loglog_slope(arm_H(), 0.0, sg)
        assert sl is not None, "no slope fitted for eps=0, sigma={}".format(sg)
        assert -0.65 < sl[0] < -0.35, \
            "expected slope ~ -0.5 at eps=0, sigma={}, got {}".format(sg, sl[0])
    # 逐点单调：N 增大 ⟹ 残差不增
    for sg in (0.10, 0.40):
        prev = None
        for N in R.N_AVG_GRID:
            got = [r["lin_resid_rms_median"] for r in arm_H()
                   if r["eps"] == 0.0 and r["sigma_obs"] == sg and r["N"] == N]
            assert len(got) == 1
            if prev is not None:
                assert got[0] <= prev * 1.15, \
                    "residual should not grow with N at eps=0"
            prev = got[0]


def test_H_averaging_plateaus_at_systematic_floor_when_violated():
    """★ R39：eps>0 时残差收敛到系统地板（∝ eps），不降到 0。

    ⟹ 这给出了操作规程：**先平均再判残差**，否则观测噪声会被读成 ¬A_D。
    """
    eps = 0.05
    # 系统地板 = 臂 E2 标定的灵敏度系数 * eps
    kappa = None
    for blk in arm_E2():
        vals = [r["resid_over_eps"] for r in blk["rows"] if r["eps"] > 0]
        kappa = sum(vals) / len(vals)
        break
    assert kappa is not None
    floor = kappa * eps
    for sg in (0.10, 0.40):
        sl = R.fit_loglog_slope(arm_H(), eps, sg)
        assert sl is not None
        assert abs(sl[0]) < 0.25, \
            "expected slope ~0 at eps={}, sigma={}, got {}".format(eps, sg, sl[0])
        big = [r["lin_resid_rms_median"] for r in arm_H()
               if r["eps"] == eps and r["sigma_obs"] == sg
               and r["N"] == max(R.N_AVG_GRID)]
        assert len(big) == 1
        assert big[0] > 0.5 * floor, \
            "residual should converge to systematic floor {} , got {}".format(
                floor, big[0])
    # sigma_obs = 0 时残差与 N 无关
    for N in R.N_AVG_GRID:
        got = [r["lin_resid_rms_median"] for r in arm_H()
               if r["eps"] == eps and r["sigma_obs"] == 0.0 and r["N"] == N]
        assert len(got) == 1
        assert abs(got[0] - floor) < 1e-3, \
            "sigma=0: residual should equal systematic floor"


# ------------------------------------------------------------------ 25 源码守卫

def _code_lines(text):
    """剔除 docstring 与 # 注释，只留可执行代码行。

    必须这么做：脚本的 docstring 里**记录了**这些陷阱（含关键字），
    若不剔除，守卫会把自己的文档当成违规代码（R29 实测到的假阳性）。
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


def _src_path():
    return os.path.join(os.path.dirname(__file__), os.pardir, "experiments",
                        "r24_ad_refutation",
                        "r24_ad_refutation_signature.py")


def test_source_does_not_use_builtin_hash():
    """防 T1 回退：不能用内置 hash() 做种子（str hash 按进程随机化）。"""
    with open(_src_path(), encoding="utf-8") as f:
        text = f.read()
    bad = [ln for ln in _code_lines(text)
           if re.search(r"\bhash\s*\(", ln) and "zlib" not in ln]
    assert not bad, "builtin hash() used: {}".format(bad)
    assert "zlib.crc32" in text


def test_source_records_traps_T12_to_T16():
    """防回退：脚本 docstring 须登记本轮陷阱 T12–T16。"""
    with open(_src_path(), encoding="utf-8") as f:
        text = f.read()
    for t in ("T12", "T13", "T14", "T15", "T16"):
        assert t in text, "trap {} not documented in source".format(t)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
