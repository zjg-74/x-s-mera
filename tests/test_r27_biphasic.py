# -*- coding: utf-8 -*-
"""
R27 回归测试 —— 认知过程的双相形状：X 单调 / M 单峰（补篇六）

设计原则（沿用 R18/R19/R21）：**只锁结构性预言，不锁具体数值。**
峰位 60/120/180 是对数斜坡的产物，换调度就会变，锁死毫无意义；
但下面这些是结构性的，必须恒真：

  1. 无噪声下，U 与 E 的峰位落在解析预测点 k*(theta) 上（用平台比对）
  2. 无噪声下，U 与 E 均单峰
  3. 峰序由 rho_res vs rho_exp 决定（推论 6.4）——A: U 先，B: E 先
  4. 中期相签名：A 为 (+,-) 钻研相；B 为 (-,+) 空洞相；两者都起于 (+,+) 终于 (-,-)
  5. 对照臂（r 恒定）无内点峰，且平均斜率符号正确
  6. 单峰性在**未平滑**时几乎检不出，平滑后恢复（B-B27，防回退）
  7. 峰位估计量在噪声下**均值稳定**（防 T3 回退：用平台下界会系统性下偏到 ~42）
  8. 源码不得用内置 hash() 做种子（防 T1 回退，那会破坏可复现性）

第 7 条是本文件最重要的防回退项：它锁的是**估计量的无偏性**，
不是某个数字。改用平台下界会让它立刻失败。
"""

import os
import re
import sys

import pytest

_EXP = os.path.join(os.path.dirname(__file__), os.pardir, "experiments", "r27_biphasic")
sys.path.insert(0, os.path.abspath(_EXP))

import r27_biphasic as R  # noqa: E402


CFG_A = R.CONFIGS[0]   # theta_E=1, theta_U=3, rho_res=3 > rho_exp=1
CFG_B = R.CONFIGS[1]   # theta_E=3, theta_U=1/3, rho_res=1/9 < rho_exp=1


# ------------------------------------------------------------------ 1 & 2


@pytest.mark.parametrize("cfg", [CFG_A, CFG_B])
def test_peak_locus_matches_analytic_prediction(cfg):
    pl = R.peak_locus(cfg, obs_sd=0.0, smooth_w=1)
    for tag in ("E", "U"):
        assert pl["pred_in_plateau_" + tag] is True, (
            "config {} / {}: analytic k*={} not inside observed plateau {}".format(
                cfg["name"], tag, pl["pred_kstar_" + tag], pl["plateau_" + tag]))


@pytest.mark.parametrize("cfg", [CFG_A, CFG_B])
def test_unimodal_noiseless(cfg):
    pl = R.peak_locus(cfg, obs_sd=0.0, smooth_w=1)
    for tag in ("E", "U"):
        assert pl["unimodal_rate_" + tag] == 1.0, (
            "config {} / {}: unimodal rate {} != 1.0".format(
                cfg["name"], tag, pl["unimodal_rate_" + tag]))


# ------------------------------------------------------------------ 3


@pytest.mark.parametrize("cfg,expected", [(CFG_A, "U_first"), (CFG_B, "E_first")])
def test_peak_order_is_decided_by_rho_comparison(cfg, expected):
    """推论 6.4：峰序 <=> rho_res vs rho_exp。"""
    po = R.peak_order(cfg, obs_sd=0.15)
    assert po["pred_order"] == expected
    assert po["rho_res_gt_rho_exp"] == (expected == "U_first")
    assert po["obs_order_acc"] >= 0.95, (
        "peak order accuracy {} too low".format(po["obs_order_acc"]))


# ------------------------------------------------------------------ 4


def _mid_combo(cfg):
    seq = R.phase_sequence(cfg)
    combos = [s["combo"] for s in seq if s["combo"] != "(+,0)"
              and s["combo"] != "(0,+)" and s["combo"] != "(0,-)"
              and s["combo"] != "(-,0)"]
    return combos


@pytest.mark.parametrize("cfg,expected_mid", [(CFG_A, "(+,-)"), (CFG_B, "(-,+)")])
def test_phase_sequence_has_predicted_mid_phase(cfg, expected_mid):
    combos = _mid_combo(cfg)
    assert combos[0] == "(+,+)", "must start in (+,+)"
    assert combos[-1] == "(-,-)", "must end in (-,-)"
    assert expected_mid in combos, (
        "mid phase {} absent from {}".format(expected_mid, combos))


def test_mid_phases_are_mutually_exclusive():
    """定理 5.2' 的互斥性在时间域的对应：A 只含 (+,-)，B 只含 (-,+)。"""
    a, b = _mid_combo(CFG_A), _mid_combo(CFG_B)
    assert "(+,-)" in a and "(-,+)" not in a
    assert "(-,+)" in b and "(+,-)" not in b


# ------------------------------------------------------------------ 5


@pytest.mark.parametrize("cfg", [CFG_A, CFG_B])
def test_control_arms_have_no_interior_peak(cfg):
    ca = R.control_arms(cfg, obs_sd=0.0)
    arms = ca["arms"]
    # 无消解相：U 单调上升，无峰
    assert arms["no_resolution_r=R_HI"]["U_slope_positive_rate"] == 1.0
    assert arms["no_resolution_r=R_HI"]["U_unimodal_rate"] <= 0.02
    # 无扩张相：U 单调下降，无峰
    assert arms["no_expansion_r=R_LO"]["U_slope_positive_rate"] == 0.0
    assert arms["no_expansion_r=R_LO"]["U_unimodal_rate"] <= 0.02


# ------------------------------------------------------------------ 6 （B-B27）


def test_unimodality_requires_coarse_graining():
    """防回退：单峰性在未平滑时几乎检不出，平滑后必须恢复。

    若有人把容差写死、或改用不平滑序列做结论，这条会失败。
    """
    sr = R.smoothing_recovery(CFG_A, obs_sd=0.15)
    rates = {row["window"]: row["unimodal_rate"] for row in sr["rows"]}
    assert rates[1] < 0.10, "raw-step unimodality should be near zero, got {}".format(rates[1])
    for w in (5, 11, 21, 41):
        assert rates[w] >= 0.95, "window {} recovery too low: {}".format(w, rates[w])


# ------------------------------------------------------------------ 7 （防 T3 回退）


def test_peak_estimator_is_unbiased_under_noise():
    """防回退：峰位估计量的均值必须在噪声下保持稳定。

    用"平台下界"做估计量时，A 组会系统性下偏到 ~42（真值 60）——
    这条会立刻失败。它锁的是**估计量的无偏性**，不是某个数值。
    """
    nr = R.noise_robustness(CFG_A)
    true_k = R.predicted_kstar(CFG_A["theta_U"])
    for row in nr["rows"]:
        assert abs(row["peak_mean_U"] - true_k) <= 3.0, (
            "obs_sd={}: peak estimate {} too far from analytic {}".format(
                row["obs_sd"], row["peak_mean_U"], true_k))
    # 方差必须随噪声单调增大（否则说明估计量没有真的在被噪声影响）
    sds = [row["peak_sd_U"] for row in nr["rows"]]
    assert sds == sorted(sds), "peak sd should grow with noise, got {}".format(sds)


# ------------------------------------------------------------------ 8 （防 T1 回退）


def _code_lines(text):
    """剔除模块 docstring 与 # 注释，只留可执行代码行。

    必须这么做：脚本的 docstring 里**记录了** T1 这个陷阱（含 "hash()" 字样），
    若不做剔除，守卫会把自己的文档当成违规代码，产生假阳性。
    """
    out = []
    in_doc = False
    for ln in text.splitlines():
        if ln.count('"""') % 2 == 1:
            in_doc = not in_doc
            continue
        if in_doc:
            continue
        ln = ln.split("#", 1)[0] if not ln.lstrip().startswith("#") else ""
        if ln.strip():
            out.append(ln)
    return out


def test_source_does_not_use_builtin_hash_for_seeding():
    """防回退：内置 hash() 对 str 按进程随机化，会破坏可复现性。"""
    src = os.path.join(os.path.dirname(__file__), os.pardir,
                       "experiments", "r27_biphasic", "r27_biphasic.py")
    with open(src, encoding="utf-8") as f:
        text = f.read()
    bad = [ln for ln in _code_lines(text)
           if re.search(r"\bhash\s*\(", ln) and "zlib" not in ln]
    assert not bad, "builtin hash() used for seeding (non-reproducible): {}".format(bad)
    assert "zlib.crc32" in text, "zlib_off must use crc32"
