#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R33 —— C 级闭合（R45）结构性测试（锁结构，不锁数字）

对应：`experiments/r33_clevel_closure/`、`16-v3.2-self-consistent-theory.md`
§11 命题 6.3 与 §12.3 R45 判据。

锁结构：
  * 四臂可读出 C 级量（e±, u±, theta, rho, kappa_E），恒等式成立
  * 自洽判定 T1–T4 跨 (theta_E, rho, 调度族) 全部通过（残差为浮点噪声）
  * 负对照（b 记错）必须破坏 T2 —— 判据有牙齿
  * rho=1（beta=0）的退化情形同样自洽

不锁：任何具体残差数值、任何具体回归系数。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r33_clevel_closure"))

import r33_clevel_closure as R  # noqa: E402


def test_module_imports():
    assert hasattr(R, "make_cfg")
    assert hasattr(R, "check_consistency")
    assert hasattr(R, "run_all")


def test_arm_readout_recovers_coefficients():
    """臂 1/2 直接读出的四系数与恒等式必须成立（C 级量可读）。"""
    for th in R.THETA_E_GRID:
        for rho in R.RHO_GRID:
            cfg = R.make_cfg(th, rho)
            d = R.arm_derived(cfg)
            assert abs(d["e_plus"] - cfg["e_plus"]) < 1e-12
            assert abs(d["e_minus"] - cfg["e_minus"]) < 1e-12
            assert abs(d["u_plus"] - cfg["u_plus"]) < 1e-12
            assert abs(d["u_minus"] - cfg["u_minus"]) < 1e-12
            assert abs(d["rho"] - d["theta_U"] / d["theta_E"]) < 1e-12
            assert abs(d["rho"] - d["rho_res"] / d["rho_exp"]) < 1e-12
            assert abs(d["kappa_E"] - 1.0 * d["e_minus"]) < 1e-12


def test_r45_consistency_holds_across_families():
    """T1–T4 跨调度族全部通过（R45 判据）。"""
    res = R.run_all()
    assert res["summary"]["all_pass"], res["summary"]
    for row in res["rows"]:
        assert row["verdict"] == "pass", row
        assert row["rel_resid"] < 1e-9, row


def test_r45_independent_of_theta_E_and_rho():
    """结构预言：自洽不依赖参数点（不是某一组参数的巧合）。"""
    res = R.run_all()
    for th in R.THETA_E_GRID:
        sub = [r for r in res["rows"] if abs(r["arm_values"]["theta_E"] - th) < 1e-12]
        assert sub and all(r["verdict"] == "pass" for r in sub)


def test_negative_control_breaks_T2():
    """b 记错必须破坏 T2（判据有牙齿），且不能破坏 T1（斜率不受 b 影响）。
    rho=1 退化点豁免：beta 预测恒为 0，负对照在该点无牙齿。"""
    res = R.run_all()
    assert res["summary"]["negative_control_breaks_T2"]
    for row in res["negative_control"]:
        assert row["T1_rho_exp"] is True, row
        if abs(row["arm_values"]["rho"] - 1.0) > 1e-12:
            assert row["T2_beta"] is False, row


def test_rho_equals_one_degenerate_but_consistent():
    """rho=1 时 beta 预测值为 0，自洽判据仍应通过（0 对 0）。"""
    cfg = R.make_cfg(0.6, 1.0)
    r = R.check_consistency(cfg, "log_ramp")
    assert r["verdict"] == "pass", r
    assert abs(r["beta_pred_from_arms"]) < 1e-12
