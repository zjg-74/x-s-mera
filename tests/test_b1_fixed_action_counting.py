#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B1 —— 固定动作计数（命题 6.3 / R33）的结构性测试（锁结构，不锁数字）

对应：`experiments/b1_fixed_action_counting/`、
`docs/xs-framework/15-v3.2-gauge-fixed-framework.md` §7.3。

性质（必须随测试保留）：本脚本是**方法 / 分析流程准备**，在合成 A_D 数据上验证
「四臂测量 → 单位系数 → 闭合 → 与回归自洽」的整条链路。**不是真人实验**，
不对真实认知过程断言任何系数。

锁的结构预言：
  * 四臂（固定展开 1 次 / 固定消解 1 次）恢复四个单位系数 (e₊,e₋,u₊,u₋)
  * 闭合恒等式：θ_E=e₋/e₊、θ_U=u₋/u₊、ρ=ρ_res/ρ_exp=θ_U/θ_E、κ_E=e₋
  * 仿射律：U = ρ_exp·E + β₀·t，β₀ = ρ_exp·e₋ − u₋ = κ_E·ρ_exp·(1−ρ)
  * 回归读数与四臂闭合读数自洽（R33 判据）
  * β₀ 符号随相位翻转：ρ>1 为负、ρ<1 为正、ρ=1 近零（有牙齿：三档不同）
  * 把 u₋ 记错 ⟹ 自洽性被破坏（判据有牙齿，不是一律通过）

不锁：任何具体残差量级（只锁"小于判阈"）、任何具体回归系数值。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "b1_fixed_action_counting"))

import b1_fixed_action_counting as R  # noqa: E402


def test_module_imports():
    for fn in ("CONFIGS", "four_arm_measure", "main_trajectory",
               "multi_linreg", "step_delta", "main"):
        assert hasattr(R, fn), fn


def test_four_arm_recovers_unit_coefficients():
    for name, cfg in R.CONFIGS.items():
        arms = R.four_arm_measure(cfg)
        for key in ("e_plus", "e_minus", "u_plus", "u_minus"):
            assert abs(arms[key] - cfg[key]) < 1e-12, (name, key)


def test_closure_identities_hold():
    for name, cfg in R.CONFIGS.items():
        a = R.four_arm_measure(cfg)
        theta_E = a["e_minus"] / a["e_plus"]
        theta_U = a["u_minus"] / a["u_plus"]
        rho_exp = a["u_plus"] / a["e_plus"]
        rho_res = a["u_minus"] / a["e_minus"]
        assert abs(theta_U / theta_E - rho_res / rho_exp) < 1e-12, name
        assert abs(rho_res / rho_exp - cfg["u_minus"] * cfg["e_plus"]
                   / (cfg["e_minus"] * cfg["u_plus"])) < 1e-12, name
        # κ_E 惯例 = e₋（与 a13/r33 一致）；β₀ = κ_E·ρ_exp·(1−ρ)
        kappa_E = a["e_minus"]
        beta0 = rho_exp * kappa_E - a["u_minus"]
        rho = rho_res / rho_exp
        assert abs(beta0 - kappa_E * rho_exp * (1.0 - rho)) < 1e-9, name


def test_affine_law_recovered_by_regression():
    """U = ρ_exp·E + β₀·t：回归斜率 = ρ_exp，时间项 = β₀ = ρ_exp·e₋ − u₋。"""
    for name, cfg in R.CONFIGS.items():
        a = R.four_arm_measure(cfg)
        rho_exp_arm = a["u_plus"] / a["e_plus"]
        beta0_arm = rho_exp_arm * a["e_minus"] - a["u_minus"]
        ks, Es, Us = R.main_trajectory(cfg)
        b = R.multi_linreg([Es, ks], Us)
        assert abs(b[1] - rho_exp_arm) / abs(rho_exp_arm) < 1e-6, (name, b[1])
        assert abs(b[2] - beta0_arm) < 1e-6, (name, b[2], beta0_arm)


def test_beta0_sign_tracks_phase():
    """β₀ = κ_E·ρ_exp·(1−ρ)：三相位给出不同符号（有牙齿）。"""
    signs = {}
    for name, cfg in R.CONFIGS.items():
        a = R.four_arm_measure(cfg)
        rho_exp = a["u_plus"] / a["e_plus"]
        rho = (a["u_minus"] / a["e_minus"]) / rho_exp
        beta0 = rho_exp * a["e_minus"] - a["u_minus"]
        signs[name] = (0 if abs(beta0) < 1e-9 else (1 if beta0 > 0 else -1), rho)
    assert signs["drill"][0] == -1 and signs["drill"][1] > 1.0
    assert signs["converge"][0] == 1 and signs["converge"][1] < 1.0
    assert signs["boundary"][0] == 0 and abs(signs["boundary"][1] - 1.0) < 1e-9


def test_self_consistency_has_teeth():
    """把 u₋ 记错（或 e₊ 记错）必须破坏自洽 —— 判据不是一律通过。"""
    cfg = R.CONFIGS["drill"]
    a = R.four_arm_measure(cfg)
    ks, Es, Us = R.main_trajectory(cfg)
    b = R.multi_linreg([Es, ks], Us)

    rho_exp_arm = a["u_plus"] / a["e_plus"]
    beta0_arm = rho_exp_arm * a["e_minus"] - a["u_minus"]
    clean = (abs(rho_exp_arm - b[1]) / abs(b[1]) < 1e-6
             and abs(beta0_arm - b[2]) < 1e-6)
    assert clean

    # 负对照 1：u₋ 记成 2 倍（破坏 β₀，不动 ρ_exp）
    bad_u = dict(a, u_minus=2.0 * a["u_minus"])
    beta0_bad = (bad_u["u_plus"] / bad_u["e_plus"]) * bad_u["e_minus"] - bad_u["u_minus"]
    assert abs(beta0_bad - b[2]) >= 1e-6

    # 负对照 2：e₊ 记成 2 倍（破坏 ρ_exp）
    bad_e = dict(a, e_plus=2.0 * a["e_plus"])
    rho_bad = bad_e["u_plus"] / bad_e["e_plus"]
    assert abs(rho_bad - b[1]) / abs(b[1]) >= 1e-6


def test_main_writes_results(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "HERE", str(tmp_path))
    R.main()
    import json
    with open(tmp_path / "b1_results.json", "r", encoding="utf-8") as fh:
        res = json.load(fh)
    for name in R.CONFIGS:
        sc = res["configs"][name]["self_consistency"]
        assert sc["rho_exp_ok"] is True, (name, sc)
        assert sc["beta0_ok"] is True, (name, sc)
    assert "非真人实验" in res["meta"]["task"] or "方法准备" in res["meta"]["task"]
