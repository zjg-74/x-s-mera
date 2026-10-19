# -*- coding: utf-8 -*-
"""
R19 回归测试 —— 定理 5.2'（U 与 Psi 的桥）与推论 5.3（震荡的生产性）

设计原则（沿用 R18）：**只锁结构性预言，不锁具体数值。**
数值随 e_pm / u_pm / 网格密度变化，锁死会让任何参数探索都失败；
但「哪几种符号组合可达」是结构性的，必须恒真。

本测试锁的五件事：
  1. rho_res > rho_exp  ->  可出现钻研型 (+,-)，**不可**出现空洞型 (-,+)
  2. rho_res < rho_exp  ->  反之
  3. 两者相等          ->  两种反号组合**都**不可达
  4. 违反 A_D 线性性   ->  两种反号组合**并存**（证明证否条件非空）
  5. 推论 5.3          ->  震荡期 U 的走向由 rho_res / rho_exp 决定

另加一条防回退：非线性配置**不得**被判为"符合推论 5.3"
（推论依赖 A_D 的线性性，写 True 会是假的确认）。
"""

import os
import sys

import pytest

_EXP = os.path.join(os.path.dirname(__file__), os.pardir, "experiments", "r19_bridge")
sys.path.insert(0, os.path.abspath(_EXP))

import r19_bridge as R  # noqa: E402


@pytest.fixture(scope="module")
def configs():
    return {c["name"]: c for c in R.CONFIGS}


def _run(name, configs, n_rep=25):
    return R.run_config(configs[name], n_rep=n_rep)


def test_rho_res_greater_excludes_hollow(configs):
    """rho_res > rho_exp：可钻研 (+,-)，不可空洞 (-,+)。"""
    got = set(_run("I_rho_res>rho_exp", configs)["实测出现组合"])
    assert "(+,-)" in got, "钻研型 (+,-) 应可达"
    assert "(-,+)" not in got, "空洞型 (-,+) 在 rho_res > rho_exp 下不可达（定理 5.2'）"


def test_rho_res_less_excludes_studious(configs):
    """rho_res < rho_exp：可空洞 (-,+)，不可钻研 (+,-)。"""
    got = set(_run("II_rho_res<rho_exp", configs)["实测出现组合"])
    assert "(-,+)" in got, "空洞型 (-,+) 应可达"
    assert "(+,-)" not in got, "钻研型 (+,-) 在 rho_res < rho_exp 下不可达（定理 5.2'）"


def test_degenerate_excludes_both_reverse_combinations(configs):
    """rho_res == rho_exp：theta_E == theta_U，两种反号都不可达。"""
    got = set(_run("III_degenerate", configs)["实测出现组合"])
    assert "(+,-)" not in got and "(-,+)" not in got
    assert "(+,+)" in got and "(-,-)" in got


def test_nonlinear_violation_breaks_exclusivity(configs):
    """
    违反 A_D 线性性 -> 两种反号组合**并存**。

    这一条锁的是**证否条件的判别力**：若"两种反号并存"在原理上不可能发生，
    那么定理 5.2' 的证否条件（R24）就是一句空话。本测试确保该签名确实可出现。
    """
    got = set(_run("IV_nonlinear_A_D_violated", configs)["实测出现组合"])
    assert "(+,-)" in got and "(-,+)" in got, \
        "违反 A_D 后两种反号组合应并存；若不并存，则 R24 的证否条件是空的"


def test_oscillation_U_direction_follows_efficiency_ratio(configs):
    """推论 5.3：震荡期 U 的走向由 rho_res / rho_exp 决定。"""
    for name in ("I_rho_res>rho_exp", "II_rho_res<rho_exp", "III_degenerate"):
        probes = R.oscillation_probe(configs[name], n_rep=25)
        assert probes, f"{name} 应有至少一个 Psi=0 点"
        for p in probes:
            assert p["符合推论5.3"] is True, f"{name} 在 r={p['r_set_to_Psi0']} 违反推论 5.3"
            ratio = p["rho_res/rho_exp"]
            if ratio > 1.0:
                assert p["震荡期U净变化均值"] < 0, "rho_res>rho_exp：震荡期 U 应下降（酝酿）"
            elif ratio < 1.0:
                assert p["震荡期U净变化均值"] > 0, "rho_res<rho_exp：震荡期 U 应上升（反刍）"


def test_nonlinear_is_marked_not_applicable(configs):
    """
    防回退：非线性配置**不得**被判为"符合推论 5.3"。

    推论 5.3 的推导依赖 A_D 的线性性；对违反 A_D 的配置返回 True
    会是假的确认（两个 Psi=0 点给出互相矛盾的 U 走向）。
    """
    for p in R.oscillation_probe(configs["IV_nonlinear_A_D_violated"], n_rep=25):
        assert p["符合推论5.3"] is None
        assert "不适用" in p["推论5.3适用性"]


def test_psi_zero_points_of_nonlinear_config(configs):
    """非线性配置的 Psi=0 应有两个点（这正是它能同时出现两种反号的原因）。"""
    pts = R.psi_zero_points(configs["IV_nonlinear_A_D_violated"])
    assert len(pts) == 2, f"非线性配置应有两个 Psi=0 点，实得 {pts}"
    assert pts[0] < pts[1]
    assert len(R.psi_zero_points(configs["I_rho_res>rho_exp"])) == 1
