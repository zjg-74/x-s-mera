#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A13 —— 第一例对应探索（换工作决策）的结构性测试（锁结构，不锁数字）

对应：`experiments/a13_first_case/`、`18-v3.3-formal-framework-charter.md` §6.1 A13、
边界地图 `27` §5「三类探索各一例」。

锁的结构预言：
  * A_D 恒等式：θ_U/θ_E = ρ_res/ρ_exp = ρ，u₊ = ρ_exp·e₊
  * U 存量非负（夹紧）—— 未决量不能为负，这是辖域条件
  * D-DRILL 的两条叙事给出**互斥类别**（钻研 vs 空洞），且峰序方向相反
    —— 这正是 D-DRILL 的 external_necessity：无框架只能给程度描述
  * D-CONV 的两条叙事给出**相反方向**（收敛 vs 发散），且各自与其 U 终值/峰高比一致
  * D-IDENT 的四级读法：A 级 ρ_exp 自足、B 级 sign(1−ρ)、C 级需外部、D 级伪目标
  * 判定程序本身有牙齿：U 恒 0 ⟹ 收敛；U 恒为峰值 ⟹ 发散（不是一律给同一答案）

不锁：峰位索引 kU/kE 的具体数值、U 终值的具体数值、收敛判阈的具体取法。

诚实边界（测试内必须保留）：本例是**演示性映射**，不是现实测量；
真实认知过程没有客观可测的 (E,U) 双序列。测试锁的是**映射内部的结构**，
不锁"现实决策确实如此"。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "a13_first_case"))

import a13_first_case as R  # noqa: E402


def _scenarioA():
    cfg = R.make_cfg(theta_E=1.0, rho=2.0)
    return cfg, R.trajectory(cfg, E0=5.0, U0=2.0)


def _scenarioB():
    cfg = R.make_cfg(theta_E=1.0, rho=0.6)
    return cfg, R.trajectory(cfg, E0=5.0, U0=20.0)


def test_module_imports():
    for fn in ("make_cfg", "trajectory", "d_drill", "d_conv", "d_ident", "main"):
        assert hasattr(R, fn), fn


def test_cfg_identities_hold():
    for rho_exp in (0.8, 1.0, 1.5):
        for rho in (0.6, 1.0, 2.0):
            c = R.make_cfg(1.0, rho, rho_exp)
            assert abs(c["theta_U"] / c["theta_E"] - rho) < 1e-12
            assert abs(c["rho_res"] / c["rho_exp"] - rho) < 1e-12
            assert abs(c["u_plus"] - rho_exp * c["e_plus"]) < 1e-12


def test_trajectory_keeps_u_nonnegative():
    for cfg, (Es, Us, _dEs, _dUs) in (_scenarioA(), _scenarioB()):
        assert min(Us) >= 0.0
        assert len(Es) == len(Us)


def test_d_drill_gives_exclusive_categories():
    """两条叙事 → 互斥类别（不是程度差异）。"""
    cfgA, (EsA, UsA, dEsA, dUsA) = _scenarioA()
    cfgB, (EsB, UsB, dEsB, dUsB) = _scenarioB()
    A = R.d_drill(cfgA, EsA, UsA, dEsA, dUsA)
    B = R.d_drill(cfgB, EsB, UsB, dEsB, dUsB)
    assert A["wedge"] == "drill"
    assert B["wedge"] == "hollow"
    assert A["wedge"] != B["wedge"]
    assert A["peak_order"] != B["peak_order"], (A, B)
    assert A["peak_order"] == "U_first"
    assert B["peak_order"] == "E_first"


def test_d_drill_wedge_matches_sign_rho():
    for rho, want in ((2.0, "drill"), (0.6, "hollow"), (1.0, "none")):
        cfg = R.make_cfg(1.0, rho)
        Es, Us, dEs, dUs = R.trajectory(cfg, 5.0, 2.0)
        assert R.d_drill(cfg, Es, Us, dEs, dUs)["wedge"] == want


def test_d_conv_gives_opposite_directions():
    cfgA, (_EA, UsA, _1, _2) = _scenarioA()
    cfgB, (_EB, UsB, _3, _4) = _scenarioB()
    cA = R.d_conv(UsA)
    cB = R.d_conv(UsB)
    assert "收敛" in cA["verdict"]
    assert "发散" in cB["verdict"]
    # 方向必须与" U 终值相对峰高"的读数自洽，而不是靠叙事标签给的
    assert cA["u_end"] < 0.10 * cA["umax"] + 1e-12
    assert cB["u_end"] > 0.10 * cB["umax"]


def test_d_conv_has_teeth():
    """判定程序有牙齿：不同 U 形态给出不同读法，不是一律收敛。"""
    flat = R.d_conv([10.0] * 50)
    assert "发散" in flat["verdict"]
    zero = R.d_conv([0.0] * 50)
    assert "收敛" in zero["verdict"]
    decay = R.d_conv([10.0 * (0.85 ** i) for i in range(20)])
    assert "收敛" in decay["verdict"]


def test_d_ident_four_levels():
    cfgA, _ = _scenarioA()
    cfgB, _ = _scenarioB()
    for cfg, sgn in ((cfgA, -1), (cfgB, 1)):
        d = R.d_ident(cfg)
        assert abs(d["A_readable"]["rho_exp"] - cfg["rho_exp"]) < 1e-12
        assert d["B_readable"]["sign(1-rho)"] == sgn
        assert d["C_needs_external"]
        assert d["D_pseudo_target"]
    # A 级自足、C 级需外部 —— 两者的可读性质相反
    assert R.d_ident(cfgA)["A_readable"] and R.d_ident(cfgA)["C_needs_external"]


def test_main_writes_contrast_record(tmp_path, monkeypatch):
    """对照记录落盘：两个场景 + 无框架人工判断 + 框架判定，三者齐全。"""
    monkeypatch.setattr(R, "HERE", str(tmp_path))
    assert R.main() == 0
    import json
    with open(tmp_path / "a13_contrast.json", "r", encoding="utf-8") as fh:
        res = json.load(fh)
    assert set(res["scenarios"]) == {"A_drill", "B_hollow"}
    for key in ("A_drill", "B_hollow"):
        sc = res["scenarios"][key]
        assert sc["no_framework_human"]
        assert set(sc["framework"]) == {"D_DRILL", "D_CONV", "D_IDENT"}
        # 诚实边界必须随结果一起落盘
        assert "演示" in sc["framework"]["D_CONV"]["honest_boundary"] or \
               ("外部" in sc["framework"]["D_CONV"]["honest_boundary"])
    assert res["scenarios"]["A_drill"]["framework"]["D_DRILL"]["wedge"] == "drill"
    assert res["scenarios"]["B_hollow"]["framework"]["D_DRILL"]["wedge"] == "hollow"
