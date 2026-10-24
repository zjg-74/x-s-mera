#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A5 —— generics 一般性检验的结构性测试（锁结构，不锁数字）

对应：`experiments/a5_generics/`（a5_genericity.py + proof.md）、
18-v3.3-formal-framework-charter.md §1.4 的 genericity 要求（A5）。

锁的结构预言：
  * genericity 的函数性形式：R = f∘Inv —— 同 Inv 组内输出恒同（不要求 f 单射）
  * 规则 1/2 的 Inv = sign(ρ−1)：跨 θ_E / ρ_exp / 调度族 / 初值，输出不变
  * 规则 2 的峰序方向由 sign(ρ−1) 决定，且 ρ>1 与 ρ<1 两档给出**相反**的读法
    （判定有牙齿：不是"总是同一答案"）
  * 规则 2 的仿射截距符号 ⟺ ρ>1（非退化点）
  * 规则 3 的四级标记：A/B 自足、C 需外部、D 伪目标
  * 规则 4 的 OLS 斜率 = ρ_exp（连续量，容差内成立）
  * 函数性检查器本身有牙齿：人造违反必须被抓出
  * 完整网格（720 实例）跑通且 all_ok

不锁：任何具体残差/分散度数值量级（只锁"小于容差"），任何具体峰位索引。
ρ=1 是退化点（β₀=0），截距符号无定义 —— 逐条豁免，与脚本一致。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "a5_generics"))

import a5_genericity as R  # noqa: E402


def _instances(rho_values):
    """小规模实例集：跨调度族 × θ_E × ρ_exp × 初值（区别于 ρ 的其余自由度）。"""
    insts = []
    for name in ("log_ramp", "linear"):
        for th in (0.5, 2.0):
            for rho_exp in (0.8, 1.5):
                for rho in rho_values:
                    cfg = R.make_cfg(th, rho, rho_exp)
                    Es, Us, dEs, dUs, _ss = R.trajectory(cfg, name, 5.0, 20.0)
                    insts.append({"cfg": cfg, "sched": name,
                                  "Es": Es, "Us": Us, "dEs": dEs, "dUs": dUs})
    return insts


def test_module_imports():
    for fn in ("make_cfg", "trajectory", "sign_rho", "rule_wedge",
               "rule_peak_order", "rule_identifiability", "rule_collinear_rho",
               "check_functionality", "main"):
        assert hasattr(R, fn), fn


def test_sign_rho_trichotomy():
    assert R.sign_rho(2.0) == 1
    assert R.sign_rho(0.5) == -1
    assert R.sign_rho(1.0) == 0


def test_cfg_identities_hold():
    """结构前提：θ_U/θ_E = ρ_res/ρ_exp = ρ，u₊ = ρ_exp·e₊（A_D 规范固定）。"""
    for rho_exp in (0.8, 1.0, 1.5):
        for rho in (1.0 / 3.0, 1.0, 2.0):
            c = R.make_cfg(1.0, rho, rho_exp)
            assert abs(c["theta_U"] / c["theta_E"] - rho) < 1e-12
            assert abs(c["rho_res"] / c["rho_exp"] - rho) < 1e-12
            assert abs(c["u_plus"] - rho_exp * c["e_plus"]) < 1e-12


def test_rule1_wedge_is_function_of_sign_only():
    """规则 1：Inv=sign(ρ−1)；同组内输出恒同，且三档齐全。"""
    rows = [{"inv": R.rule_wedge(i["cfg"])[0], "out": R.rule_wedge(i["cfg"])[1]}
            for i in _instances((1.0 / 3.0, 1.0, 2.0))]
    ok, ngroups, viol = R.check_functionality("wedge", rows)
    assert ok, viol
    assert ngroups == 3
    outs = {r["out"] for r in rows}
    assert outs == {"drill", "hollow", "none"}


def test_rule1_wedge_values_match_sign():
    for rho, want in ((2.0, "drill"), (0.5, "hollow"), (1.0, "none")):
        cfg = R.make_cfg(1.0, rho)
        assert R.rule_wedge(cfg)[1] == want


def test_rule2_peak_order_is_function_of_sign_and_has_teeth():
    """规则 2：同 sign 组内输出恒同；ρ>1 与 ρ<1 两档读法相反（判定有牙齿）。"""
    rows = []
    for i in _instances((1.0 / 3.0, 2.0)):
        inv, out, extra = R.rule_peak_order(i["cfg"], i["Es"], i["Us"],
                                            i["dEs"], i["dUs"])
        rows.append({"inv": inv, "out": out, "extra": extra})
    ok, _ng, viol = R.check_functionality("peak_order", rows)
    assert ok, viol

    by = {}
    for r in rows:
        by.setdefault(r["inv"][1], set()).add(r["out"])
    assert by[1] == {"U_first"}, by
    assert by[-1] == {"E_first"}, by
    assert by[1] != by[-1]          # 有牙齿：两档给出互斥类别，不是同一答案


def test_rule2_intercept_sign_iff_rho_gt_one():
    """截距符号 ⟺ ρ>1（ρ=1 退化点豁免，与脚本一致）。"""
    rows = []
    for i in _instances((1.0 / 3.0, 1.0, 2.0)):
        inv, _out, extra = R.rule_peak_order(i["cfg"], i["Es"], i["Us"],
                                             i["dEs"], i["dUs"])
        rows.append({"inv": inv, "extra": extra})
    nondeg = [r for r in rows if r["inv"][1] != 0]
    assert nondeg
    for r in nondeg:
        assert r["extra"]["intercept_neg"] == (r["inv"][1] > 0)


def test_rule3_identifiability_levels():
    """规则 3：A/B 自足、C 需外部、D 伪目标；Inv=(ρ_exp, sign(ρ−1))。"""
    rows = [{"inv": R.rule_identifiability(i["cfg"])[0],
             "out": R.rule_identifiability(i["cfg"])[1]}
            for i in _instances((1.0 / 3.0, 1.0, 2.0))]
    ok, _ng, viol = R.check_functionality("identifiability", rows)
    assert ok, viol
    for r in rows:
        o = r["out"]
        assert abs(o["A_slope_rho_exp"] - r["inv"][1]) < 1e-12
        assert o["B_sign(1-rho)"] == -r["inv"][2]
        assert o["levelA_self_sufficient"] is True
        assert o["levelB_self_sufficient"] is True
        assert o["levelC_needs_external"] is True
        assert o["levelD_pseudo_target"] is True


def test_rule4_collinear_slope_recovers_rho_exp():
    """规则 4：OLS(dE,dU) 斜率 = ρ_exp（连续量，容差 1e-9；组内分散 < 1e-9）。"""
    rows = []
    for i in _instances((1.0 / 3.0, 2.0)):
        inv, slope = R.rule_collinear_rho(i["cfg"], i["dEs"], i["dUs"])
        rows.append((inv[1], slope))
    by = {}
    for key, slope in rows:
        by.setdefault(key, []).append(slope)
    for key, slopes in by.items():
        assert max(slopes) - min(slopes) < 1e-9
        for s in slopes:
            assert abs(s - key) / key < 1e-9


def test_functionality_check_has_teeth():
    """检查器本身有牙齿：同 Inv 但输出不同必须被抓出，而不是一律返回通过。"""
    good = [{"inv": ("sgn", 1), "out": "drill"}, {"inv": ("sgn", 1), "out": "drill"}]
    ok, ngroups, viol = R.check_functionality("good", good)
    assert ok and ngroups == 1 and viol == []

    bad = [{"inv": ("sgn", 1), "out": "drill"},
           {"inv": ("sgn", 1), "out": "hollow"},
           {"inv": ("sgn", -1), "out": "hollow"}]
    ok, ngroups, viol = R.check_functionality("bad", bad)
    assert not ok
    assert ngroups == 2
    assert len(viol) == 1


def test_full_grid_genericity_holds(tmp_path, monkeypatch):
    """完整 720 实例网格跑通：四条规则的函数性全部成立（走脚本主流程）。"""
    monkeypatch.chdir(tmp_path)
    assert R.main() == 0
    import json
    with open(tmp_path / "a5_results.json", "r", encoding="utf-8") as fh:
        res = json.load(fh)
    assert res["overall"]["all_ok"] is True
    assert res["overall"]["n_instances"] > 0
    for key in ("wedge", "peak_order", "identifiability", "collinear_rho"):
        assert res["rules"][key]["genericity"] is True, key
        assert res["rules"][key]["n_inv_groups"] >= 2, key
    assert res["rules"]["peak_order"]["intercept_sign_crosscheck"] is True
