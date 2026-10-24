#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B2 —— R31 桥检验（A_D 一致性守卫）的结构性测试（锁结构，不锁数字）

对应：`experiments/b2_bridge_test/`、
`docs/xs-framework/15-v3.2-gauge-fixed-framework.md` §7.2（R43′）与 §7.4（R31）。

**R31 与 R43′ 必须分清（测试内保留这条边界）**：
  * R31（本脚本）= 检验 β₀ = −ρ_exp·ΔE(k*_U) 的**重言形式**，作用是
    **A_D 一致性守卫**——它不是对桥的检验，只暴露 A_D 一致性破坏。
  * R43′（真桥检验）需独立于 ΔE 的第三通道 Ψ，**不在本脚本范围内**。

锁的结构预言：
  * 内点峰探测：单峰序列给出峰索引；单调序列给出 None（有牙齿）
  * 纯 A_D 族峰点检验 ≈0（无噪声语义下的离散化余量内）
  * 注入 A_D⁺ 边界层（g(U)=s·U）后峰点检验显著 > 阈值
  * 正常族误报率 ≈0（<0.05）；最强违反族误报率 →1.0（≥0.9）
  * 回归在纯 A_D 上恢复 ρ_exp（= u₊/e₊）与 β₀（= ρ_exp·e₋ − u₋）

不锁：任何具体相对违反量级（只锁"低于/高于判阈"）、任何具体峰位。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "b2_bridge_test"))

import b2_bridge_test as R  # noqa: E402

THRESH = 0.02


def test_module_imports():
    for fn in ("CFG", "gen_trajectory", "multi_linreg",
               "find_inner_peak", "peak_test", "main"):
        assert hasattr(R, fn), fn


def test_find_inner_peak_has_teeth():
    up = [1.0, 3.0, 9.0, 4.0, 2.0]
    assert R.find_inner_peak(up) == 2
    mono = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert R.find_inner_peak(mono) is None
    down = [5.0, 4.0, 3.0, 2.0, 1.0]
    assert R.find_inner_peak(down) is None


def test_regression_recovers_ad_law_on_clean_family():
    """纯 A_D（g=0）：回归读数 = A_D 闭合读数。"""
    ks, Es, Us = R.gen_trajectory(R.CFG, 0.5, g_strength=0.0)
    b = R.multi_linreg([Es, ks], Us)
    rho_exp = R.CFG["u_plus"] / R.CFG["e_plus"]
    beta0 = rho_exp * R.CFG["e_minus"] - R.CFG["u_minus"]
    assert abs(b[1] - rho_exp) / abs(rho_exp) < 1e-6
    assert abs(b[2] - beta0) < 1e-6


def test_peak_test_near_zero_on_clean_ad():
    for pivot in (0.4, 0.6):
        ks, Es, Us = R.gen_trajectory(R.CFG, pivot, g_strength=0.0)
        b = R.multi_linreg([Es, ks], Us)
        kstar = R.find_inner_peak(Us)
        assert kstar is not None, pivot
        assert R.peak_test(Es, Us, b[1], b[2], kstar) < THRESH, pivot


def test_peak_test_fires_on_adplus_violation():
    fired = 0
    for pivot in (0.4, 0.5, 0.6, 0.7):
        ks, Es, Us = R.gen_trajectory(R.CFG, pivot, g_strength=1.2)
        b = R.multi_linreg([Es, ks], Us)
        kstar = R.find_inner_peak(Us)
        if kstar is None:
            continue
        if R.peak_test(Es, Us, b[1], b[2], kstar) > THRESH:
            fired += 1
    assert fired >= 3, fired


def test_normal_family_no_false_alarm():
    fails = 0
    total = 0
    for pivot in (0.4, 0.5, 0.6, 0.7):
        ks, Es, Us = R.gen_trajectory(R.CFG, pivot, g_strength=0.0)
        b = R.multi_linreg([Es, ks], Us)
        kstar = R.find_inner_peak(Us)
        if kstar is None:
            continue
        total += 1
        if R.peak_test(Es, Us, b[1], b[2], kstar) > THRESH:
            fails += 1
    assert total >= 4
    assert fails / total < 0.05


def test_violating_family_detected():
    worst = 0.0
    for g in (0.4, 0.8, 1.2):
        fails = 0
        total = 0
        for pivot in (0.4, 0.5, 0.6, 0.7):
            ks, Es, Us = R.gen_trajectory(R.CFG, pivot, g_strength=g)
            b = R.multi_linreg([Es, ks], Us)
            kstar = R.find_inner_peak(Us)
            if kstar is None:
                continue
            total += 1
            if R.peak_test(Es, Us, b[1], b[2], kstar) > THRESH:
                fails += 1
        if total:
            worst = max(worst, fails / total)
    assert worst >= 0.9, worst


def test_main_reports_family_rates(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "HERE", str(tmp_path))
    R.main()
    import json
    with open(tmp_path / "b2_results.json", "r", encoding="utf-8") as fh:
        res = json.load(fh)
    assert res["families"]["normal"]["fail_rate"] < 0.05
    adplus = [k for k in res["families"] if k.startswith("adplus")]
    assert adplus
    assert max(res["families"][k]["fail_rate"] for k in adplus) >= 0.9
    # 边界声明必须随结果落盘：R31 是重言形式，不是 R43′ 真桥检验
    assert "R43" in res["meta"]["note"] or "重言" in res["meta"]["note"]
