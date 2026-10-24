#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R42 —— D-SCALE 对 η 的参数敏感性扫描（礁石 R-042）结构性测试（锁结构，不锁数字）

对应：`experiments/r42_dscale_eta/`、`docs/xs-framework/reef_shturl.json` R-042、
`27-v3.4-boundary-map.md` §5（R-042 行）。

判定表（章程 §4.4，固定两观测尺度 τ_fine < τ* < τ_coarse）：

    R(τ_fine)≈0 ∧ R(τ_coarse)≈0  ⟹ 同一尺度耦合
    R(τ_fine)>0 ∧ R(τ_coarse)≈0  ⟹ 跨尺度弱扰动
    否则                        ⟹ 结构性扰动

锁的结构预言：
  * 判定表的三个分支各自可达（不是只用到其中一支）
  * A 同尺度：两尺度残差都 ≈0；B 快振：细尺度>0、粗尺度≈0（被平均掉）；
    C 慢振：粗尺度也 >0（粗窗口覆盖不足整周期）
  * 判定读数只依赖固定的 (R_fine, R_coarse)，与 η 无关（η 只决定 τ* 落点）
  * τ*(η) 随 η 单调不增；η 足够大时 τ* 退到最小窗口
  * 反常读数（细≈0、粗>0）落入结构性分支（判定表完备，不返回空）

不锁：任何具体残差量级（只锁"≈0 或 >0"）、τ* 的具体数值、
以及"合理 η 区间"（τ* 是否落在 (τ_fine, τ_coarse) 内是**登记的失效区间**，
本测试不声称失效区间非空）。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r42_dscale_eta"))

import r42_dscale_eta_sensitivity as R  # noqa: E402


def _residuals(kind):
    dE, dU = R.make_trajectory(kind)
    Rf = R.ols_resid_rms(R.block_mean(dE, R.TAU_FINE), R.block_mean(dU, R.TAU_FINE))
    Rc = R.ols_resid_rms(R.block_mean(dE, R.TAU_COARSE), R.block_mean(dU, R.TAU_COARSE))
    return dU, Rf, Rc


def test_module_imports():
    for fn in ("make_trajectory", "block_mean", "ols_resid_rms",
               "tau_star", "dscale_verdict", "main"):
        assert hasattr(R, fn), fn


def test_observation_protocol_is_ordered():
    assert R.TAU_FINE < R.TAU_COARSE


def test_verdict_table_three_branches_reachable():
    assert R.dscale_verdict(0.0, 0.0) == "same_scale_coupling"
    assert R.dscale_verdict(0.5, 0.0) == "cross_scale_weak_perturbation"
    assert R.dscale_verdict(0.5, 0.4) == "structural_perturbation"


def test_verdict_table_is_total_has_teeth():
    """判定表完备：任何 (Rf,Rc) 都落到三支之一；交换两尺度会改变读数。"""
    pairs = [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5), (1e-12, 1e-12)]
    for rf, rc in pairs:
        assert R.dscale_verdict(rf, rc) in (
            "same_scale_coupling", "cross_scale_weak_perturbation",
            "structural_perturbation")
    # 反常读数（细尺度干净、粗尺度不干净）不被静默归为同尺度
    assert R.dscale_verdict(0.0, 0.5) != "same_scale_coupling"
    assert R.dscale_verdict(0.0, 0.5) != "cross_scale_weak_perturbation"


def test_same_scale_scene_clean_at_both_scales():
    _dU, Rf, Rc = _residuals("same_scale")
    assert Rf <= R.EPS_ZERO and Rc <= R.EPS_ZERO
    assert R.dscale_verdict(Rf, Rc) == "same_scale_coupling"


def test_fast_oscillation_is_averaged_out_at_coarse_only():
    dU, Rf, Rc = _residuals("cross_scale_weak")
    assert Rf > R.EPS_ZERO
    assert Rc <= R.EPS_ZERO
    assert R.dscale_verdict(Rf, Rc) == "cross_scale_weak_perturbation"
    # 机制：粗窗口覆盖整倍数个快周期（τ_coarse 是快周期的整数倍）
    assert R.TAU_COARSE % R.PERIOD_FAST == 0


def test_slow_oscillation_survives_coarse_window():
    dU, Rf, Rc = _residuals("structural")
    assert Rf > R.EPS_ZERO
    assert Rc > R.EPS_ZERO
    assert R.dscale_verdict(Rf, Rc) == "structural_perturbation"
    # 机制：粗窗口不足慢周期的整数倍
    assert R.TAU_COARSE < R.PERIOD_SLOW


def test_tau_star_monotone_nonincreasing_in_eta():
    """η 放宽 ⟹ 分离尺度 τ* 不增（扫描从细到粗，首次满足即返回）。"""
    for kind in ("same_scale", "cross_scale_weak", "structural"):
        _dE, dU = R.make_trajectory(kind)
        seq = [R.tau_star(dU, eta) for eta in R.ETA_GRID]
        assert all(b <= a for a, b in zip(seq, seq[1:])), (kind, seq)


def test_tau_star_saturates_at_minimum_window_for_large_eta():
    dU = R.make_trajectory("cross_scale_weak")[1]
    assert R.tau_star(dU, R.ETA_GRID[-1]) == 2


def test_verdict_readout_is_independent_of_eta():
    """判定读数与 η 无关：η 只决定 τ* 落点（失效区间），不进入读数。"""
    verdicts = {}
    for kind in ("same_scale", "cross_scale_weak", "structural"):
        _dU, Rf, Rc = _residuals(kind)
        base = R.dscale_verdict(Rf, Rc)
        dU = R.make_trajectory(kind)[1]
        for eta in R.ETA_GRID:
            _t = R.tau_star(dU, eta)          # η 只影响 τ*，不影响读数
            assert R.dscale_verdict(Rf, Rc) == base
        verdicts[kind] = base
    # 三个场景分别落到三个不同分支（判定有区分力，不是一律同答）
    assert len(set(verdicts.values())) == 3, verdicts


def test_main_writes_scan_record(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    R.main()
    import json
    with open(tmp_path / "r42_results.json", "r", encoding="utf-8") as fh:
        res = json.load(fh)
    for kind in ("same_scale", "cross_scale_weak", "structural"):
        s = res["stability"][kind]
        assert s["stable"] is True, (kind, s)
    assert res["meta"]["tau_fine"] < res["meta"]["tau_coarse"]
