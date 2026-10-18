"""
R4（A8/D17 分量归类判据）的单元测试
====================================

锚点：LE 内核 F2 残差的解析结构 r_i = −(ratio−1)·gain·δ_i − base·|gain|；
判据顺序 K1→K2→K3；复合现象先分解后归类。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.r4_component_attribution import (
    collect_f2_ensemble, decompose_residual, k1_state_correlation,
    k2_repair_test, k3_filter_test, classify_le_f2, ATTRIBUTION_CRITERIA,
)
from src.le_coupling_kernel import LECouplingKernel


def test_criteria_preregistered():
    assert ATTRIBUTION_CRITERIA["K1_corr_threshold"] == 0.5
    assert ATTRIBUTION_CRITERIA["K2_elimination_min"] == 0.90
    assert "顺序" in ATTRIBUTION_CRITERIA["note"]


def test_conservation_when_c2_ok():
    """C2 通过时残差恒为零（守恒）——分解的零假设对照。"""
    trace = collect_f2_ensemble(n_runs=2, n_steps=40, seed=0,
                                c2_compatibility=True)
    # C2 通过时：C1 过的步守恒（残差=0），C1 失败的步通量被禁（残差=0）
    assert np.allclose(trace["residual"], 0.0)


def test_decomposition_recovers_true_coefficients():
    """复合场景：α、β 应精确恢复解析真值（确定性残差，R²=1）。"""
    res = classify_le_f2(n_steps=120, seed=0)
    for row in res["rows"]:
        assert row["alpha_mismatch"] == pytest.approx(res["true_alpha"], abs=1e-6)
        assert row["beta_base"] == pytest.approx(res["true_beta"], abs=1e-6)
        assert row["r_squared"] == pytest.approx(1.0, abs=1e-6)


def test_k1_mismatch_state_correlated():
    """K1：失配分量与状态差强相关 → 接口域。"""
    res = classify_le_f2(n_steps=120, seed=0)
    for row in res["rows"]:
        assert row["mismatch_state_correlated"] is True
        assert abs(row["mismatch_rho"]) >= ATTRIBUTION_CRITERIA["K1_corr_threshold"]
        assert "接口域" in row["mismatch_domain"]


def test_k2_repair_eliminates_residual():
    """K2：修复接口（C2 过）后残差消除率 ≥ 90%。"""
    rp = k2_repair_test(n_steps=120, seed=0)
    assert rp["repair_eliminates"] is True
    assert rp["residual_norm_f2"] > 0.05
    assert rp["residual_norm_repaired"] == pytest.approx(0.0, abs=1e-12)


def test_k3_base_sign_invariant():
    """K3：底噪分量在 δ 翻号两子样本中均值同号（偏置型）→ 载体域。"""
    res = classify_le_f2(n_steps=120, seed=0)
    for row in res["rows"]:
        assert row["base_sign_invariant"] is True
        assert "载体域" in row["base_domain"]


def test_pure_mismatch_control():
    """纯失配对照（底噪=0）：β 应恢复为 0，α 不变。"""
    res = classify_le_f2(n_steps=120, seed=0, ghost_base_ratio=0.0)
    row = res["rows"][0]
    assert row["alpha_mismatch"] == pytest.approx(-0.6, abs=1e-6)
    assert row["beta_base"] == pytest.approx(0.0, abs=1e-6)


def test_pure_base_control():
    """纯底噪对照（失配比=1）：α 应恢复为 0，β 不变。"""
    res = classify_le_f2(n_steps=120, seed=0, ghost_gain_ratio=1.0)
    row = res["rows"][0]
    assert row["alpha_mismatch"] == pytest.approx(0.0, abs=1e-6)
    assert row["beta_base"] == pytest.approx(-0.25, abs=1e-6)


def test_k1_uncorrelated_detected():
    """K1 的反向判别力：与 δ 无关的分量不应判为状态相关。"""
    rng = np.random.default_rng(0)
    delta = rng.normal(size=500)
    noise = rng.normal(size=500)          # 与 δ 独立
    r = k1_state_correlation(noise, delta)
    assert r["state_correlated"] is False
    assert abs(r["rho"]) < 0.2


def test_k3_insufficient_split():
    """δ 单号（不变号）时 K3 应报样本不足，而不是误判。"""
    base = -np.abs(np.linspace(0.1, 1.0, 100))
    delta = np.abs(np.linspace(0.1, 1.0, 100))     # 恒正
    r = k3_filter_test(base, delta)
    assert r["sign_invariant"] is False
    assert r.get("reason") == "insufficient_split"


def test_main_mode_conserved_in_f2():
    """F2 下主模态（分量 0）仍守恒——分解只针对副模态。"""
    trace = collect_f2_ensemble(n_runs=2, n_steps=40, seed=1)
    assert np.allclose(trace["residual"][trace["active"]][:, 0], 0.0)


def test_reproducible():
    a = classify_le_f2(n_steps=60, seed=5)["rows"][0]["alpha_mismatch"]
    b = classify_le_f2(n_steps=60, seed=5)["rows"][0]["alpha_mismatch"]
    assert a == pytest.approx(b)
