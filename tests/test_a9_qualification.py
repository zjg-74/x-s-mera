"""R5 单元测试：Wright–Fisher 仿真核、d 构念、判据②③①、资格裁决。

最重要的两个正确性锚点：
  1. 中性情形 P_fix ≈ p0（仿真的统计正确性）；
  2. 选择情形 P_fix 落在扩散近似解析值 (1-e^{-2Ns·p0})/(1-e^{-2Ns}) 上
     （这是群体遗传学的教科书公式——若仿真错，后面全部判读作废）。
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.a9_qualification import (  # noqa: E402
    DS, P0, REPS, TH_D, TH_EFF, T_GENS,
    closure_leak, coarse_predict, d_overlap, frozen_transfer,
    intervention_response, p_fix, qualify, wright_fisher,
)


# ---------------------------------------------------------------- 仿真核
def test_neutral_pfix_equals_p0():
    x = wright_fisher(P0, 0.0, 200, T_GENS, 2000, np.random.default_rng(0))
    assert abs(p_fix(x) - P0) < 0.04


def test_pfix_matches_diffusion_formula():
    """教科书解析值：P_fix = (1 − e^{−2Ns·p0}) / (1 − e^{−2Ns})。"""
    for N, s in ((50, 0.02), (200, 0.02), (100, 0.05)):
        x = wright_fisher(P0, s, N, T_GENS, 3000, np.random.default_rng(N))
        ns, p0 = N * s, P0
        exact = (1.0 - math.exp(-2 * ns * p0)) / (1.0 - math.exp(-2 * ns))
        assert abs(p_fix(x) - exact) < 0.04, f"N={N}, s={s}: {p_fix(x):.3f} vs {exact:.3f}"


def test_pfix_monotone_in_s():
    a = p_fix(wright_fisher(P0, 0.01, 100, T_GENS, 1500, np.random.default_rng(1)))
    b = p_fix(wright_fisher(P0, 0.05, 100, T_GENS, 1500, np.random.default_rng(1)))
    assert b > a + 0.15


def test_migration_pulls_toward_pool():
    x_lo = wright_fisher(P0, 0.02, 2000, T_GENS, 1000,
                         np.random.default_rng(2), m=0.10, p_ext=0.2)
    x_hi = wright_fisher(P0, 0.02, 2000, T_GENS, 1000,
                         np.random.default_rng(2), m=0.10, p_ext=0.8)
    assert x_hi.mean() > x_lo.mean() + 0.3


# ---------------------------------------------------------------- d 构念
def test_d_overlap_properties():
    rng = np.random.default_rng(3)
    a = rng.normal(0.5, 0.05, 4000)
    assert d_overlap(a, a) == pytest.approx(1.0)
    assert d_overlap(np.full(500, 0.01), np.full(500, 0.99)) == pytest.approx(0.0)
    b = rng.normal(0.5, 0.05, 4000)
    assert 0.0 < d_overlap(a, b) <= 1.0
    # 对称
    assert d_overlap(a, b) == pytest.approx(d_overlap(b, a))


# ---------------------------------------------------------------- 判据②
def test_statistical_reading_passes_everywhere():
    for N in (50, 2000, 20000):
        real = wright_fisher(P0, 0.02, N, T_GENS, 1000, np.random.default_rng(N))
        pred = coarse_predict(P0, 0.02, N, T_GENS, 1000, np.random.default_rng(9000 + N))
        assert d_overlap(pred, real) >= TH_D, f"statistical reading fails at N={N}"


def test_intervention_insensitive_at_small_Ns():
    r = intervention_response(50, seed=50, reps=1000)
    assert r["effect"] < TH_EFF, f"Ns=1 should be insensitive, eff={r['effect']}"


def test_intervention_sensitive_at_large_Ns():
    r = intervention_response(20000, seed=20000, reps=1000)
    assert r["effect"] > 0.8, f"Ns=400 should respond strongly, eff={r['effect']}"


def test_deterministic_reading_fails_at_small_Ns():
    """原降级理由的读法：点预测容差带通过率在 Ns=1 处接近 0。"""
    rng = np.random.default_rng(11)
    real = wright_fisher(P0, 0.02, 50, T_GENS, 1000, rng)
    p_det = P0
    for _ in range(T_GENS):
        p_det = p_det * 1.02 / (1.0 + 0.02 * p_det)
    assert float(np.mean(np.abs(real - p_det) < 0.05)) < 0.05


# ---------------------------------------------------------------- 判据③①
def test_frozen_transfer_fails_then_refit_restores():
    rows = frozen_transfer(seed=3, reps=800)
    worst = min(r["d_frozen"] for r in rows)
    best = max(r["d_refit"] for r in rows)
    assert worst < 0.2, f"frozen transfer should collapse, worst d={worst}"
    assert best > 0.9, f"re-instantiation should restore d, best={best}"


def test_closure_boundary():
    cl = closure_leak(N=2000, reps=800, seed=5)
    by_m = {c["m"]: c for c in cl}
    assert by_m[0.0]["closed"] is True
    assert by_m[0.10]["closed"] is False


# ---------------------------------------------------------------- 裁决
def test_qualify_verdict_structure():
    rows = [{"N": 50, "Ns": 1, "passes": False},
            {"N": 200, "Ns": 4, "passes": True},
            {"N": 20000, "Ns": 400, "passes": True}]
    t3 = [{"N": 50, "passes": False}, {"N": 200, "passes": False},
          {"N": 20000, "passes": True}]
    cl = [{"m": 0.0, "closed": True}, {"m": 0.1, "closed": False}]
    v = qualify(rows, t3, cl)
    kinds = {x["verdict"] for x in v}
    assert kinds == {"insensitive", "context_bound", "macro_in_context",
                     "boundary_leak"}
    leak_row = [x for x in v if x["verdict"] == "boundary_leak"][0]
    assert leak_row["C1"].startswith("✗(m≥0.1")


def test_thresholds_are_preregistered():
    assert TH_D == 0.80 and TH_EFF == 0.20 and DS == 0.02
    assert P0 == 0.10 and REPS == 2000
