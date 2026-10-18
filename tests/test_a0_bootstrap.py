"""
R7 方案甲（回归收敛）的单元测试
================================

锚点：解析极限 d_∞ = b/(1−g)（无噪声）、Var = σ²/(1−g²)（有噪声）。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.a0_bootstrap_convergence import (
    ProjectionChain, run_chain, fixed_point_iteration, classify_system,
    SYSTEM_MODELS,
)


def test_analytic_residual_bias():
    """无噪声时尾段张力应收敛到 b/(1−g)。"""
    g, b = 0.80, 0.05
    pc = ProjectionChain(gain=g, layer_error=b, noise_std=0.0, seed=0)
    r = pc.run(200)
    assert r["tail_mean"] == pytest.approx(b / (1 - g), rel=0.05)


def test_analytic_residual_variance():
    """有噪声时张力 d=|ε| 的稳态统计应对照解析式。

    口径注意：ε 是零均值高斯 AR(1)，Var(ε)=σ²/(1−g²)；
    而 d=|ε| 服从半正态分布，E[d]=σ_d·√(2/π)、E[d²]=Var(ε)。
    因此 tail_mean 对照的是半正态均值，二阶矩对照的才是方差本身。
    """
    g, s = 0.70, 0.10
    pc = ProjectionChain(gain=g, layer_error=0.0, noise_std=s, seed=0)
    r = pc.run(2000)
    var_theo = s ** 2 / (1 - g ** 2)
    mean_theo = np.sqrt(var_theo) * np.sqrt(2.0 / np.pi)   # 半正态 E|ε|
    assert abs(r["tail_mean"] - mean_theo) / mean_theo < 0.10
    assert abs(r["tail_sq_mean"] - var_theo) / var_theo < 0.10


def test_shrinkage_monotone_in_gain():
    """三个读法的核心事实：g 越接近 1，尾段张力越大（单调）。"""
    res = run_chain([0.70, 0.90, 1.00, 1.03], n_layers=60, layer_error=0.02,
                    noise_std=0.10, seeds=[0, 1, 2])
    tails = [r["tail_mean"] for r in res["rows"]]
    assert tails == sorted(tails), f"非单调：{tails}"


def test_verdict_classification():
    """C-1/C-2/C-3 三档判定边界正确。"""
    res = run_chain([0.70, 1.00, 1.03], n_layers=40, layer_error=0.02,
                    noise_std=0.10, seeds=[0])
    v = {r["gain"]: r["verdict"] for r in res["rows"]}
    assert "C-1" in v[0.70]
    assert "C-2" in v[1.00]
    assert "C-3" in v[1.03]


def test_fixed_point_ratio_equals_lipschitz():
    """自我回归的观测收敛比应等于算子的局部 Lipschitz 常数。"""
    for g in (0.60, 0.85):
        r = fixed_point_iteration(lambda x, g=g: (0.3 + g * x, g), x0=0.0, n_iters=400)
        assert r["converged"] is True
        assert r["observed_ratio_median"] == pytest.approx(g, abs=0.02)


def test_fixed_point_expansive_diverges():
    """扩张算子（g>1）不收敛。"""
    r = fixed_point_iteration(lambda x: (0.3 + 1.03 * x, 1.03), x0=0.0, n_iters=300)
    assert r["converged"] is False
    assert r["final_drift"] > 1.0


def test_four_class_classification():
    """R7 四类系统各自的收缩性判定与 R7 的定性结论一致。"""
    c1 = classify_system("class1_syntax", n_layers=100)
    c3 = classify_system("class3_llm", n_layers=100)
    c4 = classify_system("class4_evol", n_layers=100)
    assert "C-2" in c1["class"]        # 保距：不衰减
    assert "C-3" in c3["class"]        # 自举：发散
    assert "C-1" in c4["class"]        # 选择：强收缩


def test_strong_shrinkage_but_wrong_fixed_point():
    """方案甲成立也不能辩护 A0：类4 强收缩，但不动点是"适应峰"。

    这与 R7 的定性结论（丙为主、乙为辅）一致：收敛 ≠ 回归 A0。
    """
    c4 = classify_system("class4_evol", n_layers=200)
    assert "C-1" in c4["class"]
    # 若它真的回归 A0（偏差为 0），尾段应≈0；实际留下常驻偏置
    assert c4["tail_mean"] > 0.05


def test_reproducible():
    a = ProjectionChain(gain=0.9, layer_error=0.02, noise_std=0.1, seed=7).run(80)["d"]
    b = ProjectionChain(gain=0.9, layer_error=0.02, noise_std=0.1, seed=7).run(80)["d"]
    assert np.allclose(a, b)
