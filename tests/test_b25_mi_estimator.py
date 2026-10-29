"""B-B25 互信息估计器的测试。

覆盖：Miller-Madow 熵、条件熵、条件互信息三场景（有偏倚/无偏倚/φ 效应）。
"""

import os
import sys

import numpy as np
import pytest

EXP_DIR = os.path.join(os.path.dirname(__file__), "..", "experiments",
                       "b25_mi_estimator")
sys.path.insert(0, EXP_DIR)

from b25_mi_estimator import (  # noqa: E402
    gen_biased,
    gen_phi_effect,
    gen_unbiased,
)
from mi_estimator import (  # noqa: E402
    MIEstimator,
    cond_mutual_information,
    conditional_entropy,
    entropy_mm,
)


# ---------------- Miller-Madow 熵 ----------------

def test_entropy_singleton_zero():
    assert entropy_mm(np.array([3, 3, 3, 3])) == pytest.approx(0.0, abs=1e-9)


def test_entropy_uniform_bit():
    # 均匀 4 类 → log2(4)=2 bit（MM 修正项极小）
    n = 2000
    series = np.tile([0, 1, 2, 3], n // 4)
    assert entropy_mm(series) == pytest.approx(2.0, abs=0.01)


def test_conditional_entropy_known():
    # 按 G 分组：G=0 时 M 均匀2类(熵≈1)，G=1 时 M 单点(熵0) → H≈0.5
    # 大样本令 Miller-Madow 修正可忽略
    M = np.array([0, 1] * 500 + [2] * 1000)
    G = np.array([0] * 1000 + [1] * 1000)
    assert conditional_entropy(M, G) == pytest.approx(0.5, abs=0.05)


# ---------------- 条件互信息三场景 ----------------

def test_mi_biased_by_omega():
    """场景 A：Ω 影响 M 分布 → I(M;Ω|φ) > 0。"""
    M, Omega, phi = gen_biased(n_per=2000)
    imi = cond_mutual_information(M, Omega, phi)
    assert imi > 0.05, f"偏倚应可检出，I={imi:.4f}"


def test_mi_omega_independent():
    """场景 B：Ω 与 M 无关 → I(M;Ω|φ) ≈ 0。"""
    M, Omega, phi = gen_unbiased(n_per=2000)
    imi = cond_mutual_information(M, Omega, phi)
    assert abs(imi) < 0.02, f"应无偏倚，I={imi:.4f}"


def test_mi_phi_effect_omega_irrelevant():
    """场景 C：仅 φ 影响 M、Ω 无关 → 条件化后 I ≈ 0。"""
    M, Omega, phi = gen_phi_effect(n_per=1500)
    imi = cond_mutual_information(M, Omega, phi)
    assert abs(imi) < 0.02, f"条件化后应无偏倚，I={imi:.4f}"


def test_report_components_consistent():
    """报告分量满足 I = H(M|φ) − H(M|φ,Ω)。"""
    M, Omega, phi = gen_biased()
    r = MIEstimator(M, Omega, phi).report()
    assert r["I_M_Omega_phi"] == pytest.approx(
        r["H_M_phi"] - r["H_M_phi_Omega"], abs=1e-9)
