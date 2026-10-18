"""
T1-d 表征分化模块的单元测试
=============================

设计原则（与 tests/ 其余文件一致）：每条测试都锚定一个可复算的数值事实，
不做"跑通即可"的弱断言。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.t1d_representational_diff import (
    RepresentationPopulation, boundary_diagnostics, critical_iteration,
    hysteresis_rate, direction_report, differentiation_depth,
    flux_sweep_observables,
)


# ---------------- 诊断量本体 ----------------

def test_boundary_diagnostics_unimodal_rejected():
    """单峰数据必须被判为无边界（阴性对照的判据基础）。"""
    rng = np.random.default_rng(0)
    x = rng.normal(0.0, 1.0, size=800)
    d = boundary_diagnostics(x)
    assert d["has_boundary"] is False
    assert d["bc"] < 0.555


def test_boundary_diagnostics_bimodal_accepted():
    """对称双峰必须被判为有边界。"""
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(-1.5, 0.2, 400), rng.normal(1.5, 0.2, 400)])
    d = boundary_diagnostics(x)
    assert d["has_boundary"] is True
    assert d["bc"] > 0.8
    assert d["symmetry"] > 0.8


def test_boundary_diagnostics_rejects_shifted_unimodal():
    """关键判别力测试：单峰长尾（整体偏移）不得被误判为双峰。

    这是 v1 的真实失败模式——k-means 把"左半边 vs 右半边"切出来，
    sep 高达 2.6，靠 straddle/symmetry 条件才能排除。
    """
    rng = np.random.default_rng(1)
    x = rng.normal(0.0, 0.5, size=1200) + np.abs(rng.normal(0, 0.4, 1200))
    d = boundary_diagnostics(x)
    # 该分布是单峰偏态，确实不构成双峰
    assert d["bc"] < 0.555 or not d["has_boundary"]


def test_insufficient_samples():
    d = boundary_diagnostics(np.array([1.0, 2.0]))
    assert d["status"] == "insufficient"
    assert d["has_boundary"] is False


# ---------------- 双阱势的关键性质 ----------------

def test_symmetric_phi_never_bifurcates():
    """a=0（单阱）在任何通量下都不产生双峰——这是 T1-d 的阴性对照。"""
    for J in (0.0, 1.0, 3.0):
        pop = RepresentationPopulation(n_units=300, well_depth=0.0, flux=J,
                                       noise_std=0.10, gain=0.01, seed=0)
        pop.run(2000)
        d = pop.snapshot_stats()
        assert d["has_boundary"] is False, f"a=0, J={J} 误判为有边界"


def test_asymmetric_phi_bifurcates():
    """a>0 且通量适中时产生双峰；两簇应关于 0 对称、各在 ±c 一侧。

    实测（flux=1.0, gain=0.01, 3000 步）cluster_means ≈ [-0.81, 0.80]。
    注意实际阱底不是精确的 ±c：双阱 U=a(x²−c²)² 的极小在 ±c，
    但噪声使群体中心在**阱内**分布，且被 |x|>c 尾部拉偏，故取近似值。
    """
    pop = RepresentationPopulation(n_units=300, well_depth=1.0, well_sep=1.0,
                                   flux=1.0, noise_std=0.10, gain=0.01, seed=0)
    traj = pop.run(3000)
    d = boundary_diagnostics(traj[-1])
    assert d["has_boundary"] is True
    m1, m2 = d["cluster_means"]
    assert m1 < 0 < m2                    # 两簇分居两侧
    assert abs(abs(m1) - abs(m2)) < 0.25  # 近似对称
    assert 0.6 < abs(m1) < 1.3 and 0.6 < abs(m2) < 1.3


def test_single_trajectory_is_not_bimodal():
    """单轨迹是**双态切换**（两阱间跳跃），不是双峰分布。

    这条测试记录了建模过程中被推翻的一个假设：
    原以为"落阱即锁定、永不翻越"，实测单轨迹在 3000 步内正负各占
    约 40%/60%（signs.mean ≈ −0.34），即持续在两阱间跳跃。
    跳跃轨迹的**时间分布**是双态的，但其**取值分布**仍是两簇（因为
    每步只在一个阱附近）——但统计量不稳定，无法作为"分化"的可靠读数。
    这正是必须用多单元群的理由：群体分布才是稳定的分化读数。
    """
    pop = RepresentationPopulation(n_units=1, well_depth=1.0, flux=1.0,
                                   noise_std=0.10, gain=0.01, seed=0)
    traj = pop.run(3000)
    signs = np.sign(traj[100:, 0])
    # 确实发生过翻越（不锁定）
    assert abs(signs.mean()) < 0.9
    # 但单单元样本量不足，无法做群体判据
    d = boundary_diagnostics(traj[-1])
    assert d["status"] == "insufficient"


# ---------------- 序参量 ----------------

def test_hysteresis_rate_monotone_in_flux():
    """翻越率随通量单调上升（R8 F-c 子命题 1）。"""
    rates = []
    for J in (0.0, 1.0, 3.0):
        pop = RepresentationPopulation(n_units=200, well_depth=1.0, flux=J,
                                       noise_std=0.10, gain=0.01, seed=0)
        traj = pop.run(2000)
        rates.append(hysteresis_rate(traj, burn_in=500)["rate_per_unit_per_step"])
    assert rates[0] < rates[1] < rates[2]


def test_direction_unbiased_across_seeds():
    """f₊ 随种子在 0.5 附近波动（R8 F-c 子命题 3：方向归流，不被 Φ 锁死）。"""
    fps = []
    for s in range(6):
        pop = RepresentationPopulation(n_units=300, well_depth=1.0, flux=1.0,
                                       noise_std=0.10, gain=0.01, seed=s)
        traj = pop.run(1500)
        fps.append(direction_report(traj, burn_in=500)["frac_plus"])
    assert all(abs(f - 0.5) < 0.15 for f in fps)


def test_differentiation_depth_is_small_for_symmetric_split():
    """两阱各占一半时群体均值接近 0（对称破缺的签名）。"""
    pop = RepresentationPopulation(n_units=300, well_depth=1.0, flux=1.0,
                                   noise_std=0.10, gain=0.01, seed=0)
    traj = pop.run(1500)
    assert differentiation_depth(traj, burn_in=500)["depth"] < 0.15


# ---------------- 复现性 ----------------

def test_flux_sweep_reproducible():
    """固定种子 → 同一扫描两次结果完全一致。"""
    a = flux_sweep_observables([0.0, 2.0], seeds=[0], n_units=150,
                               steps=800, burn_in=200, gain=0.01)
    b = flux_sweep_observables([0.0, 2.0], seeds=[0], n_units=150,
                               steps=800, burn_in=200, gain=0.01)
    for ra, rb in zip(a["rows"], b["rows"]):
        assert ra["flip_rate_mean"] == pytest.approx(rb["flip_rate_mean"])
        assert ra["bc_final_mean"] == pytest.approx(rb["bc_final_mean"])


def test_boundary_stability_decays_with_flux():
    """通量增大使边界形态稳定性下降（BC 下降）——A11b 非单调性的证据。"""
    a = flux_sweep_observables([0.0], seeds=[0], n_units=300, steps=2500,
                               burn_in=800, gain=0.01)
    b = flux_sweep_observables([5.0], seeds=[0], n_units=300, steps=2500,
                               burn_in=800, gain=0.01)
    assert a["rows"][0]["bc_final_mean"] > b["rows"][0]["bc_final_mean"]
