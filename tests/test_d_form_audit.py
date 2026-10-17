"""Unit tests for src/d_form_audit.py (R2 support module)."""
import numpy as np

from src.d_form_audit import (
    synthesize, core_width, z_kurtosis, tail_mass, bimodality_coefficient,
    robust_sigma, spearman,
)


def test_synthesize_core_scale_follows_tau():
    """核心宽度应随 tau 以 1/tau 缩放（固定 coverage/paths 时）。"""
    w1 = core_width(synthesize(tau=1.0, coverage=0.99, n_paths=1,
                               a_unmodeled=100.0, n=40000, seed=1)["r"])
    w2 = core_width(synthesize(tau=10.0, coverage=0.99, n_paths=1,
                               a_unmodeled=100.0, n=40000, seed=1)["r"])
    ratio = w1 / w2
    assert 5.0 < ratio < 20.0, f"core width should shrink ~10x, ratio={ratio:.2f}"


def test_synthesize_tail_mass_follows_coverage():
    """未建模余项占比应随 coverage 下降（绝对标度阈值，远大于核心）。"""
    lo = tail_mass(synthesize(tau=10.0, coverage=0.5, n_paths=1,
                              a_unmodeled=1.0, n=40000, seed=2)["r"], tol=0.5)
    hi = tail_mass(synthesize(tau=10.0, coverage=0.95, n_paths=1,
                              a_unmodeled=1.0, n=40000, seed=2)["r"], tol=0.5)
    assert lo > hi + 0.2, f"tail mass must drop with coverage: {lo:.3f} vs {hi:.3f}"


def test_route_multiplicity_creates_modes():
    """n_paths>=2 应产生可见模态（BC 升高或双峰形态）。"""
    single = synthesize(tau=10.0, coverage=0.99, n_paths=1,
                        a_unmodeled=1.0, n=40000, seed=3)
    dual = synthesize(tau=10.0, coverage=0.99, n_paths=2,
                      a_unmodeled=1.0, a_route=1.0, n=40000, seed=3)
    # 双路径下分布的质量应明显离开 0（|r| 中位数更大）
    assert np.median(np.abs(dual["r"])) > 5 * np.median(np.abs(single["r"])) - 1e-9


def test_robust_sigma_recovers_core_scale():
    """robust_sigma 应接近 1/tau（核心主导时）。"""
    m = synthesize(tau=50.0, coverage=0.999, n_paths=1,
                   a_unmodeled=100.0, n=60000, seed=4)
    s = robust_sigma(m["r"])
    assert 0.5 / 50.0 < s < 2.0 / 50.0, f"robust sigma {s:.4f} vs 1/tau {1/50:.4f}"


def test_z_kurtosis_positive_for_heavy_tail():
    """单路径 + 低覆盖度 -> 重尾 -> 超额峰度显著为正。"""
    m = synthesize(tau=10.0, coverage=0.7, n_paths=1,
                   a_unmodeled=3.0, n=60000, seed=5)
    assert z_kurtosis(m["r"]) > 1.0


def test_bimodality_coefficient_sane():
    x = np.concatenate([np.random.default_rng(0).normal(-3, .5, 500),
                        np.random.default_rng(1).normal(3, .5, 500)])
    bc = bimodality_coefficient(x)
    assert bc > 0.555, f"clear bimodal sample should exceed 0.555, got {bc:.3f}"


def test_spearman_recovers_monotone():
    x = np.linspace(0, 1, 50)
    assert spearman(x, x ** 3) > 0.99
    assert abs(spearman(x, -x) + 1.0) < 1e-6


def test_overclaim_ratio_matches_design():
    """承诺兑现比 = 核心宽度/申报精度，正态核心下应约 1.35（IQR/sigma）。"""
    m = synthesize(tau=100.0, coverage=0.999, n_paths=1,
                   a_unmodeled=1000.0, n=60000, seed=6)
    ratio = core_width(m["r"]) * 100.0     # core / (1/tau)
    assert 1.0 < ratio < 1.8, f"IQR/sigma for normal core should be ~1.35, got {ratio:.2f}"


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and getattr(fn, "__name__", "").startswith("test_"):
            fn()
            print(f"[ok] {fn.__name__}")
    print("All d_form_audit tests passed.")
