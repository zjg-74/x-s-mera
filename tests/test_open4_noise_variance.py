#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""OPEN-4 结构性测试：有噪声下共线性法 vs 峰距法的鲁棒性预言（锁结构，不锁数字）。"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, "..", "experiments", "open4_noise_variance")
sys.path.insert(0, EXP)

from open4_noise_variance import (  # noqa: E402
    K, SIGMA_REL, NOISE_MODES, N_TRAIL,
    make_cfg, trajectory, add_noise, estimate, rng_for,
    stats,
)


def _run(rho, sched, sigma, mode):
    cfg = make_cfg(1.0, rho)
    rho_true = cfg["rho_exp"]
    Es, Us = trajectory(cfg, sched)
    std_E = (sum((v - sum(Es) / len(Es)) ** 2 for v in Es) / len(Es)) ** 0.5
    std_U = (sum((v - sum(Us) / len(Us)) ** 2 for v in Us) / len(Us)) ** 0.5
    rng = rng_for("open4t:%s:%s:%s:%s" % (sched, rho, sigma, mode))
    ec, ep, ok = [], [], 0
    for _ in range(N_TRAIL):
        sigE = sigma * std_E if mode in ("A_only", "both") else 0.0
        sigU = sigma * std_U if mode in ("B_only", "both") else 0.0
        E_obs = add_noise(Es, sigE, rng)
        U_obs = add_noise(Us, sigU, rng)
        r = estimate(E_obs, U_obs, rho_true, rho)
        ec.append(r["err_coll"]); ep.append(r["err_peak"])
        ok += 1 if r["sign_ok"] else 0
    return stats(ec)[0], stats(ep)[0], ok / float(N_TRAIL)


def test_no_noise_collinearity_exact():
    """无噪声：共线性法斜率误差 ~0（与 R30 一致）。"""
    err, _, _ = _run(2.0, "log_ramp", 0.0, "both")
    assert err < 1e-10, err


def test_a_noise_degrades_collinearity():
    """A 级（E）噪声破坏共线性斜率（差分放大）。"""
    e0, _, _ = _run(2.0, "linear", 0.0, "A_only")
    e1, _, _ = _run(2.0, "linear", 0.05, "A_only")
    assert e0 < 1e-6, e0
    assert e1 > 0.3, e1          # 大噪声下斜率判读明显退化


def test_b_noise_keeps_collinearity():
    """B 级（U）噪声下共线性斜率相对稳健。"""
    e1, _, _ = _run(2.0, "linear", 0.05, "B_only")
    assert e1 < 0.1, e1          # U 是回归 y，加方差不放大


def test_a_worse_than_b_for_collinearity():
    """同噪声档下，A(E) 对共线性法的破坏 > B(U)。"""
    ea, _, _ = _run(2.0, "linear", 0.05, "A_only")
    eb, _, _ = _run(2.0, "linear", 0.05, "B_only")
    assert ea > 5.0 * eb, (ea, eb)


def test_peak_order_small_noise_stable():
    """峰序（B 级判读）在小噪声下保持。"""
    _, _, ok = _run(2.0, "log_ramp", 0.001, "both")
    assert ok == 1.0, ok


def test_peak_distance_inherent_bias_persists():
    """峰距法固有偏差：无噪声下非对数调度已有显著偏差（§8.2）。"""
    _, ep, _ = _run(2.0, "const_blocks", 0.0, "both")
    assert ep > 0.05, ep          # 阶跃块无噪声偏差 ~20%
