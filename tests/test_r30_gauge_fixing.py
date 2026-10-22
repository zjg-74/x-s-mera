#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R30 —— A_D 规范固定的结构性测试（锁结构，不锁数字）

对应：`experiments/r30_gauge_fixing/`、审计文档 `docs/xs-framework/14-v3.2-source-refactor-audit.md`。

核心断言：补篇五 §1.2 的 A_D 参数化含冗余自由度；规范固定后等价于一条仿射关系

    dU/dt = rho_exp * (dE/dt) + kappa_E * rho_exp * (1 - rho)

本文件只锁**结构性**预言：
  * 恒等式逐样本成立（残差为浮点噪声，不是某个特定数值）
  * λ 规范变换下 (s, rho_exp, rho, kappa_E) 不变，而 (E, U) 逐点不变
  * 由仿射关系反解 rho_exp **不依赖调度族**；峰距法**依赖**（相对大小关系，不锁具体百分比）
  * 峰序只需截距符号，不需 kappa_E

不锁：任何具体百分比、任何具体残差数值、任何具体峰位。
"""

import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r30_gauge_fixing"))

import r30_gauge_fixing_identity as R  # noqa: E402

# ------------------------------------------------------------------ 基本构造


def test_module_imports():
    assert hasattr(R, "trajectory")
    assert hasattr(R, "make_cfg")
    assert hasattr(R, "sched")


def test_gauge_invariant_definitions_consistent():
    """rho = rho_res/rho_exp = theta_U/theta_E 必须恒成立（恒等式的前提）。"""
    for th in R.THETA_E_GRID:
        for rho in R.RHO_GRID:
            cfg = R.make_cfg(th, rho)
            assert abs(cfg["rho"] - cfg["rho_res"] / cfg["rho_exp"]) < 1e-12
            assert abs(cfg["rho"] - cfg["theta_U"] / cfg["theta_E"]) < 1e-12


# ------------------------------------------------------------------ V1 恒等式


def test_affine_identity_holds_pathwise():
    """(*) dU/dt = rho_exp*dE/dt + kappa_E*rho_exp*(1-rho)，逐样本。"""
    for r in R.v1_identity():
        assert r["max_abs_resid"] < 1e-12, r


def test_affine_identity_residual_is_float_noise_not_model_error():
    """残差应远小于 |dU/dt| 本身的量级（相对 < 1e-12），否则就是模型错。"""
    for r in R.v1_identity():
        assert r["scale"] > 0.0
        assert r["max_abs_resid"] / r["scale"] < 1e-12, r


def test_affine_identity_independent_of_theta_E_and_rho():
    """残差量级不应随 theta_E / rho 系统性增长（结构性，不是数值巧合）。"""
    rows = R.v1_identity()
    for th in R.THETA_E_GRID:
        sub = [r["max_abs_resid"] for r in rows if r["theta_E"] == th]
        assert max(sub) < 1e-12


def test_affine_identity_independent_of_schedule_family():
    """跨调度族都成立 ⟹ 不是某一族的巧合。"""
    rows = R.v1_identity()
    for name in R.SCHED_FAMILIES:
        sub = [r["max_abs_resid"] for r in rows if r["sched"] == name]
        assert max(sub) < 1e-12


def test_intercept_matches_closed_form():
    """截距必须等于 kappa_E * rho_exp * (1 - rho)，而不是拟合出来的某数。"""
    for th in (0.5, 1.0, 2.0):
        for rho in R.RHO_GRID:
            cfg = R.make_cfg(th, rho)
            Es, Us, dEs, dUs, ss = R.trajectory(cfg, "log_ramp")
            a, b = R._ols(dEs, dUs)
            assert abs(b - cfg["rho_exp"]) < 1e-12
            assert abs(a - cfg["kappa_E"] * cfg["rho_exp"] * (1.0 - rho)) < 1e-12


# ------------------------------------------------------------------ V2 规范不变性


def test_gauge_transform_leaves_trajectory_pointwise_invariant():
    """λ 变换后 (E, U) 逐点不变——这是定理 6.12 的内容，不是新主张。"""
    for r in R.v2_gauge():
        assert r["max_abs_E_diff"] == 0.0, r
        assert r["max_abs_U_diff"] == 0.0, r


def test_gauge_transform_changes_theta_E_but_not_invariants():
    """θ_E 必须真的变了（否则测试无意义），而不变量不能变。"""
    rows = R.v2_gauge()
    thetas = [r["theta_E_transformed"] for r in rows]
    assert max(thetas) / min(thetas) > 8.0          # 至少跨一个数量级
    for r in rows:
        assert r["max_abs_s_diff"] == 0.0, r
    assert len(set(round(r["rho_exp"], 12) for r in rows)) == 1
    assert len(set(round(r["rho"], 12) for r in rows)) == 1
    assert len(set(round(r["kappa_E"], 12) for r in rows)) == 1


def test_normalized_input_ratio_s_is_gauge_invariant():
    """s := r/theta_E 是规范不变量（这是规范固定的定义性要求）。"""
    cfg = R.make_cfg(1.0, 2.0)
    _, _, _, _, s0 = R.trajectory(cfg, "log_ramp", lam=1.0)
    for lam in R.LAM_GRID:
        _, _, _, _, sL = R.trajectory(cfg, "log_ramp", lam=lam)
        assert max(abs(a - b) for a, b in zip(sL, s0)) == 0.0


# ------------------------------------------------------------------ V3 调度无关性


def test_collinearity_estimator_is_schedule_independent():
    """共线性反解 rho_exp 的误差跨调度族恒定（全部极小）。"""
    rows = R.v3_sched_independence()
    for name in R.SCHED_FAMILIES:
        sub = [r["rel_err_collinearity"] for r in rows if r["sched"] == name]
        assert max(sub) < 1e-10, (name, max(sub))


def test_peak_distance_estimator_is_schedule_dependent():
    """峰距法必须在非 log_ramp 族上显著劣化（这是 B-B34 的结构复现）。"""
    rows = R.v3_sched_independence()
    med = {}
    for name in R.SCHED_FAMILIES:
        sub = sorted(r["rel_err_peakdistance"] for r in rows
                     if r["sched"] == name and r["rho"] != 1.0)
        med[name] = sub[len(sub) // 2]
    assert med["log_ramp"] < 0.05                  # 推论 6.5 的适用族：准
    for name in ("linear", "sigmoid"):
        assert med[name] > 0.10, (name, med[name])  # 非适用族：显著偏


def test_collinearity_strictly_beats_peak_distance_off_family():
    """在非适用调度族上，共线性法必须严格优于峰距法（这是换可观测量正当性的根据）。"""
    rows = R.v3_sched_independence()
    for name in ("linear", "sigmoid"):
        sub = [r for r in rows if r["sched"] == name and r["rho"] != 1.0]
        worst_col = max(r["rel_err_collinearity"] for r in sub)
        best_peak = min(r["rel_err_peakdistance"] for r in sub)
        assert worst_col < best_peak, (name, worst_col, best_peak)


def test_peak_order_only_needs_intercept_sign():
    """峰序（rho>1 ⟺ U 峰在前）可由截距符号读出，不需要 kappa_E。"""
    for r in R.v3_sched_independence():
        if r["rho"] == 1.0:
            continue
        assert r["peak_order_sign_ok"], r


def test_peak_order_sign_schedule_independent():
    """峰序判据跨调度族都对 ⟹ 它不依赖 κ_E，与理论一致。"""
    rows = R.v3_sched_independence()
    for name in R.SCHED_FAMILIES:
        sub = [r for r in rows if r["sched"] == name and r["rho"] != 1.0]
        assert all(r["peak_order_sign_ok"] for r in sub), name


# ------------------------------------------------------------------ 参数计数（结构性）


def test_gauge_fixed_parameterization_has_no_redundancy():
    """规范固定后：3 个无量纲不变量 + 1 个单位；λ 变换不改变其中任何一个。

    这是「B-B23 的主体属于规范冗余而非真不可辨识」的直接检验。
    """
    cfg = R.make_cfg(1.0, 2.0)
    base = (cfg["rho_exp"], cfg["rho"], cfg["kappa_E"])
    for lam in R.LAM_GRID:
        # λ 变换 = (r, e_+, u_+) -> (λr, e_+/λ, u_+/λ)
        e_p = cfg["e_plus"] / lam
        u_p = cfg["u_plus"] / lam
        e_m = cfg["e_minus"]
        u_m = cfg["u_minus"]
        th_E_L = e_m / e_p
        rho_exp_L = u_p / e_p
        rho_L = (u_m / e_m) / rho_exp_L
        kappa_L = 1.0 * e_m            # b == 1
        assert abs(rho_exp_L - base[0]) < 1e-12
        assert abs(rho_L - base[1]) < 1e-12
        assert abs(kappa_L - base[2]) < 1e-12
        # 但 θ_E 确实变了
        assert abs(th_E_L - lam * cfg["theta_E"]) < 1e-12


def test_schedule_families_are_distinct():
    """四个调度族必须真的不同，否则 V3 的结论无从谈起。"""
    mid = [R.sched(n, R.K // 2) for n in R.SCHED_FAMILIES]
    assert len(set(round(v, 9) for v in mid)) >= 3


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
