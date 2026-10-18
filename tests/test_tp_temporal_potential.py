"""
ΔT（时间势差）/ TP 谓词操作化的单元测试
==========================================

锚点：Adler 精确解 Ω_drift = √(Δω²−4K²)（解锁区）/ 0（锁定区）；
预注册判据 TP_CRITERIA（P1–P5）。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.tp_temporal_potential import (
    measure_natural_frequency, delta_T, adler_drift, coupled_phase_experiment,
    dt_zero_channel_sweep, predicate_independence_test, TP_CRITERIA,
)


def test_criteria_preregistered():
    assert "预注册" in TP_CRITERIA["note"]
    assert set(TP_CRITERIA) >= {"P1_freq_recovery", "P2_monotonicity",
                                "P3_adler_consistency", "P4_independence",
                                "P5_negative_control"}


def test_p1_frequency_recoverable():
    """P1：自由运转相位斜率恢复内禀频率，误差 < 1%。"""
    for w in (0.8, 1.3, 2.1):
        m = measure_natural_frequency(w)
        assert m["rel_err"] < 0.01, f"ω={w} 恢复误差 {m['rel_err']:.4%}"


def test_delta_T_definition():
    """ΔT 的定义性质：对称、同频为零、有界 <1。"""
    assert delta_T(1.0, 1.0) == 0.0
    assert delta_T(1.0, 2.0) == pytest.approx(delta_T(2.0, 1.0))
    assert delta_T(1.0, 2.0) == pytest.approx(1.0 / 3.0)
    assert 0.0 <= delta_T(0.1, 5.0) < 1.0


def test_adler_closed_form():
    """Adler 精确解的边界：阈值内锁定，阈值外平方根 onset。"""
    assert adler_drift(0.05, 0.15) == 0.0            # |Δω|=0.05 < 2K=0.30
    assert adler_drift(0.5, 0.10) == pytest.approx(np.sqrt(0.25 - 0.04))
    assert adler_drift(0.0, 0.0) == 0.0


def test_p3_adler_consistency_unlocked():
    """P3：解锁区实测漂移率与 Adler 精确解一致（误差 < 15%）。"""
    for rr in (1.5, 2.0):
        r = coupled_phase_experiment(1.0, rr, coupling=0.05, seed=0)
        assert not r["adler_locked"]
        err = abs(r["drift_rate"] - r["adler_predict"]) / r["adler_predict"]
        assert err < 0.15, f"ratio={rr} err={err:.3f}"


def test_locking_below_threshold():
    """锁定区：|Δω| ≤ 2K 时漂移消失（势差被接口吸收）。"""
    r = coupled_phase_experiment(1.0, 1.05, coupling=0.15, seed=0)
    assert r["adler_locked"] is True
    assert r["locked"] is True
    assert r["drift_rate"] < 0.01


def test_p5_no_coupling_no_locking():
    """P5 阴性对照：K=0 时漂移 = Δω（误差 < 5%），任意小势差都外显。"""
    r = coupled_phase_experiment(1.0, 1.5, coupling=0.0, seed=0)
    assert abs(r["drift_rate"] - 0.5) / 0.5 < 0.05
    assert r["locked"] is False
    # 无耦合时无锁定区：Δω 只要大于漂移分辨率下限（lock_tol）就必外显
    r2 = coupled_phase_experiment(1.0, 1.05, coupling=0.0, seed=0)
    assert r2["locked"] is False
    assert r2["drift_rate"] == pytest.approx(0.05, rel=0.1)


def test_p2_sweep_monotonicity_and_locking():
    """P2：扫描中漂移对 ΔT 单调（容限内），且锁定只出现在最小 ΔT 端。"""
    sw = dt_zero_channel_sweep()
    assert sw["monotone_with_tolerance"] is True
    assert sw["adler_max_rel_err"] < 0.15
    assert sw["n_locked"] >= 1
    # 锁定行必须是 ΔT 最小的行
    rows = sw["rows"]
    locked_idx = [i for i, r in enumerate(rows) if r["locked"]]
    min_dt_idx = int(np.argmin([r["delta_T"] for r in rows]))
    assert locked_idx == [min_dt_idx]


def test_p4_predicate_independence():
    """P4：固定状态、改变频率比 → d 不变、ΔT 变（谓词独立）。"""
    pi = predicate_independence_test()
    assert pi["d_cv"] < 0.05
    assert pi["dt_cv"] > 0.30
    assert pi["independent"] is True


def test_reproducible():
    a = coupled_phase_experiment(1.0, 1.5, seed=7)["drift_rate"]
    b = coupled_phase_experiment(1.0, 1.5, seed=7)["drift_rate"]
    assert a == pytest.approx(b)
