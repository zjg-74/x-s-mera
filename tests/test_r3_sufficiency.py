"""
R3-v2（C1+C2 充分性补案）的单元测试
====================================

锚点：预注册判据 SUFFICIENCY_CRITERIA——单次合规率 100% ∧ 序列成功率 0% → 强反例。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.r3_sufficiency_probe import (
    S1_cumulative_damage, S2_competitive_displacement, S3_temporal_mismatch,
    S4_mediator_contamination, S5_potential_depletion, run_full_sweep,
    SUFFICIENCY_CRITERIA,
)


def test_criteria_preregistered():
    """判据必须先于仿真写定。"""
    assert "100%" in SUFFICIENCY_CRITERIA["strong_counterexample"]
    assert "预注册" in SUFFICIENCY_CRITERIA["note"]


def test_s1_strong_counterexample():
    """S1 累积损伤：单次全合规、序列全失败。"""
    r = S1_cumulative_damage(n_trials=100)
    assert r["single_pass_rate"] == 1.0
    assert r["series_success"] == 0.0
    assert r["verdict"].startswith("强反例")
    # 损伤确实越过了阈值（不是判据设置错误）
    assert r["peak_damage_mean"] >= 0.60


def test_s2_strong_counterexample():
    """S2 竞争置换：挤出率 > 收益率 → 净存量下降。"""
    r = S2_competitive_displacement(n_trials=100)
    assert r["series_success"] == 0.0
    assert r["final_stock_mean"] < 1.0
    assert r["verdict"].startswith("强反例")


def test_s3_incremental_reading_strong():
    """S3 时序错配：增量读法下单次全过、末态越界（强反例）。"""
    r = S3_temporal_mismatch(n_trials=100)
    assert r["single_pass_rate"] == 1.0
    assert r["series_success"] == 0.0
    assert r["final_drift_mean"] > r["tolerance"]
    assert r["verdict"].startswith("强反例")


def test_s3_instant_reading_alarms():
    """S3 的对照事实：同一构造按瞬时读法会在中途报警。

    这是 S3 的核心产出——C2 的读法必须显式化为瞬时读法（绝对基准）。
    """
    r = S3_temporal_mismatch(n_trials=100)
    assert r["instant_reading_alarms"] is True


def test_s3_small_drift_no_counterexample():
    """S3 的参数敏感性：漂移足够小时累积不越界 → 无反例（构造不是恒假）。"""
    r = S3_temporal_mismatch(drift_rate=0.004, n_trials=100)
    assert r["series_success"] == 1.0
    assert r["verdict"] == "无反例"


def test_s4_strong_counterexample():
    """S4 中介污染：失配分量逐段放大，总污染越界。"""
    r = S4_mediator_contamination(n_trials=100)
    assert r["series_success"] == 0.0
    assert r["total_contamination_mean"] > r["tolerance"]
    assert r["verdict"].startswith("强反例")


def test_s5_strong_counterexample():
    """S5 势耗竭：抽取率 > 再生率 → Φ 跌破下限。"""
    r = S5_potential_depletion(n_trials=100)
    assert r["series_success"] == 0.0
    assert r["final_phi_mean"] < 0.20 + 0.05   # 允许噪声余量
    assert r["verdict"].startswith("强反例")


def test_full_sweep_all_strong():
    """全谱扫描：五类构造应全部判为强反例。"""
    res = run_full_sweep(n_trials=100)
    assert res["n_strong"] == 5
    assert res["n_invalid"] == 0
    assert "准入判据" in res["summary"]


def test_invalid_construction_detected():
    """若单次判据自身不过，判定应报'构造不合法'而非反例。"""
    # drift 大到单步越界 → 增量读法下单次合规率 < 100%
    r = S3_temporal_mismatch(drift_rate=0.05, step_tol=0.02, n_trials=50)
    assert r["single_pass_rate"] == 0.0
    assert r["verdict"].startswith("构造不合法")


def test_reproducible():
    a = S5_potential_depletion(n_trials=50, seed=9)["final_phi_mean"]
    b = S5_potential_depletion(n_trials=50, seed=9)["final_phi_mean"]
    assert a == pytest.approx(b)
