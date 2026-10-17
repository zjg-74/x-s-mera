"""R6 单元测试：A10 可观测判据表、粒度扫描、三联排除。

重点不是"判据能不能触发"，而是三件更容易出错的事：
  1. 预注册表有没有被偷偷调参（阈值快照）；
  2. 标度律斜率在已知过程上是否给出正确符号（白噪声 −0.5 / 1f ≈ 0）；
  3. 三联排除的**顺序**是否 load-bearing（同一条信号在两种顺序下裁决不同）。
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.a10_observables import (  # noqa: E402
    BETA_MAX, COMMIT_MAX, CRITERIA, DELTA_ABS, DELTA_QUOTED,
    block_average, classify, evaluate, evaluate_persistent,
    granularity_scan, repair_encapsulated, window_stats,
)


# ---------------------------------------------------------------- 基础算子
def test_block_average_shape_and_value():
    x = np.arange(10.0)
    assert np.allclose(block_average(x, 2), [0.5, 2.5, 4.5, 6.5, 8.5])
    # 不能整除时截断，不补零
    assert block_average(np.arange(7.0), 3).size == 2
    assert np.isclose(block_average(np.arange(7.0), 3)[0], 1.0)


def test_declared_scale_is_exogenous_constant():
    """Δ_quoted 必须是写死的常数，不能由 D 估计出来。"""
    assert DELTA_QUOTED == 0.01
    assert np.isclose(DELTA_ABS, 3.0 * DELTA_QUOTED)


def test_registry_is_frozen():
    """预注册表快照：阈值若被改动，此处必须同时改，否则视为违规调参。"""
    assert [c.code for c in CRITERIA] == ["P2", "P3", "P4", "P5"]
    snap = {c.code: (c.statistic, c.op, round(c.threshold, 6), c.gate) for c in CRITERIA}
    assert snap["P2"] == ("bc", "ge", round(5 / 9, 6), "C1")
    assert snap["P3"] == ("tail_abs", "ge", 0.05, "C1")
    assert snap["P4"] == ("commit", "ge", COMMIT_MAX, "C1")
    assert snap["P5"] == ("beta", "ge", BETA_MAX, "C1")


# ---------------------------------------------------------------- 标度律
def _white(n, seed=0):
    return np.random.default_rng(seed).normal(0, 0.01, n)


def _one_over_f(n, seed=0):
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, d=1.0)
    f[0] = f[1]
    band = f <= 0.45          # 不设下截止：下截止本身就是个特征尺度
    amp = np.zeros_like(f)
    amp[band] = 1.0 / np.sqrt(f[band])
    x = np.fft.irfft(amp * np.exp(1j * rng.uniform(0, 2 * np.pi, f.size)), n)
    return x / x.std() * 0.01


def test_white_noise_slope_is_minus_half():
    sc = granularity_scan(_white(8000), max_level=7)
    assert -0.62 < sc["beta"] < -0.38, f"white noise beta should be ~-0.5, got {sc['beta']}"


def test_one_over_f_slope_is_flat():
    sc = granularity_scan(_one_over_f(8000), max_level=7)
    assert sc["beta"] > BETA_MAX, f"1/f should be scale-free, got beta={sc['beta']}"


def test_scan_amplitude_decreases_for_white_increases_not_for_1f():
    w = granularity_scan(_white(8000), max_level=6)
    p = granularity_scan(_one_over_f(8000), max_level=6)
    assert w["a"][-1] < 0.35 * w["a"][0]
    assert p["a"][-1] > 0.5 * p["a"][0]


# ---------------------------------------------------------------- 判据裁决
def test_gate_closed_makes_criteria_not_evaluable():
    st = window_stats(_white(500))
    ev = evaluate(st, c1_fraction=0.0)
    assert ev["gate_open"] is False
    assert all(v is None for v in ev["fires"].values())


def test_location_invariance_blind_spot():
    """P2–P5 全为位置不变量：整体平移不得改变任何裁决（登记为盲区 G2）。"""
    x = _white(2000, seed=5)
    st0, st1 = window_stats(x), window_stats(x + 0.06)
    ev0, ev1 = evaluate(st0), evaluate(st1)
    assert ev0["fires"] == ev1["fires"]
    assert abs(st0["median"] - st1["median"] + 0.06) < 1e-9


def test_commit_ratio_scales_with_spread():
    tight = window_stats(np.random.default_rng(1).normal(0, DELTA_QUOTED, 3000))
    wide = window_stats(np.random.default_rng(1).normal(0, 5 * DELTA_QUOTED, 3000))
    assert wide["commit"] > 3 * tight["commit"]


def test_persistence_needs_two_of_three():
    rng = np.random.default_rng(3)
    # 三段等长，保证窗口边界与构造段对齐（否则窗口跨段会稀释信号）
    ok = rng.normal(0.9, DELTA_QUOTED, 2000)
    bad = rng.normal(0.9, 6 * DELTA_QUOTED, 2000)
    series = np.concatenate([ok, bad, ok])          # 只有中间窗口坏
    c1 = np.ones_like(series, bool)
    p = evaluate_persistent(series, c1, n_windows=3, min_fire=2)

    series2 = np.concatenate([ok, bad, bad])        # 两个窗口坏
    p2 = evaluate_persistent(series2, c1, n_windows=3, min_fire=2)
    assert p["fires"]["P4"] is False
    assert p2["fires"]["P4"] is True


def test_repair_check_ignores_beta():
    """修复窗口太短，估不出可靠 beta；不能因为算不准就判失败。"""
    x = np.random.default_rng(7).normal(0.9, DELTA_QUOTED, 200)
    r = repair_encapsulated(x)
    assert r["encapsulated"] is True
    assert r["fires"]["P5"] is None


# ---------------------------------------------------------------- 三联排除
def _persist(any_fire=True, gate=True):
    return {"fires": {"P2": any_fire}, "counts": {}, "gate_open": gate,
            "c1_fraction": 1.0 if gate else 0.0, "windows": [],
            "any_fire": any_fire and gate}


def test_classify_gate_excludes_f1():
    v = classify(_persist(any_fire=True, gate=False), {"encapsulated": False},
                 {"beta": -0.5})
    assert v["verdict"] == "gate"


def test_classify_no_signal_is_ok():
    v = classify(_persist(any_fire=False), {"encapsulated": True}, {"beta": -0.5})
    assert v["verdict"] == "ok"


def test_classify_reef_resolution_collapse():
    reef = classify(_persist(), {"encapsulated": True}, {"beta": 0.02})
    res = classify(_persist(), {"encapsulated": False}, {"beta": -0.5})
    col = classify(_persist(), {"encapsulated": False}, {"beta": 0.02})
    assert reef["verdict"] == "reef"
    assert res["verdict"] == "resolution"
    assert col["verdict"] == "collapse"


def test_ladder_order_is_load_bearing():
    """接口可修 + 标度律衰减：先修接口判礁石，先看标度律判分辨率伪影。

    这不是排版偏好——系统性接口伪影在标度律上可能与有限尺度噪声同形。
    """
    p, rep, scan = _persist(), {"encapsulated": True}, {"beta": -0.55}
    a = classify(p, rep, scan, order="repair_first")["verdict"]
    b = classify(p, rep, scan, order="granularity_first")["verdict"]
    assert a == "reef" and b == "resolution" and a != b


def test_scan_requires_two_decades():
    """粒度范围不足时 beta 不可用，而不是给个数硬判。"""
    sc = granularity_scan(_white(60), max_level=2)
    assert np.isfinite(sc["beta"]) or np.isnan(sc["beta"])
    assert len(sc["b"]) == 3
