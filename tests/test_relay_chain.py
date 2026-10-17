"""
Relay chain — unit tests
========================

检验命题 2′ 中继部分的三件事：
1. 串联中继的复合律是调和/短板律，不是乘积律
2. F1 段阻断整链；F2 段主模态照传但注入 ghost 污染
3. 中继数下限（min_relay_count）存在且随端距/弱端尺度增大；足够的
   中继链可从 F1 吸收态中恢复耦合

运行：pytest tests/test_relay_chain.py -v
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from src.relay_chain import (  # noqa: E402
    leg_overlap, min_relay_count, relay_bridge_demo, simulate_transport)

Z_L = np.array([1.30, 0.20, 0.00])
Z_E = np.array([0.21, 0.25, 0.00])


# ---------------------------------------------------------------- 复合律
def test_single_leg_transmission_equals_overlap():
    r = simulate_transport([0.9])
    assert r["T_measured"] == pytest.approx(0.9, abs=1e-3)


def test_series_follows_harmonic_not_product():
    r = simulate_transport([0.9, 0.9])
    assert r["T_measured"] == pytest.approx(r["T_harmonic"], abs=1e-3)
    assert r["T_harmonic"] == pytest.approx(0.45, abs=1e-6)
    assert abs(r["T_measured"] - r["T_product"]) > 0.3, "乘积律（0.81）必须被否定"


def test_series_unequal_legs_weak_leg_dominates():
    r = simulate_transport([0.95, 0.7])
    assert r["T_measured"] == pytest.approx(0.403, abs=2e-3)
    assert r["T_measured"] < min(0.95, 0.7), "串联传输受最弱段支配"


def test_product_law_falsified_across_lengths():
    d = 0.9
    for n in (2, 3, 4):
        r = simulate_transport([d] * n)
        assert r["T_measured"] == pytest.approx(d / n, abs=5e-3)      # 调和
        assert abs(r["T_measured"] - d ** n) > 0.1                    # 乘积远离


# ---------------------------------------------------------------- 失效模式
def test_F1_leg_blocks_chain():
    r = simulate_transport([0.9, 0.1])          # 0.1 < 0.25 → C1 失败
    assert r["T_measured"] < 1e-6
    assert r["T_harmonic"] == 0.0


def test_F2_leg_transmits_primary_but_contaminates():
    r_ok = simulate_transport([0.9, 0.9])
    r_f2 = simulate_transport([0.9, 0.9], leg_c2=[True, False])
    assert r_f2["T_measured"] == pytest.approx(r_ok["T_measured"], abs=1e-3), \
        "F2 段不应阻断主模态传输"
    assert r_f2["contamination"] > 0.05, "F2 段必须注入 ghost 污染"
    assert r_ok["contamination"] < 1e-6
    assert r_f2["contamination"] == pytest.approx(0.25 / 2.6, abs=0.01), \
        "稳态污染由 ghost 底噪项与两侧副模态电导分配决定：0.25·c / (c + 1.6c) = 0.25/2.6"


# ---------------------------------------------------------------- 中继数下限与恢复
def test_min_relay_count_for_collapsed_pair():
    assert min_relay_count(Z_L, Z_E) == 3


def test_min_relay_count_grows_with_gap_and_weak_end():
    close = min_relay_count(np.array([0.80, 0.22, 0.0]), Z_E)
    far = min_relay_count(Z_L, Z_E)
    assert close == 1
    assert far == 3
    assert far > close


def test_relay_recovery_succeeds_with_sufficient_relays():
    res = relay_bridge_demo(use_relay=True, n_relay=3)
    shrink = 1.0 - res["gap"][-1] / res["gap"][0]
    assert shrink > 0.7, f"足够中继应大幅收缩端距，实际收缩 {shrink:.2f}"
    assert res["d_end"][-1] > 0.9


def test_insufficient_relay_fails_like_no_relay():
    res_none = relay_bridge_demo(use_relay=False)
    res_n1 = relay_bridge_demo(use_relay=True, n_relay=1)
    assert res_n1["gap"][-1] > res_n1["gap"][0], "n=1 不足：端距应继续扩大"
    assert res_none["gap"][-1] > res_none["gap"][0]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"[pass] {fn.__name__}")
    print(f"\nAll {len(fns)} relay-chain tests passed.")
