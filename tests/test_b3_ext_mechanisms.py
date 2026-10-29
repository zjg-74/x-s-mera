"""B3-ext：LE 内核四类终态机制扩展的测试。

覆盖：默认内核向后兼容（稳态）、四类代表机制分类、机制分岔扫描、非法机制键。
"""

import os
import sys

import numpy as np
import pytest

SRC = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, SRC)

from le_coupling_kernel import LECouplingKernel  # noqa: E402

T = 2000
DT = 0.1
WARM = 600
OSC = {"oscillation_gain": 1.0, "half_override": 0.05,
       "uptake_override": 2.5, "env_r": 2.5,
       "consumer_mortality": 0.5, "conversion": 0.8}


def run(mechs, t=T, dt=DT):
    k = LECouplingKernel(mechanisms=mechs, seed=20261029)
    for _ in range(int(t / dt)):
        k.step(k.t, dt)
    return [s[0] for s in k.history['L_state']]


def classify(bio, warm=WARM):
    a = np.asarray(bio[warm:], dtype=float)
    if max(a[-100:]) < 0.05:
        return "extinction"
    if max(a[-150:]) - min(a[-150:]) > 0.10:
        return "oscillation"
    if float(np.mean(a[-150:])) < 0.5:
        return "lock"
    return "steady"


# ---------------- 默认内核向后兼容 ----------------

def test_default_kernel_still_monotone_stable():
    """无机制：保持 v3.4-le.1 单调稳定（四类不可达，稳态）。"""
    bio = run(None)
    assert classify(bio) == "steady"
    # 末期无持续波动（唯一稳态）
    assert max(bio[WARM:]) - min(bio[WARM:]) < 0.02


# ---------------- 四类代表机制 ----------------

def test_steady_reachable():
    assert classify(run({})) == "steady"


def test_extinction_reachable():
    assert classify(run({"death_threshold": 0.90, "death_rate": 1.5})) \
        == "extinction"


def test_lock_reachable():
    assert classify(run({"allee_threshold": 1.00, "allee_strength": 1.2})) \
        == "lock"


def test_oscillation_reachable():
    assert classify(run(OSC)) == "oscillation"


# ---------------- 机制分岔 ----------------

def test_extinction_threshold_bifurcation():
    assert classify(run({"death_threshold": 0.75, "death_rate": 1.5})) \
        == "steady"
    assert classify(run({"death_threshold": 0.90, "death_rate": 1.5})) \
        == "extinction"


def test_lock_threshold_bifurcation():
    assert classify(run({"allee_threshold": 0.90, "allee_strength": 1.2})) \
        == "steady"
    assert classify(run({"allee_threshold": 1.10, "allee_strength": 1.2})) \
        == "lock"


def test_oscillation_hopf_bifurcation():
    conv05 = dict(OSC, conversion=0.5)
    assert classify(run(conv05)) == "steady"
    assert classify(run(OSC)) == "oscillation"


# ---------------- 校验 ----------------

def test_unknown_mechanism_key_rejected():
    with pytest.raises(ValueError):
        LECouplingKernel(mechanisms={"bogus_key": 1})
