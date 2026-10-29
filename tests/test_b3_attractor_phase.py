"""B3 四类终态相图测绘的测试。

覆盖：分类器四类判别的正确性 + 双轨网格 smoke（LE 单调稳定 / eco 四类覆盖）。
"""

import json
import os
import sys

import numpy as np
import pytest

EXPERIMENT_DIR = os.path.join(os.path.dirname(__file__), "..", "experiments",
                              "b3_attractor_phase")
SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, EXPERIMENT_DIR)
sys.path.insert(0, SRC_DIR)

from b3_attractor_phase import (  # noqa: E402
    EcoDynamics,
    classify_tail,
    run_eco_trajectory,
    run_le_grid,
)


# ---------------- 分类器单元测试 ----------------

def test_classify_steady():
    # 收敛到非零不动点
    series = np.concatenate([np.linspace(0.5, 0.7, 200), np.full(600, 0.7)])
    assert classify_tail(series) == "steady_state"


def test_classify_extinction():
    series = np.concatenate([np.linspace(0.5, 0.02, 400), np.full(200, 0.01)])
    assert classify_tail(series) == "global_extinction"


def test_classify_oscillation():
    # 正弦持续振荡 → 多峰
    t = np.linspace(0, 20 * 2 * np.pi, 800)
    series = 0.5 + 0.2 * np.sin(t)
    assert classify_tail(series) == "periodic_oscillation"


def test_classify_lock():
    # 低幅非周期波动：峰谷幅 < OSC_PTP（不判振荡），标准差不足以判稳态 → 锁定
    rng = np.random.default_rng(7)
    series = 0.05 + 0.005 * rng.standard_normal(600)
    assert classify_tail(series) == "ecological_lock"


# ---------------- 双轨网格 smoke ----------------

def test_part_A_LE_monostable():
    """LE 内核耦合网格全部收敛稳态（单调稳定结论的网格证据）。"""
    res = run_le_grid(np.linspace(0.3, 3.0, 4), np.linspace(0.05, 0.6, 4),
                      nsteps=800, dt=0.05, seed=20261029)
    assert res["n_grid"] == 16
    assert all(d["cls"] == "steady_state" for d in res["details"]), (
        "LE 内核应单调稳定：全网格稳态")


def test_part_B_four_attractors():
    """扩展生态原型应能覆盖四类终态。"""
    counts = {}
    # 参数含确定振荡点（s=1.0, 高 a → 极限环）；s=0.05 → 灭绝
    for s in [0.05, 0.6, 1.0, 1.6]:
        for a in [0.5, 3.0, 6.0, 8.0]:
            r = run_eco_trajectory(float(s), float(a), nsteps=6000,
                                   dt=0.02, seed=20261029)
            counts[r["cls"]] = counts.get(r["cls"], 0) + 1
    expected = {"steady_state", "periodic_oscillation",
                "ecological_lock", "global_extinction"}
    assert expected.issubset(set(counts.keys())), (
        f"四类终态应全部覆盖，实际 {counts}")


def test_eco_model_bounds():
    """eco 模型状态保持非负且有界。"""
    m = EcoDynamics(s=1.0, a=3.0)
    for _ in range(200):
        r, c = m.step(0.02)
        assert r >= 0.0 and c >= 0.0
    assert r < 10.0 and c < 10.0


# ---------------- results.json 存在性 ----------------

def test_results_json_present():
    """标准运行应产出 results.json（本测试假定已先运行脚本）。"""
    p = os.path.join(EXPERIMENT_DIR, "results.json")
    if os.path.exists(p):
        data = json.load(open(p, encoding="utf-8"))
        assert "part_A_LE_kernel" in data and "part_B_eco_prototype" in data
