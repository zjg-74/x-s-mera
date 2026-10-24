#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R34 —— rho=1 顶点容差带（OPEN-1）结构性测试（锁结构，不锁数字）

对应：`experiments/r34_vertex_band/`、`27-v3.4-boundary-map.md` §3.4 OPEN-1。

容差带定义：B 级判读禁区 := {|beta_hat| <= z*SE(beta_hat)}（z=1.96，约定性）。

锁结构：
  * 判读三值输出且规则一致：|beta_hat| <= z*SE -> vertex_band，否则按符号二分
  * 无噪声：SE=0，禁区恰为 rho=1（顶点），其余点符号必须正确
  * 有噪声：顶点处留空率高于远端；带宽随噪声增大（rho=1.02 处留空率不减）
  * 判读规范不变（G_lambda 下 (E,U) 逐点不变 -> verdict 不变）

不锁：任何具体留空率数值、任何具体 beta_hat。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r34_vertex_band"))

import r34_vertex_band as R  # noqa: E402


def test_module_imports():
    assert hasattr(R, "verdict")
    assert hasattr(R, "fit_beta_with_se")
    assert hasattr(R, "run_all")


def test_verdict_rule_is_three_valued_and_consistent():
    """规则本身：|beta| <= z*SE -> 禁区；否则 beta<0 -> mid_drill。"""
    assert R.verdict(0.0, 0.0) == "vertex_band"
    assert R.verdict(-1.0, 0.0) == "mid_drill"
    assert R.verdict(1.0, 0.0) == "mid_hollow"
    assert R.verdict(0.05, 1.0) == "vertex_band"        # z*SE = 1.96 > 0.05
    assert R.verdict(-2.5, 1.0) == "mid_drill"
    assert R.verdict(2.5, 1.0) == "mid_hollow"


def test_noiseless_vertex_exact_and_signs_correct():
    """无噪声：SE=0，禁区恰为 rho=1；其余点判读符号必须正确。"""
    res = res_run()
    assert res["summary"]["noiseless_sign_correct"], res["v1_noiseless"]
    assert res["summary"]["noiseless_vertex_exact"], res["v1_noiseless"]
    for row in res["v1_noiseless"]:
        if row["rho"] == 1.0:
            assert row["verdict"] == "vertex_band", row
        else:
            assert row["verdict"] == row["true"], row


def test_band_widens_with_noise():
    """结构预言：噪声增大，靠近顶点（rho=1.02）的留空率不减。"""
    res = res_run()
    assert res["summary"]["band_widens_with_noise"]


def test_vertex_withheld_more_than_far_field():
    """结构预言：有噪声时顶点留空率高于远端（rho=0.90）。"""
    for key in ("sigma_0.05", "sigma_0.20"):
        v2 = res_run()["v2_noisy"][key]["per_rho"]
        near = next(p["withheld_rate"] for p in v2 if p["rho"] == 1.02)
        far = next(p["withheld_rate"] for p in v2 if p["rho"] == 0.90)
        assert near >= far, (key, near, far)


def test_gauge_invariance():
    """U5：G_lambda 下轨迹逐点不变、判读不变。"""
    res = res_run()
    assert res["summary"]["gauge_invariant"], res["v3_gauge"]
    g = res["v3_gauge"]
    assert g["path_pointwise_unchanged"], g
    assert g["verdict_a"] == g["verdict_b"], g


def res_run():
    """缓存 run_all，避免测试间重复计算。"""
    if not hasattr(res_run, "_cache"):
        res_run._cache = R.run_all()
    return res_run._cache
