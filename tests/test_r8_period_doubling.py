#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R8 —— M1 机制工程验证（周期倍增级联 → 混沌）结构性测试（锁结构，不锁数字）

对应：experiments/r8_gain_stepchain/（r8_period_doubling.py + r8_results.json + report.md）、
XS v3.1 补篇一 §5 R8（M1 机制）、27-v3.4-boundary-map.md §4.2 B3。

锁的结构预言（与脚本 main() 的 V1–V4 等价，用底层函数复现，不调用 main 以避免写文件）：
  * V1  至少两次周期倍增（级联 1→2→4→8 中至少 1→2、2→4 两段存在）
  * V2  混沌区存在（λ>0 的 g 占一段带宽）
  * V3  混沌阈值 g_chaos ∈ [3.4, 3.7]（logistic 解析界 3.5699 附近）
  * V4  势不进状态导数（g 仅作乘性系数，iterate 中无 g 的梯度项）

不锁：具体分岔点数字、具体 Lyapunov 数值、具体周期带边界（只锁"级联与混沌出现"）。

诚实边界（代码层可复现的数学事实，作为'有牙齿'独立校验，不依赖网格）：
  * g=4.0 必混沌（λ>0 且无非平凡短周期）；
  * g=2.5 必周期 1。
"""

import os
import sys
import math

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r8_gain_stepchain"))

import r8_period_doubling as R  # noqa: E402


def test_iterate_is_logistic_form():
    """iterate 是 logistic 型 x_{k+1}=g·x_k(1−x_k)，g 仅乘性（不进状态导数）。"""
    seq = R.iterate(3.2, 0.3, 5)
    assert len(seq) == 5
    v = 0.3
    for s in seq:
        assert abs(s - v) < 1e-12
        v = 3.2 * v * (1.0 - v)


def test_iterate_bounded_in_unit_interval():
    """logistic 在 g∈[2.5,4]、x0∈(0,1) 下有界：所有迭代值 ∈ [0,1)。"""
    for g in (2.5, 3.0, 3.5699, 4.0):
        for x in R.iterate(g, 0.3, 500):
            assert 0.0 <= x < 1.0


def test_lyapunov_positive_at_g4():
    """g=4.0 是数学上必混沌：λ>0；稳态无非平凡短周期（detect_period 返回 None）。"""
    assert R.lyapunov(4.0, 0.3) > R.LYAP_TOL
    seq = R.iterate(4.0, 0.3, 4000)[2000:]
    assert R.detect_period(seq) is None or R.detect_period(seq) > R.MAX_PERIOD


def test_period_one_at_g25():
    """g=2.5 落在周期 1 窗口：detect_period 返回 1。"""
    seq = R.iterate(2.5, 0.3, 4000)[2000:]
    assert R.detect_period(seq) == 1


def _coarse_scan(n_g=400):
    """复现 main() 的扫描逻辑（g 网格 → 周期 + Lyapunov），不写文件。"""
    gs = [R.G_LO + (R.G_HI - R.G_LO) * i / (n_g - 1) for i in range(n_g)]
    periods, lams = [], []
    for g in gs:
        seq = R.iterate(g, R.X0, R.BURN_IN + R.SEQ_LEN)[R.BURN_IN:]
        periods.append(R.detect_period(seq))
        lams.append(R.lyapunov(g, R.X0))
    return gs, periods, lams


def test_v1_period_doubling_cascade():
    """V1：出现完整的 2 倍周期倍增链 1→2→4→8（至少三段中的两段，实际四段齐备）。"""
    gs, periods, lams = _coarse_scan()
    unique = sorted({int(p) for p in periods if p is not None})
    for p in (1, 2, 4, 8):
        assert p in unique, (p, unique)
    # 至少两次 p→2p 跨越
    doubling = sum(1 for p in unique if 2 * p in unique)
    assert doubling >= 2


def test_v2_chaos_band_exists():
    """V2：混沌区存在 —— λ>0 的 g 形成一段带宽（计数足够多）。"""
    gs, periods, lams = _coarse_scan()
    chaos_idx = [i for i, v in enumerate(lams) if v > R.LYAP_TOL]
    assert len(chaos_idx) >= 20


def test_v3_chaos_threshold_near_logistic():
    """V3：混沌阈值（首个 λ>0 的 g）落在 logistic 解析界附近 [3.4, 3.7]。"""
    gs, periods, lams = _coarse_scan()
    chaos_idx = [i for i, v in enumerate(lams) if v > R.LYAP_TOL]
    g_chaos = gs[chaos_idx[0]]
    assert 3.4 <= g_chaos <= 3.7


def test_m1_mechanism_supported():
    """M1 机制被工程验证：V1∧V2∧V3∧V4 全真（V4 由 iterate 形式锁，见首条测试）。"""
    gs, periods, lams = _coarse_scan()
    unique = sorted({int(p) for p in periods if p is not None})
    doubling = sum(1 for p in unique if 2 * p in unique)
    chaos_idx = [i for i, v in enumerate(lams) if v > R.LYAP_TOL]
    g_chaos = gs[chaos_idx[0]]
    v1 = doubling >= 2
    v2 = len(chaos_idx) >= 20
    v3 = 3.4 <= g_chaos <= 3.7
    v4 = True  # g 仅乘性（test_iterate_is_logistic_form 锁）
    assert v1 and v2 and v3 and v4
