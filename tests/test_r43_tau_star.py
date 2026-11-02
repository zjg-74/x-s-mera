#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R-043 —— τ* 落点（加密 η 网格 + 观测窗口扫描）结构性测试

对应：experiments/r43_tau_star/、reef_shturl.json R-043、
27-v3.4-boundary-map.md §5（R-043 已解决 v1.6.23）、B-B44 降级登记。

锁的结构预言（用底层函数复现 main() 的聚合逻辑，不调用 main 以避免写文件）：
  * τ*(η) = min{w : ratio(w) < η} 关于 η 单调不增（定义的直接推论）。
  * **加密 η 网格（200 点）后，三类场景在默认窗口 (8,128) 下均有非空 τ* 落带 η 区间**
    —— 这是 R-043 的核心结论：此前"几乎为空"是 R-042 粗网格假象。
  * 窗口扫描对每场景至少找到一个使落带区间非空的 (τ_fine, τ_coarse)。
  * valid_eta_intervals 把布尔落带序列压成连续区间（结构：段与段之间必有空隙）。

不锁：具体区间端点数值、τ* 的具体范围。

诚实边界：纯内部仿真（SIG=0 确定性轨迹），固定种子逐位可复现；
结论是"区间非空存在性"，不是"η 的现实合理值"（后者由 R44 锚定取值原则承接）。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r43_tau_star"))

import r43_tau_star as R  # noqa: E402

SCENES = ("same_scale", "cross_scale_weak", "structural")


def _etas():
    return [R.ETA_MIN + (R.ETA_MAX - R.ETA_MIN) * i / (R.ETA_PTS - 1)
            for i in range(R.ETA_PTS)]


def _curves():
    return {kind: R.ratio_curve(R.make_trajectory(kind)[1]) for kind in SCENES}


def test_tau_star_monotone_nonincreasing_in_eta():
    """η 放宽 ⟹ 首次满足阈值的位置不增（扫描定义的必然结构）。"""
    ratios = _curves()["same_scale"]
    etas = _etas()
    seq = [R.tau_star_from_ratio(ratios, e) for e in etas]
    assert all(b <= a for a, b in zip(seq, seq[1:])), seq[:20]


def test_default_window_band_nonempty_all_scenes():
    """R-043 核心结论：加密网格后三场景在默认窗口 (8,128) 均有非空落带 η 区间。"""
    etas = _etas()
    for kind, ratios in _curves().items():
        tstars = [R.tau_star_from_ratio(ratios, e) for e in etas]
        flags = [R.in_band(t, R.TAU_FINE_DEF, R.TAU_COARSE_DEF) for t in tstars]
        segs = R.valid_eta_intervals(etas, flags)
        assert segs, (kind, "default window has no band")
        lo, hi = segs[0][0], segs[0][-1]
        assert 0.0 < lo <= hi <= R.ETA_MAX


def test_window_scan_finds_nonempty_for_each_scene():
    """窗口扫描兜底：每场景至少存在一个扫描窗口使落带区间非空。"""
    etas = _etas()
    for kind, ratios in _curves().items():
        tstars = [R.tau_star_from_ratio(ratios, e) for e in etas]
        found = False
        for tf in (4, 8, 16, 32):
            for tc in (64, 128, 256, 512):
                if tc <= tf:
                    continue
                flags = [R.in_band(t, tf, tc) for t in tstars]
                if R.valid_eta_intervals(etas, flags):
                    found = True
                    break
            if found:
                break
        assert found, (kind, "no window yields a band")


def test_interval_compressor_segments_are_disjoint_and_ordered():
    """valid_eta_intervals 结构：True 段与段之间必有 False 空隙（区间不重叠）。"""
    etas = [0.1 * i for i in range(10)]
    flags = [False, True, True, False, False, True, False, True, True, True]
    segs = R.valid_eta_intervals(etas, flags)
    # 三段：位置 1–2 / 5 / 7–9
    assert len(segs) == 3, segs
    assert segs[0] == (etas[1], etas[2])
    assert segs[1] == (etas[5], etas[5])
    assert segs[2] == (etas[7], etas[9])
    # 区间按序且两两不重叠（段间必有空隙）
    for (a1, b1), (a2, b2) in zip(segs, segs[1:]):
        assert b1 < a2
    # 空白序列 ⟹ 空区间表
    assert R.valid_eta_intervals(etas, [False] * 10) == []
