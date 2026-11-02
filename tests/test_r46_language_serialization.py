#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R46 —— 语言化序列化算子 + 判别方向稳定性（B 线 P4 方法论）结构性测试

对应：experiments/r23_real_cognitive/r46_language_serialization.py、r46_report.md、
v1.6.35（污染 = 抹平而非翻转 + 现实参数方向保持 100%）。

锁的结构预言（调底层函数复现 main() 的聚合逻辑，不调用 main 以避免写文件）：
  * 原始（未语言化）判别：已知方向的合成轨迹**全部判对**
    （down→ruminate_down / up→ruminate_up）——ground truth 可用性。
  * 核心反证命题：语言化算子（时间细分 m × 聚焦 focus）**不翻转方向**——
    模式 A（纯细分）与模式 B（聚焦锯齿）的全部格子翻转率 = 0。
  * 现实采集窗口（m∈{2,3}, focus∈{0,0.3}）方向保持率 ≥ 0.9。
  * subdivide_keep 结构：每子步质量 ≈ 总量/m，且不注入新的相对分布偏置。

不锁：具体保持率数值（只锁判定阈值）、flat 出现在哪些病态格子。

诚实边界：合成轨迹上的方法论自检（S2），锁的是"判别器对语言化污染的稳健性结构"，
不锁"真人 think-aloud 也如此"（真人数据归约属 B 线 P4 正题）。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r23_real_cognitive"))

import r46_language_serialization as R  # noqa: E402


def _trajs():
    down = [R.make_traj(i, "ruminate_down") for i in range(R.N_TRAJ)]
    up = [R.make_traj(i, "ruminate_up") for i in range(R.N_TRAJ)]
    return [("down", t) for t in down] + [("up", t) for t in up]


def _expected(dirn):
    return "ruminate_down" if dirn == "down" else "ruminate_up"


def _grid_stats(mode):
    """复现 main() 的格子统计：{key: {down:..., up:...}}。"""
    import random
    rng = random.Random(R.SEED)
    stat = {}
    for m in R.MS:
        for foc in R.FOCUSES:
            key = "%s,m%d,f%.1f" % (mode, m, foc)
            st = {"down": {"keep": 0, "flip": 0, "flat": 0},
                  "up": {"keep": 0, "flip": 0, "flat": 0}}
            for dirn, t in _trajs():
                seq = []
                for w in t:
                    if mode == "keep":
                        seq.extend(R.subdivide_keep(w, m))
                    else:
                        seq.extend(R.subdivide_focus(w, m, foc, rng))
                v, _, _ = R.judge(seq)
                if v == _expected(dirn):
                    st[dirn]["keep"] += 1
                elif v is None or v == "flat":
                    st[dirn]["flat"] += 1
                else:
                    st[dirn]["flip"] += 1
            stat[key] = st
    return stat


def test_original_verdict_all_correct():
    """ground truth：未语言化的已知方向轨迹全部判对。"""
    for dirn, t in _trajs():
        v, net, seg = R.judge(t)
        assert v == _expected(dirn), (dirn, v)


def test_no_flip_any_cell_both_modes():
    """核心反证命题：语言化不翻转方向（污染 = 抹平 flat，非翻转）。"""
    for mode in ("keep", "focus"):
        for key, st in _grid_stats(mode).items():
            total_flip = st["down"]["flip"] + st["up"]["flip"]
            assert total_flip == 0, (mode, key, st)


def test_realistic_window_keep_rate_at_least_09():
    """现实采集窗口（m=2/3, focus=0~0.3）聚焦模式下方向保持率 ≥ 0.9。"""
    stats = _grid_stats("focus")
    for m in (2, 3):
        for foc in (0.0, 0.3):
            key = "focus,m%d,f%.1f" % (m, foc)
            st = stats[key]
            keep = (st["down"]["keep"] + st["up"]["keep"]) / (2 * R.N_TRAJ)
            assert keep >= 0.9, (key, keep)


def test_subdivide_keep_mass_and_no_injection():
    """模式 A 结构：每子步总质量 ≈ per = 总量/m；比例结构与原权重同序。"""
    w = {"A": 6, "B": 2, "C": 1, "D": 1}
    total = sum(w.values())
    for m in (1, 2, 3, 5):
        subs = R.subdivide_keep(w, m)
        assert len(subs) == m
        per = total / m
        for sub in subs:
            assert abs(sum(sub.values()) - round(per)) <= 1, (m, sub)
        # 相对结构与原权重同序（不注入新偏置）
        order_orig = sorted(w, key=w.get, reverse=True)
        for sub in subs:
            order_sub = sorted(sub, key=sub.get, reverse=True)
            assert order_sub == order_orig, (m, order_sub, order_orig)
