#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R47 —— δ-时间化检验（震荡/方向是否时间化算子的产物）结构性测试

对应：experiments/r23_real_cognitive/r47_delta_temporalization.py、r47_report.md、
v1.6.36（快思维本体不可测 Q0 / 时间化不造方向 Q1 / 结构跨算子稳健 Q2）。

锁的结构预言（调底层函数复现 main() 的聚合逻辑，不调用 main 以避免写文件）：
  * Q0：纯 δ 单点（无时间延展）不构成震荡段 ⟹ 判别器对快思维本体**不可测**
    —— 判别需要时间结构，这是"慢思维轨迹"判读的必要条件。
  * Q1：无结构的随机 δ 经时间化**不被凭空造出方向**（方向检出率峰值 ≤ 0.5）
    —— 时间化只展开已有结构，不是方向的来源。
  * Q2：带潜在结构（酝酿/反刍）的 δ 序列在现实窗口（m=2/3, focus=0~0.3）
    跨算子展开方向保持率 ≥ 0.8 ⟹ 判出的方向部分反映真实结构。
  * temporalize 结构：每子步总质量守恒 ≈ T/m。

不锁：具体保持率数值、各格三分分布细节。

诚实边界：合成 δ 事件 + 时间化算子的 S2 检验；结论是"方向非纯时间化伪影"，
不是"快思维本体有方向"（Q0 恰恰登记快思维本体不可测）。
"""

import os
import sys
from collections import Counter

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r23_real_cognitive"))

import r47_delta_temporalization as R  # noqa: E402


def _grid():
    """复现 main() 的 (structure × m × focus) 判别三分网格（固定种子流推进）。"""
    import random
    rng = random.Random(R.SEED)
    grid = {}
    for structure in R.STRUCTURES:
        for m in R.MS:
            for foc in R.FOCUSES:
                verds = []
                for _ in range(30):
                    events = R.make_events(structure, R.C, R.N, rng)
                    traj = R.temporalize(events, m, foc, R.C)
                    v, _, _ = R.judge(traj)
                    verds.append(v if v else "none")
                grid["%s,m%d,f%.1f" % (structure, m, foc)] = dict(Counter(verds))
    return grid


def test_q0_delta_untemporalized_unmeasurable():
    """Q0：纯 δ 单点（无时间延展）无震荡段 ⟹ 快思维本体不可测。"""
    import random
    rng = random.Random(R.SEED)
    n_none = 0
    for _ in range(30):
        R.make_events("random", R.C, R.N, rng)
        v, net, seg = R.judge([{"A": R.T}])   # 单点：无时间结构
        if seg is None:
            n_none += 1
    assert n_none == 30, n_none


def test_q1_temporalization_does_not_create_direction():
    """Q1：随机 δ 经时间化方向检出率峰值 ≤ 0.5（时间化不凭空造方向）。"""
    grid = _grid()
    for m in R.MS:
        for foc in R.FOCUSES:
            c = grid["random,m%d,f%.1f" % (m, foc)]
            tot = sum(c.values())
            dr = (c.get("ruminate_down", 0) + c.get("ruminate_up", 0)) / tot
            assert dr <= 0.5, (m, foc, dr, c)


def test_q2_structure_recovered_across_operators():
    """Q2：现实窗口（m=2/3, focus=0~0.3）酝酿/反刍方向保持率 ≥ 0.8。"""
    grid = _grid()
    for structure in ("ruminate_down", "ruminate_up"):
        for m in (2, 3):
            for foc in (0.0, 0.3):
                c = grid["%s,m%d,f%.1f" % (structure, m, foc)]
                tot = sum(c.values())
                keep = c.get(structure, 0) / tot
                assert keep >= 0.8, (structure, m, foc, keep, c)


def test_temporalize_mass_conservation():
    """temporalize 结构：每子步总质量 ≈ per = T/m（舍入差 ≤ 1）。"""
    events = [0.3 + 0.5 * k / 11 for k in range(12)]
    for m in R.MS:
        subs = R.temporalize(events, m, 0.3, R.C)
        assert len(subs) == len(events) * m
        per = R.T / m
        for sub in subs:
            assert abs(sum(sub.values()) - round(per)) <= 1, (m, sub)
