#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R44 —— B-B44 取值原则锚定验证（落带区间内 η 的可靠性 / 误读率）结构性测试

对应：experiments/r44_eta_principle/（r44_eta_principle.py + r44_eta_principle.json + report.md）、
27-v3.4-boundary-map.md §4.2 / §5（B-B44 闭合、R-043 落带区间、R-044 锚定）、33-v3.4-eta-scaling-principle.md。

锁的结构预言（用底层函数复现 main() 的聚合逻辑，不调用 main 以避免写文件）：
  * τ* 关于 η 单调不增（η 放宽 ⟹ 分离尺度 τ* 不增）。
  * 在 R-043 登记的落带区间内取 η，三场景各自存在非空"可靠子区间"——
    即存在 η 使 τ* 落带率（P(τ*∈(τ_fine,τ_coarse))）≥ 0.90 ⟹ 取值原则可锚定
    （"可接受错误率"语义被兑现，η 可从约定根变为有根）。

不锁：具体可靠子区间端点、具体落带率数值（只锁"区间内存在可取值的 η"这一结构结论）。

诚实边界：固定种子可复现的纯内部仿真；M 取较小值仅为加速，不影响"区间非空"结论；
B 跨尺度弱场景可靠子区间需收窄到落带区间边界以内（场景相关声明）。
"""

import os
import sys
import random

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r44_eta_principle"))

import r44_eta_principle as R  # noqa: E402


def test_tau_star_monotone_nonincreasing_in_eta():
    """η 放宽 ⟹ τ* 不增（扫描从细到粗，首次满足即返回，随 η 单调前移）。"""
    kind = "same_scale"
    _, dU = R.make_trajectory(kind, 0.0, R.SEED)
    ratios = R.ratio_curve(dU)
    etas = [R.BANDS[kind][0] + (R.BANDS[kind][1] - R.BANDS[kind][0]) * i / 39
            for i in range(40)]
    seq = [R.tau_star(ratios, eta) for eta in etas]
    assert all(b <= a for a, b in zip(seq, seq[1:])), seq


def test_reliable_subinterval_nonempty():
    """R44 核心结论：三场景下，R-043 落带区间内存在 η 使落带率 ≥ 0.90（可靠子区间非空）。"""
    sig = 0.05
    m_traj = 80
    for kind in R.BANDS:
        lo, hi = R.BANDS[kind]
        etas = [lo + (hi - lo) * i / 39 for i in range(40)]
        curves = []
        for m in range(m_traj):
            seed = R.SEED + m * 131 + int(sig * 1000) * 17 + sum(ord(c) for c in kind)
            _, dU = R.make_trajectory(kind, sig, seed)
            curves.append(R.ratio_curve(dU))
        rels = []
        for eta in etas:
            n_band = sum(1 for rc in curves
                         if R.TAU_FINE < R.tau_star(rc, eta) < R.TAU_COARSE)
            rels.append(n_band / len(curves))
        # 可靠子区间非空：区间内至少一个 η 的落带率 ≥ 0.90
        assert max(rels) >= 0.90, (kind, round(max(rels), 3))
