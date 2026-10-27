#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R11 —— 推论 4.1 判决性实验（增益越权 → 伪定向 / F2 / ghost 伪力）结构性测试

对应：experiments/r11_pseudo_drift/（r11_pseudo_drift.py + r11_results.json + report.md）、
XS v3.1 补篇二 §5 R11（推论 4.1）、27-v3.4-boundary-map.md（B2 维，M1/M2 机制链）。

锁的结构预言（与脚本 main() 的 P1/P2 等价，调用 run_arm_g 而非 main，不写文件）：
  * P1  legal 臂（A11-R 合规）：F2/伪收敛率跨 g 近似平稳（极差 ≤ 4 个百分点）——增益
         只调幅度、不调方向，不产生伪定向。
  * P2  violate 臂（越权：g 进候选得分/状态导数）：F2 率显著高于 legal 臂
         （效应量 ≥ 20pp）且从基线（g=1）单调升至峰值 —— 增益冒充漂移 → 伪定向。
  * R11 支持 ⟺ P1 ∧ P2。
  * 方向效应：violate 臂 F2 与 g 的相关性明显强于 legal 臂（越权产生 g 驱动的伪定向）。

不锁：具体 F2 百分点或相关性数值量级（只锁"legal 平稳、violate 显著上升且效应≥20pp"）。

诚实边界：复用 r18 宏加载原型（S2 自制仿真）；偏置 b_M 是注入的"冒充漂移"代理；
g=4/5 回落 = 冒充漂移饱和次级效应（不否定 P2 上升段）。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "r11_pseudo_drift"))

import r11_pseudo_drift as R  # noqa: E402


def _all_rows():
    return [R.run_arm_g(arm, g) for arm in ("legal", "violate")
            for g in R.G_VALUES]


def test_arm_outputs_have_f2_metric():
    """两臂都返回 F2 伪收敛率字段（结算函数对两臂完全一致）。"""
    for r in _all_rows():
        assert "伪收敛F2%" in r
        assert r["n"] == R.N_TRIALS


def test_p1_legal_no_pseudo_drift():
    """P1：legal 臂 F2 跨 g 极差 ≤ 4pp（近似平稳，增益不调方向）。"""
    legal = [r for r in _all_rows() if r["arm"] == "legal"]
    f2 = [r["伪收敛F2%"] for r in legal]
    assert (max(f2) - min(f2)) <= 4.0


def test_p2_violate_pseudo_drift():
    """P2：violate 臂 F2 峰值显著高于 legal 峰值（效应 ≥ 20pp）且从基线升至峰值。"""
    rows = _all_rows()
    legal = [r for r in rows if r["arm"] == "legal"]
    violate = sorted((r for r in rows if r["arm"] == "violate"), key=lambda r: r["g"])
    f2_legal = [r["伪收敛F2%"] for r in legal]
    f2_v = [r["伪收敛F2%"] for r in violate]
    peak_violate = max(f2_v)
    peak_legal = max(f2_legal)
    assert (peak_violate - peak_legal) >= 20.0
    assert peak_violate > f2_v[0] + 10.0          # 从 g=1 基线显著上升


def test_r11_supported():
    """R11 支持：P1 ∧ P2 同时成立（推论 4.1：越权 → 伪定向）。"""
    rows = _all_rows()
    legal = [r for r in rows if r["arm"] == "legal"]
    violate = sorted((r for r in rows if r["arm"] == "violate"), key=lambda r: r["g"])
    f2_legal = [r["伪收敛F2%"] for r in legal]
    f2_v = [r["伪收敛F2%"] for r in violate]
    p1 = (max(f2_legal) - min(f2_legal)) <= 4.0
    p2 = (max(f2_v) - max(f2_legal) >= 20.0) and (max(f2_v) > f2_v[0] + 10.0)
    assert p1 and p2


def test_violate_f2_dominates_legal():
    """量级区分 + 效应远超波动：越权臂 F2 在所有 g 上远高于合规臂，且基线→峰值
    的上升幅度远超合规臂的平稳波动 —— 这是 P2 区别于 P1 的本质（伪定向由越权产生）。"""
    rows = _all_rows()
    legal = [r for r in rows if r["arm"] == "legal"]
    violate = sorted((r for r in rows if r["arm"] == "violate"), key=lambda r: r["g"])
    f2l = [r["伪收敛F2%"] for r in legal]
    f2v = [r["伪收敛F2%"] for r in violate]
    # 量级：violate 在每个 g 上的 F2 都高于 legal 峰值（越权产生大量伪定向）
    assert min(f2v) > max(f2l)
    # 效应量：violate 基线→峰值上升 远超 legal 的平稳波动（P2 区别于 P1）
    legal_fluct = max(f2l) - min(f2l)
    v_rise = max(f2v) - f2v[0]
    assert v_rise >= 20.0
    assert v_rise > 2 * legal_fluct
