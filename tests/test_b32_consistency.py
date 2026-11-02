#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B-B32 —— 震荡两判据一致性（形状 vs 概率）结构性测试

对应：experiments/b32_shape_prob_consistency/b32_consistency.py、
08-v3.1-supp7 缺口 B-B32、v1.6.38（一致率 0.75 不成立 ⟹ 互补需合并覆盖）。

锁的结构预言（调底层函数与 main()，不落盘——main() 只返回 dict，写文件在 __main__ 块）：
  * **一致性不成立**：网格上一致率 < 0.9，且**双向均有反例**
    （形状判震荡/概率不判 与 形状不判/概率判震荡 两类都存在）
    ⟹ 两判据互补而非冗余，需合并覆盖（这是 B-B32 的登记结论）。
  * 判据各有牙齿：
    - 形状判据（推论 6.11）抓「局部收敛被重开」：高收敛+无重开（c=0.9,o=0）⟹ 非震荡；
      反复重开（o 高）⟹ 震荡。
    - 概率判据抓「未决熵下降未归零」：人工构造 U 下降且末端未归零 ⟹ 震荡；
      U 恒 0 ⟹ 非震荡（判据不是一律同一答案）。
  * prob_judge 分量结构：u_positive / declining / unset 三条件与总判定的合取关系。

不锁：一致率的具体数值（0.75）、不一致格子的具体 (c,o) 分布。

诚实边界：S2 合成候选竞争模型，判据操作化是建模选择；结论是框架内两判据的
逻辑关系（互补），不是现实认知过程的主张。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "b32_shape_prob_consistency"))

import b32_consistency as R  # noqa: E402


def test_agreement_fails_with_counterexamples_both_ways():
    """B-B32 登记结论：一致率 < 0.9 且双向均有反例 ⟹ 互补、需合并覆盖。"""
    res = R.main()
    s = res["summary"]
    assert s["total"] == len(res["grid"]["c_rates"]) * len(res["grid"]["o_rates"])
    assert s["agree_rate"] < 0.9, s["agree_rate"]
    d = s["disagree_cases"]
    shape_no_prob_yes = [x for x in d if not x["shape_osc"] and x["prob_osc"]]
    shape_yes_prob_no = [x for x in d if x["shape_osc"] and not x["prob_osc"]]
    assert shape_no_prob_yes, "no case: shape=no & prob=yes"
    assert shape_yes_prob_no, "no case: shape=yes & prob=no"


def test_shape_criterion_has_teeth():
    """形状判据：高收敛无重开 ⟹ 非震荡；反复重开 ⟹ 震荡。"""
    ms_conv, _ = R.simulate(0.9, 0.0, R.SEED_BASE)
    osc_conv, det_conv = R.shape_judge(ms_conv)
    assert osc_conv is False and det_conv["reopen"] == 0, (osc_conv, det_conv)

    ms_open, _ = R.simulate(0.9, 0.3, R.SEED_BASE)
    osc_open, det_open = R.shape_judge(ms_open)
    assert osc_open is True and det_open["reopen"] >= 1, (osc_open, det_open)


def test_prob_criterion_has_teeth_on_synthetic_series():
    """概率判据：人工序列上 U 下降且末端未归零 ⟹ 震荡；U 恒 0 ⟹ 非震荡。"""
    # U 单调下降但末端未归零（未决振荡）
    us = [2.0 - 0.004 * k for k in range(R.LEN)]
    osc, det = R.prob_judge(us)
    assert osc is True, det
    assert det["u_positive"] and det["declining"] and det["unset"], det

    # U 恒 0（已归零）
    us0 = [0.0] * R.LEN
    osc0, det0 = R.prob_judge(us0)
    assert osc0 is False, det0
    assert not det0["u_positive"] and not det0["unset"], det0


def test_prob_criterion_conjunctive_structure():
    """prob_judge 是三条件的合取：拆掉任一条件 ⟹ 判非震荡。"""
    us_declining_settled = [2.0 - 0.01 * k for k in range(R.LEN)]  # 降到 0 以下，末端已归零
    osc, det = R.prob_judge(us_declining_settled)
    assert osc is False and det["declining"] and not det["unset"], det

    us_flat_positive = [1.5] * R.LEN   # U>0 但无下降趋势
    osc2, det2 = R.prob_judge(us_flat_positive)
    assert osc2 is False and det2["u_positive"] and not det2["declining"], det2
