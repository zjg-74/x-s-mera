#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R8-B9 —— 势消除混沌的唯一机制归约（M1/M2/M3 消融判别）结构性测试

对应：experiments/r8_gain_stepchain/r8b_uniqueness.py、补篇一缺口 B-B9、
v1.6.37（框架内唯一性成立；跨框架全局唯一性转外部开放）。

锁的结构预言：
  * U1：增益调制（M1）在混沌窗内可把 λ 由正转负（移出混沌）——纯函数独立复算。
  * U3：合法对称破缺 = 各候选增益不对称（乘性），λ(s) 单调 ⟹ M2 合法形式 ⊆ M1
    —— 独立复算 lam2 的单调性。
  * U4：强制归宿（μ）单调压缩混沌连续寿命（M3 与拓扑正交）——独立复算。
  * U2：势进状态导数（加性偏向）被 A11-R 合规检查排除（结构性登记）。
  * 综合判据：out["uniqueness_framework_internal"] 为真（模块级执行产物）。

不锁：g_chaos 的具体数值、寿命的具体数值、λ 的具体数值。

诚实边界：X-S-MERA 离散迭代 + A11-R 框架内的唯一性；跨框架全局唯一性
非本框架可裁（登记为外部开放项，测试同样不声称全局唯一）。

实现说明：该实验是脚本式（模块级直接执行并写 JSON 到 CWD），故在导入前
把 CWD 切到 pytest 临时目录，避免污染仓库；导入产物 out 直接断言。
"""

import importlib.util
import os

import pytest

EXP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "experiments", "r8_gain_stepchain", "r8b_uniqueness.py")


@pytest.fixture(scope="module")
def m(tmp_path_factory):
    """在临时目录中执行脚本式实验一次，返回模块（out 为完整结果）。"""
    d = tmp_path_factory.mktemp("r8b_cwd")
    old = os.getcwd()
    os.chdir(str(d))
    try:
        spec = importlib.util.spec_from_file_location("r8b_uniqueness", EXP)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
    finally:
        os.chdir(old)
    return m


def test_u1_gain_modulation_leaves_chaos(m):
    """U1 独立复算：λ(3.8)>0（混沌）且 λ(3.5)≤0（调低增益移出混沌）。"""
    lam_hi = m.lam_g(3.8)
    lam_lo = m.lam_g(3.5)
    assert lam_hi > m.TOL_LAM, lam_hi
    assert lam_lo <= m.TOL_LAM, lam_lo
    assert m.out["U1"]["gain_is_multiplicative_only"] is True


def test_u2_additive_bias_channel_excluded(m):
    """U2：势进状态导数（加性项）被 A11-R 合规检查排除——结构性登记。"""
    u2 = m.out["U2"]
    assert u2["additive_bias_enters_state_derivative"] is True
    assert u2["channel_excluded_by_framework"] is True


def test_u3_symmetry_breaking_reduces_to_m1(m):
    """U3 独立复算：s 只进乘性增益时 λ(s) 单调 ⟹ M2 合法形式 ⊆ M1。"""
    s_list = [0.0, 0.02, 0.05, 0.08, 0.10, 0.15]
    vals = [m.lam2(3.8, s) for s in s_list]
    mono = (all(b > a for a, b in zip(vals, vals[1:])) or
            all(b < a for a, b in zip(vals, vals[1:])))
    assert mono, vals
    assert m.out["U3"]["s_only_in_multiplicative_gain"] is True


def test_u4_forced_absorption_compresses_lifetime_only(m):
    """U4 独立复算：μ 增大 ⟹ 混沌连续寿命单调下降（阈值不动，拓扑正交）。"""
    mus = [0.0, 0.01, 0.05, 0.2]
    life = [m.chaos_lifetime(3.8, mu) for mu in mus]
    assert all(b < a for a, b in zip(life, life[1:])), life
    assert m.out["U4"]["g_chaos_independent_of_mu"] is True


def test_verdict_framework_internal_uniqueness(m):
    """综合判据：U1∧U2∧U3∧U4 ⟹ 框架内唯一性成立（不声称跨框架全局唯一）。"""
    v = m.out["verdicts"]
    assert v["U1_M1_sufficient_compliant"] is True
    assert v["U2_illegal_channel_excluded"] is True
    assert v["U3_M2_legal_reduces_to_M1"] is True
    assert v["U4_M3_orthogonal_compresses_lifetime"] is True
    assert m.out["uniqueness_framework_internal"] is True
    # 诚实边界：结论文本不得声称全局唯一
    assert "跨框架" in m.out["conclusion"] or "框架内" in m.out["conclusion"]
