# -*- coding: utf-8 -*-
"""
R21 回归测试 —— 内约束强度的操作化（O_adm 剪枝量）

背景：R18-b 用 tau（接受阈值）与 p_opt（注入最优性）编码"内约束强度"，
两者都不成立（前者同义反复、后者无效应），已撤回。R21 改用 O_adm 剪枝量。

本测试锁的四件事：
  1. F2 率关于 Q 呈**倒 U**（非单调）—— 这是 R21 的核心发现
  2. **条件于真值宏未被剪掉**的子样本中，F2 单调下降、真收敛上升
     （两个效应必须分层，否则会把"真值被剪"的平凡效应写成发现）
  3. 真收敛率单调下降、F1 单调上升
  4. **防回退**：实验不得改动停止规则/阈值/判定函数
     （R18-b 的失败正是因为动了阈值；只准改 O_adm）

另加一条：倒 U 在高分辨率下保持（排除 sigma 量化伪像）。
"""

import os
import re
import sys

import pytest

_EXP = os.path.join(os.path.dirname(__file__), os.pardir,
                    "experiments", "r21_constraint_specificity")
sys.path.insert(0, os.path.abspath(_EXP))

import r21_constraint_specificity as R  # noqa: E402
import r21_sweep as S  # noqa: E402

N = 350  # 测试用小样本；结构性签名在小样本下已稳定


@pytest.fixture(scope="module")
def sweep():
    return [S.run_K(K, n=N) for K in S.K_GRID]


def _f2(rows):
    return [r["F2%"] for r in rows]


def test_F2_is_inverted_U(sweep):
    """核心预言：F2 率关于 Q 非单调，呈倒 U。"""
    f2 = _f2(sweep)
    mono_up = all(f2[i] <= f2[i + 1] for i in range(len(f2) - 1))
    mono_dn = all(f2[i] >= f2[i + 1] for i in range(len(f2) - 1))
    assert not (mono_up or mono_dn), f"F2 应非单调（倒 U），实得序列 {f2}"
    peak = max(range(len(f2)), key=lambda i: f2[i])
    # 峰值应在中部，不在两端
    assert 0 < peak < len(f2) - 1, f"倒 U 峰值应在内部，实得 index={peak} / {f2}"


def test_survivor_subset_F2_decreases(sweep):
    """
    非平凡层：条件于真值宏 M* **未被**剪掉，剪枝是纯收益 —— F2 应下降。

    这一条同时是我自己提出的"虚假唯一性"预言的证否证据：
    gap 增大了，但 F2 反而下降。锁死它，防止日后再被当作发现重新提出。
    """
    surv = [r["存活_F2%"] for r in sweep]
    first, last = surv[0], surv[-1]
    assert last < first, \
        f"存活子样本 F2 应下降（真值未被剪时剪枝是纯收益），实得 {first} -> {last}"


def test_survivor_true_convergence_rises(sweep):
    """存活子样本的真收敛率应上升。"""
    surv = [r["存活_真收敛%"] for r in sweep]
    assert surv[-1] > surv[0], f"存活子样本真收敛应上升，实得 {surv[0]} -> {surv[-1]}"


def test_true_convergence_decreases_and_F1_rises(sweep):
    """全样本：真收敛单调下降、F1 单调上升（允许小幅噪声，故只锁首尾）。"""
    tc = [r["真收敛%"] for r in sweep]
    f1 = [r["F1%"] for r in sweep]
    assert tc[-1] < tc[0], f"真收敛应下降，实得 {tc[0]} -> {tc[-1]}"
    assert f1[-1] > f1[0], f"F1 应上升，实得 {f1[0]} -> {f1[-1]}"


def test_star_ban_rate_matches_design(sweep):
    """P(M* 被剪) 应等于 1 - K/|O_all|（设计值），用于确认剪枝确在实施。"""
    for r in sweep[:4]:  # 强剪枝段存活样本太小，只校验前几档
        expect = round(100 * (1 - r["K"] / R.N_MACRO), 1)
        assert abs(r["M*被剪%"] - expect) < 6.0, \
            f"K={r['K']}: M*被剪率 {r['M*被剪%']} 偏离设计值 {expect} 过多"


def test_source_does_not_touch_stop_rule_or_thresholds():
    """
    防回退（本测试最重要的一条）：实验**不得**改动停止规则、接受阈值或判定函数。

    R18-b 的失败正是因为用接受阈值 tau 编码"内约束强度"——那是同义反复。
    R21 只允许改 O_adm（可用宏集）。本测试检查源码没有给不同层级分派不同
    阈值/停止条件。
    """
    src = open(os.path.join(_EXP, "r21_constraint_specificity.py"), encoding="utf-8").read()
    # 只检查代码行
    code = "\n".join(ln for ln in src.splitlines() if not ln.lstrip().startswith("#"))
    # 不得出现按 n_ban / K 分支的阈值调整
    for pat in (r"N_TH\s*\*", r"M_TH\s*\*", r"N_TH\s*\+", r"M_TH\s*\+",
                r"EPS_GAP\s*\*", r"B_BUDGET\s*\*"):
        assert not re.search(pat, code), \
            f"检测到阈值被缩放（{pat}）：内约束强度只许操作化在 O_adm 上，不许动阈值"
    # 判定函数只出现一次（对所有层级共用）
    assert code.count('outcome = "真收敛"') == 1, \
        "判定函数应唯一，对所有剪枝层级共用；出现多次说明分层判定了"


def test_inverted_U_holds_at_higher_resolution():
    """排除"倒 U 是 sigma 离散化量子的伪像"：提高分辨率后应仍非单调。"""
    import r21_highres as H
    rows = [H.run_K(K, n=250) for K in (40, 26, 20, 12, 8, 4)]
    f2 = [r["F2%"] for r in rows]
    mono_up = all(f2[i] <= f2[i + 1] for i in range(len(f2) - 1))
    mono_dn = all(f2[i] >= f2[i + 1] for i in range(len(f2) - 1))
    assert not (mono_up or mono_dn), \
        f"高分辨率下 F2 应仍非单调，实得 {f2}"
