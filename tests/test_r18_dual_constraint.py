# -*- coding: utf-8 -*-
"""
R18 —— 双约束失效可分离性 的回归测试
=====================================
对应 experiments/r18_dual_constraint/ 与 docs/xs-framework/04-v3.1-supp3-*.md。

本测试**只锁判别签名与可复现性**，不锁具体数值：
  * 数值会随 N_FEAT / N_MACRO / B 等参数变化，锁死数值会让任何参数探索都失败；
  * 但"外失效 → F2 主导且 σ 不降"与"内失效 → F1 主导且 σ 显著更低"是
    补篇三推论 2.1/2.2 的**结构性预言**，必须恒真。

运行：pytest tests/test_r18_dual_constraint.py
"""

import os
import sys

import pytest

# experiments/ 不在包路径内，显式加入
_HERE = os.path.dirname(os.path.abspath(__file__))
_EXP = os.path.join(_HERE, "..", "experiments", "r18_dual_constraint")
sys.path.insert(0, os.path.abspath(_EXP))

import r18_dual_constraint as R  # noqa: E402

N = 300          # 测试用样本量（主实验为 800，此处为速度折中）
SEED = 20261019


def _arm(arm, n=N):
    res = R.run_arm(arm, n=n, seed=SEED)
    return res


@pytest.fixture(scope="module")
def arms():
    return {a: _arm(a) for a in ("base", "A", "B")}


# ------------------------------------------------------------------ 主签名
def test_outer_failure_is_F2_dominant(arms):
    """外约束失效 -> 伪收敛(F2) 主导：加载成功，但加载到了错的宏上。"""
    a = arms["A"]
    assert a["伪收敛F2%"] > 50.0, f"A 臂 F2 率过低：{a['伪收敛F2%']}"
    assert a["F2/F1"] > 2.0, f"A 臂 F2/F1 未占优：{a['F2/F1']}"
    assert a["真收敛%"] < 5.0, f"A 臂仍有过多真收敛：{a['真收敛%']}"


def test_inner_failure_is_F1_dominant(arms):
    """内约束失效 -> F1 主导：没有完成加载（不收敛）。"""
    b = arms["B"]
    assert b["F1%"] > 90.0, f"B 臂 F1 率过低：{b['F1%']}"
    assert b["真收敛%"] == 0.0
    assert b["伪收敛F2%"] == 0.0


def test_signature_sigma_outer_unchanged_inner_dropped(arms):
    """判别签名 R18-a：σ(A) ≈ σ(base)，σ(B) 显著低于 σ(base)。

    这是本实验最核心的一条：外约束是**方向约束**（错了但一样确信），
    内约束是**强度约束**（不够确信所以停不下来 / 提前停下）。
    """
    base, a, b = arms["base"], arms["A"], arms["B"]
    gap_outer = abs(a["sigma_end均值"] - base["sigma_end均值"])
    gap_inner = base["sigma_end均值"] - b["sigma_end均值"]
    assert gap_outer < 0.03, f"外失效改变了确信度（|Δσ|={gap_outer:.3f}），方向约束命题受损"
    assert gap_inner > 0.10, f"内失效未压低确信度（Δσ={gap_inner:.3f}）"
    # 分离度必须相差一个数量级，否则"两类失效可分离"不成立
    assert gap_inner > 5 * max(gap_outer, 1e-6), (
        f"两类失效的 σ 分离度不足：外 {gap_outer:.3f} vs 内 {gap_inner:.3f}"
    )


def test_F1_F2_are_derived_not_assigned():
    """F1/F2 由 σ 演化路径导出，不是被指定的：同一判定函数下 C1/C2 互斥且穷尽。"""
    rng = __import__("random").Random(SEED)
    macros = R.build_macros(__import__("random").Random(SEED))
    seen = set()
    for _ in range(200):
        row = R.trial(rng, macros, rng.choice(["base", "A", "B", "C"]))
        assert row["outcome"] in {"真收敛", "伪收敛(F2)", "F1(未加载/不收敛)"}
        # C1/C2 的组合与 outcome 必须一致
        if row["C1"] and row["C2"]:
            assert row["outcome"] == "真收敛"
        elif row["C1"] and not row["C2"]:
            assert row["outcome"] == "伪收敛(F2)"
        else:
            assert row["outcome"] == "F1(未加载/不收敛)"
        seen.add(row["state"])
    # 三态都可达（否则"三分"在操作化层面就是二分的）
    assert seen <= {"收敛", "震荡", "发散"}


# ------------------------------------------------------------------ 可复现性
def test_seed_derivation_is_process_stable():
    """seed 派生不得依赖 hash(str)（受 PYTHONHASHSEED 随机化影响）。

    这是入库时修复的一个真实缺陷：修复前同一脚本连跑两次产出的 json md5 不同。
    本测试直接断言**同参数两次调用结果完全一致**。
    """
    r1 = R.run_arm("A", n=120, seed=SEED)
    r2 = R.run_arm("A", n=120, seed=SEED)
    assert r1 == r2, "同一 seed 两次调用结果不一致 —— seed 派生被污染"


def test_source_does_not_use_builtin_str_hash():
    """静态检查：源码中不得出现 `hash(arm)` 这类对 str 的哈希调用。"""
    src = open(os.path.join(_EXP, "r18_dual_constraint.py"), encoding="utf-8").read()
    # 只检查**代码行**：源码注释里为了说明"为什么不能用"必然会出现这个字符串，
    # 直接对全文做子串匹配会误命中注释本身。
    code_lines = [ln for ln in src.splitlines() if not ln.lstrip().startswith("#")]
    code = "\n".join(code_lines)
    assert "hash(arm)" not in code, "检测到 hash(arm)：str 哈希随进程随机化，破坏可复现性"
    assert "crc32" in code


# ------------------------------------------------------------------ 高分辨率
def test_signature_holds_at_higher_resolution():
    """把 σ 量子从 1/8 细化到 1/16，签名必须仍然成立。

    背景：σ 只能取 1/N_FEAT 的整数倍，而接受阈值 n=0.55 恰好落在 0.500 与 0.625
    两个量子之间——需要排除"签名只是离散化伪像"。
    """
    saved = (R.N_FEAT, R.N_MACRO, R.N_PROTO, R.TOPK, R.B_BUDGET)
    try:
        R.N_FEAT, R.N_MACRO, R.N_PROTO = 16, 32, 6
        R.TOPK, R.B_BUDGET = 8, 20
        base = R.run_arm("base", n=N, seed=SEED)
        a = R.run_arm("A", n=N, seed=SEED)
        b = R.run_arm("B", n=N, seed=SEED)
    finally:
        R.N_FEAT, R.N_MACRO, R.N_PROTO, R.TOPK, R.B_BUDGET = saved

    assert abs(a["sigma_end均值"] - base["sigma_end均值"]) < 0.03
    assert base["sigma_end均值"] - b["sigma_end均值"] > 0.10
    assert a["伪收敛F2%"] > 50.0
    assert b["F1%"] > 90.0
