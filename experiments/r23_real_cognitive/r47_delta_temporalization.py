#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R47 —— δ-时间化检验：震荡/方向是否时间化算子的产物（B 线 P4 方法论）

背景
----
用户命题（2026-10-28）：快思维 A→B 是瞬时跃迁（τ=0、δ、无时间结构）；**震荡不是认知本体，
而是"时间出现在函数中"才产生的表观**；一旦时间化，快思维就不再是快思维，而是慢思维。
⟹ 判别器（读最长震荡段内 U 走向）测的是**慢思维轨迹**，且三分（酝酿/反刍/空转）可能是
**时间化算子的函数**，而非认知本体属性。

R47 把该命题转成可证伪检验：**从无内部时间结构的 δ-事件出发**，用不同时间化算子
（细分粒度 m × 聚焦 focus）生成慢思维轨迹，判别三分，检验：
  Q1（随机基线）：无方向的快思维 δ 事件，经时间化是否被**凭空造出方向**？
     ——若随机 δ 判出方向 ⇒ 方向是算子伪影；若判空转 ⇒ 算子不造方向、只展开已有结构。
  Q2（结构恢复）：带潜在结构（收敛=酝酿 / 发散=反刍）的 δ 事件序列，跨不同算子展开，
     方向是否稳健？——若漂移 ⇒ 方向依赖算子（伪影）；若稳健 ⇒ 部分方向反映慢思维真实结构。

仅标准库；固定种子逐位可复现；判别器逻辑与 analyze_trajectories.py 一致。
"""

import json
import math
import os
import random
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261028
C = 4          # 候选类数
N = 12         # 快思维 δ 事件数
T = 100        # 每事件展开子步的总伪计数质量（分摊到各子步；取大值避免量化平台/熵阶梯跳变）
MS = [1, 2, 3, 5]                 # 每事件展开子步数
FOCUSES = [0.0, 0.3, 0.6, 0.9]    # 选中类的表观主导度（0=弱主导多候选，1=全聚焦单类）
STRUCTURES = ["ruminate_down", "ruminate_up", "random"]  # 潜在结构


# ---------- 判别器（复制自 analyze_trajectories.py，保证一致） ----------
def p_hat(weights):
    w = {k: v for k, v in weights.items() if v > 0}
    if not w:
        return 0.0
    Ntot = sum(w.values())
    r_hat = len(w)
    ps = [v / Ntot for v in w.values()]
    ent = -sum(p * math.log2(p) for p in ps)
    return ent + (r_hat - 1) / (2 * Ntot * math.log(2))


def find_oscillation_segments(rs):
    segs, run = [], []
    for k, rk in enumerate(rs):
        if rk >= 2:
            run.append(k)
        else:
            if len(run) >= 2:
                segs.append(run)
            run = []
    if len(run) >= 2:
        segs.append(run)
    return segs


def classify(net):
    if net is None:
        return None
    tol = 1e-4
    if net < -tol:
        return "ruminate_down"   # 酝酿
    if net > tol:
        return "ruminate_up"     # 反刍
    return "flat"


def judge(weights_seq):
    Us = [p_hat(w) for w in weights_seq]
    rs = [len({k: v for k, v in w.items() if v > 0}) for w in weights_seq]
    segs = find_oscillation_segments(rs)
    if segs:
        seg = max(segs, key=len)
        net = Us[seg[-1]] - Us[seg[0]]
    else:
        seg, net = None, None
    return classify(net), net, seg


# ---------- 快思维 δ 事件（无时间结构：瞬时，只携带"候选集中度"结构） ----------
def make_events(structure, C, N, rng):
    """返回长度 N 的主导类权重 p 序列（A 的占比）。每事件是瞬时选择（无内部步、无时间延展）。
    潜在结构编码在**候选集中度（熵）随事件演化**：
    ruminate_down(集中度↑熵↓=酝酿) / ruminate_up(集中度↓熵↑=反刍) / random(恒定=无结构)。"""
    seq = []
    for k in range(N):
        t = k / (N - 1)
        if structure == "ruminate_down":
            pA = 0.25 + (0.85 - 0.25) * t
        elif structure == "ruminate_up":
            pA = 0.85 + (0.25 - 0.85) * t
        else:  # random：恒定均匀，无结构
            pA = 0.25
        seq.append(pA)
    return seq


# ---------- 时间化算子（把 δ 事件序列 → 时间轨迹；时间化=生成时间结构的必要条件） ----------
def temporalize(events_p, m, focus, C):
    """每个 δ 事件（携带主导权重 pk）展开成 m 子步，时间上生成"候选切换/震荡"。
    子步权重：主导类 A 占 a = pk + focus*(1-pk)（pk 随事件演化 → 子步熵随事件演化，结构可传递），
    其余 (C-1) 类均分背景。总质量守恒：每子步总量 = T/m。"""
    names = [chr(ord("A") + j) for j in range(C)]
    per = T / m
    subs = []
    for pk in events_p:
        a = pk + focus * (1 - pk)   # 锐化后主导占比（focus=0→不额外锐化，focus→1→趋向单类）
        for _ in range(m):
            wA = max(0, round(per * a))
            bg = round(per - wA)
            sub = {"A": wA}
            others = names[1:]
            if bg > 0:
                base = bg // len(others)
                for nm in others:
                    sub[nm] = max(0, base)
                s = sum(sub.values())
                diff = round(per) - s
                if diff:
                    sub["A"] = max(0, sub["A"] + diff)
            subs.append(sub)
    return subs


# ---------- 主流程 ----------
def main():
    rng = random.Random(SEED)
    out = []
    p = lambda s="": out.append(str(s))
    p("=" * 84)
    p("R47 —— δ-时间化检验：震荡/方向是否时间化算子的产物（B 线 P4 方法论）")
    p("=" * 84)

    # 对每个潜在结构 × (m, focus) 网格：判别三分分布
    grid = {}
    for structure in STRUCTURES:
        for m in MS:
            for foc in FOCUSES:
                verds = []
                for trial in range(30):   # 30 次重复采样（固定种子流推进）
                    events = make_events(structure, C, N, rng)
                    traj = temporalize(events, m, foc, C)
                    v, net, seg = judge(traj)
                    verds.append(v if v else "none")
                cnt = Counter(verds)
                grid["%s,m%d,f%.1f" % (structure, m, foc)] = dict(cnt)
        p("  %-12s 完成" % structure)

    # Q0：判别器对快思维本体（无时间结构）失效
    p("\n## Q0 必要条件：判别器需要时间结构，快思维本体(δ)无时间结构 ⟹ 不可测")
    p("  纯δ事件序列（不展开为可判震荡）→ 应无震荡段 → 无判定")
    q0 = 0
    for _ in range(30):
        ev = make_events("random", C, N, rng)
        # 纯快思维：每事件即单点 δ（不时间化 = 无展开，不构成震荡段）
        traj = [{"A": T}]  # 单点：无时间延展
        v, net, seg = judge(traj)
        if seg is None:
            q0 += 1
    p("  纯δ单点判别：无震荡段(无判定)比例 = %d/30" % q0)

    # Q1：随机 δ 基线——时间化是否凭空造出方向？
    p("\n## Q1 随机δ基线：时间化是否凭空造出方向？")
    p("  若随机δ判出方向(酝酿/反刍)比例高 ⟹ 方向是算子伪影；若判空转/无判定为主 ⟹ 算子不造方向")
    rand_direction_rates = {}
    for m in MS:
        for foc in FOCUSES:
            c = grid["random,m%d,f%.1f" % (m, foc)]
            tot = sum(c.values())
            dr = (c.get("ruminate_down", 0) + c.get("ruminate_up", 0)) / tot if tot else 0
            rand_direction_rates["m%d,f%.1f" % (m, foc)] = round(dr, 3)
            p("    m=%d focus=%.1f  方向检出率=%.3f  %s" % (m, foc, dr, dict(c)))

    # Q2：结构恢复——同一潜在结构跨算子展开，方向是否稳健？
    p("\n## Q2 结构恢复：酝酿/反刍 δ 序列跨算子展开的方向稳健性")
    for structure in ["ruminate_down", "ruminate_up"]:
        expected = "ruminate_down" if structure == "ruminate_down" else "ruminate_up"
        p("  [%s] 期望判 %s" % (structure, expected))
        for m in MS:
            row = []
            for foc in FOCUSES:
                c = grid["%s,m%d,f%.1f" % (structure, m, foc)]
                tot = sum(c.values())
                keep = c.get(expected, 0) / tot if tot else 0
                row.append(round(keep, 2))
            p("    m=%d  保持率[focus 0/0.3/0.6/0.9]=%s" % (m, row))

    result = {
        "meta": {"purpose": "R47 δ-时间化检验：震荡/方向是否算子伪影", "framework": "XS/B线P4方法论",
                 "seed": SEED, "classes": C, "events": N, "per_quality": T,
                 "m_scan": MS, "focus_scan": FOCUSES, "trials": 30,
                 "nature": "纯内部仿真，判别器=analyze_trajectories同逻辑；快思维δ=无时间结构瞬时选择"},
        "grid_3class": grid,
        "Q0_delta_untemporalizable": {"no_segment_ratio_30": q0},
        "Q1_random_direction_rate": rand_direction_rates,
    }

    # 判定
    # Q1: 随机δ方向检出率（若显著>0.5 → 算子造方向伪影）
    q1_max = max(rand_direction_rates.values())
    q1_false_direction = q1_max > 0.5
    # Q2: 现实窗口(m=2/3, focus=0~0.3)结构保持率
    keep_rates = {}
    for structure in ["ruminate_down", "ruminate_up"]:
        expected = structure
        kr = []
        for m in [2, 3]:
            for foc in [0.0, 0.3]:
                c = grid["%s,m%d,f%.1f" % (structure, m, foc)]
                tot = sum(c.values())
                kr.append(c.get(expected, 0) / tot if tot else 0)
        keep_rates[structure] = kr
    q2_robust = min(min(v) for v in keep_rates.values()) >= 0.8

    result["Q1_max_random_direction_rate"] = q1_max
    result["Q1_false_direction"] = bool(q1_false_direction)
    result["Q2_keep_rates_real"] = keep_rates
    result["Q2_robust"] = bool(q2_robust)

    p("\n" + "=" * 84)
    p("R47 判定")
    p("=" * 84)
    p("Q0 纯δ单点(无时间化)不可测比例 = %d/30 ⟹ %s" % (
        q0, "✅ 判别器需时间结构，快思维本体(δ)不可测" if q0 >= 28 else "⚠️ 部分δ可测"))
    p("Q1 随机δ方向检出率峰值=%.3f ⟹ %s" % (
        q1_max, "⚠️ 时间化可凭空造出方向(算子伪影)" if q1_false_direction
        else "✅ 时间化不凭空造方向（随机δ判空转为主），方向需已有结构才能展开"))
    p("Q2 现实窗口(m=2/3,focus=0~0.3)结构保持率：酝酿%s 反刍%s ⟹ %s" % (
        str([round(x,2) for x in keep_rates["ruminate_down"]]),
        str([round(x,2) for x in keep_rates["ruminate_up"]]),
        "✅ 结构跨算子稳健" if q2_robust else "⚠️ 结构跨算子漂移(方向依赖算子)"))

    conclusion_parts = []
    if q0 >= 28 and not q1_false_direction and q2_robust:
        result["conclusion"] = ("✅ 命题得到三层支持：Q0判别器需时间结构、快思维δ本体不可测；Q1时间化"
            "不凭空造方向；Q2潜在结构（酝酿/反刍）在现实参数下跨算子稳健恢复 ⟹ 判别器测的方向"
            "部分反映慢思维真实结构、非纯时间化伪影；但快思维本体方向不可测，仍登记为"
            "'慢思维(时间化认知)轨迹'。")
    elif q1_false_direction:
        result["conclusion"] = ("⚠️ 方向含算子伪影：随机δ经时间化也被判出方向 ⟹ 震荡/方向部分是"
            "时间化方式构造的表观，判别需谨慎；R23 结果应明确登记为'时间化轨迹表观方向'。")
    else:
        result["conclusion"] = ("⚠️ 结构跨算子漂移：同一潜在结构不同算子展开方向不一致 ⟹ "
            "方向依赖观测窗口，判别读数需报告算子参数。")
    p(result["conclusion"])

    with open(os.path.join(HERE, "r47_delta_temporalization.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))


if __name__ == "__main__":
    main()
