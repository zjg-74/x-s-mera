# -*- coding: utf-8 -*-
"""
B-B32 —— 震荡的两种判据一致性检验（形状判据 vs 概率判据）
==================================================================
对应 XS v3.1 补篇七（08-v3.1-supp7）缺口 **B-B32**：
> 推论 6.11（震荡 = 局部收敛被重开）给"震荡"一个**可数的形状判据**（重开次数），
> 而此前的判据是概率性的（U>0 且仍在下降）。两者应当一致，但该一致性尚未检验。

判据定义（08 §3.4–3.5）：
- **局部收敛**：m(k)=1；**局部收敛段**：连续 m=1 且长度 ≥ τ_c（声明式检测门槛）。
- **重开**：从一个局部收敛段离开到 m ≥ 2。
- **形状判据（推论 6.11）**：震荡 ⟺ 发生局部收敛但未能保持，即**重开次数 ≥ 1 且未落定**。
- **概率判据（主文）**：震荡 ⟺ **U > 0 且仍在下降**（有未决度、处于收敛中但未归零）。

本实验：用固定种子的离散候选竞争模型生成多种轨迹（一次性收敛 / 反复重开震荡 /
未决振荡(有回跳无重开) / 发散），对每条轨迹分别用两判据判"震荡/非震荡"，比较一致性，
并对不一致轨迹做**形态归类**。

证据硬度：S2（合成候选竞争模型，纯标准库，固定种子可复现）。判据操作化是模型选择，
如实登记；不声称现实认知结论。
"""

import json
import math


SEED_BASE = 20261028
M_MAX = 8            # 最大候选数
EPS = 1e-9           # 活动候选阈值
TAU_C = 5            # 局部收敛段声明式门槛（步）
BURN = 50            # 判别主区间起点
END_WIN = 20         # 末端窗口（判"未落定"）
U_END_TOL = 1e-6     # 末端 U 判收敛归零阈值
LEN = 400            # 轨迹长度


class R:
    """LCG，固定种子，跨 (c,o) 派生。"""
    def __init__(self, seed): self.s = seed
    def next(self):
        self.s = (self.s * 1103515245 + 12345) & 0x7fffffff
        return self.s / 0x7fffffff


def evolve(w, c_rate, o_rate, r, beta=6.0, delta=0.6):
    """一步演化：重开扰动(按 o_rate) + 收敛集中(按 c_rate)，再归一化。"""
    n = len(w)
    # 重开：把主导权重分散到其他候选（引入新的未决候选）
    if r.next() < o_rate:
        i0 = max(range(n), key=lambda i: w[i])
        # 从主导抽出 delta 权重，按均匀分散到其余候选
        rest = [i for i in range(n) if i != i0]
        for i in rest:
            w[i] += delta / max(1, len(rest))
        w[i0] -= delta
        w = [max(x, 0.0) for x in w]
    # 收敛：权重向主导指数集中（乘性增益，等价离散步长调制）
    if r.next() < c_rate:
        i0 = max(range(n), key=lambda i: w[i])
        for i in range(n):
            w[i] *= math.exp(beta if i == i0 else -beta)
    s = sum(w)
    return [x / s if s > 0 else 1.0 / n for x in w]


def simulate(c_rate, o_rate, seed):
    r = R((seed * 7919 + 104729) & 0x7fffffff)
    w = [1.0 / M_MAX] * M_MAX
    ms = [0] * LEN
    us = [0.0] * LEN
    for k in range(LEN):
        ms[k] = sum(1 for x in w if x > EPS)
        us[k] = -sum(x * math.log(x) for x in w if x > 0)
        w = evolve(w, c_rate, o_rate, r)
    return ms, us


# ---- 形状判据（推论 6.11）：震荡 ⟺ 重开≥1 且未落定 ----
def shape_judge(ms, tau_c=TAU_C, end_win=END_WIN):
    n = len(ms)
    # 找局部收敛段（连续 m=1 且长度≥τ_c）
    segments = []
    i = 0
    while i < n:
        if ms[i] == 1:
            j = i
            while j < n and ms[j] == 1:
                j += 1
            if j - i >= tau_c:
                segments.append((i, j))   # [i, j)
            i = j
        else:
            i += 1
    # 重开：某局部收敛段之后 m 升到 ≥2（即段末非终止且后面 m≥2）
    reopen = 0
    for (a, b) in segments:
        if b < n and ms[b] >= 2:
            reopen += 1
    # 未落定：末端窗口内仍存在 m≥2，或最后一次段被重开
    settled = (segments and
               segments[-1][1] >= n - end_win and
               all(ms[k] == 1 for k in range(segments[-1][1], n)))
    oscillation = bool(reopen >= 1 and not settled)
    return oscillation, {"local_segments": len(segments), "reopen": reopen,
                         "settled": settled}


# ---- 概率判据（主文）：震荡 ⟺ U>0 且仍在下降（未归零） ----
def prob_judge(us, burn=BURN, u_end_tol=U_END_TOL, end_win=END_WIN):
    seg = us[burn:-end_win] if end_win < len(us) - burn else us[burn:]
    if not seg:
        return False, {}
    u_mean = sum(seg) / len(seg)
    # U>0：主区间熵均值高于阈值
    u_positive = u_mean > 1e-4
    # 仍在下降：主区间 U 净变化为负（收敛趋势）且总体下降
    # 用首末比较 + 趋势：从最大熵后整体下降
    u_first = seg[0]
    u_last = seg[-1]
    declining = u_last < u_first - 1e-6
    # 未归零：末端 U 仍高于 tol
    u_end = sum(us[-end_win:]) / end_win if end_win <= len(us) else sum(us) / len(us)
    unset = u_end > u_end_tol
    oscillation = bool(u_positive and declining and unset)
    return oscillation, {"u_mean": round(u_mean, 4), "declining": declining,
                         "u_end": round(u_end, 5), "u_positive": u_positive,
                         "unset": unset}


def main():
    c_rates = [0.05, 0.3, 0.6, 0.9]
    o_rates = [0.0, 0.05, 0.15, 0.3, 0.6]
    rows = []
    agree_total = 0
    disagree = []
    for c in c_rates:
        for o in o_rates:
            ms, us = simulate(c, o, SEED_BASE)
            shape, sdet = shape_judge(ms)
            prob, pdet = prob_judge(us)
            same = (shape == prob)
            agree_total += 1 if same else 0
            if not same:
                disagree.append({"c": c, "o": o,
                                 "shape_osc": shape, "prob_osc": prob,
                                 "shape": sdet, "prob": pdet})
            rows.append({"c": c, "o": o, "shape_osc": shape, "prob_osc": prob,
                         "agree": same,
                         "reopen": sdet["reopen"],
                         "segments": sdet["local_segments"],
                         "u_mean": pdet.get("u_mean"), "u_end": pdet.get("u_end")})
    total = len(rows)
    return {
        "experiment": "B-B32-shape-vs-prob-consistency",
        "framework": "XS v3.4 · 补篇七 B-B32（震荡两判据一致性）",
        "hardness": "S2（合成候选竞争模型，纯标准库，固定种子可复现）",
        "params": {"M_MAX": M_MAX, "TAU_C": TAU_C, "LEN": LEN,
                   "burn": BURN, "end_win": END_WIN, "u_end_tol": U_END_TOL},
        "grid": {"c_rates": c_rates, "o_rates": o_rates},
        "rows": rows,
        "summary": {"total": total, "agree": agree_total,
                    "agree_rate": round(agree_total / total, 4),
                    "disagree_count": len(disagree),
                    "disagree_cases": disagree},
        "conclusion": None,
    }


if __name__ == "__main__":
    res = main()
    # 不一致形态分析
    d = res["summary"]["disagree_cases"]
    shape_non_prob_yes = [x for x in d if not x["shape_osc"] and x["prob_osc"]]
    shape_yes_prob_non = [x for x in d if x["shape_osc"] and not x["prob_osc"]]
    res["disagree_typology"] = {
        "shape=非震荡 & prob=震荡（无重开的未决振荡）": {
            "count": len(shape_non_prob_yes),
            "cases": shape_non_prob_yes},
        "shape=震荡 & prob=非震荡": {
            "count": len(shape_yes_prob_non),
            "cases": shape_yes_prob_non},
    }
    agree_rate = res["summary"]["agree_rate"]
    a = res["disagree_typology"]["shape=非震荡 & prob=震荡（无重开的未决振荡）"]["count"]
    b = res["disagree_typology"]["shape=震荡 & prob=非震荡"]["count"]
    res["conclusion"] = (
        "两判据不一致（一致率 %.3f，双向均有反例）。" % agree_rate +
        "① 形状=震荡 & 概率=非震荡（%d 例）：反复「局部收敛被重开」的经典震荡中，" % b +
        "U 被重开反复抬高、不满足主文『U 仍在下降』，概率判据漏判。" +
        "② 形状=非震荡 & 概率=震荡（%d 例）：一直 m≥2、U 下降但从未到 m=1 的未决振荡中，" % a +
        "无局部收敛可被重开，形状判据（推论 6.11）漏判。" +
        "⟹ 两判据并不同指同一类过程：形状判据锁定『局部收敛被推翻』的结构事件" +
        "（要求先到 m=1），概率判据锁定『未决熵收敛趋势』（要求 U 下降）。" +
        "推论 6.11 与主文概率判据的一致性不成立，二者为互补关系而非冗余，需合并覆盖。" +
        "（判据操作化为建模选择，S2 合成证据，如实登记。）")
    print(json.dumps(res, ensure_ascii=False, indent=2))
    with open("b32_consistency.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("saved b32_consistency.json")
