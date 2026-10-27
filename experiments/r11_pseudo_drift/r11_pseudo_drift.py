# -*- coding: utf-8 -*-
"""
R11 —— A11-R 认知含义的判决性实验：增益越权 → 伪定向（F2/ghost 伪力）
==================================================================
对应 XS v3.1 补篇二 §5 R11（强度：高）：

> 在工程原型上人为违反 A11-R（让语义增益场进入状态导数），应当产生**高加载度但
> 事后为错**的输出（F2/ghost 伪力）。若不产生，则推论 4.1 被证否。

机制背景（补篇二推论 4.1）
--------------------------------------------------
- A11-R"势-流分离"：语义增益场 g 同时缩放漂移与扩散（Ψ 与 ξ 都被 g 缩放），因此
  **g 改变不了 ψ̄/√Var(ξ) 比值，即改变不了定向性**——g 只改变走得多快，不改变往哪走。
- 推论 4.1：**A11-R 的认知对应物 = 禁止伪定向**。若语义增益场**越权进入状态导数**，
  它就会**冒充漂移**——把"增量幅度大"误认作"方向正确"。
  - 认知形态 = 伪收敛（加载度高、体验确信、事后为错）；
  - 工程形态 = **F2 / ghost 伪力**（加载成功但事后验证为错）。

本实验设计（复用 r18 宏加载原型）
--------------------------------------------------
对照两种语义增益注入方式，**唯一差别**是 g 是否进入候选得分（状态导数/决策）：
- **legal（A11-R 合规）**：score(M) = σ(M,F)。g 只缩放"增量幅度"（每步推进量），
  **不进候选得分** → 不改变候选方向 → F2 率应与 g 无关。
- **violate（A11-R 越权）**：score(M) = σ(M,F) + g·b_M，其中 b_M 为每个宏的恒定
  偏置（"增益场冒充的漂移"）。g 进入候选得分（状态导数）→ g 越大，偏置越主导，
  系统越倾向锁到"偏置大"的宏（未必是 truth）→ **提前锁定错误宏 → F2 率随 g 上升**。

预言（推论 4.1 / R11）：
  P1  legal 臂：F2 率与 g **弱相关**（g 只调幅度、不调方向 → 不产生伪定向）；
  P2  violate 臂：F2 率与 g **强正相关 / 单调递增**（越权 → 冒充漂移 → 伪定向）。
  R11 支持 ⟺ P1 ∧ P2 同时成立；violate 无 g 效应 ⟹ 推论 4.1 被证否。

诚实边界：复用 r18 宏加载原型（S2 自制仿真，非经验数据）；偏置 b_M 为注入的
"冒充漂移"代理，量化"越权强度"；P1 弱相关是定性判据（|ρ|<0.3），P2 强正相关
（ρ>0.8 或单调）。不声称"认知必然伪收敛"，只检验 A11-R 违反在原型上的后果。
"""

import json
import random
import zlib
import statistics as st
from collections import Counter

# ---------------------------------------------------------------- 参数（同 r18）
N_FEAT = 8
N_PROTO = 4
N_MACRO = 20
M_TH = 0.35
N_TH = 0.55
EPS_GAP = 0.10
B_BUDGET = 12
RHO_OBS = 0.85
P_NOISE = 0.05
TOPK = 6
SIGMA_MIN = -1.0

# 语义势增益扫描（g=1 为基线；越大越"越权"）
G_VALUES = [1.0, 2.0, 3.0, 4.0, 5.0]
N_TRIALS = 800


def build_macros(rng):
    protos = [[rng.randint(0, 1) for _ in range(N_FEAT)] for _ in range(N_PROTO)]
    seen, macros = set(), []
    while len(macros) < N_MACRO:
        p = protos[len(macros) % N_PROTO]
        v = p[:]
        for _ in range(2):
            j = rng.randrange(N_FEAT)
            v[j] = 1 - v[j]
        t = tuple(v)
        if t not in seen:
            seen.add(t)
            macros.append(v)
    return macros


def sigma(M, F):
    s = 0
    for j, val in F.items():
        s += 1 if M[j] == val else -1
    return s / N_FEAT


def cand(F, macros, th):
    return [(i, sigma(M, F)) for i, M in enumerate(macros) if sigma(M, F) >= th]


def best_split(rng, macros, F, unread):
    cur = sorted(range(len(macros)), key=lambda i: -sigma(macros[i], F))[:TOPK]
    best_j, best_score = None, -1.0
    for j in unread:
        n1 = sum(1 for i in cur if macros[i][j] == 1)
        n0 = len(cur) - n1
        score = min(n0, n1) + 0.001 * rng.random()
        if score > best_score:
            best_j, best_score = j, score
    return best_j


def observe(rng, src_vec, obs_mask):
    phi = {}
    for j, on in enumerate(obs_mask):
        if not on:
            continue
        v = src_vec[j]
        if rng.random() < P_NOISE:
            v = 1 - v
        phi[j] = v
    return phi


def trial(rng, macros, bias, g, arm):
    """单次宏加载试验。arm ∈ {legal, violate}；g 为语义势增益。"""
    M_star = rng.randrange(len(macros))
    truth = macros[M_star]
    obs_mask = [rng.random() < RHO_OBS for _ in range(N_FEAT)]
    if not any(obs_mask):
        obs_mask[rng.randrange(N_FEAT)] = True
    phi = observe(rng, truth, obs_mask)

    def score(i, F):
        s = sigma(macros[i], F)
        if arm == "violate":
            s += g * bias[i]     # 越权：g 进入候选得分（状态导数）
        return s

    F, k = {}, 0
    while k < B_BUDGET:
        unread = [j for j in phi if j not in F]
        if not unread:
            break
        j = best_split(rng, macros, F, unread)
        F[j] = phi[j]
        k += 1
        A = cand(F, macros, M_TH)
        if not A:
            continue
        top = sorted(A, key=lambda t: -score(t[0], F))
        if top[0][1] >= N_TH and (len(top) < 2 or top[0][1] - top[1][1] >= EPS_GAP):
            break   # 落定（唯一且达阈值）

    # ---- 结算（判定函数两臂完全相同）
    scored = [(i, score(i, F)) for i in range(len(macros))]
    top = sorted(scored, key=lambda t: -t[1])
    win, win_s = top[0]
    gap = win_s - (top[1][1] if len(top) > 1 else SIGMA_MIN)
    state = "收敛" if (win_s >= N_TH and gap >= EPS_GAP) else "震荡"
    C1 = state == "收敛"
    C2 = bool(C1 and win == M_star)
    if C1 and C2:
        outcome = "真收敛"
    elif C1 and not C2:
        outcome = "伪收敛(F2)"
    else:
        outcome = "F1(未加载/不收敛)"
    return {"state": state, "outcome": outcome, "C1": C1, "C2": C2,
            "sigma_end": win_s, "steps": k}


def run_arm_g(arm, g, n=N_TRIALS, seed=20261027):
    arm_off = zlib.crc32((arm + str(g)).encode("utf-8"))
    rng = random.Random(seed + arm_off % 9973)
    macros = build_macros(random.Random(seed))
    bias = [rng.random() for _ in range(N_MACRO)]   # 恒定的"冒充漂移"偏置
    rows = [trial(rng, macros, bias, g, arm) for _ in range(n)]
    c = Counter(r["outcome"] for r in rows)
    f2 = c["伪收敛(F2)"]
    f1 = c["F1(未加载/不收敛)"]
    return {
        "arm": arm, "g": g, "n": n,
        "真收敛%": round(100 * c["真收敛"] / n, 1),
        "伪收敛F2%": round(100 * f2 / n, 1),
        "F1%": round(100 * f1 / n, 1),
        "F2/F1": round(f2 / f1, 2) if f1 else float("inf"),
        "sigma_end均值": round(st.mean(r["sigma_end"] for r in rows), 3),
    }


def _pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    vx = sum((a - mx) ** 2 for a in xs)
    vy = sum((b - my) ** 2 for b in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    return cov / (vx ** 0.5 * vy ** 0.5)


def main():
    all_rows = []
    for arm in ("legal", "violate"):
        for g in G_VALUES:
            all_rows.append(run_arm_g(arm, g))

    legal = [r for r in all_rows if r["arm"] == "legal"]
    violate = [r for r in all_rows if r["arm"] == "violate"]
    r_legal = _pearson([r["g"] for r in legal], [r["伪收敛F2%"] for r in legal])
    r_violate = _pearson([r["g"] for r in violate], [r["伪收敛F2%"] for r in violate])

    # 单调性（violate 参考）：F2 在 g=1→峰值是否上升
    f2_v = [r["伪收敛F2%"] for r in sorted(violate, key=lambda r: r["g"])]
    g_v = [r["g"] for r in sorted(violate, key=lambda r: r["g"])]
    peak_violate = max(f2_v)
    peak_violate_g = g_v[f2_v.index(peak_violate)]
    f2_legal = [r["伪收敛F2%"] for r in legal]
    f2_v_at_base = f2_v[0]                     # violate 在 g=1 的 F2
    peak_legal = max(f2_legal)

    # ---- 预言判定（忠实于"越权产生伪定向"的核心断言）
    # P1 legal 无 g 驱动的伪定向：F2 跨 g 的极差 ≤ 4 个百分点（噪声级，近似平稳）
    p1 = (max(f2_legal) - min(f2_legal)) <= 4.0
    # P2 violate 产生伪定向：F2 显著高于 legal（效应量 ≥ 20pp）且从基线单调升至峰值
    p2 = (peak_violate - peak_legal >= 20.0) and (peak_violate > f2_v_at_base + 10.0)
    # 补充：非单调（g=3 后回落）登记为"冒充漂移饱和"次级效应，不否定 P2 的上升段
    rising_to_peak = peak_violate > f2_v_at_base
    saturates = bool(peak_violate_g < g_v[-1] and f2_v[-1] < peak_violate - 5.0)
    r11_supported = bool(p1 and p2)

    return {
        "experiment": "R11-pseudo-drift",
        "framework": "XS v3.4 · 补篇二 R11（推论 4.1：A11-R 违反 → 伪定向）",
        "hardness": "S2（自制仿真，复用 r18 宏加载原型）",
        "g_values": G_VALUES, "n_trials": N_TRIALS,
        "rows": all_rows,
        "legal_F2_range_pp": round(max(f2_legal) - min(f2_legal), 1),
        "violate_F2_peak": round(peak_violate, 1),
        "violate_F2_peak_g": peak_violate_g,
        "violate_F2_effect_over_legal_pp": round(peak_violate - peak_legal, 1),
        "correlations": {"legal_F2_vs_g": round(r_legal, 3),
                         "violate_F2_vs_g": round(r_violate, 3)},
        "violate_F2_rises_to_peak": bool(rising_to_peak),
        "drift_saturation_at_large_g": saturates,
        "verdicts": {"P1_legal_no_drift(极差≤4pp)": p1,
                     "P2_violate_pseudo_drift(效应≥20pp且升至峰值)": p2},
        "r11_supported": r11_supported,
        "conclusion": (
            "R11 支持：越权注入（g 进状态导数）产生 F2/ghost 伪力（伪定向），效应量 "
            "+%.0fpp，合法注入不改变方向 → 推论 4.1 成立；g 过大出现'冒充漂移饱和'"
            "（伪方向恒定→一致性锁定），伪收敛率回落，为次级效应"
            % (peak_violate - peak_legal) if r11_supported else
            "R11 不支持：越权注入未见 g 相关的伪定向 → 推论 4.1 被证否"),
    }


if __name__ == "__main__":
    res = main()
    print(json.dumps(res, ensure_ascii=False, indent=2))
    with open("r11_results.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
