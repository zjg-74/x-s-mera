# -*- coding: utf-8 -*-
"""
R21 —— 内约束强度的重新操作化：措辞特异性（取代已撤回的 R18-b）
====================================================================
对应 XS v3.1 补篇五（B-B21）与补篇三 定义 3.2（Q = 操作空间剪枝量）。

为什么撤回 R18-b
--------------------------------------------------------------------
R18-b 用两个代理量编码"内约束强度"，拆解后发现都不成立：
  * 接受阈值 tau    —— 同义反复（tau 低于阈值当然不收敛）
  * 注入最优性 p_opt —— 在预算宽裕时完全无效应（68/70/68/70/70）
两者都不是"内约束强度"的操作化。

新的操作化：措辞特异性 = 允许操作集 O_adm 的剪枝量
--------------------------------------------------------------------
内约束是一条**禁止性指令**。它的强度不是"要求有多严"，
而是**它排除了多少种操作**：

    W0  「对这句话做任何事」              |O_adm| = 20   Q = 0
    W1  「理解这句话」                     |O_adm| = 15   Q = 0.415
    W2  「理解，且不要按字面义理解」         |O_adm| = 10   Q = 1.000
    W3  「理解，且不要字面义 / 不要引申义」   |O_adm| =  5   Q = 2.000

    Q := log2(|O_all| / |O_adm|)           —— 补篇三 定义 3.2

关键设计约束（避免重蹈 R18-b 的同义反复）
--------------------------------------------------------------------
**不改停止规则、不改接受阈值、不改判定函数。**
只改 O_adm（可用宏集）。这样测到的任何效应都来自"可操作对象变少"，
而不是"判定变松/变紧"。四臂的判定函数完全相同。

两种剪枝模式（对照）
--------------------------------------------------------------------
  cluster：按语义簇剪枝——禁止整簇宏，模拟"不要按 X 义理解"（措辞式）
  random ：随机剪掉同等数量的宏（数量对照）
若两者效应相同 -> 效应只是"候选变少"的平凡数量效应；
若不同       -> "措辞的结构性"（禁止的是一整类理解方式）有独立作用。

非平凡预言（核心）
--------------------------------------------------------------------
剪枝会把"真值宏 M*"一起剪掉，这本身必然抬高 F2 —— 这是**设计上的平凡效应**，
不能算发现。因此必须分层：

  【P1 平凡层】全样本 F2 率随剪枝单调上升（P(M* 被剪) = n_ban / N_PROTO）

  【P2 非平凡层】**条件于 M* 未被剪掉**的子样本中：
      剪枝制造"虚假唯一性"—— 唯一性间隙 gap 增大，不是因为证据更充分，
      而是因为第二名候选被**人为排除**了。
      预言：该子样本的 gap 单调增大、收敛判定率上升、而 F2 率**不下降**
      （甚至上升）。若 F2 率下降到 0，则"虚假唯一性"主张被证否。

  【P3 端点层】|O_adm| 极小时 A_k = ∅ -> 发散率上升
      （补篇二 推论 3.1 的分叉：A_k ≠ ∅ -> 伪收敛；A_k = ∅ -> 发散）
"""

import json
import math
import random
import statistics as st
import zlib
from collections import Counter

# ------------------------------------------------------------------ 参数（与 R18 主实验对齐）
N_FEAT = 8
N_PROTO = 4       # 语义簇数（= 原型数）
N_MACRO = 20
M_TH = 0.35
N_TH = 0.55
EPS_GAP = 0.10
B_BUDGET = 12
RHO_OBS = 0.85
P_NOISE = 0.05
TOPK = 6
SIGMA_MIN = -1.0

N_TRIAL = 1200


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


def cluster_of(i):
    """宏 i 所属的语义簇（build_macros 按 i % N_PROTO 派生）。"""
    return i % N_PROTO


def sigma(M, F):
    s = 0
    for j, val in F.items():
        s += 1 if M[j] == val else -1
    return s / N_FEAT


def cand(F, macros, th, allowed):
    return [(i, sigma(macros[i], F)) for i in allowed if sigma(macros[i], F) >= th]


def best_split(rng, macros, F, unread, allowed):
    cur = sorted(allowed, key=lambda i: -sigma(macros[i], F))[:TOPK]
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


def trial(rng, macros, n_ban, mode="cluster"):
    """
    n_ban：禁止的簇数 / 禁止的宏比例（措辞特异性层级）
    mode ：cluster = 按簇剪枝（措辞式）；random = 随机剪枝（数量对照）
    """
    M_star = rng.randrange(len(macros))
    truth = macros[M_star]

    obs_mask = [rng.random() < RHO_OBS for _ in range(N_FEAT)]
    if not any(obs_mask):
        obs_mask[rng.randrange(N_FEAT)] = True
    phi = observe(rng, truth, obs_mask)

    # ---- 内约束的措辞特异性：剪掉一部分操作（宏）
    if n_ban <= 0:
        banned = set()
    elif mode == "cluster":
        # 禁止整簇：模拟「不要按 X 义理解」
        ban_clusters = rng.sample(range(N_PROTO), n_ban)
        banned = {i for i in range(len(macros)) if cluster_of(i) in ban_clusters}
    else:
        # 随机剪掉同等数量（每簇 5 个，故 n_ban 簇 = 5*n_ban 个宏）
        banned = set(rng.sample(range(len(macros)), min(N_MACRO - 1, 5 * n_ban)))
    allowed = [i for i in range(len(macros)) if i not in banned]
    star_banned = M_star in banned

    # ---- 迭代（停止规则与 R18 主实验的「理解」臂**完全相同**，未改 tau）
    F, k = {}, 0
    while k < B_BUDGET:
        unread = [j for j in phi if j not in F]
        if not unread:
            break
        j = best_split(rng, macros, F, unread, allowed)
        F[j] = phi[j]
        k += 1
        A = cand(F, macros, M_TH, allowed)
        if not A:
            continue
        top = sorted(A, key=lambda t: -t[1])
        gap = top[0][1] - (top[1][1] if len(top) > 1 else SIGMA_MIN)
        if top[0][1] >= N_TH and gap >= EPS_GAP:
            break

    # ---- 结算（判定函数对所有层级完全相同）
    A_m = cand(F, macros, M_TH, allowed)
    if not A_m:
        state, win, win_s, gap = "发散", None, 0.0, 0.0
    else:
        top = sorted(A_m, key=lambda t: -t[1])
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

    return {
        "outcome": outcome, "state": state,
        "sigma_end": win_s if A_m else max((sigma(macros[i], F) for i in allowed), default=0.0),
        "gap": round(gap, 4), "steps": k,
        "n_cand": len(A_m), "star_banned": star_banned,
        "n_allowed": len(allowed),
    }


def run(n_ban, mode="cluster", n=N_TRIAL, seed=20261019):
    off = zlib.crc32(f"{mode}:{n_ban}".encode("utf-8")) % 9973
    rng = random.Random(seed + off)
    macros = build_macros(random.Random(seed))
    rows = [trial(rng, macros, n_ban, mode) for _ in range(n)]

    alive = [r for r in rows if not r["star_banned"]]   # 条件子样本：真值宏未被剪

    def summarize(rs, tag):
        if not rs:
            return {"tag": tag, "n": 0}
        c = Counter(r["outcome"] for r in rs)
        return {
            "tag": tag, "n": len(rs),
            "真收敛%": round(100 * c["真收敛"] / len(rs), 1),
            "伪收敛F2%": round(100 * c["伪收敛(F2)"] / len(rs), 1),
            "F1%": round(100 * c["F1(未加载/不收敛)"] / len(rs), 1),
            "gap均值": round(st.mean(r["gap"] for r in rs), 4),
            "步数均值": round(st.mean(r["steps"] for r in rs), 2),
            "sigma_end均值": round(st.mean(r["sigma_end"] for r in rs), 3),
            "候选数均值": round(st.mean(r["n_cand"] for r in rs), 2),
        }

    n_all = rows[0]["n_allowed"]
    return {
        "n_ban": n_ban, "mode": mode,
        "|O_adm|": n_all,
        "Q=log2(20/|O_adm|)": round(math.log2(N_MACRO / n_all), 3) if n_all else None,
        "全样本": summarize(rows, "全样本"),
        "真值未被剪子样本": summarize(alive, "M*未被剪"),
        "M*被剪%": round(100 * sum(r["star_banned"] for r in rows) / len(rows), 1),
    }


if __name__ == "__main__":
    out = {"cluster剪枝": [run(k, "cluster") for k in (0, 1, 2, 3)],
           "random剪枝对照": [run(k, "random") for k in (0, 1, 2, 3)]}
    print(json.dumps(out, ensure_ascii=False, indent=2))
    with open("r21_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
