# -*- coding: utf-8 -*-
"""
R18 —— 双约束失效的可分离性检验
================================================
对应 XS v3.1 补篇三：
  推论 2.1  外约束失效 -> F2（加载成功但事后为错）
            内约束失效 -> F1（加载失败 / 不收敛）
  推论 2.2  F1 <-> 内约束失效；F2 <-> 外约束失效
  定义 3.1/3.2  sigma（加载度）为中介变量
  定理 4.1  T_conv ~ max(0, Psi) / I_eff
  推论 4.4  最优注入策略 = argmax I(M; f | F_k)

本实验额外检验一条**新判别签名**：
  R18-a  两条件的**终局加载度 sigma 应显著不同**——
         外失效：sigma 照常爬升至高位（高确信、错结果）
         内失效：sigma 停滞在低位（低确信、无结果）
         若两条件 sigma 无差异，则"sigma 为中介"的机制被证否，
         推论 2.1 退化为一张无机制的对应表。

设计要点（避免循环论证）
--------------------------------------------------
两个条件**只改各自那一个约束**，不改落定规则、不改阈值、不改判定函数：
  * 外约束（物理约束）失效 = **对象被替换**：系统诚实读取，但读到的物理形态
    已被换成另一个宏的形态（phi -> phi'）。系统照常努力填槽。
  * 内约束（理解约束）失效 = **任务被替换**：物理形态不变、读取照常诚实，
    但任务由"理解"改为"续写"。由此**停止条件**改变：
        理解 -> 停止条件 = 唯一且达到接受阈值 n（需要消歧）
        续写 -> 停止条件 = 存在任一宏达到可用阈值 m 即可（不要求唯一）
    F1/F2 不是被指定的，而是由 sigma 的演化路径**导出**的。

四臂
--------------------------------------------------
  base  正常 phi + 理解        -> 预期：真收敛主导
  A     篡改 phi + 理解        -> 预期：F2（伪收敛）主导，终局 sigma 高
  B     正常 phi + 续写        -> 预期：F1 主导，终局 sigma 低
  C     篡改 phi + 续写        -> 预期：F1 主导（内约束的早停截断错误累积）
"""

import json
import math
import random
import zlib
import statistics as st
from collections import Counter

# ------------------------------------------------------------------ 参数
N_FEAT = 8        # 句子物理形态的特征维度
N_PROTO = 4       # 宏库原型数（产生歧义簇）
N_MACRO = 20      # 宏库规模
M_TH = 0.35       # m：可用阈值（够用即可，不构成接受）
N_TH = 0.55       # n：接受阈值（落定）
EPS_GAP = 0.10    # eps：唯一性间隙
B_BUDGET = 12     # B：迭代预算
RHO_OBS = 0.85    # 物理形态的观测率（省略 = 1 - RHO_OBS）
P_NOISE = 0.05    # 载体噪声
TOPK = 6          # 最优注入时考察的头部候选数

SIGMA_MIN = -1.0


def build_macros(rng):
    """宏库：4 个原型，每个原型派生若干宏（扰动 2 个特征）-> 簇内歧义，簇间可分。"""
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
    """加载度：已填槽中一致记 +1、冲突记 -1（补篇一 §3.7：被证伪必须下降）。"""
    s = 0
    for j, val in F.items():
        s += 1 if M[j] == val else -1
    return s / N_FEAT


def cand(F, macros, th):
    return [(i, sigma(M, F)) for i, M in enumerate(macros) if sigma(M, F) >= th]


def best_split(rng, macros, F, unread):
    """最优注入策略（推论 4.4）：选能把头部候选集最接近二分的那个特征。"""
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
    """按物理形态取观测：src_vec 的掩码部分 + 载体噪声。"""
    phi = {}
    for j, on in enumerate(obs_mask):
        if not on:
            continue
        v = src_vec[j]
        if rng.random() < P_NOISE:
            v = 1 - v
        phi[j] = v
    return phi


def trial(rng, macros, arm):
    M_star = rng.randrange(len(macros))
    truth = macros[M_star]

    # ---- 物理形态 phi（外约束的载体）
    obs_mask = [rng.random() < RHO_OBS for _ in range(N_FEAT)]
    if not any(obs_mask):
        obs_mask[rng.randrange(N_FEAT)] = True

    outer_ok = arm not in ("A", "C")
    if outer_ok:
        phi = observe(rng, truth, obs_mask)
    else:
        # 外约束失效：对象被替换 —— 换成另一个宏的物理形态
        other = rng.choice([i for i in range(len(macros)) if i != M_star])
        obs_mask = [rng.random() < RHO_OBS for _ in range(N_FEAT)]
        if not any(obs_mask):
            obs_mask[rng.randrange(N_FEAT)] = True
        phi = observe(rng, macros[other], obs_mask)

    task = "generate" if arm in ("B", "C") else "understand"

    # ---- 迭代
    F, k = {}, 0
    while k < B_BUDGET:
        unread = [j for j in phi if j not in F]
        if not unread:
            break
        j = best_split(rng, macros, F, unread) if task == "understand" else rng.choice(unread)
        F[j] = phi[j]
        k += 1

        A = cand(F, macros, M_TH)
        if not A:
            continue
        top = sorted(A, key=lambda t: -t[1])
        best_s = top[0][1]
        if task == "understand":
            # 停止条件：唯一且达到接受阈值
            gap = top[0][1] - (top[1][1] if len(top) > 1 else SIGMA_MIN)
            if top[0][1] >= N_TH and gap >= EPS_GAP:
                break
        else:
            # 停止条件：够用即可
            if best_s >= M_TH:
                break

    # ---- 结算（判定函数对四臂完全相同）
    A_m = cand(F, macros, M_TH)
    if not A_m:
        state = "发散"
        win, win_s, gap = None, 0.0, 0.0
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
        "state": state,
        "outcome": outcome,
        "C1": C1,
        "C2": C2,
        "sigma_end": win_s if A_m else max((sigma(M, F) for M in macros), default=0.0),
        "steps": k,
        "n_cand": len(A_m),
    }


def run_arm(arm, n=800, seed=20261019):
    # 注意：不得用 hash(arm)——str 哈希受 PYTHONHASHSEED 随机化影响，
    # 会导致同一 seed 在不同进程产出不同结果（入库前实测：两次运行 json 的
    # md5 不一致）。改用 crc32，保证跨进程、跨平台可复现。
    arm_off = zlib.crc32(arm.encode("utf-8"))
    rng = random.Random(seed + arm_off % 9973)
    macros = build_macros(random.Random(seed))
    rows = [trial(rng, macros, arm) for _ in range(n)]
    c = Counter(r["outcome"] for r in rows)
    sig = [r["sigma_end"] for r in rows]
    steps = [r["steps"] for r in rows]
    stc = Counter(r["state"] for r in rows)
    f2 = c["伪收敛(F2)"]
    f1 = c["F1(未加载/不收敛)"]
    return {
        "arm": arm,
        "n": n,
        "真收敛%": round(100 * c["真收敛"] / n, 1),
        "伪收敛F2%": round(100 * f2 / n, 1),
        "F1%": round(100 * f1 / n, 1),
        "F2/F1": round(f2 / f1, 2) if f1 else float("inf"),
        "sigma_end均值": round(st.mean(sig), 3),
        "sigma_end标准差": round(st.pstdev(sig), 3),
        "步数均值": round(st.mean(steps), 2),
        "终局态分布": dict(stc),
    }


if __name__ == "__main__":
    res = [run_arm(a) for a in ("base", "A", "B", "C")]
    print(json.dumps(res, ensure_ascii=False, indent=2))
    with open("r18_results.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
