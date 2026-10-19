# -*- coding: utf-8 -*-
"""
R18 补充实验：剂量-反应 与 参数鲁棒性
================================================
r18_dual_constraint.py 只做了二值对照（有/无约束），无法排除"分离只是二值伪像"。
本脚本把两个约束都做成**连续变量**，看分离是否随剂量单调演化：

  R18-b  内约束强度 kappa ∈ [0,1] 的剂量-反应
         kappa 同时控制两件事（都来自"内约束 = 必须理解"这一条）：
           (a) 注入策略的最优性：kappa=1 用最优消歧注入（推论 4.4），kappa=0 随机取
           (b) 停止标准：tau = m + kappa*(n-m)，gap 要求 = kappa*eps
             理解要求"唯一且达标"，续写只要求"够用"——kappa 在两者间插值
         预言：收敛率与终局 sigma 随 kappa 单调上升；F1 率随 kappa 单调下降。

  R18-c  外约束失效程度 p_swap ∈ [0,1] 的剂量-反应
         p_swap = 物理形态中被替换成"另一宏形态"的特征比例（对象被替换的程度）
         预言（本实验最关键的一条）：
           F2 率随 p_swap 单调上升，而**终局 sigma 保持平坦**。
         若成立，则"外约束是方向约束、不是强度约束"得到支持：
         外约束失效让系统照常到达同样高的确信度，只是到达了错的地方。

  鲁棒性  在若干组 (N_MACRO, RHO_OBS, P_NOISE, B) 下重复 A/B 分离检验，
         看 σ(A) ≈ σ(base) 且 σ(B) < σ(base) 这一签名是否稳定。
"""

import json
import random
import statistics as st
from collections import Counter

CFG = dict(
    N_FEAT=8, N_PROTO=4, N_MACRO=20,
    M_TH=0.35, N_TH=0.55, EPS_GAP=0.10,
    B_BUDGET=12, RHO_OBS=0.85, P_NOISE=0.05, TOPK=6,
)


def build_macros(rng, cfg):
    nf, nproto, nm = cfg["N_FEAT"], cfg["N_PROTO"], cfg["N_MACRO"]
    protos = [[rng.randint(0, 1) for _ in range(nf)] for _ in range(nproto)]
    seen, macros = set(), []
    while len(macros) < nm:
        p = protos[len(macros) % nproto]
        v = p[:]
        for _ in range(2):
            j = rng.randrange(nf)
            v[j] = 1 - v[j]
        t = tuple(v)
        if t not in seen:
            seen.add(t)
            macros.append(v)
    return macros


def sigma(M, F, nf):
    s = 0
    for j, val in F.items():
        s += 1 if M[j] == val else -1
    return s / nf


def cand(F, macros, th, nf):
    return [(i, sigma(M, F, nf)) for i, M in enumerate(macros) if sigma(M, F, nf) >= th]


def best_split(rng, macros, F, unread, nf, topk):
    cur = sorted(range(len(macros)), key=lambda i: -sigma(macros[i], F, nf))[:topk]
    bj, bs = None, -1.0
    for j in unread:
        n1 = sum(1 for i in cur if macros[i][j] == 1)
        n0 = len(cur) - n1
        s = min(n0, n1) + 0.001 * rng.random()
        if s > bs:
            bj, bs = j, s
    return bj


def observe(rng, src, mask, p_noise):
    phi = {}
    for j, on in enumerate(mask):
        if not on:
            continue
        v = src[j]
        if rng.random() < p_noise:
            v = 1 - v
        phi[j] = v
    return phi


def trial(rng, macros, cfg, kappa=1.0, p_swap=0.0):
    nf = cfg["N_FEAT"]
    n = len(macros)
    M_star = rng.randrange(n)
    truth = macros[M_star]

    mask = [rng.random() < cfg["RHO_OBS"] for _ in range(nf)]
    if not any(mask):
        mask[rng.randrange(nf)] = True
    phi = observe(rng, truth, mask, cfg["P_NOISE"])

    if p_swap > 0:
        other = rng.choice([i for i in range(n) if i != M_star])
        for j in list(phi):
            if rng.random() < p_swap:
                v = macros[other][j]
                if rng.random() < cfg["P_NOISE"]:
                    v = 1 - v
                phi[j] = v

    tau = cfg["M_TH"] + kappa * (cfg["N_TH"] - cfg["M_TH"])
    gap_req = kappa * cfg["EPS_GAP"]

    F, k = {}, 0
    while k < cfg["B_BUDGET"]:
        unread = [j for j in phi if j not in F]
        if not unread:
            break
        if rng.random() < kappa:
            j = best_split(rng, macros, F, unread, nf, cfg["TOPK"])
        else:
            j = rng.choice(unread)
        F[j] = phi[j]
        k += 1

        A = cand(F, macros, cfg["M_TH"], nf)
        if not A:
            continue
        top = sorted(A, key=lambda t: -t[1])
        gap = top[0][1] - (top[1][1] if len(top) > 1 else -1.0)
        if top[0][1] >= tau and gap >= gap_req:
            break

    A_m = cand(F, macros, cfg["M_TH"], nf)
    if not A_m:
        state, win, win_s = "发散", None, max(sigma(M, F, nf) for M in macros)
    else:
        top = sorted(A_m, key=lambda t: -t[1])
        win, win_s = top[0]
        gap = win_s - (top[1][1] if len(top) > 1 else -1.0)
        state = "收敛" if (win_s >= cfg["N_TH"] and gap >= cfg["EPS_GAP"]) else "震荡"

    C1 = state == "收敛"
    C2 = bool(C1 and win == M_star)
    outcome = "真收敛" if (C1 and C2) else ("伪收敛(F2)" if C1 else "F1")
    return dict(outcome=outcome, state=state, sigma_end=win_s, steps=k)


def summarize(rows):
    n = len(rows)
    c = Counter(r["outcome"] for r in rows)
    sig = [r["sigma_end"] for r in rows]
    f2, f1 = c["伪收敛(F2)"], c["F1"]
    return {
        "真收敛%": round(100 * c["真收敛"] / n, 1),
        "F2%": round(100 * f2 / n, 1),
        "F1%": round(100 * f1 / n, 1),
        "sigma_end": round(st.mean(sig), 3),
        "steps": round(st.mean([r["steps"] for r in rows]), 2),
    }


def sweep(cfg, n=600, seed=7):
    rng = random.Random(seed)
    macros = build_macros(random.Random(seed), cfg)
    out = {}

    out["R18b_内约束强度kappa"] = {}
    for kk in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        rows = [trial(rng, macros, cfg, kappa=kk, p_swap=0.0) for _ in range(n)]
        out["R18b_内约束强度kappa"][f"kappa={kk}"] = summarize(rows)

    out["R18c_外约束失效程度p_swap"] = {}
    for ps in (0.0, 0.25, 0.5, 0.75, 1.0):
        rows = [trial(rng, macros, cfg, kappa=1.0, p_swap=ps) for _ in range(n)]
        out["R18c_外约束失效程度p_swap"][f"p_swap={ps}"] = summarize(rows)
    return out


def robustness(n=500, seed=11):
    grid = [
        dict(N_MACRO=20, RHO_OBS=0.85, P_NOISE=0.05, B_BUDGET=12),
        dict(N_MACRO=40, RHO_OBS=0.85, P_NOISE=0.05, B_BUDGET=12),
        dict(N_MACRO=20, RHO_OBS=0.65, P_NOISE=0.05, B_BUDGET=12),
        dict(N_MACRO=20, RHO_OBS=0.85, P_NOISE=0.15, B_BUDGET=12),
        dict(N_MACRO=20, RHO_OBS=0.85, P_NOISE=0.05, B_BUDGET=20),
        dict(N_MACRO=60, RHO_OBS=0.75, P_NOISE=0.10, B_BUDGET=16),
    ]
    res = {}
    for i, upd in enumerate(grid):
        cfg = dict(CFG)
        cfg.update(upd)
        rng = random.Random(seed + i)
        macros = build_macros(random.Random(seed + i), cfg)
        base = summarize([trial(rng, macros, cfg, kappa=1.0, p_swap=0.0) for _ in range(n)])
        armA = summarize([trial(rng, macros, cfg, kappa=1.0, p_swap=1.0) for _ in range(n)])
        armB = summarize([trial(rng, macros, cfg, kappa=0.0, p_swap=0.0) for _ in range(n)])
        res[str(upd)] = {
            "base": base,
            "A_外失效": armA,
            "B_内失效": armB,
            "签名成立": bool(
                armA["F2%"] > base["F2%"] and armB["F1%"] > base["F1%"]
                and armB["sigma_end"] < base["sigma_end"] < armA["sigma_end"] + 0.02
                and abs(armA["sigma_end"] - base["sigma_end"]) < 0.05
            ),
        }
    return res


if __name__ == "__main__":
    rep = {"主参数下扫描": sweep(CFG), "鲁棒性网格": robustness()}
    print(json.dumps(rep, ensure_ascii=False, indent=2))
    with open("r18_sweep.json", "w", encoding="utf-8") as f:
        json.dump(rep, f, ensure_ascii=False, indent=2)
