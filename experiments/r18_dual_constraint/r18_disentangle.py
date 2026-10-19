# -*- coding: utf-8 -*-
"""
R18 拆解实验：把 kappa 里被我混在一起的两个子机制分开
========================================================
r18_sweep.py 的 R18-b 把"内约束强度 kappa"同时施加在两处：
    (a) 注入策略的最优性（是否按推论 4.4 做最优消歧注入）
    (b) 停止标准（tau 与 gap 要求）
结果**不是单调的**，而是阶跃的。这可能是真发现，也可能只是两个机制的
阈值被我叠在了一起。本脚本把它们分开扫，以判定属于哪一种。

  D1  kappa 细扫 0.0 -> 1.0（步长 0.1），确认阶跃形状
  D2  只动注入最优性 p_opt（tau=n, gap=eps 固定）
  D3  只动停止标准 tau（注入恒为最优）
"""

import json
import random
import statistics as st
from collections import Counter

from r18_sweep import CFG, build_macros, sigma, cand, best_split, observe, summarize


def trial(rng, macros, cfg, p_opt=1.0, tau=None, gap_req=None, p_swap=0.0):
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

    if tau is None:
        tau = cfg["N_TH"]
    if gap_req is None:
        gap_req = cfg["EPS_GAP"]

    F, k = {}, 0
    while k < cfg["B_BUDGET"]:
        unread = [j for j in phi if j not in F]
        if not unread:
            break
        if rng.random() < p_opt:
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


if __name__ == "__main__":
    n = 600
    rng = random.Random(23)
    macros = build_macros(random.Random(23), CFG)
    m, nn, eps = CFG["M_TH"], CFG["N_TH"], CFG["EPS_GAP"]
    out = {}

    out["D1_kappa细扫"] = {}
    for i in range(11):
        kk = i / 10
        tau = m + kk * (nn - m)
        gr = kk * eps
        rows = [trial(rng, macros, CFG, p_opt=kk, tau=tau, gap_req=gr) for _ in range(n)]
        out["D1_kappa细扫"][f"{kk:.1f}"] = summarize(rows)

    out["D2_只动注入最优性"] = {}
    for po in (0.0, 0.25, 0.5, 0.75, 1.0):
        rows = [trial(rng, macros, CFG, p_opt=po, tau=nn, gap_req=eps) for _ in range(n)]
        out["D2_只动注入最优性"][f"p_opt={po}"] = summarize(rows)

    out["D3_只动停止标准tau"] = {}
    for f in (0.0, 0.25, 0.5, 0.75, 1.0):
        tau = m + f * (nn - m)
        rows = [trial(rng, macros, CFG, p_opt=1.0, tau=tau, gap_req=eps) for _ in range(n)]
        out["D3_只动停止标准tau"][f"tau={tau:.3f}"] = summarize(rows)

    print(json.dumps(out, ensure_ascii=False, indent=2))
    with open("r18_disentangle.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
