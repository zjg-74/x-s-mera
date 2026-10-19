# -*- coding: utf-8 -*-
"""
R21 细粒度扫描（与 r21_constraint_specificity 同分辨率，但剪枝粒度更细）
========================================================================
主脚本 r21_constraint_specificity.py 按**簇**剪枝，只有 4 档（n_ban = 0..3）。
4 个点不足以支撑「F2 非单调」的主张——3 个点也可能因噪声呈现假的非单调。
本脚本改为按**宏数量** K 剪枝，扫 10 档（K = 20..2），n = 1500/档。

结果：F2 序列 9.9 -> 15.3 -> 16.7 -> 20.4 -> 23.9 -> 23.9 -> 22.7 -> 19.3 -> 18.7 -> 11.2，
峰值在 Q ~ 1.0（K = 10~12），倒 U 明确。
"""

import json
import math
import random
import statistics as st
import zlib

import r21_constraint_specificity as R

N_TRIAL = 1500
K_GRID = (20, 18, 16, 14, 12, 10, 8, 6, 4, 2)


def trial_K(rng, macros, K):
    """剪枝按**宏数量** K：allowed = 随机 K 个宏。"""
    M_star = rng.randrange(len(macros))
    truth = macros[M_star]
    obs_mask = [rng.random() < R.RHO_OBS for _ in range(R.N_FEAT)]
    if not any(obs_mask):
        obs_mask[rng.randrange(R.N_FEAT)] = True
    phi = R.observe(rng, truth, obs_mask)

    allowed = sorted(rng.sample(range(len(macros)), K))
    banned = set(range(len(macros))) - set(allowed)
    star_banned = M_star in banned

    F, k = {}, 0
    while k < R.B_BUDGET:
        unread = [j for j in phi if j not in F]
        if not unread:
            break
        j = R.best_split(rng, macros, F, unread, allowed)
        F[j] = phi[j]
        k += 1
        A = R.cand(F, macros, R.M_TH, allowed)
        if not A:
            continue
        top = sorted(A, key=lambda t: -t[1])
        gap = top[0][1] - (top[1][1] if len(top) > 1 else R.SIGMA_MIN)
        if top[0][1] >= R.N_TH and gap >= R.EPS_GAP:
            break

    A_m = R.cand(F, macros, R.M_TH, allowed)
    if not A_m:
        state, win, win_s, gap = "发散", None, 0.0, 0.0
    else:
        top = sorted(A_m, key=lambda t: -t[1])
        win, win_s = top[0]
        gap = win_s - (top[1][1] if len(top) > 1 else R.SIGMA_MIN)
        state = "收敛" if (win_s >= R.N_TH and gap >= R.EPS_GAP) else "震荡"

    C1 = state == "收敛"
    C2 = bool(C1 and win == M_star)
    outcome = "真收敛" if (C1 and C2) else ("伪收敛(F2)" if C1 else "F1(未加载/不收敛)")
    return {
        "outcome": outcome, "gap": round(gap, 4), "steps": k,
        "sigma_end": win_s if A_m else max((R.sigma(macros[i], F) for i in allowed), default=0.0),
        "n_cand": len(A_m), "star_banned": star_banned,
    }


def run_K(K, n=N_TRIAL, seed=20261019):
    rng = random.Random(seed + zlib.crc32(f"K{K}".encode("utf-8")) % 9973)
    macros = R.build_macros(random.Random(seed))
    rows = [trial_K(rng, macros, K) for _ in range(n)]
    alive = [r for r in rows if not r["star_banned"]]

    def pct(rs, key):
        return round(100 * sum(1 for r in rs if r["outcome"] == key) / len(rs), 1) if rs else None

    return {
        "K": K, "Q": round(math.log2(R.N_MACRO / K), 3),
        "F2%": pct(rows, "伪收敛(F2)"), "真收敛%": pct(rows, "真收敛"),
        "F1%": pct(rows, "F1(未加载/不收敛)"),
        "gap": round(st.mean(r["gap"] for r in rows), 4),
        "步数": round(st.mean(r["steps"] for r in rows), 2),
        "M*被剪%": round(100 * sum(r["star_banned"] for r in rows) / n, 1),
        "存活子样本n": len(alive),
        "存活_F2%": pct(alive, "伪收敛(F2)"), "存活_真收敛%": pct(alive, "真收敛"),
        "存活_gap": round(st.mean(r["gap"] for r in alive), 4) if alive else None,
    }


if __name__ == "__main__":
    rows = [run_K(K) for K in K_GRID]
    f2 = [r["F2%"] for r in rows]
    peak_i = max(range(len(f2)), key=lambda i: f2[i])
    mono = all(f2[i] <= f2[i + 1] for i in range(len(f2) - 1)) or \
           all(f2[i] >= f2[i + 1] for i in range(len(f2) - 1))
    out = {"N_TRIAL": N_TRIAL, "扫描": rows,
           "F2序列": f2, "峰值K": rows[peak_i]["K"], "峰值Q": rows[peak_i]["Q"],
           "单调": mono, "倒U": not mono}
    print(json.dumps(out, ensure_ascii=False, indent=2))
    with open("r21_sweep.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
