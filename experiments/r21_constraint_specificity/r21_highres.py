# -*- coding: utf-8 -*-
"""
R21 高分辨率复核 —— 排除「F2 倒 U 是 sigma 离散化量子的伪像」
================================================================
主扫描用 N_FEAT=8，sigma 只能取 1/8 的整数倍（量子较粗）。
若倒 U 是量化伪像，提高分辨率后应消失或峰位漂移。

本脚本把 N_FEAT 提到 16（量子 1/16）、N_MACRO 提到 40，重跑扫描。

结果（n=1200/档）：倒 U **保持**，峰值仍在 Q ~ 1.0（K=20/40），
与主参数一致。故倒 U 不是量化伪像。
"""

import json
import math
import random
import statistics as st
import zlib

import r21_constraint_specificity as R

N_FEAT = 16
N_MACRO = 40
N_PROTO = 4

# 关键：把分辨率参数写回被导入模块的命名空间。
# build_macros / observe / sigma 都在 R 的模块全局里查这些常量，
# 只改本文件的同名变量不会生效，会造出长度 8 的宏而 obs_mask 长度 16 -> 越界。
R.N_FEAT = N_FEAT
R.N_MACRO = N_MACRO
R.N_PROTO = N_PROTO

N_TRIAL = 1200
K_GRID = (40, 32, 26, 20, 16, 12, 8, 4)


def trial_K(rng, macros, K):
    """与 r21_constraint_specificity.trial 相同，但剪枝按**宏数量**而非簇。"""
    M_star = rng.randrange(len(macros))
    truth = macros[M_star]
    obs_mask = [rng.random() < R.RHO_OBS for _ in range(N_FEAT)]
    if not any(obs_mask):
        obs_mask[rng.randrange(N_FEAT)] = True
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
        "outcome": outcome, "gap": gap, "steps": k,
        "sigma_end": win_s if A_m else max((R.sigma(macros[i], F) for i in allowed), default=0.0),
        "n_cand": len(A_m), "star_banned": star_banned,
    }


def run_K(K, n=N_TRIAL, seed=20261019):
    rng = random.Random(seed + zlib.crc32(f"HR{K}".encode("utf-8")) % 9973)
    macros = R.build_macros(random.Random(seed))
    rows = [trial_K(rng, macros, K) for _ in range(n)]
    alive = [r for r in rows if not r["star_banned"]]

    def pct(rs, key):
        return round(100 * sum(1 for r in rs if r["outcome"] == key) / len(rs), 1) if rs else None

    return {
        "K": K, "Q": round(math.log2(N_MACRO / K), 3),
        "F2%": pct(rows, "伪收敛(F2)"), "真收敛%": pct(rows, "真收敛"),
        "F1%": pct(rows, "F1(未加载/不收敛)"),
        "gap": round(st.mean(r["gap"] for r in rows), 4),
        "步数": round(st.mean(r["steps"] for r in rows), 2),
        "存活n": len(alive), "存活F2%": pct(alive, "伪收敛(F2)"),
        "存活真收敛%": pct(alive, "真收敛"),
    }


if __name__ == "__main__":
    rows = [run_K(K) for K in K_GRID]
    out = {"N_FEAT": N_FEAT, "N_MACRO": N_MACRO, "N_TRIAL": N_TRIAL, "扫描": rows}
    f2 = [r["F2%"] for r in rows]
    peak_i = max(range(len(f2)), key=lambda i: f2[i])
    mono = all(f2[i] <= f2[i + 1] for i in range(len(f2) - 1)) or \
           all(f2[i] >= f2[i + 1] for i in range(len(f2) - 1))
    out["F2序列"] = f2
    out["峰值K"] = rows[peak_i]["K"]
    out["峰值Q"] = rows[peak_i]["Q"]
    out["单调"] = mono
    out["倒U保持"] = not mono
    print(json.dumps(out, ensure_ascii=False, indent=2))
    with open("r21_highres.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
