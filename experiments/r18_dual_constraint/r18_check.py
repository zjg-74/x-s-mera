# -*- coding: utf-8 -*-
"""
R18 自检：两个可疑结果的证伪/证实
====================================
可疑点 1  D1/D3 的"阶跃阈值"是不是 sigma 离散化伪像？
          N_FEAT=8 -> sigma 只能取 1/8=0.125 的整数倍；而接受阈值 n=0.55
          恰好落在 0.500 与 0.625 两个量子之间。于是任何 tau>0.50 都迫使
          系统一路涨到 0.625，从而必然越过 n。
          检验：把 N_FEAT 提到 16（量子 0.0625），重扫 kappa，看阶跃是否消失。

可疑点 2  D2 显示"最优注入策略无效"是不是因为预算不受约束？
          主参数下 B=12 而平均只需 ~5.6 步，资源远远不紧张，
          于是注入顺序无关紧要——反正最后都读完了。
          推论 4.4 预言的是**每一步的消歧效率**，只有当 B 接近所需步数时
          才应显出差别。
          检验：收紧 B，看 p_opt=1 是否开始优于 p_opt=0。
"""

import json
import random

from r18_sweep import build_macros, summarize
from r18_disentangle import trial

BASE = dict(
    N_FEAT=8, N_PROTO=4, N_MACRO=20,
    M_TH=0.35, N_TH=0.55, EPS_GAP=0.10,
    B_BUDGET=12, RHO_OBS=0.85, P_NOISE=0.05, TOPK=6,
)


def cfg16(**over):
    c = dict(BASE)
    c.update(N_FEAT=16, N_PROTO=6, N_MACRO=40, B_BUDGET=22, TOPK=8)
    c.update(over)
    return c


def kappa_sweep(cfg, n=600, seed=31):
    rng = random.Random(seed)
    macros = build_macros(random.Random(seed), cfg)
    m, nn, eps = cfg["M_TH"], cfg["N_TH"], cfg["EPS_GAP"]
    rows = {}
    for i in range(0, 11, 2):
        kk = i / 10
        tau = m + kk * (nn - m)
        r = [trial(rng, macros, cfg, p_opt=kk, tau=tau, gap_req=kk * eps) for _ in range(n)]
        rows[f"kappa={kk:.1f}"] = summarize(r)
    return rows


def budget_sweep(cfg, budgets, n=600, seed=37):
    rng = random.Random(seed)
    macros = build_macros(random.Random(seed), cfg)
    out = {}
    for b in budgets:
        c = dict(cfg)
        c["B_BUDGET"] = b
        r_opt = [trial(rng, macros, c, p_opt=1.0, tau=c["N_TH"], gap_req=c["EPS_GAP"]) for _ in range(n)]
        r_rnd = [trial(rng, macros, c, p_opt=0.0, tau=c["N_TH"], gap_req=c["EPS_GAP"]) for _ in range(n)]
        so, sr = summarize(r_opt), summarize(r_rnd)
        out[f"B={b}"] = {
            "最优注入": so,
            "随机注入": sr,
            "真收敛差值(最优-随机)": round(so["真收敛%"] - sr["真收敛%"], 1),
            "步数差值": round(so["steps"] - sr["steps"], 2),
        }
    return out


def ab_separation(cfg, n=600, seed=41):
    rng = random.Random(seed)
    macros = build_macros(random.Random(seed), cfg)
    base = summarize([trial(rng, macros, cfg, p_opt=1.0, tau=cfg["N_TH"],
                            gap_req=cfg["EPS_GAP"], p_swap=0.0) for _ in range(n)])
    armA = summarize([trial(rng, macros, cfg, p_opt=1.0, tau=cfg["N_TH"],
                            gap_req=cfg["EPS_GAP"], p_swap=1.0) for _ in range(n)])
    armB = summarize([trial(rng, macros, cfg, p_opt=0.0, tau=cfg["M_TH"],
                            gap_req=0.0, p_swap=0.0) for _ in range(n)])
    return {
        "base": base, "A_外失效": armA, "B_内失效": armB,
        "σ(A)≈σ(base)": abs(armA["sigma_end"] - base["sigma_end"]) < 0.05,
        "σ(B)<σ(base)": armB["sigma_end"] < base["sigma_end"] - 0.05,
        "A以F2为主": armA["F2%"] > armA["F1%"],
        "B以F1为主": armB["F1%"] > 80,
    }


if __name__ == "__main__":
    c16 = cfg16()
    rep = {
        "可疑点1_提高分辨率后重扫kappa(N_FEAT=16)": kappa_sweep(c16),
        "可疑点1_对照_低分辨率(N_FEAT=8)": kappa_sweep(dict(BASE)),
        "可疑点2_收紧预算(N_FEAT=16)": budget_sweep(c16, budgets=[8, 10, 12, 14, 18, 22]),
        "主结果在高分辨率下是否仍成立": ab_separation(c16),
    }
    print(json.dumps(rep, ensure_ascii=False, indent=2))
    with open("r18_check.json", "w", encoding="utf-8") as f:
        json.dump(rep, f, ensure_ascii=False, indent=2)
