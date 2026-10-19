# -*- coding: utf-8 -*-
"""
R18 高分辨率复核：N_FEAT=16 / N_MACRO=40 下的 base / A / B 三臂对照
=====================================================================
补做本实验的原因（如实记录，勿删）：

  初版实验报告曾写"高分辨率（N_FEAT=16, N_MACRO=40）复现 base σ = A σ = 0.581
  精确相等"。这是**引用错误**——0.581 取自 r18_check.json 的 κ 扫描
  （N_FEAT=16、κ=1.0 那一行的 σ），**不是 base/A 的对照值**；
  r18_sweep.json 的鲁棒性网格虽然含 N_MACRO=40，但 N_FEAT 仍为 8。

  本脚本补做真正的高分辨率 base/A/B 对照，用于判定"σ(A) ≈ σ(base)"这条
  判别签名在细化 σ 量子（1/8 = 0.125 -> 1/16 = 0.0625）之后是否仍成立。

用法：python r18_highres.py  ->  r18_highres.json
"""

import json
import importlib
import statistics as st
from collections import Counter

import r18_dual_constraint as R

# 高分辨率设置
R.N_FEAT = 16          # sigma 量子 0.125 -> 0.0625
R.N_MACRO = 40
R.N_PROTO = 6
R.TOPK = 8
R.B_BUDGET = 20        # 维度翻倍，预算相应放宽，避免被预算截断主导

N_TRIAL = 600
SEED = 20261019


def run(arm, n=N_TRIAL, seed=SEED):
    rng = __import__("random").Random(seed + R.zlib.crc32(arm.encode("utf-8")) % 9973)
    macros = R.build_macros(__import__("random").Random(seed))
    rows = [R.trial(rng, macros, arm) for _ in range(n)]
    c = Counter(r["outcome"] for r in rows)
    sig = [r["sigma_end"] for r in rows]
    return {
        "arm": arm,
        "n": n,
        "真收敛%": round(100 * c["真收敛"] / n, 1),
        "伪收敛F2%": round(100 * c["伪收敛(F2)"] / n, 1),
        "F1%": round(100 * c["F1(未加载/不收敛)"] / n, 1),
        "sigma_end均值": round(st.mean(sig), 3),
        "sigma_end标准差": round(st.pstdev(sig), 3),
        "步数均值": round(st.mean([r["steps"] for r in rows]), 2),
    }


if __name__ == "__main__":
    res = [run(a) for a in ("base", "A", "B")]
    res.append({
        "判别签名": "σ(A) ≈ σ(base) 且 σ(B) < σ(base)",
        "|Δσ(A-base)|": round(abs(res[1]["sigma_end均值"] - res[0]["sigma_end均值"]), 4),
        "σ(base) - σ(B)": round(res[0]["sigma_end均值"] - res[2]["sigma_end均值"], 4),
        "签名成立": bool(
            abs(res[1]["sigma_end均值"] - res[0]["sigma_end均值"]) <= 0.02
            and res[2]["sigma_end均值"] < res[0]["sigma_end均值"] - 0.05
        ),
    })
    print(json.dumps(res, ensure_ascii=False, indent=2))
    with open("r18_highres.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
