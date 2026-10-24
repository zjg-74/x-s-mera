#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R-043 —— τ* 落点：η 的合理取值区间（礁石 R-043 / B-B44）

背景
----
R-042（experiments/r42_dscale_eta）固定观测窗口 τ_fine=8、τ_coarse=128，在 η 网格
[0.005,0.01,0.02,0.05,0.1,...,1.0] 上检查 τ*（由 Var(ξ_τ)/ψ̄²<η 反解）是否落在
(τ_fine, τ_coarse) 内，发现 A/C 场景全空、B 仅 η=0.05 一格。⟹ 判定读数对 η 不敏感
**不等于** η 的合理取值区间已定（B-B44 未解除）。

R-043 问题：是否存在**非空**的 η 区间使 τ* ∈ (τ_fine, τ_coarse)？若有，区间边界如何随场景变化？
纯内部可解：加密 η/τ 网格 + 重选观测窗口，不需外部数据。

方法
----
1. 加密 η 网格（0.001–1.0，200 点线性）——R-042 的 η 网格过粗（0.02→0.05 跳变），
   可能错过 τ* 落带区间。
2. 对三类场景（A 同尺度 / B 跨尺度弱 / C 结构性，复用 R-042 内核），
   预计算 ratio(w)=Var(块均值_dU,w)/ψ̄²，再对每 η 反解 τ*(η)=min{w:ratio(w)<η}。
3. 固定窗口 (8,128)：统计 τ* 落带的 η 区间（连续段）。
4. 若某场景在固定窗口下无落带，扫描观测窗口 (τ_fine, τ_coarse)
   （τ_fine∈{4,8,16}，τ_coarse∈{64,128,256,512}），找使落带 η 区间非空且较宽的窗口。

仅标准库；种子固定，逐位可复现。
"""

import json
import math
import random
import statistics
from pathlib import Path

SEED = 20261024
N = 1024
RHO_EXP = 1.3
C_CONST = 0.4
SIG = 0.0
AMP_FAST = 0.55
PERIOD_FAST = 16
AMP_SLOW = 0.5
PERIOD_SLOW = N // 4

TAU_FINE_DEF = 8
TAU_COARSE_DEF = 128

ETA_MIN, ETA_MAX, ETA_PTS = 0.001, 1.0, 200


def rng_for(off):
    return random.Random(SEED + sum(ord(c) for c in off))


def make_trajectory(kind: str):
    """复用 R-042 内核：生成 (dE, dU) 增量序列。"""
    rng = rng_for(kind)
    dE, dU = [], []
    for k in range(N):
        base = 1.0 + 0.4 * math.sin(2 * math.pi * k / (N / 3.0))
        du = RHO_EXP * base + C_CONST
        if kind == "same_scale":
            de = base
        elif kind == "cross_scale_weak":
            de = base
            du += AMP_FAST * math.sin(2 * math.pi * k / PERIOD_FAST)
        elif kind == "structural":
            de = base
            du += AMP_SLOW * math.sin(2 * math.pi * k / PERIOD_SLOW)
        else:
            raise ValueError(kind)
        du += rng.gauss(0.0, SIG)
        dE.append(de)
        dU.append(du)
    return dE, dU


def block_mean(seq, w: int):
    n = len(seq)
    out = []
    i = 0
    while i < n:
        j = min(i + w, n)
        out.append(sum(seq[i:j]) / (j - i))
        i = j
    return out


def ratio_curve(dU):
    """ratio(w)=Var(块均值_dU,w)/ψ̄²，w∈[2, N//2]。ψ̄ 对 A_D 轨迹非零。"""
    psi_bar = sum(dU) / len(dU)
    denom = psi_bar * psi_bar
    out = []
    for w in range(2, len(dU) // 2 + 1):
        meanw = block_mean(dU, w)
        var = statistics.variance(meanw) if len(meanw) > 1 else 0.0
        ratio = var / denom if abs(psi_bar) > 1e-12 else float("inf")
        out.append(ratio)
    return out


def tau_star_from_ratio(ratios, eta):
    """τ*(η)=min{w: ratio(w)<η}；w 从 2 起。"""
    for idx, r in enumerate(ratios, start=2):
        if r < eta:
            return idx
    return len(ratios) + 1


def in_band(tstar, tf, tc):
    return tf < tstar < tc


def valid_eta_intervals(etas, flags):
    """把布尔落带序列压成连续 η 区间 [lo,hi]。"""
    segs = []
    run = []
    for e, f in zip(etas, flags):
        if f:
            run.append(e)
        else:
            if run:
                segs.append((run[0], run[-1]))
            run = []
    if run:
        segs.append((run[0], run[-1]))
    return segs


def main():
    out = []
    p = lambda s="": out.append(str(s))
    p("=" * 84)
    p("R-043 —— τ* 落点：加密 η 网格 + 观测窗口扫描，找 τ* 落带的非空 η 区间")
    p("=" * 84)

    scenes = {
        "A 同尺度": "same_scale",
        "B 跨尺度弱扰动": "cross_scale_weak",
        "C 结构性扰动": "structural",
    }
    etas = [ETA_MIN + (ETA_MAX - ETA_MIN) * i / (ETA_PTS - 1) for i in range(ETA_PTS)]

    summary = {}
    window_scan = {}
    for label, kind in scenes.items():
        dE, dU = make_trajectory(kind)
        ratios = ratio_curve(dU)
        tstars = [tau_star_from_ratio(ratios, e) for e in etas]
        flags_def = [in_band(t, TAU_FINE_DEF, TAU_COARSE_DEF) for t in tstars]
        segs_def = valid_eta_intervals(etas, flags_def)
        span_def = (segs_def[0][0], segs_def[0][-1]) if segs_def else None
        summary[kind] = {
            "label": label,
            "eta_intervals_in_default_window": segs_def,
            "default_window": (TAU_FINE_DEF, TAU_COARSE_DEF),
            "tau_star_range": (min(tstars), max(tstars)),
        }
        p("")
        p(f"--- 场景 {label} ({kind}) · 固定窗口 ({TAU_FINE_DEF},{TAU_COARSE_DEF}) ---")
        p(f"  τ* 范围: [{min(tstars)}, {max(tstars)}]  ψ̄={sum(dU)/len(dU):.4f}")
        if span_def:
            p(f"  ✅ τ* 落带的 η 连续区间: [{span_def[0]:.4f}, {span_def[1]:.4f}]  (宽 {(span_def[1]-span_def[0]):.4f})")
        else:
            p("  ⚠️ 固定窗口下 τ* 无落带 η 区间")

        # 窗口扫描：找使落带区间非空且最宽的 (τ_fine, τ_coarse)
        best = None
        for tf in (4, 8, 16, 32):
            for tc in (64, 128, 256, 512):
                if tc <= tf:
                    continue
                flags = [in_band(t, tf, tc) for t in tstars]
                segs = valid_eta_intervals(etas, flags)
                if segs:
                    lo, hi = segs[0][0], segs[0][-1]
                    width = hi - lo
                    if best is None or width > best[0]:
                        best = (width, tf, tc, lo, hi)
        window_scan[kind] = best
        if best:
            w, tf, tc, lo, hi = best
            p(f"  🎯 最佳观测窗口 ({tf},{tc}): τ* 落带 η 区间 [{lo:.4f}, {hi:.4f}] 宽 {w:.4f}")
        else:
            p("  ❌ 所有扫描窗口下均无 τ* 落带 η 区间")

    p("")
    p("=" * 84)
    p("R-043 判定")
    p("=" * 84)
    ok = 0
    for kind, s in summary.items():
        if s["eta_intervals_in_default_window"]:
            ok += 1
    # 判定：默认窗口下三类场景是否都有非空落带区间
    if ok == len(scenes):
        conclusion = ("三类场景在默认窗口下均有非空 τ* 落带 η 区间 ⟹ B-B44 可取经验值范围；"
                      "登记各场景 η 合理区间，B-B44 降级为「已给出区间、待锚定取值」")
    else:
        conclusion = ("默认窗口下仅部分场景有落带 η 区间 ⟹ 需重选观测窗口；"
                      "窗口扫描给出每场景的最佳窗口与 η 区间（见 window_scan），"
                      "B-B44 保持开放但范围收窄")
    p(conclusion)

    result = {
        "meta": {
            "purpose": "R-043 τ* 落点：η 合理取值区间（加密网格 + 窗口扫描）",
            "framework": "XS v3.4 / reef_shturl.json R-043 / B-B44",
            "seed": SEED, "N": N, "eta_grid_points": ETA_PTS,
            "eta_min": ETA_MIN, "eta_max": ETA_MAX,
            "default_window": [TAU_FINE_DEF, TAU_COARSE_DEF],
            "nature": "纯内部仿真（无需外部数据），固定种子逐位可复现",
        },
        "summary": summary,
        "window_scan": {k: (v if v is None else list(v)) for k, v in window_scan.items()},
        "conclusion": conclusion,
    }
    with open(str(Path(__file__).parent / "r43_results.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))


if __name__ == "__main__":
    main()
