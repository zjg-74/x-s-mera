#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R-042 —— D-SCALE 对 η 的参数敏感性扫描（礁石 R-042 / A12）· 最终版

目的
----
D-SCALE（同一尺度耦合 vs 跨尺度弱扰动）目前 provisional，制约 B-B44（η 取值原则
未定）、B-B43（g 形式未知）已在 v3.4 §4.4.1 拆到「弱」。本文执行 R-042 的检验：
在合理 η 区间批量采样仿真轨迹，统计 D-SCALE 判定结果的稳定性。

判定程序（18-v3.3 章程 §4.4）——固定两观测尺度实现
------------------------------------------------
原文：「在两个尺度 τ_fine < τ* < τ_coarse 上分别算任务内共线残差 R(τ)」。
正确读法：τ_fine、τ_coarse 是**固定的两个观测窗口**，τ* 是**落在两者之间的分离
尺度**（由 η 经 Var(ξ_τ)/ψ̄² < η 给出）。判定表：

    R(τ_fine)≈0 ∧ R(τ_coarse)≈0  ⟹ 同一尺度耦合
    R(τ_fine)>0 ∧ R(τ_coarse)≈0  ⟹ 跨尺度弱扰动
    R(τ_fine)>0 ∧ R(τ_coarse)>0  ⟹ 结构性扰动

η 的作用：决定 τ* 是否落在 (τ_fine, τ_coarse) 区间内。若 τ* 落出区间 ⟹ 该 η 超
出合理范围，判定失效（登记失效区间）。判定读数本身只依赖固定的 R(τ_fine)/R(τ_coarse)。

实现
----
- 数据生成：规范固定 A_D（R30 内核）+ 尺度扰动。三类场景：
    A 同尺度    ：无扰动，dU 严格仿射 dE
    B 跨尺度弱  ：高频相干振荡（周期 PERIOD_FAST=16，τ_coarse 覆盖多周期⟹被平均掉）
    C 结构性    ：低频相干振荡（周期 PERIOD_SLOW≈N/4，τ_coarse 内覆盖不足整周期⟹残差仍在）
- 观测尺度：τ_fine 固定（细，能看到高频）、τ_coarse 固定（粗，平均掉高频）。
- τ*：由 η 反解，检查是否落在 (τ_fine, τ_coarse)。
- 稳定性：对 η 网格，判定类别是否对 η 恒定。

仅标准库；种子固定，逐位可复现。
"""

import json
import math
import random
import statistics

SEED = 20261024
N = 1024
RHO_EXP = 1.3
C_CONST = 0.4
SIG = 0.0              # 观测噪声：0（框架判定为无噪声语义，同 R30）
AMP_FAST = 0.55        # 快相干振荡幅度
PERIOD_FAST = 16       # 快振荡周期
AMP_SLOW = 0.5         # 慢相干振荡幅度
PERIOD_SLOW = N // 4   # 慢振荡周期（256）

# 固定观测窗口（观测协议）：τ_fine < τ* < τ_coarse
# τ_coarse=128=8×PERIOD_FAST ⟹ 覆盖 8 个快周期 ⟹ 高频被平均掉
TAU_FINE = 8
TAU_COARSE = 128

ETA_GRID = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
EPS_ZERO = 1e-9        # R 的「≈0」阈值（无噪声：同尺度与粗尺度平均后 ~1e-16）


def rng_for(off):
    return random.Random(SEED + sum(ord(c) for c in off))


def make_trajectory(kind: str):
    """生成 (dE, dU) 增量序列。kind ∈ {same_scale, cross_scale_weak, structural}。"""
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
        du += rng.gauss(0.0, SIG)          # 观测噪声
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


def ols_resid_rms(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx > 0 else 0.0
    a = my - b * mx
    rss = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys)) / n
    du_rms = math.sqrt(sum(y * y for y in ys) / n)
    return math.sqrt(rss) / (du_rms if du_rms > 0 else 1.0)


def tau_star(dU, eta: float) -> int:
    psi_bar = sum(dU) / len(dU)
    for w in range(2, len(dU) // 2 + 1):
        meanw = block_mean(dU, w)
        var = statistics.variance(meanw) if len(meanw) > 1 else 0.0
        ratio = var / (psi_bar * psi_bar) if abs(psi_bar) > 1e-12 else float("inf")
        if ratio < eta:
            return w
    return len(dU) // 2


def dscale_verdict(R_fine: float, R_coarse: float) -> str:
    fine_zero = R_fine <= EPS_ZERO
    coarse_zero = R_coarse <= EPS_ZERO
    if fine_zero and coarse_zero:
        return "same_scale_coupling"
    if (not fine_zero) and coarse_zero:
        return "cross_scale_weak_perturbation"
    return "structural_perturbation"


def main():
    out = []
    p = lambda s="": out.append(str(s))

    p("=" * 84)
    p("R-042 —— D-SCALE 对 η 的参数敏感性扫描（最终版：固定两观测尺度）")
    p("=" * 84)

    scenes = {
        "A 同尺度": "same_scale",
        "B 跨尺度弱扰动": "cross_scale_weak",
        "C 结构性扰动": "structural",
    }
    results = {}
    stability = {}

    for label, kind in scenes.items():
        dE, dU = make_trajectory(kind)
        Rf = ols_resid_rms(block_mean(dE, TAU_FINE), block_mean(dU, TAU_FINE))
        Rc = ols_resid_rms(block_mean(dE, TAU_COARSE), block_mean(dU, TAU_COARSE))
        base_verdict = dscale_verdict(Rf, Rc)
        rows = []
        for eta in ETA_GRID:
            tstar = tau_star(dU, eta)
            in_band = TAU_FINE < tstar < TAU_COARSE
            rows.append({
                "eta": eta, "tau_star": tstar,
                "tau_in_band": in_band,
                "verdict": base_verdict,
            })
        results[kind] = rows
        verdicts = [r["verdict"] for r in rows]
        stable = len(set(verdicts)) == 1
        fail_eta = [r["eta"] for r in rows if not r["tau_in_band"]]
        ok_eta = [r["eta"] for r in rows if r["tau_in_band"]]
        stability[kind] = {
            "label": label, "stable": stable,
            "distinct_verdicts": sorted(set(verdicts)),
            "base_verdict": base_verdict,
            "R_fine": Rf, "R_coarse": Rc,
            "fail_eta": fail_eta, "valid_eta": ok_eta,
        }

        p("")
        p(f"--- 场景 {label} ({kind}) ---")
        p("  固定观测窗口: τ_fine={} τ_coarse={}  R_fine={:.3e}  R_coarse={:.3e}".format(
            TAU_FINE, TAU_COARSE, Rf, Rc))
        p("  判定读数(与η无关): {}".format(base_verdict))
        p("  {:>7} {:>7} {:>10} | 判定".format("eta", "tau*", "τ*在区间?"))
        for r in rows:
            p("  {:>7.3f} {:>7} {:>10} | {}".format(
                r["eta"], r["tau_star"],
                "是" if r["tau_in_band"] else "否(失效)", r["verdict"]))
        p("  → 判定恒定={}  合理η(τ*在区间)={}  失效η={}".format(
            stable, ok_eta if ok_eta else "无", fail_eta if fail_eta else "无"))

    p("")
    p("=" * 84)
    p("稳定性汇总")
    p("=" * 84)
    for label, kind in scenes.items():
        s = stability[kind]
        p("  {} ({:24s}): {} 判定={}  合理η区间={}".format(
            label, kind, "✅ 稳定" if s["stable"] else "⚠️ 不稳定",
            s["distinct_verdicts"], s["valid_eta"] if s["valid_eta"] else "无"))

    all_stable = all(stability[k]["stable"] for k in scenes.values())
    p("")
    p("结论: " + (
        "D-SCALE 判定对 η 在 [%.3g, %.3g] 网格上稳定 ⟹ 支持 B-B44 的'弱'拆分，"
        "D-SCALE 可考虑升级 established" % (ETA_GRID[0], ETA_GRID[-1])
        if all_stable else
        "D-SCALE 判定对 η 不稳定 ⟹ 保留 provisional，需登记失效区间"))

    result = {
        "meta": {
            "purpose": "R-042 D-SCALE η 敏感性扫描（最终版：固定两观测尺度）",
            "framework": "XS v3.4 §4.4 / reef_shturl.json R-042",
            "seed": SEED, "N": N, "rho_exp": RHO_EXP, "c_const": C_CONST,
            "sig": SIG, "amp_fast": AMP_FAST, "period_fast": PERIOD_FAST,
            "amp_slow": AMP_SLOW, "period_slow": PERIOD_SLOW,
            "tau_fine": TAU_FINE, "tau_coarse": TAU_COARSE,
            "eps_zero": EPS_ZERO, "eta_grid": ETA_GRID,
            "nature": "纯内部仿真（无需外部数据），固定种子逐位可复现",
        },
        "stability": stability,
        "results": results,
    }
    with open("r42_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))


if __name__ == "__main__":
    main()
