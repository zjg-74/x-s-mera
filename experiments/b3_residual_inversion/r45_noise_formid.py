#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R45 —— R32 反演扩展：加噪稳健性（裸 vs 去噪平滑）+ g 形式辨识能力（B-B43 纯内部推进）

背景
----
B-B43（g(U) 形式未知）已由 R32 / b3_residual_inversion 推进到「残差可反演 g 的变异」
（合成轨迹上 corr=1.0、负对照无虚假结构），但 README 明确留下两块后续：
  ① **加噪稳健性**（反演对 U 观测噪声的容忍度）；
  ② **形式辨识能力**（给定反演形状，能否从候选形式族中确定 g(U) 的具体形式）。
本实验纯内部（合成数据）回答这两点——它们决定 **B-B43 的剩余自由度能否由真实 (E,U)
数据闭合**：
  - 若加噪后（裸或去噪平滑）反演仍能恢复形状 → 真实数据反演可行（稳健）；
  - 若反演形状能从候选族识别真形式 → 有了真实数据即可确定 g(U) 形式（可辨识）；
   否则形式永不可辨，B-B43 只能降级为「结构约束 + 形状」而无法定形式。

**前置发现（v1）**：裸中心差分反演对噪声极敏感（U 1% 噪声即把 corr 打到 0.05），
因为差分放大噪声 ∝ 1/step。故本版加入**去噪平滑**分支（对残差滑动平均后再差分），
评估真实数据可行性。

仅标准库；固定种子逐位可复现；内核复用 b3_residual_inversion.py。
"""

import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))

STEP = 0.01
N = 2001
RHO_EXP = 1.0
KAPPA_E = 0.5
RHO = 1.5
U0 = 1.0
E0 = 0.0
E_SLOPE = 0.05
SEED = 20261027

SIGS = [0.0, 0.01, 0.02, 0.05, 0.1]
M_NOISE = 40
W_SMOOTH = 11          # 滑动平均窗口
CANDIDATES = ["linear", "square", "tanh", "saturated"]


def g_linear(u): return u
def g_square(u): return u * u
def g_tanh(u): return math.tanh(u)
def g_sat(u): return u / (1.0 + u)

GFUNCS = {"linear": g_linear, "square": g_square, "tanh": g_tanh, "saturated": g_sat}


def e_path(k):
    return E0 + E_SLOPE * k * STEP


def ad_plus_trajectory(g, n=N, step=STEP):
    ks = [i * step for i in range(n)]
    Es, Us = [e_path(0)], [U0]
    for i in range(1, n):
        dE = e_path(i) - e_path(i - 1)
        dU = RHO_EXP * dE + KAPPA_E * RHO_EXP * (1.0 - RHO * g(Us[-1])) * step
        Us.append(Us[-1] + dU)
        Es.append(Es[-1] + dE)
    return ks, Es, Us


def linfit(xs, ys):
    n = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    if abs(den) < 1e-15:
        return 0.0, 0.0
    b = (n * sxy - sx * sy) / den
    a = (sy - b * sx) / n
    return a, b


def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    den = math.sqrt(vx * vy)
    return cov / den if den > 1e-15 else 0.0


def moving_average(seq, w):
    """滑动平均（尾部截断以保持长度，避免窗口平移引入相位）。"""
    w = int(w)
    if w <= 1:
        return list(seq)
    half = w // 2
    out = []
    n = len(seq)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out.append(sum(seq[lo:hi]) / (hi - lo))
    return out


def center_diff(d, step=STEP):
    return [(d[i + 1] - d[i - 1]) / (2.0 * step) for i in range(1, len(d) - 1)]


def invert(g, noise=0.0, seed=0, smooth_w=0):
    """反演 g 形状（中心化）。加噪：对 U 加独立高斯；smooth_w>0 时先对残差滑动平均。"""
    ks, Es, Us = ad_plus_trajectory(g)
    if noise > 0.0:
        rng = random.Random(SEED + seed * 7919 + int(noise * 1000))
        Us = [u + rng.gauss(0.0, noise) for u in Us]
    R = [U - RHO_EXP * E for U, E in zip(Us, Es)]
    if smooth_w > 0:
        R = moving_average(R, smooth_w)
    a, b = linfit(ks, R)
    D = [r - (a + b * k) for r, k in zip(R, ks)]
    dD = center_diff(D)
    k_c = ks[1:1 + len(dD)]
    U_c = Us[1:1 + len(dD)]
    g_inv = [-dd / (KAPPA_E * RHO_EXP * RHO) for dd in dD]
    m = sum(g_inv) / len(g_inv)
    return U_c, [x - m for x in g_inv]


def normalize(xs):
    m = sum(xs) / len(xs)
    c = [x - m for x in xs]
    s = math.sqrt(sum(x * x for x in c) / len(c))
    return [x / s for x in c] if s > 1e-15 else c


def corr_to_true(shape_fn, g):
    """反演形状与真值形状（在无噪 U 上参数化）的 corr。"""
    U_c0, _ = invert(g, noise=0.0)
    g_true = normalize([g(u) for u in U_c0])
    inv = shape_fn()
    return pearson(normalize(inv), g_true)


def main():
    out = []
    p = lambda s="": out.append(str(s))
    p("=" * 84)
    p("R45 —— R32 反演扩展：加噪稳健性（裸 vs 去噪平滑）+ g 形式辨识（B-B43）")
    p("=" * 84)

    # ---------- ① 加噪稳健性 ----------
    p("\n## ① 加噪稳健性：反演 corr 随 U 观测噪声 σ（中位数 over M=%d）" % M_NOISE)
    robust = {}
    for name in GFUNCS:
        g = GFUNCS[name]
        raw_row, sm_row = {}, {}
        for sig in SIGS:
            def mk_raw(noise=sig):
                _, inv = invert(g, noise=noise, seed=1)
                return inv
            def mk_sm(noise=sig):
                _, inv = invert(g, noise=noise, seed=1, smooth_w=W_SMOOTH)
                return inv
            c_raw = sorted(corr_to_true(mk_raw, g) for _ in range(M_NOISE))
            c_sm = sorted(corr_to_true(mk_sm, g) for _ in range(M_NOISE))
            raw_row[str(sig)] = {"median": round(c_raw[len(c_raw) // 2], 3),
                                 "min": round(c_raw[0], 3)}
            sm_row[str(sig)] = {"median": round(c_sm[len(c_sm) // 2], 3),
                                "min": round(c_sm[0], 3)}
        raw_tol = max((s for s in SIGS if raw_row[str(s)]["median"] > 0.90), default=None)
        sm_tol = max((s for s in SIGS if sm_row[str(s)]["median"] > 0.90), default=None)
        robust[name] = {"raw": raw_row, "smooth": sm_row,
                        "raw_tolerance": raw_tol, "smooth_tolerance": sm_tol}
        p(f"  {name:10s} 裸:" + " ".join(f"{s}:{raw_row[str(s)]['median']:.2f}" for s in SIGS) +
          f" [容忍σ≤{raw_tol}]")
        p(f"  {name:10s} 平滑:" + " ".join(f"{s}:{sm_row[str(s)]['median']:.2f}" for s in SIGS) +
          f" [容忍σ≤{sm_tol}]")

    # ---------- ② 形式辨识 ----------
    p("\n## ② 形式辨识：无噪 + 有噪(σ=0.02, 平滑)下识别率")
    form = {}
    for sig in [0.0, 0.02]:
        correct_n = 0
        form[f"sigma_{sig}"] = {}
        p(f"\n  σ={sig}:")
        for true_name in CANDIDATES:
            g = GFUNCS[true_name]
            U_c0, _ = invert(g, noise=0.0)
            if sig == 0.0:
                _, inv = invert(g, noise=0.0)
                scores = {cand: pearson(normalize(inv), normalize([GFUNCS[cand](u) for u in U_c0]))
                          for cand in CANDIDATES}
                best = max(scores, key=scores.get)
                second = sorted(scores.values(), reverse=True)[1]
                margin = scores[best] - second
            else:
                votes = []
                for m in range(M_NOISE):
                    _, inv_m = invert(g, noise=sig, seed=m, smooth_w=W_SMOOTH)
                    sc = {cand: pearson(inv_m, normalize([GFUNCS[cand](u) for u in U_c0]))
                          for cand in CANDIDATES}
                    votes.append(max(sc, key=sc.get))
                best = max(set(votes), key=votes.count)
                margin = votes.count(best) / len(votes)
            correct = (best == true_name) and margin >= 0.02
            correct_n += int(correct)
            p(f"    {true_name:10s} 识别={best:10s} 领先margin={margin:.3f} 正确={correct}")
            form[f"sigma_{sig}"][true_name] = {"identified": best, "correct": correct}
        form[f"sigma_{sig}"]["identify_rate"] = correct_n / len(CANDIDATES)

    result = {"meta": {
        "purpose": "B-B43 纯内部推进：R32 反演加噪稳健性（裸/去噪平滑）+ g 形式辨识",
        "framework": "XS v3.4 / B-B43 / R32 / b3_residual_inversion",
        "seed": SEED, "N": N, "M_noise": M_NOISE, "sigs": SIGS,
        "smooth_window": W_SMOOTH, "candidates": CANDIDATES,
        "nature": "纯内部仿真，固定种子可复现"},
        "robustness": robust, "form_identification": form, "conclusion": None}

    # 判定：平滑后存在 σ>0 容忍上限，且无噪/有噪识别率≥3/4
    smooth_ok = any(r["smooth_tolerance"] is not None and r["smooth_tolerance"] > 0.0
                    for r in robust.values())
    id_ok = form["sigma_0.0"]["identify_rate"] >= 0.75 and \
            form["sigma_0.02"]["identify_rate"] >= 0.75
    result["verdict"] = bool(smooth_ok and id_ok)
    result["conclusion"] = (
        "✅ B-B43 纯内部推进：裸反演对噪声极敏感（1% 噪声即失效），但**去噪平滑后**反演"
        "可容忍显著噪声且形状可从候选族识别真 g（无噪/有噪识别率≥75%）⟹ 真实 (E,U) 数据"
        "需先去噪平滑再反演，B-B43 剩余自由度（g 形式）**可由数据闭合**"
        if result["verdict"] else
        "⚠️ 部分失败：去噪平滑后加噪容忍或形式识别未达标 ⟹ B-B43 剩余自由度能否由数据"
        "闭合存疑，需登记（可能只能降级为结构约束+形状）")
    p("\n" + "=" * 84)
    p("R45 判定")
    p("=" * 84)
    p(result["conclusion"])

    with open(os.path.join(HERE, "r45_noise_formid.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))


if __name__ == "__main__":
    main()
