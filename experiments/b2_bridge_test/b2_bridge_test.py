#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B2 —— 桥检验（R31：A_D 一致性守卫）：方法 + 脚本准备

对应 18-v3.3-formal-framework-charter.md §6.2 B2（缺口 B-B19 / R43′）。
理论锚点：docs/xs-framework/15-v3.2-gauge-fixed-framework.md §7.2（R43′）与 §7.4 R31。

**背景（R43′ 与 R31 的关系，必须分清）**：
  - **R43′（真桥检验）**：桥 B-B19 写成 `dU/dt = ρ_exp·Ψ + β₀`，其可证伪形式为
    `β₀ = −ρ_exp·Ψ(k*_U)`，其中 **Ψ 必须来自独立于 ΔE 的第三通道**（补篇四独立标定程序）。
    若用 ΔE 代替 Ψ，则**退化为重言**（U 内点峰蕴含 ΔU(k*)=0，由 ΔU=ρ_exp·ΔE+β₀ 即得）。
  - **R31（本脚本实现）**：检验 `β₀ = −ρ_exp·ΔE(k*_U)`（重言形式）。它的作用是
    **A_D 一致性守卫**——不是对桥的检验，而是检验"数据是否满足 A_D + 内点峰自洽"：
    在 A_D 上应≈0 误报，在故意违反 A_D 的注入族上误报率应→1.0（R31 判据原文）。

**性质（诚实声明）**：
  1. 本脚本是**方法 / 检验流程准备**：在合成数据上验证 R31 峰点检验的判别能力
     （正常族≈0、违反族→1.0）。不是对现实认知数据的检验。
  2. **R43′（独立标定 Ψ 的真桥检验）不在本脚本范围内**：它需补篇四的独立标定程序
     （经验分布熵 + Miller–Madow；Q=I(M;C_in|φ)），需要真实认知过程数据，无法在合成
     A_D 数据上忠实模拟。留待 B 组获得真实数据时实施。
  3. R31 是重言形式的守卫，**不能**被误读为"桥成立"的证据——它只暴露 A_D 一致性破坏。

产出：A_D 正常族 vs 注入违反族的多实例统计，峰点检验误报率。
结果写入 b2_results.json。仅标准库；种子固定，机械可复核。
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

SEED = 20261024
DT = 0.01
K = 3000

# 钻研型配置（ρ_res=3 > ρ_exp=1，θ_U=3）：r 穿过 θ_U 时 U 有内点峰
CFG = {"e_plus": 1.0, "e_minus": 1.0, "u_plus": 1.0, "u_minus": 3.0}

# A_D⁺ 边界层参数（违反族用）：ΔU += κ_E·ρ_exp·(1 − ρ·g_strength·U)·dt
KAPPA_E = 0.5
RHO = 1.5
RHO_EXP = 1.0


def r_curve(k, K, pivot, r_hi=9.0, r_lo=0.02):
    """分段投入曲线（同 a13）：前期高、后期低；pivot 控制峰位。"""
    t = k / float(K - 1)
    r_mid = math.sqrt(r_hi * r_lo)
    if t <= pivot:
        return math.exp(math.log(r_hi) + (math.log(r_mid) - math.log(r_hi)) * t / pivot)
    return math.exp(math.log(r_mid) + (math.log(r_lo) - math.log(r_mid)) * (t - pivot) / (1 - pivot))


def gen_trajectory(cfg, pivot, g_strength=0.0, u0=0.5, k=K, dt=DT):
    """生成 (E, U)。g_strength>0 时注入 A_D⁺ 边界层违反 g(U)=g_strength·U（违反纯 A_D）。"""
    Es, Us = [0.0], [u0]
    for i in range(1, k):
        a = r_curve(i, k, pivot)
        dE = dt * (a * cfg["e_plus"] - cfg["e_minus"])
        dU = dt * (a * cfg["u_plus"] - cfg["u_minus"])
        if g_strength != 0.0:
            dU += KAPPA_E * RHO_EXP * (1.0 - RHO * g_strength * Us[-1]) * dt
        Es.append(Es[-1] + dE)
        Us.append(Us[-1] + dU)
    ks = [i * dt for i in range(k)]
    return ks, Es, Us


# ---------- 统计工具（仅标准库） ----------
def multi_linreg(cols, y):
    n = len(y); m = len(cols) + 1
    X = [[1.0] + [c[i] for c in cols] for i in range(n)]
    XtX = [[0.0] * m for _ in range(m)]; Xty = [0.0] * m
    for i in range(n):
        xi = X[i]
        for r in range(m):
            Xty[r] += xi[r] * y[i]
            for c in range(m):
                XtX[r][c] += xi[r] * xi[c]
    aug = [row[:] + [Xty[r]] for r, row in enumerate(XtX)]
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(aug[r][col]))
        aug[col], aug[piv] = aug[piv], aug[col]
        for r in range(m):
            if r != col:
                fac = aug[r][col] / aug[col][col]
                for c in range(col, m + 1):
                    aug[r][c] -= fac * aug[col][c]
    return [aug[r][m] / aug[r][r] for r in range(m)]


def find_inner_peak(Us):
    """内点峰：内部点 U[i] ≥ 两侧且取全局最大。返回索引；无内点峰返回 None。"""
    n = len(Us)
    best, bi = None, None
    for i in range(1, n - 1):
        if Us[i] >= Us[i - 1] and Us[i] >= Us[i + 1]:
            if best is None or Us[i] > best:
                best, bi = Us[i], i
    return bi


def peak_test(Es, Us, rho_exp, beta0, kstar):
    """峰点检验：A_D 预测的每步 ΔU = ρ_exp·ΔE_c + β₀·dt，峰点处应 ≈0（ΔU_c≈0）。

    返回相对偏离 S/(|ρ_exp·ΔE_c|+|β₀·dt|)。A_D 正常族应≈0，违反族应明显 >0。
    """
    dE_c = (Es[kstar + 1] - Es[kstar - 1]) / 2.0   # 每 dt 的 dE（中心差分）
    S = abs(rho_exp * dE_c + beta0 * DT)           # A_D 预测的每步 dU
    scale = abs(rho_exp * dE_c) + abs(beta0 * DT) + 1e-12
    return S / scale


# ---------- 主流程 ----------
def main():
    pivots = [0.4, 0.5, 0.6, 0.7]
    g_levels = [0.0, 0.4, 0.8, 1.2]   # A_D⁺ 边界层强度；0=纯 A_D（正常族）

    report = {"meta": {
        "task": "B2/R31 桥检验（A_D 一致性守卫）——方法准备，非经验验证",
        "model": "A_D（15-v3.2 §7.2/§7.4 R31）",
        "dt": DT, "K": K,
        "note": "R31 为重言形式（用 ΔE），是 A_D 一致性守卫，非 R43′ 独立 Ψ 真桥检验；"
                "违反族 = A_D⁺ 边界层注入（g(U)=g_strength·U）",
        "engine": "仅标准库；确定性（pivot 网格）",
    }, "families": {}}

    for g_strength in g_levels:
        trials = []
        fails = 0
        for pivot in pivots:
            ks, Es, Us = gen_trajectory(CFG, pivot, g_strength=g_strength)
            b = multi_linreg([Es, ks], Us)
            rho_exp, beta0 = b[1], b[2]
            kstar = find_inner_peak(Us)
            if kstar is None:
                trials.append({"pivot": pivot, "status": "no_inner_peak", "violation": None})
                continue
            rel_viol = peak_test(Es, Us, rho_exp, beta0, kstar)
            is_fail = rel_viol > 0.02
            if is_fail:
                fails += 1
            trials.append({"pivot": pivot, "kstar": kstar, "rel_violation": rel_viol, "fail": is_fail})
        rate = fails / len(pivots)
        report["families"][("normal" if g_strength == 0.0 else f"adplus_g{g_strength}")] = {
            "g_strength": g_strength, "trials": trials, "fail_rate": rate}
        print(f"g={g_strength:5.2f}  fail_rate={rate:.3f}  (误报{'≈0' if g_strength==0 else '→1.0 目标'})")
        for t in trials:
            if "rel_violation" in t:
                print(f"    pivot={t['pivot']:.1f} k*={t['kstar']} rel_viol={t['rel_violation']:.4f} fail={t['fail']!s}")

    out_path = os.path.join(HERE, "b2_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print(f"\nwrote: {out_path}")

    print("\n== 判据（R31）==")
    normal = report["families"]["normal"]["fail_rate"]
    violate = max(report["families"][k]["fail_rate"] for k in report["families"] if k.startswith("adplus"))
    ok_normal = normal < 0.05
    ok_violate = violate >= 0.9
    print(f"A_D 正常族误报率 ≈0 (<0.05): {normal:.3f} -> {'PASS' if ok_normal else 'FAIL'}")
    print(f"违反族误报率 →1.0 (≥0.9): {violate:.3f} -> {'PASS' if ok_violate else 'FAIL'}")


if __name__ == "__main__":
    main()
