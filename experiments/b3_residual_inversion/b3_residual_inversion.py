#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B3 —— 残差反演 g（R32）：把「残差是噪声」变成可检验假设（方法 + 脚本准备）

对应 18-v3.3-formal-framework-charter.md §6.2 B3（缺口 B-B43 / R32）：
    「残差反演 g：把『残差是噪声』变成可检验假设。」

理论锚点：docs/xs-framework/15-v3.2-gauge-fixed-framework.md §5.7（推论 R43）。
在 A_D⁺（带边界层 g(U)，g(0)=0、g 单调增）下，U 方程为：
    dU/dk = ρ_exp·dE/dk + κ_E·ρ_exp·(1 − ρ·g(U))
⟹ U(k) − ρ_exp·E(k) = α + βk − κ_E·ρ_exp·ρ·∫₀ᵏ (g(U(t)) − ḡ) dt
⟹ 推论 R43：共线残差 ≡ 0 ⟺ g(U(k)) 沿路径为常数；残差的形状即 g 沿路径变化的反演。
即：残差不是噪声，是信号——这是把「残差是噪声」变成可检验假设的方法基础。

**性质（诚实声明）**：
  1. 本脚本是**方法验证**：在**已知 g 的合成轨迹**上检验反演器是否恢复真值（R32 原文判据
     「在已知 g 的合成轨迹上，残差反演值与真值一致」），不是对现实数据的经验验证。
  2. 反演恢复 g 的**变异 / 形状**（相对路径均值），**绝对水平不可辨识**——常数平移被
     ∫(g−ḡ) 吸收，这是规范自由度（B-B23）在本方法上的体现，不是缺陷，是边界。
  3. 「真实认知残差确实是 g 信号」**不由此声称**——那属于 B 组 B4（跨任务/跨人群验证）。

产出：对若干单调增、g(0)=0 的 g 函数做 A_D⁺ 合成 → 残差反演 → 与真值对比；
      负对照（g 恒为 0 或常数）确认反演不产生虚假结构。
结果写入 b3_results.json。仅标准库；无随机（确定性合成），机械可复核。
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- 确定性合成参数 ----------
STEP = 0.01          # 离散步长（时间/序号步）
N = 2001             # 采样点数（含端点）
RHO_EXP = 1.0        # 映射系数 ρ_exp
KAPPA_E = 0.5        # 尺度系数 κ_E
RHO = 1.5            # 钻研相（ρ>1），g 边界层更显著
U0 = 1.0             # U 初值（保持 U 远离 0，B-B43 辖域）
E0 = 0.0             # E 初值
E_SLOPE = 0.05       # E 线性驱动斜率


# ---------- g 函数族（单调增、g(0)=0） ----------
def g_linear(u):
    return u

def g_square(u):
    return u * u

def g_tanh(u):
    return math.tanh(u)

def g_sat(u):
    # 有界饱和型：U/(1+U)，g(0)=0、单调增、趋近 1
    return u / (1.0 + u)

GFUNCS = {
    "linear":  g_linear,
    "square":  g_square,
    "tanh":    g_tanh,
    "saturated": g_sat,
}

# ---------- 负对照 ----------
def g_zero(u):
    return 0.0

def g_const(u):
    return 0.5  # 常数：变异为零


# ---------- A_D⁺ 合成轨迹 ----------
def e_path(k):
    return E0 + E_SLOPE * k * STEP


def ad_plus_trajectory(g, n=N, step=STEP):
    """欧拉推进 A_D⁺：dU = ρ_exp·dE + κ_E·ρ_exp·(1−ρ·g(U))·dk。U 不夹紧（保持远离 0）。"""
    ks = [i * step for i in range(n)]
    Es, Us = [e_path(0)], [U0]
    for i in range(1, n):
        dE = e_path(i) - e_path(i - 1)
        dU = RHO_EXP * dE + KAPPA_E * RHO_EXP * (1.0 - RHO * g(Us[-1])) * step
        Us.append(Us[-1] + dU)
        Es.append(Es[-1] + dE)
    return ks, Es, Us


# ---------- 统计工具（仅标准库） ----------
def linfit(xs, ys):
    """最小二乘直线 a + b·x。"""
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


def center_diff(d, step=STEP):
    """内部点中心差分 d'(i) = (d[i+1]-d[i-1])/(2·step)，返回 (导数列表, 有效起点索引)。"""
    out = []
    for i in range(1, len(d) - 1):
        out.append((d[i + 1] - d[i - 1]) / (2.0 * step))
    return out, 1


# ---------- 反演主体 ----------
def invert_g(g, n=N, step=STEP):
    """对给定 g 做完整反演，返回诊断。"""
    ks, Es, Us = ad_plus_trajectory(g, n=n, step=step)
    R = [U - RHO_EXP * E for U, E in zip(Us, Es)]

    # 1) 对残差去仿射趋势（拟合 a+bk，吸收 A_D 基线 + g 积分的线性部分）
    a, b = linfit(ks, R)
    D = [r - (a + b * k) for r, k in zip(R, ks)]

    # 2) 中心差分 → 反演 g（残差等式：dD/dk = −κ_E·ρ_exp·ρ·(g(U)−ḡ)）
    dD, start = center_diff(D, step)
    k_c = ks[start:start + len(dD)]
    U_c = Us[start:start + len(dD)]
    g_inv = [-dd / (KAPPA_E * RHO_EXP * RHO) for dd in dD]

    # 3) 真值（沿路径，参数化）
    g_true = [g(u) for u in U_c]
    m_true = sum(g_true) / len(g_true)
    m_inv = sum(g_inv) / len(g_inv)
    g_true_c = [x - m_true for x in g_true]
    g_inv_c = [x - m_inv for x in g_inv]

    corr = pearson(g_inv_c, g_true_c)
    std_true = math.sqrt(sum(x * x for x in g_true_c) / len(g_true_c))
    std_inv = math.sqrt(sum(x * x for x in g_inv_c) / len(g_inv_c))
    scale = (std_inv / std_true) if std_true > 1e-15 else 0.0
    # 归一化 RMSE：反演缩放（除以 scale）后与真值变异比的 RMSE（再除以真值 std）
    inv_scaled = [x / scale for x in g_inv_c] if scale > 1e-12 else [0.0] * len(g_inv_c)
    rmse = math.sqrt(sum((x - y) ** 2 for x, y in zip(inv_scaled, g_true_c)) / len(g_true_c))

    # 结构显著度：反演变异的 std（相对于 g 真值变异），负对照应 ≈0
    return {
        "g": None,
        "corr": corr,
        "scale_std": scale,
        "norm_rmse": rmse / (std_true if std_true > 1e-15 else 1.0),
        "g_inv_std": std_inv,
        "g_true_std": std_true,
        "u_range": (min(U_c), max(U_c)),
        "n_points": len(g_true_c),
    }


# ---------- 主流程 ----------
def main():
    results = {"meta": {
        "task": "B3/R32 残差反演 g（方法准备，非经验验证）",
        "model": "A_D⁺（15-v3.2-gauge-fixed §5.7 推论 R43）",
        "rho_exp": RHO_EXP, "kappa_E": KAPPA_E, "rho": RHO,
        "step": STEP, "n": N,
        "engine": "仅标准库；确定性合成（无随机）",
    }, "cases": {}, "negative_controls": {}}

    for name, g in GFUNCS.items():
        diag = invert_g(g)
        diag["g"] = name
        results["cases"][name] = diag
        print(f"[case] g={name:10s} corr={diag['corr']:.6f} "
              f"scale={diag['scale_std']:.4f} norm_rmse={diag['norm_rmse']:.4f} "
              f"u_range=({diag['u_range'][0]:.3f},{diag['u_range'][1]:.3f})")

    # 负对照：g≡0（A_D 退化）与 g≡const（无变异）
    for name, g in [("g_zero", g_zero), ("g_const", g_const)]:
        diag = invert_g(g)
        results["negative_controls"][name] = diag
        print(f"[neg ] g={name:10s} corr={diag['corr']:.4f} g_inv_std={diag['g_inv_std']:.2e} "
              f"g_true_std={diag['g_true_std']:.2e}")

    out_path = os.path.join(HERE, "b3_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\nwrote: {out_path}")

    # 判据汇总
    print("\n== 判据 ==")
    for name, d in results["cases"].items():
        ok = d["corr"] > 0.99 and d["norm_rmse"] < 0.05
        print(f"{name:10s} corr>0.99:{d['corr']>0.99!s} norm_rmse<0.05:{d['norm_rmse']<0.05!s} -> {'PASS' if ok else 'FAIL'}")
    for name, d in results["negative_controls"].items():
        ok = d["g_inv_std"] < 0.02
        print(f"{name:10s} 无虚假结构(g_inv_std<0.02):{d['g_inv_std']<0.02!s} -> {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
