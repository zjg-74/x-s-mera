#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R30 —— A_D 的规范固定：恒等式、规范不变性、与调度无关的可观测量

本脚本是「工程补丁源头归因审计」的支撑证据，不是新实验臂。

背景
----
补篇五 §1.2 的 A_D 写作

    dE = a*e_+ - b*e_-          dU = a*u_+ - b*u_-        r := a/b

这里 (a, b, e_+, e_-) 四个量只通过 a*e_+ 与 b*e_- 两个组合进入动力学，
因此存在**参数化冗余**。补篇九定理 6.12 已把这个冗余显式化为 λ-规范自由度

    G_lambda :  r -> lambda*r ,  e_+ -> e_+/lambda ,  u_+ -> u_+/lambda

它保持 (E(k), U(k)) 逐点不变。⟹ 四个系数的绝对值不是物理量。

本轮的源头修改（规范固定）
--------------------------
引入**归一化投入比** s := r/theta_E（无量纲）与**单一时间尺度** kappa_E := b*e_-，
则 A_D 等价于

    dE/dt = kappa_E * (s - 1)
    dU/dt = kappa_E * rho_exp * (s - rho)          rho := rho_res/rho_exp = theta_U/theta_E

消去 s ⟹ **一条仿射关系**

    dU/dt = rho_exp * (dE/dt)  +  kappa_E * rho_exp * (1 - rho)          (*)

内禀自由参数从「4 系数 + r 过程」降到 **3 个无量纲数 (s, rho_exp, rho) + 1 个时间单位 kappa_E**，
且 (s, rho_exp, rho, kappa_E) 全部 λ-规范不变 ⟹ **规范自由度在定义层被消灭**。

(*) 是本审计的核心断言。它一次性蕴含：

  定理 6.10（U - rho_exp*E 仿射守恒）        = (*) 的积分形式
  定理 6.11（增量云共线、残差恒为 0）        = (*) 本身
  推论 6.4（峰序 rho>1 ⟺ U 峰在前）          = (*) 截距的符号
  推论 6.5（峰距 ∝ |ln rho|）                = (*) + 调度假设（**唯一需要调度的一环**）

本脚本验证三件事
----------------
V1  (*) 逐样本精确成立（残差 ~1e-15）
V2  λ 规范变换下 (s, rho_exp, rho, kappa_E) 不变，且 (E, U) 逐点不变
V3  由 (*) 反解 rho_exp **不依赖调度族**；而峰距法（推论 6.5）**依赖**
    ⟹ 直接量化 B-B34 的 36.8%–76.5% 偏差，并证明换用共线性可消除之

方法学陷阱（本轮新增，T17–T18）
--------------------------------
T17. **「不可辨识」有两种，必须分开**：(a) 规范冗余造成的**伪**不可辨识，
     源头改参数化即可消灭；(b) 真不可辨识（如 kappa_E 需要外部单位）。
     把 (a) 当成 (b) 会徒劳地去找外部标定。B-B23 属于 (a) 而非 (b)。
T18. **调度依赖性的来源可定位到具体参数**：ρ 的**符号**（峰序）不需要 kappa_E，
     ρ 的**数值**（峰距）需要。⟹ 「峰距法调度依赖」不是公式的缺陷，
     而是它试图读取一个本就需要额外标定的量。

仅标准库；种子固定，逐位可复现。
"""

import json
import math
import random
import zlib

# ------------------------------------------------------------------ 参数

K = 241
DT = 0.01
R_HI = 9.0
R_LO = 1.0 / 9.0
LN_HI = math.log(R_HI)
LN_LO = math.log(R_LO)
L_SPAN = abs(LN_LO - LN_HI)          # 4.394449...
E_START = 5.0
U_START = 30.0                        # 抬高，避免 U 触底（承 R29 先例）
SEED = 20261022

SCHED_FAMILIES = ("log_ramp", "linear", "sigmoid", "const_blocks")
LAM_GRID = (0.25, 0.5, 1.0, 2.0, 4.0)
RHO_GRID = (1.0 / 3.0, 0.5, 1.0, 2.0, 3.0)
THETA_E_GRID = (0.5, 1.0, 2.0)


def rng_for(off):
    return random.Random(SEED + zlib.crc32(off.encode("utf-8")))


def make_cfg(theta_E, rho, e_plus=1.0, u_plus=1.0):
    """以 (theta_E, rho, e_+, u_+) 参数化；rho = rho_res/rho_exp。"""
    e_minus = theta_E * e_plus
    rho_exp = u_plus / e_plus
    rho_res = rho * rho_exp
    u_minus = rho_res * e_minus
    theta_U = u_minus / u_plus
    return {
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": theta_E, "theta_U": theta_U,
        "rho_exp": rho_exp, "rho_res": rho_res, "rho": rho,
        "kappa_E": e_minus,               # b == 1
    }


def sched(name, k, KK=None):
    """r 调度族。KK 用于归一化。"""
    KK = KK or K
    t = k / float(KK - 1) if KK > 1 else 0.0
    if name == "log_ramp":
        return math.exp(LN_HI + (LN_LO - LN_HI) * t)
    if name == "linear":
        return R_HI + (R_LO - R_HI) * t
    if name == "sigmoid":
        # 在 log 空间做 logistic：两端与 log_ramp 一致，中部陡峭
        z = 1.0 / (1.0 + math.exp(-10.0 * (t - 0.5)))
        z0 = 1.0 / (1.0 + math.exp(5.0))
        z1 = 1.0 / (1.0 + math.exp(-5.0))
        return math.exp(LN_HI + (LN_LO - LN_HI) * (z - z0) / (z1 - z0))
    if name == "const_blocks":
        nb = 8
        i = min(int(t * nb), nb - 1)
        tt = (i + 0.5) / nb
        return math.exp(LN_HI + (LN_LO - LN_HI) * tt)
    raise ValueError(name)


def trajectory(cfg, sched_name, KK=None, lam=1.0):
    """b == 1，a = r(k)。lam 为规范变换因子（作用于 r, e_+, u_+）。"""
    KK = KK or K
    e_p = cfg["e_plus"] / lam
    u_p = cfg["u_plus"] / lam
    e_m = cfg["e_minus"]
    u_m = cfg["u_minus"]
    Es = [E_START]
    Us = [U_START]
    dEs = []
    dUs = []
    ss = []
    for k in range(KK):
        r = lam * sched(sched_name, k, KK)
        dE = DT * (r * e_p - e_m)
        dU = DT * (r * u_p - u_m)
        dEs.append(dE / DT)
        dUs.append(dU / DT)
        ss.append(r / (e_m / e_p))          # s := r / theta_E（变换后的 theta_E）
        Es.append(Es[-1] + dE)
        Us.append(Us[-1] + dU)
    return Es, Us, dEs, dUs, ss


def argmax_plateau(vals):
    """平顶取最左（承 R24/T3 修正）。"""
    m = max(vals)
    for i, v in enumerate(vals):
        if v >= m - 1e-15:
            return i
    return 0


def _ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx > 0 else 0.0
    return my - b * mx, b


# ------------------------------------------------------------------ V1

def v1_identity():
    """(*) dU/dt = rho_exp*(dE/dt) + kappa_E*rho_exp*(1-rho)  逐样本成立？"""
    rows = []
    for th in THETA_E_GRID:
        for rho in RHO_GRID:
            for name in SCHED_FAMILIES:
                cfg = make_cfg(th, rho)
                Es, Us, dEs, dUs, ss = trajectory(cfg, name)
                pred_c = cfg["kappa_E"] * cfg["rho_exp"] * (1.0 - cfg["rho"])
                worst = 0.0
                for de, du in zip(dEs, dUs):
                    resid = du - cfg["rho_exp"] * de - pred_c
                    worst = max(worst, abs(resid))
                rows.append({
                    "theta_E": th, "rho": rho, "sched": name,
                    "max_abs_resid": worst,
                    "scale": max(abs(v) for v in dUs),
                })
    return rows


# ------------------------------------------------------------------ V2

def v2_gauge():
    """λ 规范变换下 (s, rho_exp, rho, kappa_E) 不变 + (E,U) 逐点不变？"""
    rows = []
    th, rho = 1.0, 2.0
    cfg = make_cfg(th, rho)
    E0, U0, dEs0, dUs0, s0 = trajectory(cfg, "log_ramp", lam=1.0)
    for lam in LAM_GRID:
        EL, UL, dEsL, dUsL, sL = trajectory(cfg, "log_ramp", lam=lam)
        dE = max(abs(a - b) for a, b in zip(EL, E0))
        dU = max(abs(a - b) for a, b in zip(UL, U0))
        ds = max(abs(a - b) for a, b in zip(sL, s0))
        rows.append({
            "lambda": lam,
            "theta_E_transformed": (cfg["e_minus"]) / (cfg["e_plus"] / lam),
            "max_abs_E_diff": dE,
            "max_abs_U_diff": dU,
            "max_abs_s_diff": ds,
            "rho_exp": cfg["rho_exp"],
            "rho": cfg["rho"],
            "kappa_E": cfg["kappa_E"],
        })
    return rows


# ------------------------------------------------------------------ V3

def v3_sched_independence():
    """共线性反解 rho_exp 是否调度无关？峰距法是否调度依赖？

    两条路线：
      A 共线性（本轮提出）：OLS  dU/dt ~ dE/dt  ⟹ 斜率 = rho_exp
      B 峰距（推论 6.5）  ：|ln rho|_hat = |k*_E - k*_U| * L / (K-1)
    """
    rows = []
    for th in THETA_E_GRID:
        for rho in RHO_GRID:
            cfg = make_cfg(th, rho)
            lnrho_true = abs(math.log(rho))
            for name in SCHED_FAMILIES:
                Es, Us, dEs, dUs, ss = trajectory(cfg, name)
                # A：共线性
                a, b = _ols(dEs, dUs)
                rho_exp_hat = b
                err_A = abs(rho_exp_hat - cfg["rho_exp"]) / cfg["rho_exp"]
                # 截距符号 ⟹ 峰序（rho>1 <=> 截距<0）
                sign_ok = (a < 0) == (cfg["rho"] > 1.0)
                # B：峰距
                kE = argmax_plateau(Es)
                kU = argmax_plateau(Us)
                dk = abs(kE - kU)
                lnrho_hat = dk * L_SPAN / (K - 1)
                err_B = (abs(lnrho_hat - lnrho_true) / lnrho_true
                         if lnrho_true > 1e-12 else float("nan"))
                rows.append({
                    "theta_E": th, "rho": rho, "sched": name,
                    "rho_exp_hat": rho_exp_hat,
                    "rel_err_collinearity": err_A,
                    "intercept": a,
                    "peak_order_sign_ok": sign_ok,
                    "kE": kE, "kU": kU, "dk": dk,
                    "ln_rho_true": lnrho_true,
                    "ln_rho_hat": lnrho_hat,
                    "rel_err_peakdistance": err_B,
                })
    return rows


# ------------------------------------------------------------------ 输出

def main():
    out = []

    def p(s=""):
        out.append(str(s))

    p("=" * 78)
    p("R30 —— A_D 规范固定的三项验证")
    p("=" * 78)

    p("")
    p("--- V1  恒等式 (*)  dU/dt = rho_exp*(dE/dt) + kappa_E*rho_exp*(1-rho) ---")
    V1 = v1_identity()
    p("  {:>8} {:>6} {:>13} | {:>14} {:>14}".format(
        "theta_E", "rho", "sched", "max|resid|", "scale|dU/dt|"))
    for r in V1:
        p("  {:>8} {:>6} {:>13} | {:>14.3e} {:>14.4f}".format(
            r["theta_E"], r["rho"], r["sched"],
            r["max_abs_resid"], r["scale"]))
    w1 = max(r["max_abs_resid"] for r in V1)
    p("  逐样本最大残差 = {:.3e}".format(w1))
    p("  预言 V1：残差为浮点噪声（~1e-15），与 theta_E / rho / 调度族均无关")

    p("")
    p("--- V2  规范不变性：lambda 变换下 (s, rho_exp, rho, kappa_E) 与 (E,U) ---")
    V2 = v2_gauge()
    p("  {:>7} {:>10} | {:>14} {:>14} {:>14} | {:>9} {:>6} {:>9}".format(
        "lambda", "theta_E'", "max|dE|", "max|dU|", "max|ds|",
        "rho_exp", "rho", "kappa_E"))
    for r in V2:
        p("  {:>7} {:>10.4f} | {:>14.3e} {:>14.3e} {:>14.3e} | {:>9.4f} {:>6.2f} {:>9.4f}".format(
            r["lambda"], r["theta_E_transformed"],
            r["max_abs_E_diff"], r["max_abs_U_diff"], r["max_abs_s_diff"],
            r["rho_exp"], r["rho"], r["kappa_E"]))
    p("  预言 V2：theta_E 变了 lam 倍，而 (E,U,s,rho_exp,rho,kappa_E) 全不变")

    p("")
    p("--- V3  调度无关性：共线性法 vs 峰距法 ---")
    V3 = v3_sched_independence()
    p("  {:>8} {:>6} {:>13} | {:>16} | {:>10} {:>7} | {:>16} {:>7}".format(
        "theta_E", "rho", "sched", "err(共线性)", "截距符号", "峰序",
        "err(峰距)", "dk"))
    for r in V3:
        p("  {:>8} {:>6} {:>13} | {:>16.3e} | {:>10} {:>7} | {:>16.4f} {:>7}".format(
            r["theta_E"], r["rho"], r["sched"],
            r["rel_err_collinearity"],
            "ok" if r["peak_order_sign_ok"] else "WRONG",
            "-" if r["rho"] == 1.0 else ("+" if r["rho"] > 1 else "-"),
            r["rel_err_peakdistance"], r["dk"]))
    p("")
    p("  按调度族汇总（峰距法相对误差的中位数 / 共线性法的最大误差）：")
    for name in SCHED_FAMILIES:
        sub = [r for r in V3 if r["sched"] == name and r["rho"] != 1.0]
        if not sub:
            continue
        eb = sorted(r["rel_err_peakdistance"] for r in sub)
        med_b = eb[len(eb) // 2]
        max_a = max(r["rel_err_collinearity"] for r in sub)
        p("    {:>13} | 峰距 err 中位 {:>8.2%} | 共线性 err 最大 {:>10.3e}".format(
            name, med_b, max_a))
    p("  预言 V3：共线性误差跨调度族恒定（~1e-15）；峰距误差只在 log_ramp 下小")

    p("")
    p("=" * 78)

    result = {
        "meta": {
            "K": K, "DT": DT, "R_HI": R_HI, "R_LO": R_LO,
            "L_SPAN": L_SPAN, "E_START": E_START, "U_START": U_START,
            "SEED": SEED,
            "SCHED_FAMILIES": list(SCHED_FAMILIES),
            "LAM_GRID": list(LAM_GRID),
            "nature": "规范固定恒等式的代数/数值验证，不是实证验证",
        },
        "V1_identity": V1,
        "V2_gauge": V2,
        "V3_sched_independence": V3,
    }
    with open("r30_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))


if __name__ == "__main__":
    main()
