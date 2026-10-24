#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""OPEN-4 —— A/B 级噪声方差（纯仿真）：有噪声下「共线性法 vs 峰距法」的鲁棒性

背景
----
16-v3.2 §8.2 的无噪声原型确认（= R30 的 V3）显示：共线性法（OLS dU/dt ~ dE/dt 读 ρ_exp）
跨调度族误差 ~1e-15，而峰距法（推论 6.5：|ln ρ|_hat = |k*_E - k*_U|·L/(K-1)）在非对数调度下
偏 20%–76%。但该比较是**无噪声**的；§8.2 边界与 §13.2 明确登记：
**"有噪声时两种方法的方差对比尚未做"**（即 OPEN-4）。

本实验补上这一格：给 A_D 轨迹加 A/B 级**观测噪声**，重跑两种判读，比较方差与鲁棒性。
- A 级噪声：加在 E 观测上（ρ_exp 斜率读取所依赖的可观测量）
- B 级噪声：加在 U 观测上（截距符号 / 峰序读取所依赖的可观测量）

方法学纪律（承 R30）
--------------------
- 仅标准库；固定种子；确定性（R30 的 rng_for 模式）；机械可复核。
- 合成 A_D 数据，非经验验证；结论只回答"有噪声下哪种判读更稳"，不构成对现实断言。

判读方法（复用 R30 的 V3）
--------------------------
- 共线性法：从带噪观测 (E_obs, U_obs) 做一阶差分得 (dE, dU)，OLS 斜率 ⟹ ρ_exp_hat；
  截距符号 ⟹ 峰序（ρ>1 ⟺ 截距<0）。
- 峰距法：从带噪 E_obs, U_obs 取峰位，dk = |k*_E - k*_U|，|ln ρ|_hat = dk·L/(K-1)。

噪声尺度
--------
σ_A = σ_rel · std(E)  （A 级，相对 E 信号幅度）
σ_B = σ_rel · std(U)  （B 级，相对 U 信号幅度）
σ_rel ∈ {0, 0.001, 0.01, 0.05}；0 = 无噪声基准对照。
"""

import json
import math
import random
import zlib
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ 参数（承 R30）

K = 241
DT = 0.01
R_HI = 9.0
R_LO = 1.0 / 9.0
LN_HI = math.log(R_HI)
LN_LO = math.log(R_LO)
L_SPAN = abs(LN_LO - LN_HI)
E_START = 5.0
U_START = 30.0
SEED = 20261024

SCHED_FAMILIES = ("log_ramp", "linear", "sigmoid", "const_blocks")
RHO_GRID = (0.5, 2.0)          # 避开 ρ=1 顶点带（R34 已解）
SIGMA_REL = (0.0, 0.001, 0.01, 0.05)
NOISE_MODES = ("A_only", "B_only", "both")   # A=加在E观测; B=加在U观测
N_TRAIL = 200


def rng_for(off):
    return random.Random(SEED + zlib.crc32(off.encode("utf-8")))


# 以下四函数源自 experiments/r30_gauge_fixing/r30_gauge_fixing_identity.py（复用，不修改）

def make_cfg(theta_E, rho, e_plus=1.0, u_plus=1.0):
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
        "kappa_E": e_minus,
    }


def sched(name, k, KK=None):
    KK = KK or K
    t = k / float(KK - 1) if KK > 1 else 0.0
    if name == "log_ramp":
        return math.exp(LN_HI + (LN_LO - LN_HI) * t)
    if name == "linear":
        return R_HI + (R_LO - R_HI) * t
    if name == "sigmoid":
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


def trajectory(cfg, sched_name, KK=None):
    KK = KK or K
    e_p = cfg["e_plus"]; e_m = cfg["e_minus"]
    u_p = cfg["u_plus"]; u_m = cfg["u_minus"]
    Es = [E_START]; Us = [U_START]
    for k in range(KK):
        r = sched(sched_name, k, KK)
        dE = DT * (r * e_p - e_m)
        dU = DT * (r * u_p - u_m)
        Es.append(Es[-1] + dE)
        Us.append(Us[-1] + dU)
    return Es, Us


def argmax_plateau(vals):
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


# ------------------------------------------------------------------ OPEN-4 判读

def add_noise(sig, sigma, rng):
    return [v + rng.gauss(0.0, sigma) for v in sig]


def estimate(Es, Us, rho_true, rho):
    """从 (E, U) 观测做两种判读，返回字典。"""
    # 共线性法：一阶差分 OLS
    dE = [Es[i + 1] - Es[i] for i in range(len(Es) - 1)]
    dU = [Us[i + 1] - Us[i] for i in range(len(Us) - 1)]
    a, b = _ols(dE, dU)
    rho_exp_hat = b
    err_coll = abs(rho_exp_hat - rho_true) / rho_true
    sign_ok = (a < 0) == (rho > 1.0)
    # 峰距法
    kE = argmax_plateau(Es)
    kU = argmax_plateau(Us)
    dk = abs(kE - kU)
    lnrho_hat = dk * L_SPAN / (K - 1)
    lnrho_true = abs(math.log(rho))
    err_peak = (abs(lnrho_hat - lnrho_true) / lnrho_true
                if lnrho_true > 1e-12 else float("nan"))
    return {
        "rho_exp_hat": rho_exp_hat, "err_coll": err_coll,
        "sign_ok": sign_ok, "kE": kE, "kU": kU, "dk": dk,
        "lnrho_hat": lnrho_hat, "err_peak": err_peak,
    }


def stats(xs):
    n = len(xs)
    m = sum(xs) / n
    var = sum((x - m) ** 2 for x in xs) / n
    return m, math.sqrt(var)


# ------------------------------------------------------------------ 主流程

def main():
    theta_E = 1.0
    report = {"meta": {
        "task": "OPEN-4 A/B 级噪声方差（纯仿真）——有噪声下共线性法 vs 峰距法鲁棒性",
        "nature": "合成 A_D 数据；非经验验证；补齐 16 §8.2 / §13.2 '有噪声未做' 的登记格",
        "params": {
            "K": K, "DT": DT, "theta_E": theta_E, "E_START": E_START, "U_START": U_START,
            "RHO_GRID": list(RHO_GRID), "SCHED_FAMILIES": list(SCHED_FAMILIES),
            "SIGMA_REL": list(SIGMA_REL), "N_TRAIL": N_TRAIL, "SEED": SEED,
        },
        "noise_definition": "A级=加在E观测(σ·std(E)); B级=加在U观测(σ·std(U)); 0=无噪声基准; "
                            "mode: A_only/B_only/both",
    }, "rows": []}

    lines = []
    def p(s=""):
        lines.append(str(s))

    p("=" * 100)
    p("OPEN-4 A/B 级噪声方差：共线性法 vs 峰距法（有噪声，N=%d trail/配置）" % N_TRAIL)
    p("=" * 100)

    for rho in RHO_GRID:
        cfg = make_cfg(theta_E, rho)
        rho_true = cfg["rho_exp"]
        for sched_name in SCHED_FAMILIES:
            Es, Us = trajectory(cfg, sched_name)
            std_E = (sum((v - sum(Es) / len(Es)) ** 2 for v in Es) / len(Es)) ** 0.5
            std_U = (sum((v - sum(Us) / len(Us)) ** 2 for v in Us) / len(Us)) ** 0.5
            for sigma in SIGMA_REL:
                for mode in NOISE_MODES:
                    rng = rng_for("open4:%s:%s:%s:%s" % (sched_name, rho, sigma, mode))
                    err_coll_list = []
                    err_peak_list = []
                    sign_ok_count = 0
                    for _ in range(N_TRAIL):
                        sigma_E = sigma * std_E if mode in ("A_only", "both") else 0.0
                        sigma_U = sigma * std_U if mode in ("B_only", "both") else 0.0
                        E_obs = add_noise(Es, sigma_E, rng)
                        U_obs = add_noise(Us, sigma_U, rng)
                        r = estimate(E_obs, U_obs, rho_true, rho)
                        err_coll_list.append(r["err_coll"])
                        err_peak_list.append(r["err_peak"])
                        if r["sign_ok"]:
                            sign_ok_count += 1
                    m_coll, s_coll = stats(err_coll_list)
                    m_peak, s_peak = stats(err_peak_list)
                    sign_ok_rate = sign_ok_count / float(N_TRAIL)
                    p("  rho=%.2f sched=%-12s sigma=%-6g mode=%-7s | "
                      "err_coll %.3e(std %.3e) | err_peak %.3e(std %.3e) | sign_ok %.2f"
                      % (rho, sched_name, sigma, mode,
                         m_coll, s_coll, m_peak, s_peak, sign_ok_rate))
                    report["rows"].append({
                        "rho": rho, "sched": sched_name, "sigma_rel": sigma, "mode": mode,
                        "std_E": std_E, "std_U": std_U,
                        "err_coll_mean": m_coll, "err_coll_std": s_coll,
                        "err_peak_mean": m_peak, "err_peak_std": s_peak,
                        "sign_ok_rate": sign_ok_rate,
                    })

    p("")
    p("=" * 100)
    p("汇总（跨 rho×sched 中位）：共线性法 err_coll vs 峰距法 err_peak，按 mode")
    p("=" * 100)
    for mode in NOISE_MODES:
        p("  mode=%-7s" % mode)
        for sigma in SIGMA_REL:
            sub = [r for r in report["rows"] if r["sigma_rel"] == sigma and r["mode"] == mode]
            if not sub:
                continue
            mc = sorted(r["err_coll_mean"] for r in sub)
            mp = sorted(r["err_peak_mean"] for r in sub)
            sc = sorted(r["err_coll_std"] for r in sub)
            sp = sorted(r["err_peak_std"] for r in sub)
            p("    sigma=%-6g | 共线性 err 中位 %10.3e (std %10.3e) | "
              "峰距 err 中位 %10.3e (std %10.3e)"
              % (sigma, mc[len(mc) // 2], sc[len(sc) // 2], mp[len(mp) // 2], sp[len(sp) // 2]))

    out_path = os.path.join(HERE, "open4_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)

    txt = "\n".join(lines)
    print(txt)
    print("\nwrote:", out_path)


if __name__ == "__main__":
    main()
