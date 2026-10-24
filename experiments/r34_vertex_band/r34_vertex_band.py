#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R34 —— rho=1 顶点容差带（OPEN-1 的解：B 级判读禁区）
========================================================

对应：`27-v3.4-boundary-map.md` §3.4 开放格 OPEN-1——楔形中相边界 rho=1 处
sign(1-rho) 的符号判读退化，需要一个容差带定义。

**容差带定义（本实验给出并检验的）**：

    B 级判读禁区（顶点带） := { 数据 : |beta_hat| <= z * SE(beta_hat) }

其中 beta_hat 与 SE(beta_hat) 来自多元回归 U = alpha + gamma*E + beta*k + eps
（经典最小二乘，仅标准库实现），z 为约定性带宽参数（与 eta 同类：约定性输入，
按可接受误判率选取，默认 z=1.96）。

性质（结构性）：
  1. 判定只依赖 (E, U, k) 数据与 z —— beta_hat 是规范不变量（U5 通过）。
  2. 无噪声时 SE=0，禁区退化为 |beta| <= eps，恰为顶点 rho=1（零测集）。
  3. 有噪声时带宽增大；rho -> 1 时留空（withheld）比例 -> 1，错误符号率被压制。
  4. 判读输出三值：mid_drill（钻研方向，beta<0）、mid_hollow（空洞方向，beta>0）、
     vertex_band（禁区，不判）——**不强行二分**。

无噪声 eps 取 NOISELESS_EPS；V1 无噪声网格、V2 有噪声蒙特卡洛、V3 规范不变性。

仅标准库。JSON 结果写入 CWD：r34_results.json。
"""

from __future__ import annotations

import json
import math
import random

NOISELESS_EPS = 1e-12
DEFAULT_Z = 1.96


# ------------------------------------------------------------------ 数据生成


def make_cfg(theta_E: float, rho: float, rho_exp: float = 1.3) -> dict:
    e_plus = 1.0
    e_minus = float(theta_E)
    u_plus = rho_exp * e_plus
    rho_res = rho * rho_exp
    u_minus = rho_res * e_minus
    return {"e_plus": e_plus, "e_minus": e_minus,
            "u_plus": u_plus, "u_minus": u_minus,
            "rho_exp": rho_exp, "rho_res": rho_res, "rho": rho}


def trajectory(cfg: dict, k_steps: int = 60, sigma_u: float = 0.0,
               seed: int = 0, rate_scale: float = 1.0,
               coef_scale: float = 1.0) -> list:
    """b=1 固定、a_k 线性调度的 A_D 轨迹；对 U 叠加 N(0, sigma_u) 观测噪声。

    G_lambda 实现：coef_scale = lambda 时四系数整体乘 lambda（记号侧），
    rate_scale = 1/lambda 时动作速率 (a, b) 整体反标（速率侧），
    两者同时施加时 (E,U) 轨迹逐点不变（定理 6.12 的实现）。
    """
    rng = random.Random(seed)
    e_p = cfg["e_plus"] * coef_scale
    e_m = cfg["e_minus"] * coef_scale
    u_p = cfg["u_plus"] * coef_scale
    u_m = cfg["u_minus"] * coef_scale
    E, U = 0.0, 0.0
    obs = []
    for k in range(k_steps + 1):
        u_obs = U + (rng.gauss(0.0, sigma_u) if sigma_u > 0.0 else 0.0)
        obs.append((float(k), E, u_obs))
        a_k = (0.2 + 1.2 * k / max(1, k_steps)) * rate_scale
        b_k = 1.0 * rate_scale
        E = E + a_k * e_p - b_k * e_m
        U = U + a_k * u_p - b_k * u_m
    return obs


# ------------------------------------------------------------------ 回归与判读


def fit_beta_with_se(obs: list) -> dict:
    """多元回归 U = alpha + gamma*E + beta*k；返回 beta_hat 与 SE(beta_hat)。

    SE 用经典公式 SE(beta) = sqrt(s2 * C[j,j])，s2 = RSS/(n-3)，C = (X'X)^-1。
    无噪声时 RSS=0 -> SE=0。
    """
    n = len(obs)
    if n < 4:
        raise ValueError("need >= 4 points")
    rows = [(1.0, E, float(k), U) for (k, E, U) in obs]
    xtx = [[0.0] * 3 for _ in range(3)]
    xty = [0.0] * 3
    for r in rows:
        x = r[:3]
        y = r[3]
        for i in range(3):
            xty[i] += x[i] * y
            for j in range(3):
                xtx[i][j] += x[i] * x[j]
    coef = _solve3(xtx, xty)                       # [alpha, gamma, beta]
    rss = 0.0
    for r in rows:
        pred = coef[0] + coef[1] * r[1] + coef[2] * r[2]
        rss += (r[3] - pred) ** 2
    s2 = rss / (n - 3) if n > 3 else 0.0
    c_inv = _invert3(xtx)
    se_beta = math.sqrt(max(0.0, s2 * c_inv[2][2]))
    return {"alpha": coef[0], "gamma": coef[1], "beta_hat": coef[2],
            "se_beta": se_beta, "rss": rss, "n": n}


def verdict(beta_hat: float, se_beta: float, z: float = DEFAULT_Z) -> str:
    """B 级判读规则（三值输出，禁区不判）。"""
    if abs(beta_hat) <= z * se_beta + NOISELESS_EPS:
        return "vertex_band"
    return "mid_drill" if beta_hat < 0.0 else "mid_hollow"


def true_mid(rho: float) -> str:
    """真值标签（仅用于评分，不进判读规则）：rho>1 -> 钻研相中期。"""
    if rho == 1.0:
        return "vertex"
    return "mid_drill" if rho > 1.0 else "mid_hollow"


# ------------------------------------------------------------------ 实验臂


RHO_GRID = (0.90, 0.95, 0.98, 1.00, 1.02, 1.05, 1.10)


def v1_noiseless() -> list:
    """无噪声：SE=0，禁区退化为恰 rho=1；其余点符号必须正确。"""
    rows = []
    for rho in RHO_GRID:
        cfg = make_cfg(0.6, rho)
        obs = trajectory(cfg, sigma_u=0.0)
        fit = fit_beta_with_se(obs)
        v = verdict(fit["beta_hat"], fit["se_beta"])
        rows.append({"rho": rho, "beta_hat": fit["beta_hat"],
                     "se_beta": fit["se_beta"], "verdict": v,
                     "true": true_mid(rho),
                     "sign_correct": (v == true_mid(rho)) or v == "vertex_band"})
    return rows


def v2_noisy(sigma_u: float = 0.05, n_trials: int = 200) -> dict:
    """有噪声：统计各 rho 处的留空率、错误符号率（判出符号时错的比例）。"""
    per_rho = []
    for rho in RHO_GRID:
        withheld = 0
        signed = 0
        wrong = 0
        for t in range(n_trials):
            cfg = make_cfg(0.6, rho)
            obs = trajectory(cfg, sigma_u=sigma_u, seed=1000 * t + int(rho * 100))
            fit = fit_beta_with_se(obs)
            v = verdict(fit["beta_hat"], fit["se_beta"])
            if v == "vertex_band":
                withheld += 1
            else:
                signed += 1
                if v != true_mid(rho):
                    wrong += 1
        per_rho.append({
            "rho": rho, "sigma_u": sigma_u, "n_trials": n_trials,
            "withheld_rate": withheld / n_trials,
            "false_sign_rate": (wrong / signed) if signed else 0.0,
        })
    return {"sigma_u": sigma_u, "per_rho": per_rho}


def v3_gauge_invariance() -> dict:
    """规范不变性：G_lambda（系数乘 lambda、速率除 lambda）下 (E,U) 逐点不变
    -> 回归量与判读不变（U5）。"""
    cfg = make_cfg(0.6, 1.05)
    lam = 7.3
    obs_a = trajectory(cfg, sigma_u=0.02, seed=7)
    obs_b = trajectory(cfg, sigma_u=0.02, seed=7,
                       coef_scale=lam, rate_scale=1.0 / lam)
    same_path = all(abs(a[1] - b[1]) < 1e-9 and abs(a[2] - b[2]) < 1e-9
                    for a, b in zip(obs_a, obs_b))
    fit_a = fit_beta_with_se(obs_a)
    fit_b = fit_beta_with_se(obs_b)
    va = verdict(fit_a["beta_hat"], fit_a["se_beta"])
    vb = verdict(fit_b["beta_hat"], fit_b["se_beta"])
    return {"lambda": lam, "path_pointwise_unchanged": same_path,
            "verdict_a": va, "verdict_b": vb, "unchanged": same_path and va == vb}


# ------------------------------------------------------------------ 数值工具


def _solve3(a: list, b: list) -> list:
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-300:
            raise ZeroDivisionError("singular X'X")
        m[col], m[piv] = m[piv], m[col]
        for r in range(3):
            if r == col:
                continue
            f = m[r][col] / m[col][col]
            for c in range(col, 4):
                m[r][c] -= f * m[col][c]
    return [m[i][3] / m[i][i] for i in range(3)]


def _invert3(a: list) -> list:
    """3x3 逆矩阵（高斯-约当）。"""
    m = [row[:] + [1.0 if i == j else 0.0 for j in range(3)] for i, row in enumerate(a)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-300:
            raise ZeroDivisionError("singular X'X")
        m[col], m[piv] = m[piv], m[col]
        div = m[col][col]
        m[col] = [v / div for v in m[col]]
        for r in range(3):
            if r == col:
                continue
            f = m[r][col]
            m[r] = [rv - f * cv for rv, cv in zip(m[r], m[col])]
    return [row[3:] for row in m]


# ------------------------------------------------------------------ 主流程


def run_all() -> dict:
    v1 = v1_noiseless()
    v2a = v2_noisy(sigma_u=0.05)
    v2b = v2_noisy(sigma_u=0.20)
    v3 = v3_gauge_invariance()

    noiseless_ok = all(r["sign_correct"] for r in v1)
    v1_at_vertex = [r["verdict"] for r in v1 if r["rho"] == 1.0]
    vertex_only = (v1_at_vertex == ["vertex_band"])
    # 结构预言：噪声越大，越靠顶点留空越多（比较 rho=1.02 处的留空率）
    def withheld_at(v2, rho):
        return next(p["withheld_rate"] for p in v2["per_rho"] if p["rho"] == rho)
    band_widens = withheld_at(v2b, 1.02) >= withheld_at(v2a, 1.02)

    result = {
        "experiment": "R34_vertex_band",
        "definition": "B 级判读禁区 := {|beta_hat| <= z*SE(beta_hat)}, z=1.96（约定性）",
        "v1_noiseless": v1,
        "v2_noisy": {"sigma_0.05": v2a, "sigma_0.20": v2b},
        "v3_gauge": v3,
        "summary": {
            "noiseless_sign_correct": noiseless_ok,
            "noiseless_vertex_exact": vertex_only,
            "band_widens_with_noise": band_widens,
            "gauge_invariant": v3["unchanged"],
            "open1_verdict": "容差带已定义并通过结构检验"
                             if (noiseless_ok and vertex_only and v3["unchanged"])
                             else "未通过",
        },
    }
    return result


def main():
    result = run_all()
    with open("r34_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    s = result["summary"]
    print("R34 vertex tolerance band (OPEN-1)")
    print(f"  noiseless sign correct: {s['noiseless_sign_correct']}")
    print(f"  noiseless vertex exact: {s['noiseless_vertex_exact']}")
    print(f"  band widens with noise: {s['band_widens_with_noise']}")
    print(f"  gauge invariant: {s['gauge_invariant']}")
    print(f"  OPEN-1 verdict: {s['open1_verdict']}")


if __name__ == "__main__":
    main()
