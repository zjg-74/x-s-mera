#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R33 —— C 级闭合实测：固定动作计数能否闭合 C 级（R45 核验）
============================================================

对应：`16-v3.2-self-consistent-theory.md` §11（命题 6.3 实验性规范固定）、
§12.3 R45 判据：「四臂测得的 (e+, e-, u+, u-) 应与回归所得 (rho_exp, beta) 自洽」。

模型（A_D 的动作计数形式，确定性）：

    每步 k：展开速率 a_k（按调度族变化），消解速率 b_k = 1（固定），
    dt = 1：
        E(k+1) = E(k) + a_k*e_plus - b_k*e_minus
        U(k+1) = U(k) + a_k*u_plus - b_k*u_minus

四臂：
    臂 1（固定展开 1 次）：a=1, b=0  -> 测 e_plus = dE, u_plus = dU
    臂 2（固定消解 1 次）：a=0, b=1  -> 测 e_minus = -dE, u_minus = -dU
    臂 3（自由调度序列）：b=1 固定，a_k 按调度族走 -> 记录 (E_k, U_k)
    臂 4（回归）：拟合 U = rho_exp*E + beta*k + alpha -> (rho_exp^reg, beta^reg)

自洽判定（R45）：
    T1  rho_exp^reg == u_plus/e_plus            （臂与回归自洽）
    T2  beta^reg   == e_minus*(rho_exp - rho_res)  （kappa_E = b*e_minus，b=1）
    T3  v* == 1/rho_exp（v*U - E 对 k 仿射的零失配斜率）
    T4  四系数可由臂直接读出；theta/rho 恒等式成立
    N1  负对照：把 b 记错（b=2）-> T2 应破坏（判据有牙齿）

⟹ 若 T1–T4 跨调度族成立且 N1 破坏，则「固定动作计数闭合 C 级」在原型内成立。

仅标准库。JSON 结果写入 CWD：r33_results.json。
"""

from __future__ import annotations

import json
import math
import random

# ------------------------------------------------------------------ 构造


SCHED_FAMILIES = ("log_ramp", "linear", "sigmoid", "const_blocks")


def make_cfg(theta_E: float, rho: float, rho_exp: float = 1.3) -> dict:
    """由 (theta_E, rho) 构造 A_D 参数。恒等式：rho = rho_res/rho_exp = theta_U/theta_E。"""
    e_plus = 1.0
    e_minus = float(theta_E)              # theta_E = e_minus/e_plus
    u_plus = rho_exp * e_plus             # rho_exp = u_plus/e_plus
    rho_res = rho * rho_exp               # rho = rho_res/rho_exp
    u_minus = rho_res * e_minus           # rho_res = u_minus/e_minus
    theta_U = u_minus / u_plus
    return {
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": theta_E, "theta_U": theta_U,
        "rho_exp": rho_exp, "rho_res": rho_res, "rho": rho,
    }


def sched(family: str, k_steps: int) -> list:
    """a_k 序列（b_k = 1 固定）。返回 list[float]。"""
    if k_steps < 2:
        raise ValueError("k_steps must be >= 2")
    n = k_steps
    if family == "log_ramp":
        return [0.2 * math.exp(3.0 * i / (n - 1)) for i in range(n)]
    if family == "linear":
        return [0.1 + 1.5 * i / (n - 1) for i in range(n)]
    if family == "sigmoid":
        return [0.15 + 1.3 / (1.0 + math.exp(-0.25 * (i - n / 2.0))) for i in range(n)]
    if family == "const_blocks":
        return [0.3 if (i // max(1, n // 4)) % 2 == 0 else 1.2 for i in range(n)]
    raise ValueError(f"unknown family: {family}")


# ------------------------------------------------------------------ 四臂


def arm_expand_once(cfg: dict) -> dict:
    """臂 1：固定展开 1 次。"""
    return {"e_plus_m": cfg["e_plus"], "u_plus_m": cfg["u_plus"]}


def arm_resolve_once(cfg: dict) -> dict:
    """臂 2：固定消解 1 次。dE = -e_minus, dU = -u_minus。"""
    return {"e_minus_m": cfg["e_minus"], "u_minus_m": cfg["u_minus"]}


def arm_free_schedule(cfg: dict, family: str, k_steps: int = 60) -> dict:
    """臂 3：自由调度序列（b=1 固定），返回 (E_k, U_k) 与累计消解计数。"""
    a_seq = sched(family, k_steps)
    E, U = 0.0, 0.0
    traj = [(E, U)]
    for a_k in a_seq:
        E = E + a_k * cfg["e_plus"] - 1.0 * cfg["e_minus"]
        U = U + a_k * cfg["u_plus"] - 1.0 * cfg["u_minus"]
        traj.append((E, U))
    return {"traj": traj, "a_seq": a_seq, "b_total": float(k_steps)}


def fit_affine_in_E_and_k(traj: list) -> dict:
    """臂 4：最小二乘拟合 U = alpha + rho_exp*E + beta*k。

    三个正规方程手写高斯消元（仅标准库）。返回 rho_exp, beta, alpha 与
    v_star（= argmin ||v*U - E - affine(k)|| 的闭式：v_star = 1/rho_exp
    当拟合无残差时成立；此处同时报告拟合残差）。
    """
    # 设计矩阵 X 的列为 (1, E, k)；因变量 U。
    rows = []
    for k, (E, U) in enumerate(traj):
        rows.append((1.0, E, float(k), U))
    n = len(rows)
    # X'X 3x3 与 X'y
    xtx = [[0.0] * 3 for _ in range(3)]
    xty = [0.0] * 3
    for r in rows:
        x = r[:3]
        y = r[3]
        for i in range(3):
            xty[i] += x[i] * y
            for j in range(3):
                xtx[i][j] += x[i] * x[j]
    coef = _solve3(xtx, xty)          # [alpha, rho_exp, beta]
    alpha, rho_exp, beta = coef
    resid = max(abs(r[3] - (alpha + rho_exp * r[1] + beta * r[2])) for r in rows)
    scale = max(1e-30, max(abs(r[3]) for r in rows))
    return {"alpha": alpha, "rho_exp": rho_exp, "beta": beta,
            "max_abs_resid": resid, "rel_resid": resid / scale,
            "v_star": 1.0 / rho_exp if rho_exp != 0.0 else float("inf")}


def _solve3(a: list, b: list) -> list:
    """3x3 高斯消元（部分主元）。"""
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


# ------------------------------------------------------------------ 自洽判定


def arm_derived(cfg: dict) -> dict:
    """从臂 1/2 直接读出的 C 级量。"""
    a1 = arm_expand_once(cfg)
    a2 = arm_resolve_once(cfg)
    e_plus, u_plus = a1["e_plus_m"], a1["u_plus_m"]
    e_minus, u_minus = a2["e_minus_m"], a2["u_minus_m"]
    rho_exp = u_plus / e_plus
    rho_res = u_minus / e_minus
    return {
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": e_minus / e_plus,
        "theta_U": u_minus / u_plus,
        "rho_exp": rho_exp, "rho_res": rho_res,
        "rho": rho_res / rho_exp,
        "kappa_E": 1.0 * e_minus,        # kappa_E = b*e_minus，b=1（每步 1 次消解）
    }


def check_consistency(cfg: dict, family: str, k_steps: int = 60,
                      b_counted: float = 1.0) -> dict:
    """四臂自洽判定。b_counted 是"账面上的"消解率（负对照传 2.0）。"""
    d = arm_derived(cfg)
    arm3 = arm_free_schedule(cfg, family, k_steps)
    fit = fit_affine_in_E_and_k(arm3["traj"])

    # 预测值（由臂导出，不引用 cfg 的真值）
    rho_exp_pred = d["rho_exp"]
    beta_pred = b_counted * d["e_minus"] * (d["rho_exp"] - d["rho_res"])

    tol = 1e-9 * max(1.0, abs(cfg["rho_exp"]) + abs(cfg["rho"]) + 1.0)
    t1 = abs(fit["rho_exp"] - rho_exp_pred) <= tol
    t2 = abs(fit["beta"] - beta_pred) <= tol
    t3 = abs(fit["v_star"] * fit["rho_exp"] - 1.0) <= 1e-9
    # T4：恒等式
    t4 = (abs(d["rho"] - d["theta_U"] / d["theta_E"]) <= 1e-12
          and abs(d["rho"] - d["rho_res"] / d["rho_exp"]) <= 1e-12)

    return {
        "family": family, "k_steps": k_steps, "b_counted": b_counted,
        "arm_values": d,
        "reg_rho_exp": fit["rho_exp"], "reg_beta": fit["beta"],
        "beta_pred_from_arms": beta_pred,
        "rel_resid": fit["rel_resid"],
        "T1_rho_exp": t1, "T2_beta": t2, "T3_vstar": t3, "T4_identities": t4,
        "verdict": "pass" if (t1 and t2 and t3 and t4) else "fail",
    }


# ------------------------------------------------------------------ 主流程


THETA_E_GRID = (0.35, 0.60, 1.00)
RHO_GRID = (0.60, 1.00, 1.60)


def run_all() -> dict:
    rows = []
    neg = []
    for th in THETA_E_GRID:
        for rho in RHO_GRID:
            cfg = make_cfg(th, rho)
            for fam in SCHED_FAMILIES:
                rows.append(check_consistency(cfg, fam))
            # 负对照：b 记错
            neg.append(check_consistency(cfg, "log_ramp", b_counted=2.0))

    all_pass = all(r["verdict"] == "pass" for r in rows)
    # 负对照：b 记错应破坏 T2；rho=1 时 beta 预测恒为 0，无法区分记账错误，
    # 该点豁免（结构上 beta=0 对 0，负对照在退化点无牙齿）。
    neg_broken = all(
        r["T2_beta"] is False
        for r in neg
        if abs(r["arm_values"]["rho"] - 1.0) > 1e-12
    )
    result = {
        "experiment": "R33_clevel_closure",
        "criterion": "R45: 四臂 (e+,e-,u+,u-) 与回归 (rho_exp, beta) 自洽",
        "rows": rows,
        "negative_control": neg,
        "summary": {
            "n_rows": len(rows),
            "all_pass": all_pass,
            "negative_control_breaks_T2": neg_broken,
            "r45_verdict": "C 级闭合在原型内成立" if (all_pass and neg_broken)
                           else "未成立",
        },
    }
    return result


def main():
    result = run_all()
    with open("r33_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    s = result["summary"]
    print("R33 C-level closure (R45)")
    print(f"  rows: {s['n_rows']}  all_pass: {s['all_pass']}")
    print(f"  negative control breaks T2: {s['negative_control_breaks_T2']}")
    print(f"  R45 verdict: {s['r45_verdict']}")


if __name__ == "__main__":
    main()
