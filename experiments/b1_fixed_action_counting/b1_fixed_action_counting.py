#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B1 —— 固定动作计数注入（R33 / 命题 6.3）：方法 + 脚本准备

对应 18-v3.3-formal-framework-charter.md §6.2 B1（缺口 B-B45 / R45，闭合 C 级）。
理论锚点：docs/xs-framework/15-v3.2-gauge-fixed-framework.md §7.3（命题 6.3 首选方案）。

**命题 6.3（固定动作计数）**：
  固定展开 1 次（a=1, b=0）→ 测 ΔE = e₊、ΔU = u₊；
  固定消解 1 次（b=1, a=0）→ 测 ΔE = −e₋、ΔU = −u₋。
  ⟹ 四个单位系数 (e₊, e₋, u₊, u₋) 全部可测 ⟹ 一次闭合：
      θ_E = e₋/e₊、θ_U = u₋/u₊、ρ_exp = u₊/e₊、ρ_res = u₋/e₋、κ_E = e₋。

**R33 判据**：由四臂测得的 (e₊,e₋,u₊,u₋) 与轨迹回归所得 (ρ_exp, β₀) 是否自洽。

**性质（诚实声明）**：
  1. 本脚本是**方法 / 分析流程准备**：在**合成 A_D 数据**上验证「四臂测量 → 单位系数 →
     闭合 → 与回归自洽」的整个分析流程正确（R33 判据原文）。**不是真人实验**，不对真实认知
     过程断言任何系数。
  2. 真实实验需受试者做固定动作计数（真人测量），属 B 组 B1 的外部数据部分；本脚本在其
     数据到达前先冻结可执行的分析流程与判据。
  3. 「C 级闭合（R45）真正达成」**只在真人实验数据上才算**，本脚本不声称。

产出：对若干真实参数配置，合成 A_D 轨迹 + 四臂注入 → 闭合系数 → 与回归自洽检验。
结果写入 b1_results.json。仅标准库；无随机（确定性合成），机械可复核。
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

DT = 0.01
K = 3000      # 主轨迹步数（r 驱动）
ARMS_STEP = 1.0  # 四臂"固定动作计数"：单次动作的欧拉推进（a=1 或 b=1，步长 1 步）


# ---------- 真实参数配置（借 R27 CONFIGS，见 tests/test_r27_biphasic.py） ----------
# 每项 = (e_plus, e_minus, u_plus, u_minus)；A_D 要求 u_plus = rho_exp * e_plus 自动满足
CONFIGS = {
    # theta_E=1, theta_U=3, rho_res=3 > rho_exp=1   （钻研型）
    "drill":  {"e_plus": 1.0, "e_minus": 1.0, "u_plus": 1.0, "u_minus": 3.0},
    # theta_E=3, theta_U=1/3, rho_res=1/9 < rho_exp=1 （收敛型）
    "converge": {"e_plus": 1.0, "e_minus": 3.0, "u_plus": 1.0, "u_minus": 1.0 / 3.0},
    # theta_E=1, theta_U=1, rho_res=1 == rho_exp=1   （边界退化）
    "boundary": {"e_plus": 1.0, "e_minus": 1.0, "u_plus": 1.0, "u_minus": 1.0},
}


def r_rampdown(k, K, r_hi=9.0, r_lo=0.02, pivot=0.6):
    """分段投入强度（同 a13_first_case）：前期高投入（展开/钻研），后期低投入（消解/收敛）。"""
    t = k / float(K - 1)
    r_mid = math.sqrt(r_hi * r_lo)
    if t <= pivot:
        return math.exp(math.log(r_hi) + (math.log(r_mid) - math.log(r_hi)) * t / pivot)
    return math.exp(math.log(r_mid) + (math.log(r_lo) - math.log(r_mid)) * (t - pivot) / (1 - pivot))


# ---------- A_D 动力学（展开 a、消解 b，单位系数 e_±/u_±） ----------
def step_delta(cfg, a, b):
    """一步动作 (a, b) 产生的 (ΔE, ΔU)。"""
    dE = DT * (a * cfg["e_plus"] - b * cfg["e_minus"])
    dU = DT * (a * cfg["u_plus"] - b * cfg["u_minus"])
    return dE, dU


def main_trajectory(cfg, k=K, dt=DT):
    """r 驱动主轨迹（b=1 恒定，a=r(k)）。"""
    Es, Us = [0.0], [0.5]   # U0 远离 0（B-B43 辖域）
    for i in range(1, k):
        a = r_rampdown(i, k)
        dE, dU = step_delta(cfg, a, 1.0)
        Es.append(Es[-1] + dE)
        Us.append(Us[-1] + dU)
    ks = [i * dt for i in range(k)]
    return ks, Es, Us


# ---------- 统计工具（仅标准库） ----------
def multi_linreg(cols, y):
    """多重线性回归 y = b0 + Σ bj·colj，正规方程。cols: list of lists（同长）。"""
    n = len(y)
    m = len(cols) + 1
    X = [[1.0] + [c[i] for c in cols] for i in range(n)]   # 设计矩阵（含截距列）
    # X^T X 与 X^T y
    XtX = [[0.0] * m for _ in range(m)]
    Xty = [0.0] * m
    for i in range(n):
        xi = X[i]
        for r in range(m):
            Xty[r] += xi[r] * y[i]
            for c in range(m):
                XtX[r][c] += xi[r] * xi[c]
    # 高斯消元解 XtX · b = Xty
    aug = [row[:] + [Xty[r]] for r, row in enumerate(XtX)]
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(aug[r][col]))
        aug[col], aug[piv] = aug[piv], aug[col]
        if abs(aug[col][col]) < 1e-15:
            raise ValueError("singular design")
        for r in range(m):
            if r != col:
                fac = aug[r][col] / aug[col][col]
                for c in range(col, m + 1):
                    aug[r][c] -= fac * aug[col][c]
    b = [aug[r][m] / aug[r][r] for r in range(m)]
    return b  # [b0, b_E, b_k]


# ---------- 四臂注入（固定动作计数） ----------
def four_arm_measure(cfg, step=ARMS_STEP, dt=DT):
    """模拟四臂：臂1 纯展开 1 次，臂2 纯消解 1 次。

    用与主轨迹相同的欧拉规则（步长=step 步）测量单位系数。
    返回 (e_plus, e_minus, u_plus, u_minus)，按 A_D 定义（ΔE/ΔU 除以 step·dt 得单位速率）。
    """
    dE_exp, dU_exp = step_delta(cfg, 1.0, 0.0)   # a=1, b=0 → ΔE=e₊, ΔU=u₊
    dE_dis, dU_dis = step_delta(cfg, 0.0, 1.0)   # a=0, b=1 → ΔE=−e₋, ΔU=−u₋
    scale = step * dt
    e_plus = dE_exp / scale
    u_plus = dU_exp / scale
    e_minus = -dE_dis / scale
    u_minus = -dU_dis / scale
    return {"e_plus": e_plus, "e_minus": e_minus, "u_plus": u_plus, "u_minus": u_minus}


# ---------- 主流程 ----------
def main():
    results = {"meta": {
        "task": "B1/R33 固定动作计数注入（方法准备，非真人实验）",
        "model": "A_D（15-v3.2 §7.3 命题 6.3）",
        "dt": DT, "K": K,
        "engine": "仅标准库；确定性合成（无随机）",
    }, "configs": {}}

    for name, cfg_true in CONFIGS.items():
        # 1) 四臂测量
        arms = four_arm_measure(cfg_true)
        ep, em, up, um = arms["e_plus"], arms["e_minus"], arms["u_plus"], arms["u_minus"]
        # 2) 闭合（命题 6.3）
        theta_E = em / ep
        theta_U = um / up
        rho_exp_arm = up / ep
        rho_res_arm = um / em
        kappa_E = em                       # 惯例：κ_E = e_minus（a13 make_cfg）
        beta0_arm = rho_exp_arm * em - um  # A_D：β₀ = ρ_exp·e₋ − u₋

        # 3) 主轨迹回归 → (ρ_exp, β₀)
        ks, Es, Us = main_trajectory(cfg_true)
        b = multi_linreg([Es, ks], Us)     # U = b0 + b_E·E + b_k·k
        rho_exp_reg = b[1]
        beta0_reg = b[2]                   # β₀（每步，按 dt 归一）
        # 回归 β₀ 是"每 dt 步"的增量，与 arm 的 β₀（单位速率·dt?）对齐
        beta0_reg_scale = beta0_reg * DT   # 若 b_k 是 U 对 k(=t) 的斜率

        # 自洽检验：ρ_exp 用相对差（非零）；β₀ 用绝对差（β₀ 真值可为 0，相对差会因除零失效）
        rel_r = abs(rho_exp_arm - rho_exp_reg) / max(abs(rho_exp_reg), 1e-15)
        abs_b = abs(beta0_arm - beta0_reg)
        rel_b = abs_b / max(abs(beta0_reg), 1e-15)   # 仅报告，不用于判据

        rec = {
            "config": name,
            "true": {"e_plus": ep, "e_minus": em, "u_plus": up, "u_minus": um},
            "arms": arms,
            "closed": {
                "theta_E": theta_E, "theta_U": theta_U,
                "rho_exp_arm": rho_exp_arm, "rho_res_arm": rho_res_arm,
                "kappa_E": kappa_E, "beta0_arm": beta0_arm,
            },
            "regression": {"rho_exp_reg": rho_exp_reg, "beta0_reg": beta0_reg},
            "self_consistency": {
                "rel_rho_exp": rel_r, "abs_beta0": abs_b, "rel_beta0": rel_b,
                "rho_exp_ok": rel_r < 1e-6, "beta0_ok": abs_b < 1e-6,
            },
        }
        results["configs"][name] = rec
        print(f"[{name:9s}] arms(ρ_exp={rho_exp_arm:.4f} ρ_res={rho_res_arm:.4f} "
              f"θ_E={theta_E:.4f} θ_U={theta_U:.4f})  "
              f"reg(ρ_exp={rho_exp_reg:.6f} β₀={beta0_reg:.6f})  "
              f"自洽 rel_ρ={rel_r:.2e} abs_β={abs_b:.2e}")

    out_path = os.path.join(HERE, "b1_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\nwrote: {out_path}")

    print("\n== 判据（R33：四臂 vs 回归自洽）==")
    all_ok = True
    for name, rec in results["configs"].items():
        ok = rec["self_consistency"]["rho_exp_ok"] and rec["self_consistency"]["beta0_ok"]
        all_ok = all_ok and ok
        print(f"{name:9s} rho_exp 自洽:{rec['self_consistency']['rho_exp_ok']!s} "
              f"beta0 自洽:{rec['self_consistency']['beta0_ok']!s} -> {'PASS' if ok else 'FAIL'}")
    print("ALL:", "PASS" if all_ok else "FAIL")


if __name__ == "__main__":
    main()
