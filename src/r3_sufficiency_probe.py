"""
X-S-MERA V3.4 · R3-v2：C1+C2 的充分性补案
============================================

R3 的局限（修订登记 §12.8 第 3 条）明确写道：

  "案例选择虽按全谱覆盖原则，仍可能遗漏'C1 过 C2 过但历史失败'型案例
   （若有，将直接证伪 C1+C2 的充分性——目前案例库只检验了必要性方向）。"

本模块构造这类反例。方法与 R 系列一致：**先预注册判据，再构造全谱，再判定**。

为什么反例必须从 C1/C2 的**定义缝隙**里长出来
----------------------------------------------
C1 = d（边界重合度）超阈值；C2 = 接口算子 S-兼容。
二者都是**单次交互**的判据。而"历史成功/失败"是**序列**属性。
因此充分性的证伪窗口只有一处：**单次合规但序列失败**。
这不是逻辑漏洞，而是条款作用域的隐含假设——C1+C2 默认
"交互可重复、环境平稳、单次收益可累加"。反例就长在这三条假设的破口上。

构造出五类候选反例（全部满足 C1∧C2）：

  S1 累积损伤型：每次交互局部合规，但每次都在同一方向累积不可逆损伤
     （连续小剂量毒性 / 反复摩擦）。单次可逆，序列不可逆。
  S2 竞争置换型：交互本身成功，但挤出率高于收益（竞争排斥）。
  S3 时序错配型：C1/C2 在每步都成立，但"接口版本"随交互进程漂移，
     后段交互与前段不兼容（协议演化 / 病原体逃逸）。
  S4 中介污染型：直接交互合规，但经中介链时，中介自身的 F2 污染累积
     （用 R4 的 P4 分解：底噪分量归 D17，但失配分量归 F2，逐段叠加）。
  S5 势耗竭型：交互合规且成功，但势（可分化性潜能的载体）被耗尽，
     后续无势可用——这是 R8 A11c 的直接推论：Φ 由流生成，也会被流耗尽。

判定与结论
----------
本模块**不预设反例一定成立**。它给出五类构造的数值化版本，逐类测量
"单次合规率"与"序列成功率"，并输出：

  单次合规 100%、序列成功 0% → **充分性反例成立**（C1+C2 必须加序列条件）
  单次合规 100%、序列成功 >0 但远低于预期 → 弱反例（需登记为条件性）

无论结果如何，产出都是确定的：C1+C2 的作用域须显式声明为
**单次交互的准入判据**，而非**历史成败的充分条件**。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "S1_cumulative_damage",
    "S2_competitive_displacement",
    "S3_temporal_mismatch",
    "S4_mediator_contamination",
    "S5_potential_depletion",
    "run_full_sweep",
    "SUFFICIENCY_CRITERIA",
]

SUFFICIENCY_CRITERIA = {
    "strong_counterexample": "单次合规率 = 100% 且 序列成功率 = 0%",
    "weak_counterexample": "单次合规率 = 100% 且 序列成功率 < 20%",
    "no_counterexample": "单次合规率 < 100%（该构造不是合法交互）或 序列成功率 ≥ 20%",
    "note": "判据先于仿真写定（预注册），不得事后调整阈值。",
}


def _verdict(single_pass_rate: float, series_success: float) -> str:
    if single_pass_rate < 1.0 - 1e-9:
        return "构造不合法（单次判据未全过）"
    if series_success <= 0.0:
        return "强反例：充分性被证伪"
    if series_success < 0.20:
        return "弱反例：充分性不成立（条件性）"
    return "无反例"


# ----------------------------------------------------------------------
# S1 累积损伤型
# ----------------------------------------------------------------------

def S1_cumulative_damage(n_rounds: int = 40, dose: float = 0.06,
                         repair: float = 0.04, damage_threshold: float = 0.60,
                         n_trials: int = 200, seed: int = 0) -> dict:
    """每次交互都满足 C1∧C2（d 高、接口兼容），但损伤按剂量累积。

    动力学：Dmg_{t+1} = Dmg_t + dose − repair，下界 0。
    单次"成功"定义：本次交互造成的增量损伤 |Δ| < repair（即单次可逆）。
    序列"成功"定义：全过程 Dmg 未越过 damage_threshold。
    """
    rng = np.random.default_rng(seed)
    successes = 0
    single_flags = []
    peak = []
    for _ in range(n_trials):
        dmg = 0.0
        single_ok = True
        ok = True
        for _t in range(n_rounds):
            inc = dose + rng.normal(0.0, dose * 0.10)
            dmg = max(0.0, dmg + inc - repair)
            if inc > repair:          # 单次是否"可逆"
                single_ok = False
            if dmg >= damage_threshold:
                ok = False
        single_flags.append(single_ok)
        peak.append(dmg)
        if ok and single_ok:
            successes += 1
    return {"case": "S1 累积损伤型",
            "single_pass_rate": 1.0,      # C1/C2 恒过（构造保证）
            "single_reversible_rate": float(np.mean(single_flags)),
            "series_success": successes / n_trials,
            "peak_damage_mean": float(np.mean(peak)),
            "verdict": _verdict(1.0, successes / n_trials)}


# ----------------------------------------------------------------------
# S2 竞争置换型
# ----------------------------------------------------------------------

def S2_competitive_displacement(n_rounds: int = 30, benefit: float = 0.05,
                                displacement: float = 0.07,
                                n_trials: int = 200, seed: int = 1) -> dict:
    """交互成功且每次都有正收益，但同类挤出率高于收益 → 净增长为负。

    单次成功：本次收益 > 0（C1∧C2 恒过，故恒真）。
    序列成功：最终净存量 > 初始存量。
    """
    rng = np.random.default_rng(seed)
    ok = 0
    finals = []
    for _ in range(n_trials):
        stock = 1.0
        for _t in range(n_rounds):
            stock += benefit + rng.normal(0, 0.01) - displacement * (1.0 + 0.1 * rng.normal())
            stock = max(stock, 0.0)
        finals.append(stock)
        if stock > 1.0:
            ok += 1
    return {"case": "S2 竞争置换型",
            "single_pass_rate": 1.0,
            "series_success": ok / n_trials,
            "final_stock_mean": float(np.mean(finals)),
            "verdict": _verdict(1.0, ok / n_trials)}


# ----------------------------------------------------------------------
# S3 时序错配型
# ----------------------------------------------------------------------

def S3_temporal_mismatch(n_rounds: int = 40, drift_rate: float = 0.015,
                         tolerance: float = 0.25, step_tol: float = 0.04,
                         n_trials: int = 200, seed: int = 2) -> dict:
    """接口版本随交互进程漂移：**每一步的增量**都在可消化范围内，但累积越界。

    接口兼容性有两读法，本构造正是要暴露这个双读法的后果：

      读法 (i) 瞬时读法：|vL(τ) − vE(τ)| < tolerance。
               这是 C2 的字面读法（"这一刻接口是否匹配"）。
      读法 (ii) 增量读法：|ΔvL(τ) − ΔvE(τ)| < step_tol。
               这是工程实际读法（"这一次交互的接口变化能否被吃掉"）。

    构造：drift 取到让**每步增量**远小于 step_tol、但总量越过 tolerance 的区间。
      - 若 C2 按增量读法执行 → 单次判据全程通过 → 序列失败 = **强反例**；
      - 若 C2 按瞬时读法执行 → 中途就会报警 → 不构成反例。
    **所以 S3 的结论本身就是有条件的**：它证明的是"若 C2 按增量读法，
    则 C1+C2 不充分"，并给出即时纠正（C2 必须按瞬时读法+绝对基准判定）。
    这是本模块最有价值的一条：反例的价值不在"打败条款"，
    而在逼出条款读法必须显式化。
    """
    rng = np.random.default_rng(seed)
    ok = 0
    drifts = []
    inc_all_ok = True
    for _ in range(n_trials):
        vL = 0.0
        vE = 0.0
        for _t in range(n_rounds):
            dvL = drift_rate * (1 + 0.3 * rng.normal())
            dvE = drift_rate * 0.4 * (1 + 0.3 * rng.normal())
            vL += dvL
            vE += dvE
            if abs(dvL - dvE) >= step_tol:      # 增量读法下的单次判定
                inc_all_ok = False
        drifts.append(abs(vL - vE))
        if abs(vL - vE) < tolerance:
            ok += 1
    single = 1.0 if inc_all_ok else 0.0
    out = {"case": "S3 时序错配型（增量读法）",
           "single_pass_rate": single,
           "series_success": ok / n_trials,
           "final_drift_mean": float(np.mean(drifts)),
           "tolerance": tolerance,
           "step_tol": step_tol,
           "verdict": _verdict(single, ok / n_trials)}
    # 同一构造在"瞬时读法"下的对照：报警时机
    out["instant_reading_alarms"] = bool(float(np.mean(drifts)) >= tolerance)
    return out


# ----------------------------------------------------------------------
# S4 中介污染型（接 R4 的 P4 分解）
# ----------------------------------------------------------------------

def S4_mediator_contamination(n_legs: int = 3, mismatch: float = 0.16,
                              base_noise: float = 0.02, tol: float = 0.30,
                              n_trials: int = 300, seed: int = 3) -> dict:
    """每段中继都满足 C1∧C2，但中介自身 F2 的失配分量逐段累积。

    复用中继链的调和律结论（v1.1.0）：串联传输 T = 1/Σ(1/dᵢ)。
    污染按 R4 §13.4 P4 的分解：失配分量（状态相关）会**逐段叠加并放大**，
    底噪分量（D17）不放大。
    序列失败判据：总污染 > tol。
    """
    rng = np.random.default_rng(seed)
    ok = 0
    totals = []
    for _ in range(n_trials):
        total = 0.0
        for i in range(n_legs):
            # 每段的失配分量（状态相关，随段序放大：后段接前段已污染的输入）
            total += mismatch * (1.0 + 0.5 * i) * (1.0 + 0.1 * rng.normal())
            total += base_noise * rng.normal()      # 底噪（不放大）
        totals.append(total)
        if total <= tol:
            ok += 1
    return {"case": "S4 中介污染型",
            "single_pass_rate": 1.0,
            "series_success": ok / n_trials,
            "total_contamination_mean": float(np.mean(totals)),
            "tolerance": tol,
            "verdict": _verdict(1.0, ok / n_trials)}


# ----------------------------------------------------------------------
# S5 势耗竭型（R8 A11c 的直接推论）
# ----------------------------------------------------------------------

def S5_potential_depletion(n_rounds: int = 40, draw: float = 0.05,
                           regen: float = 0.02, phi_min: float = 0.20,
                           n_trials: int = 200, seed: int = 4) -> dict:
    """交互合规且每步成功，但"可分化性潜能"Φ 被流消耗（抽取率 > 再生率）。

    R8 的 A11c 说：Φ 可由流生成。它的对偶——Φ 也可被流耗尽——此前未被写出。
    动力学：Φ_{t+1} = Φ_t − draw + regen，下界 0。
    序列失败：Φ 低于 phi_min（此后新边界不可生成，交互失去意义）。
    """
    rng = np.random.default_rng(seed)
    ok = 0
    finals = []
    for _ in range(n_trials):
        phi = 1.0
        good = True
        for _t in range(n_rounds):
            phi = max(0.0, phi - draw + regen + 0.004 * rng.normal())
            if phi < phi_min:
                good = False
                break
        finals.append(phi)
        if good:
            ok += 1
    return {"case": "S5 势耗竭型",
            "single_pass_rate": 1.0,
            "series_success": ok / n_trials,
            "final_phi_mean": float(np.mean(finals)),
            "verdict": _verdict(1.0, ok / n_trials)}


# ----------------------------------------------------------------------
# 全谱扫描
# ----------------------------------------------------------------------

def run_full_sweep(n_trials: int = 200) -> dict:
    """跑五类构造，返回判定表与总结论。"""
    cases = [
        S1_cumulative_damage(n_trials=n_trials),
        S2_competitive_displacement(n_trials=n_trials),
        S3_temporal_mismatch(n_trials=n_trials),
        S4_mediator_contamination(n_trials=n_trials),
        S5_potential_depletion(n_trials=n_trials),
    ]
    n_strong = sum(1 for c in cases if c["verdict"].startswith("强反例"))
    n_weak = sum(1 for c in cases if c["verdict"].startswith("弱反例"))
    n_invalid = sum(1 for c in cases if c["verdict"].startswith("构造不合法"))
    if n_strong > 0:
        summary = ("充分性被证伪：存在单次全合规而序列必失败的可实现构造。"
                   "C1+C2 的作用域须限定为**单次交互的准入判据**。")
    elif n_weak > 0:
        summary = ("充分性不成立（条件性）：存在单次全合规而序列成功率极低的构造。"
                   "C1+C2 须补序列条件方能主张充分性。")
    else:
        summary = "未找到反例：五类构造的序列成功率均 ≥20%。"
    return {"cases": cases, "n_strong": n_strong, "n_weak": n_weak,
            "n_invalid": n_invalid, "summary": summary}
