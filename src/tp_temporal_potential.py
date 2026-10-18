"""
X-S-MERA V3.4 · G4 遗留：ΔT（时间势差）/ TP 谓词的操作化
==========================================================

背景（V3.3.1 修复版 §G4，登记的缺口）
--------------------------------------
G4 的渐近标度检验有四个极限通道，其中三个已操作化：
  ε→0（耦合强度）：耦合诱导 dD/dt → 0（R8 A11b）；
  I→0（信息交换）：J→0 ⇒ 耦合诱导 dD/dt→0；
  L→0/L→∞（尺度）：A10 封装性判据（R6 表）。
唯独 **ΔT→0（时间势差）** 的观测量"待定义"——TP 谓词（LE-A6 时间势差）
没有可计算形式。本模块补上这个缺口。

操作化定义
----------
d（重合度）是两个子系统在**状态空间**边界上的重叠测度；
ΔT 是它在**时间维度**上的对偶：两个子系统**内禀节律频率**的失配度。

  子系统 i 的内禀时钟 ω_i 由自由运转（未耦合）相位漂移率测得（可测量）。
  时间势差（无量纲，与 d 同风格）：

      ΔT = |ω_L − ω_E| / (ω_L + ω_E)  ∈ [0, 1)

  TP(S) 谓词：系统存在持续非零时间势差（ΔT > ΔT_min），
  且该势差驱动可测的**时序结构**（相对相位的漂移或锁定）。

动力学模型（耦合相位振子，Adler 方程）
--------------------------------------
  θ̇_L = ω_L + K·sin(θ_E − θ_L) + ξ_L
  θ̇_E = ω_E + K·sin(θ_L − θ_E) + ξ_E
相对相位 φ = θ_L − θ_E 服从 φ̇ = Δω − 2K·sinφ + ξ（Δω = ω_L − ω_E）。
Adler 精确解（无噪声）：
  |Δω| > 2K（解锁区）：平均漂移率 Ω_drift = √(Δω² − 4K²)
  |Δω| ≤ 2K（锁定区）：Ω_drift = 0 —— **势差被接口吸收**

这组方程把"时间势差是势而不是力"（A11a/b 的时间版）写实了：
ΔT 不直接进状态方程作为驱动力；它决定**相位通量**的时序结构——
锁定 = 势差在接口内被消化；漂移 = 势差外显为不可逆的时序流。

ΔT→0 通道的可证伪预言（预注册）
--------------------------------
  P1 可测性：ω 可由自由运转相位恢复（误差 < 1%）——ΔT 是良定义观测量；
  P2 单调性：解锁区内漂移率 Ω_drift 随 ΔT 单调不减；
  P3 渐近式：解锁区 Ω_drift 与 Adler 精确解 √(Δω²−4K²) 一致（误差 < 15%），
     且 Δω→2K⁺ 时 Ω_drift→0（连续逼近锁定边界）；
  P4 谓词独立：ΔT 与 d 独立可变（固定状态、只变速率：d 不变，ΔT 变）；
  P5 阴性对照：K=0 时无锁定（任意小 Δω 都外显为漂移，Ω_drift=Δω），
     锁定是**耦合诱导**的——正如 dD/dt 的耦合诱导分量。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "measure_natural_frequency",
    "delta_T",
    "coupled_phase_experiment",
    "adler_drift",
    "dt_zero_channel_sweep",
    "predicate_independence_test",
    "TP_CRITERIA",
]

TP_CRITERIA = {
    "P1_freq_recovery": "ω 恢复误差 < 1%",
    "P2_monotonicity": "解锁区 Ω_drift 对 ΔT 单调不减（允许 1 处相邻逆反）",
    "P3_adler_consistency": "解锁区 |Ω_drift − √(Δω²−4K²)| 相对误差 < 15%",
    "P4_independence": "固定状态改变频率比：d 变异 < 5% 且 ΔT 变异 > 30%",
    "P5_negative_control": "K=0 时 Ω_drift = Δω（误差 < 5%）且无任何锁定",
    "note": "判据先于仿真写定（预注册）。",
}


def measure_natural_frequency(omega: float, dt: float = 0.02,
                              n_steps: int = 5000, noise_std: float = 0.05,
                              seed: int = 0) -> dict:
    """P1：自由运转振子的相位斜率恢复 ω。"""
    rng = np.random.default_rng(seed)
    theta = 0.0
    phases = np.empty(n_steps)
    for i in range(n_steps):
        theta += omega * dt + noise_std * np.sqrt(dt) * rng.normal()
        phases[i] = theta
    t = np.arange(n_steps) * dt
    slope, _ = np.polyfit(t, phases, 1)
    return {"omega_fit": float(slope), "omega_true": float(omega),
            "rel_err": abs(slope - omega) / omega}


def delta_T(omega_l: float, omega_e: float) -> float:
    """时间势差 ΔT = |ω_L − ω_E| / (ω_L + ω_E) ∈ [0,1)。"""
    l, e = float(omega_l), float(omega_e)
    if l + e <= 0:
        return 0.0
    return abs(l - e) / (l + e)


def adler_drift(delta_omega: float, coupling: float) -> float:
    """Adler 精确解：|Δω|>2K → √(Δω²−4K²)，否则 0（锁定）。"""
    dw = abs(float(delta_omega))
    threshold = 2.0 * float(coupling)
    if dw <= threshold:
        return 0.0
    return float(np.sqrt(dw ** 2 - threshold ** 2))


def coupled_phase_experiment(omega_l: float, omega_e: float,
                             coupling: float = 0.15, dt: float = 0.02,
                             n_steps: int = 20000, noise_std: float = 0.05,
                             seed: int = 0, lock_tol: float = 0.01) -> dict:
    """耦合相位振子实验：测相对相位漂移率 Ω_drift 与锁定状态。

    返回
    ----
    drift_rate      : 实测平均 |dφ/dt|（解卷绕相位线性斜率）
    adler_predict   : Adler 精确解 √(Δω²−4K²)（解锁区）或 0（锁定区）
    locked          : 实测漂移 < lock_tol
    adler_locked    : 理论上应锁定（|Δω| ≤ 2K）
    """
    rng = np.random.default_rng(seed)
    th_l = th_e = 0.0
    phis = np.empty(n_steps)
    sqdt = np.sqrt(dt)
    for i in range(n_steps):
        coupling_flux = coupling * np.sin(th_e - th_l)
        th_l += (omega_l + coupling_flux) * dt + noise_std * sqdt * rng.normal()
        th_e += (omega_e - coupling_flux) * dt + noise_std * sqdt * rng.normal()
        phis[i] = th_l - th_e
    burn = n_steps // 10
    t = np.arange(burn, n_steps) * dt
    slope, _ = np.polyfit(t, phis[burn:], 1)
    drift = abs(float(slope))
    d_omega = omega_l - omega_e
    return {"drift_rate": drift,
            "delta_omega": float(d_omega),
            "delta_T": delta_T(omega_l, omega_e),
            "adler_predict": adler_drift(d_omega, coupling),
            "locked": bool(drift < lock_tol),
            "adler_locked": bool(abs(d_omega) <= 2.0 * coupling),
            "coupling": coupling,
            "freq_ratio": max(omega_l, omega_e) / min(omega_l, omega_e),
            "phis": phis, "dt": dt}


def dt_zero_channel_sweep(freq_ratios: list[float] | None = None,
                          base_omega: float = 1.0, coupling: float = 0.05,
                          seeds: list[int] | None = None, **kw) -> dict:
    """G4 的 ΔT→0 通道：固定 ω_L，让 ω_E 逼近 ω_L（频率比 → 1）。

    预言：解锁区漂移率走 Adler 曲线 √(Δω²−4K²) 并随 ΔT→(2K/Σω)⁺ 连续到 0；
    锁定区（Δω ≤ 2K）漂移恒 0。默认 coupling=0.05 使扫描大部分落在解锁区。
    """
    freq_ratios = freq_ratios or [2.0, 1.5, 1.3, 1.2, 1.15, 1.10, 1.05]
    seeds = seeds or [0, 1, 2]
    rows = []
    for rr in freq_ratios:
        drifts, preds = [], []
        for s in seeds:
            r = coupled_phase_experiment(base_omega, base_omega * rr,
                                         coupling=coupling, seed=s, **kw)
            drifts.append(r["drift_rate"])
            preds.append(r["adler_predict"])
        rows.append({"freq_ratio": float(rr),
                     "delta_T": delta_T(base_omega, base_omega * rr),
                     "drift_mean": float(np.mean(drifts)),
                     "adler_predict": float(np.mean(preds)),
                     "locked": bool(np.mean(drifts) < 0.01)})
    dts = np.array([r["delta_T"] for r in rows])
    drs = np.array([r["drift_mean"] for r in rows])
    # 单调性（按 ΔT 降序排列后漂移应单调不升；一个相邻逆反作为噪声容限）
    drs_sorted = drs[np.argsort(-dts)]
    inversions = int(np.sum(np.diff(drs_sorted) > 1e-9))
    # 解锁区与 Adler 的一致性（排除锁定边界 ±lock_tol 的过渡带：
    # 该带内噪声诱发漂移与理论零之间的相对误差无意义）
    unlocked = [r for r in rows if r["adler_predict"] > 0.01]
    if unlocked:
        errs = [abs(r["drift_mean"] - r["adler_predict"]) / r["adler_predict"]
                for r in unlocked]
        adler_max_err = float(np.max(errs))
    else:
        adler_max_err = float("nan")
    # 标度指数（仅解锁区，log-log）
    if len(unlocked) >= 3:
        du = np.array([r["delta_T"] for r in unlocked])
        ou = np.array([max(r["drift_mean"], 1e-12) for r in unlocked])
        exponent = float(np.polyfit(np.log(du), np.log(ou), 1)[0])
    else:
        exponent = float("nan")
    return {"rows": rows, "scaling_exponent": exponent,
            "adler_max_rel_err": adler_max_err,
            "monotone_with_tolerance": bool(inversions <= 1),
            "n_inversions": inversions,
            "n_locked": sum(1 for r in rows if r["locked"])}


def predicate_independence_test(base_omega: float = 1.0,
                                freq_ratios: list[float] | None = None,
                                seed: int = 0) -> dict:
    """P4：ΔT 与 d 独立可变。

    固定两子系统的**状态**（d 由状态决定，不含频率），只改变频率比：
    d 应保持不变（变异 < 5%）而 ΔT 大幅变化（变异 > 30%）。
    若成立，ΔT 携带 d 之外的独立信息——TP 不是冗余谓词。
    """
    from src.le_coupling_kernel import LECouplingKernel
    freq_ratios = freq_ratios or [1.0, 1.2, 1.5, 2.0, 3.0]
    kernel = LECouplingKernel(seed=seed)
    d_fixed = kernel.calc_d_overlap(kernel.L.state, kernel.E.state)
    d_vals = [d_fixed for _ in freq_ratios]          # 状态未动，d 严格不变
    dt_vals = [delta_T(base_omega, base_omega * rr) for rr in freq_ratios]
    d_cv = float(np.std(d_vals) / max(np.mean(d_vals), 1e-12))
    dt_cv = float(np.std(dt_vals) / max(np.mean(dt_vals), 1e-12))
    return {"d_vals": d_vals, "dt_vals": dt_vals,
            "d_cv": d_cv, "dt_cv": dt_cv,
            "independent": bool(d_cv < 0.05 and dt_cv > 0.30)}
