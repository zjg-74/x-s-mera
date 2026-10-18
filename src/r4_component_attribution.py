"""
X-S-MERA V3.4 · R4 遗留：A8/D17 冲突的分量归类判据（LE 内核验证）
==================================================================

背景（修订登记 §十三 / R4 结论）
--------------------------------
A8 前置条款与 D17 载体泄漏条款的**真冲突**：交集非空（交互中的异常信号，
ghost force 本身容许两种叙事）且处置不相容（D17 滤波留偏置 → 礁石洗净锁定；
F2 误用修接口 → 真通道误判非法）。根源是按**来源叙事**划管辖，而来源不可观测。

R4 的解法：管辖按**可测属性**分配，顺序 load-bearing：
  K1 状态相关性（首要）→ K2 修复检验 → K3 滤波检验
复合现象（LE 内核 F2 = ghost_gain_ratio 失配 + ghost_base_ratio 底噪）
必须**先分解后归类**，禁止整体归类。

泄漏三域与处置：
  载体域（D17，处置=消除/滤波）：底噪分量，方向偏置型，与状态差方向无关；
  接口域（A8-F2，处置=修复）：失配分量，与状态差强相关，修复接口后消失；
  语义域（A1-Σ，处置=监测）：d∈D 是信号不是缺陷。

本模块在 LECouplingKernel 的 F2 阶段实测这三步判据。

LE 内核 F2 残差的解析结构（用于对照）
--------------------------------------
守恒残差 r = flux_E + flux_L。正常耦合（C2 过）时 r ≡ 0。
F2（C1∧¬C2）时副模态（分量 i ≥ 1）：
  r_i = −(ghost_gain_ratio − 1)·gain·δ_i  −  ghost_base_ratio·|gain|
其中 δ = E.state − L.state，gain = coupling_gain·∇Φ·d。
  第一项 = **失配分量**：∝ δ_i（状态相关、随 δ 翻号）→ 接口域；
  第二项 = **底噪分量**：∝ |gain|（方向不变偏置，与 δ 翻号无关）→ 载体域。
两项在**同一个残差向量**里叠加——这正是"复合现象先分解后归类"的测试床。

判据（预注册）
--------------
  K1 状态相关性：|corr(comp, δ)| ≥ 0.5 → 状态相关（接口域候选）；
                 方向不变（δ 翻号下分量均值同号）→ 偏置型（载体域候选）。
  K2 修复检验：置 c2_compatibility=True 重跑，状态相关分量应被消除
               （残差范数下降 ≥ 90%）。
  K3 滤波检验：偏置分量在 δ>0 与 δ<0 两个子样本中**均值同号**
               （符号不变 ⇒ 可按方向偏置滤波消除，不动动力学本身）。
"""

from __future__ import annotations

import numpy as np

from src.le_coupling_kernel import LECouplingKernel

__all__ = [
    "collect_f2_trace",
    "collect_f2_ensemble",
    "decompose_residual",
    "k1_state_correlation",
    "k2_repair_test",
    "k3_filter_test",
    "classify_le_f2",
    "ATTRIBUTION_CRITERIA",
]

ATTRIBUTION_CRITERIA = {
    "K1_corr_threshold": 0.5,
    "K2_elimination_min": 0.90,
    "K3_sign_invariance": "δ 翻号两子样本均值同号（偏置型）",
    "note": "判据先于仿真写定（预注册），顺序 K1→K2→K3 不可调换。",
}


def _kernel_gain(kernel: LECouplingKernel, d: float) -> float:
    return kernel.coupling_gain * kernel.semantic_potential_grad * d


def collect_f2_trace(kernel: LECouplingKernel, n_steps: int = 400,
                     dt: float = 0.05) -> dict:
    """在 F2 阶段采集守恒残差轨迹。

    每步先按当前状态计算通量（与 kernel.step 内部同一调用序），登记
    残差 r = flux_E + flux_L、状态差 δ、增益 gain，然后推进内核。
    """
    residuals, deltas, gains, ds, actives = [], [], [], [], []
    t = kernel.t
    for _ in range(n_steps):
        d = kernel.calc_d_overlap(kernel.L.state, kernel.E.state)
        C1, C2, region, fail = kernel.judge_interaction_conditions(d)
        flux_L, flux_E = kernel.compute_coupling_flux(d, C1, C2)
        n = min(kernel.L.state.size, kernel.E.state.size)
        residuals.append(flux_E + flux_L)
        deltas.append(kernel.E.state[:n] - kernel.L.state[:n])
        gains.append(_kernel_gain(kernel, d))
        ds.append(d)
        actives.append(bool(C1))       # C1 失败的步通量被禁，残差恒零，须屏蔽
        kernel.step(t + dt, dt)
        t += dt
    return {"residual": np.asarray(residuals), "delta": np.asarray(deltas),
            "gain": np.asarray(gains), "d": np.asarray(ds),
            "active": np.asarray(actives, dtype=bool)}


def collect_f2_ensemble(n_runs: int = 8, n_steps: int = 120, dt: float = 0.05,
                        seed: int = 0, dim: int = 3, **kernel_kw) -> dict:
    """多初值集成采集：打破 gain·δ 与 |gain| 的共线性。

    单条轨迹中 δ 往往长时间保持同号且 gain 缓变，导致失配分量与底噪分量
    在回归中不可辨识（共线性）。用随机初值的多条短轨迹池化样本，
    使 δ 的符号与幅值充分变化——这是"先分解后归类"的可辨识性前提。
    """
    rng = np.random.default_rng(seed)
    parts = []
    for _ in range(n_runs):
        kw = dict(c2_compatibility=False, seed=int(rng.integers(0, 2 ** 31)))
        kw.update(kernel_kw)
        k = LECouplingKernel(dim_life=dim, dim_env=dim, **kw)
        # 随机化初值（保留分量数，范围覆盖两侧高低组合）
        k.L.state = rng.uniform(0.15, 1.20, size=dim)
        k.E.state = rng.uniform(0.15, 1.20, size=dim)
        parts.append(collect_f2_trace(k, n_steps=n_steps, dt=dt))
    out = {key: np.concatenate([p[key] for p in parts], axis=0)
           for key in ("residual", "delta", "gain", "d", "active")}
    return out


def decompose_residual(trace: dict) -> dict:
    """把残差按模态分解为失配分量（∝ gain·δ）与底噪分量（∝ |gain|）。

    对分量 i ≥ 1 做线性回归 r_i = α_i·(gain·δ_i) + β_i·|gain|。
    解析真值：α = −(ghost_gain_ratio − 1)，β = −ghost_base_ratio。
    主模态 i=0 在 F2 下守恒（r_0 ≡ 0），不参与分解。
    C1 失败的步（F1 禁通量，残差恒零而回归元非零）先行屏蔽，
    否则零残差点会污染分解（可辨识性前提）。
    """
    mask = trace["active"]
    r = trace["residual"][mask]
    delta = trace["delta"][mask]
    gain = trace["gain"][mask]
    gain_abs = np.abs(gain)
    n_modes = r.shape[1]
    out = {"n_modes": n_modes, "mask": mask, "components": []}
    for i in range(1, n_modes):
        driver = gain * delta[:, i]               # gain·δ_i（带号）
        X = np.column_stack([driver, gain_abs])
        coef, *_ = np.linalg.lstsq(X, r[:, i], rcond=None)
        pred = X @ coef
        resid = r[:, i] - pred
        ss_res = float(np.sum(resid ** 2))
        ss_tot = float(np.sum((r[:, i] - r[:, i].mean()) ** 2)) + 1e-15
        mismatch_part = coef[0] * driver          # 失配分量的重构
        base_part = coef[1] * gain_abs            # 底噪分量的重构
        out["components"].append({
            "mode": i,
            "alpha_mismatch": float(coef[0]),
            "beta_base": float(coef[1]),
            "r_squared": 1.0 - ss_res / ss_tot,
            "mismatch_part": mismatch_part,
            "base_part": base_part,
            "delta_i": delta[:, i],
        })
    return out


def k1_state_correlation(component_part: np.ndarray, delta_i: np.ndarray) -> dict:
    """K1：分量与状态差的相关性。返回是否状态相关（≥0.5 为是）。"""
    if np.std(component_part) < 1e-12 or np.std(delta_i) < 1e-12:
        rho = 0.0
    else:
        rho = float(np.corrcoef(component_part, delta_i)[0, 1])
    return {"rho": rho, "state_correlated": bool(abs(rho) >=
            ATTRIBUTION_CRITERIA["K1_corr_threshold"])}


def k2_repair_test(dim: int = 3, n_steps: int = 400, dt: float = 0.05,
                   seed: int = 0, **kernel_kw) -> dict:
    """K2：修复检验。同一参数分别跑 F2 与已修复（C2 过）两个内核。

    接口域分量的定义性特征：修复接口后该分量消失。
    返回残差范数均值（修复前/后）与消除率。
    """
    kw = dict(c2_compatibility=False, seed=seed)
    kw.update(kernel_kw)
    k_f2 = LECouplingKernel(dim_life=dim, dim_env=dim, **kw)
    tr_f2 = collect_f2_trace(k_f2, n_steps=n_steps, dt=dt)

    kw["c2_compatibility"] = True
    k_ok = LECouplingKernel(dim_life=dim, dim_env=dim, **kw)
    tr_ok = collect_f2_trace(k_ok, n_steps=n_steps, dt=dt)

    norm_f2 = float(np.linalg.norm(tr_f2["residual"], axis=1).mean())
    norm_ok = float(np.linalg.norm(tr_ok["residual"], axis=1).mean())
    elim = 1.0 - norm_ok / max(norm_f2, 1e-15)
    return {"residual_norm_f2": norm_f2, "residual_norm_repaired": norm_ok,
            "elimination_ratio": float(elim),
            "repair_eliminates": bool(elim >= ATTRIBUTION_CRITERIA["K2_elimination_min"])}


def k3_filter_test(base_part: np.ndarray, delta_i: np.ndarray) -> dict:
    """K3：滤波检验（偏置型判定）。

    底噪/载体分量的操作特征：**与状态差方向无关的方向偏置**。
    把样本按 sign(δ_i) 分两组，若分量均值同号 → 符号不变 → 偏置型，
    可按方向偏置滤波消除（D17 处置=消除），不需要修接口。
    """
    pos = base_part[delta_i > 0]
    neg = base_part[delta_i < 0]
    if pos.size < 5 or neg.size < 5:
        return {"sign_invariant": False, "reason": "insufficient_split",
                "mean_pos": float("nan"), "mean_neg": float("nan")}
    m_pos, m_neg = float(pos.mean()), float(neg.mean())
    invariant = bool(np.sign(m_pos) == np.sign(m_neg))
    return {"sign_invariant": invariant, "mean_pos": m_pos, "mean_neg": m_neg,
            "filterable_as_bias": invariant}


def classify_le_f2(n_steps: int = 400, dt: float = 0.05, seed: int = 0,
                   dim: int = 3, **kernel_kw) -> dict:
    """完整管线：F2 采集 → 分解 → K1/K2/K3 → 三域归类。

    返回每个副模态的归类判定与整体结论。处置映射：
      状态相关 ∧ 修复可消除 → 接口域（A8-F2，处置=修复）；
      方向不变偏置 ∧ 可滤波  → 载体域（D17，处置=消除）。
    """
    probe = LECouplingKernel(dim_life=dim, dim_env=dim,
                             **dict(c2_compatibility=False, seed=seed, **kernel_kw))
    trace = collect_f2_ensemble(n_runs=8, n_steps=n_steps, dt=dt, seed=seed,
                                dim=dim, **kernel_kw)
    dec = decompose_residual(trace)
    repair = k2_repair_test(dim=dim, n_steps=n_steps, dt=dt, seed=seed, **kernel_kw)

    true_alpha = -(probe.ghost_gain_ratio - 1.0)
    true_beta = -probe.ghost_base_ratio

    rows = []
    for comp in dec["components"]:
        k1_mm = k1_state_correlation(comp["mismatch_part"], comp["delta_i"])
        k3_base = k3_filter_test(comp["base_part"], comp["delta_i"])
        rows.append({
            "mode": comp["mode"],
            "alpha_mismatch": comp["alpha_mismatch"],
            "beta_base": comp["beta_base"],
            "r_squared": comp["r_squared"],
            "mismatch_state_correlated": k1_mm["state_correlated"],
            "mismatch_rho": k1_mm["rho"],
            "mismatch_domain": "接口域（A8-F2，处置=修复）"
                if k1_mm["state_correlated"] and repair["repair_eliminates"]
                else "待定",
            "base_sign_invariant": k3_base["sign_invariant"],
            "base_domain": "载体域（D17，处置=消除）"
                if k3_base["sign_invariant"] else "待定",
        })
    return {"rows": rows, "repair": repair, "trace": trace,
            "decomposition": dec,
            "true_alpha": float(true_alpha), "true_beta": float(true_beta),
            "summary": ("复合异常已分解：失配分量归接口域（修复），"
                        "底噪分量归载体域（消除）；禁止整体归类。")}
