"""R5: A9 判据②统计版 + 判据③重合度重述下，推论 3（演化无向择优）的资格审计.

被审对象
--------
推论 3【演化无向择优】（v3 定稿原文）：
    "无向变异提供素材，资源与环境约束完成筛选；演化仅趋向局部适配最优，
     不必然复杂化、不必然进化。"
原降级理由（A9-03）：内部干预后下游结果不具备**可预测的确定性**，路径依赖
→ 判据②（确定性读法）失败 → 降级为机制，不具 M 资格。

R5 的问题
---------
判据②改成统计版（预测**分布**而非点值）、判据③改成 d 重合度重述
（"各语境投影中 d 重合度超阈值"）之后，这个降级还成立吗？

操作化（全部预注册，先于任何运行）
----------------------------------
现实（细粒度真值）  : 个体级 Wright–Fisher（二项采样，精确）
M 层粗模型          : 教科书式扩散近似（同 s、同 N，高斯增量 WF）
d(x, y)             : 两分布的直方图交叠（概率质量重合度，∈[0,1]）
                      ——与框架的"边界重合度 d"同一构念

判据②统计版（预测） : d(粗模型预测分布, 现实分布) ≥ θ_d
判据②干预版（核心） : 干预响应 eff(s→s+Δs) = 1 − d(下游|s, 下游|s+Δs) ≥ θ_eff
                      ——内部干预必须**可测量地**改变下游，且幅度可预测
判据③d 重述版       : 参照语境标定的粗模型**冻结参数**后，
                      d(冻结预测, 新语境现实) ≥ θ_d 对所有语境成立
判据①闭合版         : 钳制外部状态（迁入池频率 p_ext），下游分布不得改变；
                      d(下游|p_ext=低, 下游|p_ext=高) ≥ θ_d 才算闭合
推论 3 的"序层无向"  : D_seq = P(优等位固定) − p0。中性/漂移语境 D_seq ≈ 0，
                      选择主导语境 D_seq 显著 > 0——"无向"本身是语境量。

阈值（预注册）: θ_d = 0.80，θ_eff = 0.20，Δs = 0.02，T = 200 代，
                p0 = 0.10，reps = 2000，直方图 60 bins。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "TH_D", "TH_EFF", "DS", "T_GENS", "P0", "REPS", "N_BINS",
    "wright_fisher", "coarse_predict", "d_overlap", "p_fix",
    "intervention_response", "closure_leak", "frozen_transfer", "qualify",
]

# ---------------------------------------------------------------- 预注册阈值
TH_D = 0.80          # 预测/迁移重合度阈值
TH_EFF = 0.20        # 干预响应阈值（效应量下限）
DS = 0.02            # 干预幅度（选择系数增量）
T_GENS = 200         # 演化代数
P0 = 0.10            # 初始频率
REPS = 2000          # 重复实现数（分布分辨率）
N_BINS = 60          # 直方图 bins


# ---------------------------------------------------------------- 仿真核
def wright_fisher(p0: float, s: float, N: int, T: int, reps: int,
                  rng: np.random.Generator,
                  m: float = 0.0, p_ext: float = 0.5) -> np.ndarray:
    """个体级 Wright–Fisher（基因型选择，适应度 1+s vs 1）。

    m > 0 时每代以比例 m 从外部池（频率 p_ext）迁入——这是判据①的
    边界泄漏源。返回 reps 条实现的最终频率。
    """
    p = np.full(int(reps), float(p0))
    for _ in range(int(T)):
        p_sel = p * (1.0 + s) / (1.0 + s * p)
        if m > 0.0:
            p_sel = (1.0 - m) * p_sel + m * float(p_ext)
        p = rng.binomial(int(N), np.clip(p_sel, 0.0, 1.0)) / float(N)
    return p


def coarse_predict(p0: float, s: float, N: int, T: int, reps: int,
                   rng: np.random.Generator,
                   N_model: int | None = None,
                   m: float = 0.0, p_ext: float = 0.5) -> np.ndarray:
    """M 层粗模型：教科书扩散近似（高斯增量的 Wright–Fisher）。

    与现实的唯一差别是采样分布（高斯 vs 二项）。N_model 允许"冻结参数"：
    用参照语境的 N 标定后拿到新语境去用（判据③的迁移测试）。
    """
    N_use = int(N if N_model is None else N_model)
    p = np.full(int(reps), float(p0))
    for _ in range(int(T)):
        p_det = p * (1.0 + s) / (1.0 + s * p)
        if m > 0.0:
            p_det = (1.0 - m) * p_det + m * float(p_ext)
        sd = np.sqrt(np.clip(p_det * (1.0 - p_det), 0.0, None) / N_use)
        p = np.clip(p_det + rng.normal(0.0, sd), 0.0, 1.0)
    return p


# ---------------------------------------------------------------- d 构念
def d_overlap(a: np.ndarray, b: np.ndarray, bins: int = N_BINS) -> float:
    """两个样本分布的直方图交叠 = 概率质量重合度 ∈ [0,1]。

    与内核 le_coupling_kernel 的 d（两状态的高斯核重合度）同一构念：
    1 = 完全重合，0 = 支撑不相交。
    """
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    edges = np.linspace(0.0, 1.0, bins + 1)
    ha, _ = np.histogram(a, bins=edges)
    hb, _ = np.histogram(b, bins=edges)
    ha = ha / max(ha.sum(), 1)
    hb = hb / max(hb.sum(), 1)
    return float(np.minimum(ha, hb).sum())


def p_fix(x: np.ndarray) -> float:
    """优等位固定频率 P(p_T > 0.5)。推论 3"序层方向"的度量基座。"""
    return float(np.mean(np.asarray(x, float) > 0.5))


# ---------------------------------------------------------------- 判据
def intervention_response(N: int, s: float = 0.02, ds: float = DS,
                          reps: int = REPS, seed: int = 0) -> dict:
    """判据②干预版：改 s 之后下游分布动了多少（效应量 = 1 − 重合度）。"""
    rng = np.random.default_rng(seed)
    x1 = wright_fisher(P0, s, N, T_GENS, reps, rng)
    x2 = wright_fisher(P0, s + ds, N, T_GENS, reps, rng)
    d = d_overlap(x1, x2)
    return {"N": N, "Ns": N * s, "overlap": d, "effect": 1.0 - d,
            "passes": bool(1.0 - d >= TH_EFF),
            "pfix1": p_fix(x1), "pfix2": p_fix(x2)}


def closure_leak(N: int, s: float = 0.02, m_values=(0.0, 0.02, 0.10),
                 p_ext_pair=(0.2, 0.8), reps: int = REPS,
                 seed: int = 1) -> list[dict]:
    """判据①闭合版：钳制外部池频率，下游分布不得改变。

    d(下游|p_ext=低, 下游|p_ext=高) ≥ θ_d ⇒ 边界闭合（外部状态不进入）。
    注意粗模型不知道 m：同时报 d(粗模型, 现实) 看"忽略边界条件"的代价。
    """
    out = []
    for m in m_values:
        rng = np.random.default_rng(seed + int(m * 1000))
        lo = wright_fisher(P0, s, N, T_GENS, reps, rng, m=m, p_ext=p_ext_pair[0])
        hi = wright_fisher(P0, s, N, T_GENS, reps, rng, m=m, p_ext=p_ext_pair[1])
        rng2 = np.random.default_rng(seed + 7)
        pred = coarse_predict(P0, s, N, T_GENS, reps, rng2)   # 粗模型无 m
        out.append({"m": m, "d_ext": d_overlap(lo, hi),
                    "closed": bool(d_overlap(lo, hi) >= TH_D),
                    "d_model": d_overlap(pred, lo)})
    return out


def frozen_transfer(s: float = 0.02, N_ref: int = 20000,
                    N_targets=(50, 200, 2000, 20000),
                    reps: int = REPS, seed: int = 2) -> list[dict]:
    """判据③d 重述版：参照语境标定、冻结参数、逐语境迁移。

    ③ 通过 ⇔ 每个语境 d(冻结预测, 现实) ≥ θ_d。
    同时报"按语境重标定后"的 d 作对照——区分"模型错了"与"参数没跟着走"。
    """
    out = []
    for N in N_targets:
        rng = np.random.default_rng(seed + N)
        real = wright_fisher(P0, s, N, T_GENS, reps, rng)
        frozen = coarse_predict(P0, s, N, T_GENS, reps,
                                np.random.default_rng(seed + 100 + N),
                                N_model=N_ref)
        refit = coarse_predict(P0, s, N, T_GENS, reps,
                               np.random.default_rng(seed + 200 + N))
        d_fr, d_rf = d_overlap(frozen, real), d_overlap(refit, real)
        out.append({"N": N, "Ns": N * s, "d_frozen": d_fr, "d_refit": d_rf,
                    "passes": bool(d_fr >= TH_D),
                    "pfix": p_fix(real)})
    return out


def qualify(rows_2: list[dict], rows_3: list[dict],
            closure: list[dict]) -> list[dict]:
    """按语境汇总三判据，给出 M 资格裁决（含资格条件字符串）。"""
    by_n = {r["N"]: r for r in rows_2}
    tr = {r["N"]: r for r in rows_3}
    failed = [c["m"] for c in closure if not c["closed"]]
    leak = min(failed) if failed else None
    verdicts = []
    for r in rows_2:
        N, t3 = r["N"], tr[r["N"]]
        c2 = r["passes"]
        c3 = t3["passes"]
        verdict = ("macro_in_context" if (c2 and c3) else
                   "insensitive" if not c2 else "context_bound")
        q = ("本语境内宏 M（资格条件需逐语境申报）" if verdict == "macro_in_context"
             else f"无效（干预失敏, Ns={r['Ns']:.0f}<≈5）" if verdict == "insensitive"
             else f"语境内宏 M，但参数不可跨语境迁移（N≈{N:,.0f}）")
        verdicts.append({"N": N, "Ns": r["Ns"], "C1": "✓", "C2": c2, "C3": c3,
                         "verdict": verdict, "qualifier": q})
    if leak is not None:
        verdicts.append({"N": None, "Ns": None, "C1": f"✗(m≥{leak:g})", "C2": None,
                         "C3": None, "verdict": "boundary_leak",
                         "qualifier": f"m≥{leak:g} 时外部状态进入下游，边界不闭合"})
    else:
        verdicts.append({"N": None, "Ns": None, "C1": "✓", "C2": None,
                         "C3": None, "verdict": "closed_at_tested_m",
                         "qualifier": "受测 m 范围内边界闭合"})
    return verdicts
