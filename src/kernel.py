"""认知架构内核 src/kernel.py（v1.6.40 落地）。

对应仓库登记缺口：README「src/kernel.py ⏳ 计划中：认知架构内核（快慢耦合
非对称势阱），当前不在仓库」；「认知架构仿真（流形投影 + 非对称势阱，
T1–T6 自证伪测试）❌ 尚不在本仓库」。

设计依据：`workbuddy/2026-10-28-09-09-48/认知架构内核-设计草案.md`（v0.1，
材料为 x-s-mera v1.6.37 四条线盘点）与 `kernel_proposal/kernel.py`（提案稿）。

本文件实现设计草案 §6「第一步只做这些」的第 1–4 条：
  1. 内核骨架：状态 (E, U) + 规范固定参数 (s, rho_exp, rho, kappa_E)
     + 仿射积分器 + A_D+ 的 g(U) 结构约束（B-B43 处置：不实现形式拟合）
  2. 内核级自证伪：T-affine（共线性残差）、T-gauge（G_lambda 规范不变性）
  3. D-SCALE 判别接口迁入：楔形分类三值判读（R34：{|beta_hat|<=z*SE} 禁区）
  4. 输出适配器：内核轨迹 → R23 判别器（analyze_trajectories.py）可消费格式

诚实边界（沿用仓库口径，逐条写死）：
  - 数值原型，非经验验证；参数未标定到真实认知数据。
  - A_D 是未证建模假设（B-B22），带单一驱动附款（S2）：dE 与 dU 必须由调用方
    保证来自同一 r 过程——内核无法自检，这是调用方契约（违反时峰序判据可
    以 41% 误差静默判反，B-B37）。
  - 快思维本体不可测（R47/Q0）：本内核不实现、不声称读出快思维本体；快/慢只
    作为时间化参数的两个尺度通道存在。
  - g(U) 只实现结构约束（g(0)=0、单调、渐近→1），形式选择是约定（R45）。
  - 所有读数条件于 A_D 与调用方声明的尺度 (w, h, eta)。
  - D-SCALE 判定读数对 eta 稳定、g 不进入判定流程（R-042）；z=1.96 为约定性带宽。

仅标准库。入库前须过三道 CI 门（xs_trace_gate + xs_formal_checks + pytest），
新记号先入 symbol_registry.json（本内核所用 E/U/rho_exp/rho/kappa_E/s/g 均已登记）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# 1a. 规范固定后的参数集（审计 S1：lambda 规范自由度在定义层已消灭）
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GaugeFixedParams:
    """内禀参数：3 个无量纲数 + 1 个时间单位。

    s       := r/theta_E   调度速率（规范固定后无量纲）
    rho_exp  展开速率（仿射斜率，可直接标定且调度无关，R30）
    rho     := rho_res/rho_exp 之比（符号=截距符号；数值需 kappa_E，B-B42）
    kappa_E := b*e_-       时间单位（C 级闭合由固定动作计数读出，R33/R45）
    """
    s: float
    rho_exp: float
    rho: float
    kappa_E: float

    def __post_init__(self):
        if self.rho_exp <= 0:
            raise ValueError("rho_exp 必须为正（展开是认知存在的前提，补篇三/四）")
        if self.kappa_E <= 0:
            raise ValueError("kappa_E 必须为正（时间单位）")


# ---------------------------------------------------------------------------
# 1b. A_D+ 的边界层函数 g(U)：只实现结构约束（B-B43 处置）
# ---------------------------------------------------------------------------

def g_saturated(U: float) -> float:
    """g(U) = U/(1+U)：g(0)=0、单调增、渐近→1 的默认约定形式。

    注意：这是**约定**，不是辨识结果。R45 已证 {linear,square,tanh,saturated}
    在现实噪声水平下不可区分 ⟹ 换形式不得改变内核的判定读数（R-042 确认 g
    不进入 D-SCALE 判定流程）。
    """
    if U <= 0.0:
        return 0.0
    return U / (1.0 + U)


def g_linear(U: float) -> float:
    """g(U) = min(U, 1)：线性区的结构约束形式（对照用）。"""
    if U <= 0.0:
        return 0.0
    return 1.0 if U >= 1.0 else U


# ---------------------------------------------------------------------------
# 1c. 内核状态与仿射积分器
# ---------------------------------------------------------------------------

@dataclass
class KernelState:
    E: float = 0.0
    U: float = 0.0
    k: int = 0
    history: List[Tuple[int, float, float]] = field(default_factory=list)
    # 可选：每步候选权重（外部观测代理，供 R23 判别器消费；内核不虚构）
    cand_weights: List[Optional[Dict[str, float]]] = field(default_factory=list)

    def snapshot(self) -> Tuple[int, float, float]:
        return (self.k, self.E, self.U)


def step(state: KernelState, p: GaugeFixedParams, dE: float, dt: float = 1.0,
         g: Callable[[float], float] = g_saturated,
         cand_weights: Optional[Dict[str, float]] = None) -> KernelState:
    """单步推进（A_D + 单一驱动附款 + A_D+ 状态依赖）。

    仿射定律（v3.2 主定理）：
        dU = rho_exp*dE + kappa_E*rho_exp*(1-rho)*dt
    A_D+ 修正（审计 S6）：U>=0 由截断点变稳定不动点，经 g(U) 作用于消解项。
    单一驱动附款（S2）：dE 与 dU 必须由调用方保证来自同一 r 过程——内核无法
    自检，这是调用方契约（见模块 docstring 与 B-B37）。
    """
    state.E += dE
    drift = p.kappa_E * p.rho_exp * (1.0 - p.rho) * dt
    state.U = max(0.0, state.U + p.rho_exp * dE + drift * g(state.U))
    state.k += 1
    state.history.append(state.snapshot())
    state.cand_weights.append(cand_weights)
    return state


# ---------------------------------------------------------------------------
# 2. 内核级自证伪（先于 T1–T6 实验臂入库的两条）
# ---------------------------------------------------------------------------

def t_affine(states: List[Tuple[float, float]], p: GaugeFixedParams,
             tol: float = 2.0e-15) -> float:
    """T-affine：共线性残差须逐样本 <= tol（现行基准 2.0e-15，R30）。

    输入 [(E_k, U_k), ...]，返回 max|U_k - rho_exp*E_k - (U_0 - rho_exp*E_0)|。
    超界 ⟹ A_D 在该轨迹上被证否（补篇九推论 6.13 是 A_D 的有效证否途径）。

    注意（冒烟实测）：step() 默认启用 A_D+ 的 g(U) 调制，残差会随轨迹非线性地
    远离 2e-15 基准——这是 g 起作用的如实表现，不是 bug。2e-15 基准只适用于
    纯 A_D（g 恒 1 或线性区、U 远离 0）的轨迹；辖域限定「共线性只在 U 远离 0
    时成立」是本质的（B-B43/R45）。
    """
    if len(states) < 2:
        return 0.0
    E0, U0 = states[0]
    c0 = U0 - p.rho_exp * E0
    return max(abs(U - p.rho_exp * E - c0) for E, U in states)


def t_gauge(traj_a: List[Tuple[float, float]],
            traj_b: List[Tuple[float, float]],
            tol: float = 0.0) -> float:
    """T-gauge：G_lambda 规范变换下 (E,U) 须逐点不变（定理 6.12）。

    输入同一过程在两个规范下的轨迹，返回逐点最大偏差。不成立 ⟹ 参数化仍有
    未消灭的冗余自由度（B-B23 回潮）。R30 基准：lambda 跨 16x 时偏差 0.00。
    """
    n = min(len(traj_a), len(traj_b))
    return max(
        max(abs(traj_a[i][0] - traj_b[i][0]), abs(traj_a[i][1] - traj_b[i][1]))
        for i in range(n)
    ) if n else 0.0


# ---------------------------------------------------------------------------
# 3. D-SCALE 判别接口（R34 迁入：楔形分类三值判读）
# ---------------------------------------------------------------------------

DEFAULT_Z = 1.96
NOISELESS_EPS = 1e-12


def fit_beta_with_se(obs: List[Tuple[float, float, float]]) -> dict:
    """多元回归 U = alpha + gamma*E + beta*k；返回 beta_hat 与 SE(beta_hat)。

    经典最小二乘，仅标准库（3x3 高斯消元 + 求逆，仿 R34）。
    SE 用经典公式 SE(beta) = sqrt(s2 * C[j,j])，s2 = RSS/(n-3)，C=(X'X)^-1。
    无噪声时 RSS=0 -> SE=0。
    输入 obs 为 [(k, E, U), ...]，至少 4 点。
    """
    n = len(obs)
    if n < 4:
        raise ValueError("need >= 4 points for D-SCALE regression")
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


def dscale_verdict(beta_hat: float, se_beta: float, z: float = DEFAULT_Z) -> str:
    """B 级判读规则（三值输出，禁区不判）。楔形分类：
    - |beta_hat| <= z*SE + eps  → vertex_band（禁区，不判；R34 OPEN-1 解）
    - beta_hat < 0              → mid_drill（钻研方向，U 随 k 降）
    - beta_hat > 0              → mid_hollow（空洞方向，U 随 k 升）
    """
    if abs(beta_hat) <= z * se_beta + NOISELESS_EPS:
        return "vertex_band"
    return "mid_drill" if beta_hat < 0.0 else "mid_hollow"


def dscale_classify(states: List[Tuple[float, float, float]],
                    z: float = DEFAULT_Z) -> dict:
    """D-SCALE 判别接口：对 (E,U) 轨迹做楔形三值分类。

    输入 states 为 [(k, E, U), ...]（内核 history 同构）。返回 {fit, verdict}。
    对应设计草案 §6.3：D-SCALE 判定程序作为内核的判别接口迁入。
    """
    fit = fit_beta_with_se(states)
    return {"fit": fit, "verdict": dscale_verdict(fit["beta_hat"], fit["se_beta"], z)}


# ---------------------------------------------------------------------------
# 4. 输出适配器：内核轨迹 → R23 判别器（analyze_trajectories.py）可消费格式
# ---------------------------------------------------------------------------

def to_r23_json(states: List[Tuple[int, float, float]],
                cand_weights: Optional[List[Optional[Dict[str, float]]]] = None,
                task_id: str = "kernel") -> dict:
    """把内核轨迹序列化为 R23 判别器输入格式。

    R23 判别器（experiments/r23_real_cognitive/analyze_trajectories.py）消费：
        {"trajectories": [{"id", "steps": [{"weights": {类:伪计数}}]}]}
    其程序 P-hat 由 weights 算 U、由活动类数算候选数 r。

    诚实声明：内核标准状态是 (E,U)，无候选权重分布。若调用方提供了每步
    cand_weights（外部观测代理），透传记录；否则 weights 为空、该轨迹对 R23
    判 U=0（不可判）。内核不虚构权重——直接从 (E,U) 出发的判别应走
    `dscale_classify`（D-SCALE 接口）而非 weights 路径。
    """
    cw = cand_weights if cand_weights is not None else [None] * len(states)
    steps = []
    for (k, E, U), w in zip(states, cw):
        steps.append({"k": k, "E": E, "U": U,
                      "weights": w if w is not None else {}})
    return {"trajectories": [{"id": task_id, "steps": steps}]}


def rho_res_verdict(us: List[float], rs: List[int], tol: float = 1e-4) -> Optional[str]:
    """复刻 R23 判别器的 ρ_res>ρ_exp 分类（推论 5.3 / R23 逆用）。

    输入内核轨迹的 U 序列与候选数 r 序列；取最长震荡期（连续 r>=2 段），
    判其 U 净变化：
      降 ⟹ "rho_res > rho_exp"；升 ⟹ "rho_res < rho_exp"；不变 ⟹ "="。
    这是内核连接 B 线 ρ_res>ρ_exp 判定的直接接口（P1/P4/S0 已验证通道）。
    无震荡期返回 None（不可判）。
    """
    # 找震荡期（r>=2 的连续段，长度>=2）
    segs = []
    run = []
    for k, rk in enumerate(rs):
        if rk >= 2:
            run.append(k)
        else:
            if len(run) >= 2:
                segs.append(run)
            run = []
    if len(run) >= 2:
        segs.append(run)
    if not segs:
        return None
    seg = max(segs, key=len)
    net = us[seg[-1]] - us[seg[0]]
    if net < -tol:
        return "rho_res > rho_exp"
    if net > tol:
        return "rho_res < rho_exp"
    return "rho_res = rho_exp"


# ---------------------------------------------------------------------------
# 数值工具（仿 R34，私有）
# ---------------------------------------------------------------------------

def _solve3(a: List[List[float]], b: List[float]) -> List[float]:
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


def _invert3(a: List[List[float]]) -> List[List[float]]:
    m = [row[:] + [1.0 if i == j else 0.0 for j in range(3)]
         for i, row in enumerate(a)]
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
