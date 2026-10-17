"""
X-S-MERA V3.4 · Relay chain（中继链）模块
=========================================

命题 2′ 的中继部分：F1 吸收态下，直接交互被禁止，交互可经中介投影系统
分段中继，每段仍须满足局部 C1（+C2）。

本模块把"每段的重合度"与"端到端传输"解耦，用于检验**复合律**：

    命题 2′ 原表述：总成功率 ≈ ∏(pᵢ·qᵢ)   （乘积律）
    本模块检验对象：串联中继链的实测传输率 T 究竟服从乘积律还是调和律

串联扩散链的物理：各段电导 cᵢ = gain·Φ·dᵢ，稳态端到端等效电导
    c_eff = 1 / Σ(1/cᵢ)          →  T = c_eff / c_max = 1 / Σ(1/dᵢ)
等长同质时 T = d/n（短板/调和律），而乘积律给 dⁿ。二者差异巨大
（d=0.9, n=3：调和 0.30 vs 乘积 0.73）。仿真将裁决哪一方成立。

另检验：C2 失效段（F2）在串联链中的角色——它是**阻断**还是**污染**。
"""

from __future__ import annotations

import numpy as np

from .le_coupling_kernel import LECouplingKernel

__all__ = ["RelayChain", "simulate_transport", "relay_bridge_demo", "min_relay_count"]


class RelayChain:
    """L — R₁ … Rₙ — E 的串联中继链。

    Parameters
    ----------
    overlaps : list[float] | None
        每段的接口重合度 pᵢ（**已测定的常量**，对应 MMS 实践中"重叠区性质
        一次性测定"）。给 None 时退化为按状态差实时计算。
    leg_c2 : list[bool] | None
        每段的接口算子兼容性；None = 全部兼容。
    with_self_dynamics : bool
        True 时端点叠加载体自演化（借用 LECouplingKernel 的 RHS），
        用于恢复实验；False 时纯输运，用于传输率测量。
    clamped_ends : bool
        True 时端点状态固定（测传输率）；False 时端点可动。
    """

    def __init__(self, n_relay: int = 1, dim: int = 3, gain: float = 6.0,
                 pot_grad: float = 0.40, d_threshold_c1: float = 0.25,
                 d_sigma: float = 0.45, overlaps: list[float] | None = None,
                 leg_c2: list[bool] | None = None,
                 ghost_gain_ratio: float = 1.60, ghost_base_ratio: float = 0.25,
                 with_self_dynamics: bool = False, clamped_ends: bool = False,
                 clip_negative: bool = True,
                 z_start: np.ndarray | None = None, z_end: np.ndarray | None = None,
                 seed: int | None = None):
        self.n_legs = n_relay + 1
        self.dim = dim
        self.gain = float(gain)
        self.pot_grad = float(pot_grad)
        self.d_threshold_c1 = float(d_threshold_c1)
        self.d_sigma = float(d_sigma)
        self.overlaps = list(overlaps) if overlaps is not None else None
        self.leg_c2 = list(leg_c2) if leg_c2 is not None else [True] * self.n_legs
        self.ghost_gain_ratio = float(ghost_gain_ratio)
        self.ghost_base_ratio = float(ghost_base_ratio)
        self.with_self_dynamics = bool(with_self_dynamics)
        self.clamped_ends = bool(clamped_ends)
        # 生物学语义下状态非负；纯输运测量时关闭截断，否则 ghost 副模态的
        # 负向偏移会被 max(·,0) 吞掉，污染指标失真。
        self.clip_negative = bool(clip_negative)

        self._rng = np.random.default_rng(seed)
        # 借用内核的载体 RHS（生命/环境自演化）
        self._kernel = LECouplingKernel(dim_life=dim, dim_env=dim, seed=seed)

        z0 = np.asarray(z_start, dtype=float) if z_start is not None else np.array([1.30, 0.20, 0.00][:dim])
        ze = np.asarray(z_end, dtype=float) if z_end is not None else np.array([0.21, 0.25, 0.00][:dim])
        self.z: list[np.ndarray] = [z0.copy()]
        for i in range(n_relay):
            w = (i + 1) / self.n_legs
            self.z.append(z0 * (1 - w) + ze * w)          # 中继节点取线性插值初值
        self.z.append(ze.copy())

        self.history = {k: [] for k in ("t", "legs_d", "legs_c2", "flux_primary",
                                        "ghost_magnitude", "zL", "zE")}

    # ---------------- 每段判定与通量 ----------------
    def _leg_flux(self, i: int) -> tuple[np.ndarray, float, bool]:
        """返回第 i 段的 (通量向量, |ghost|, 是否导通)。约定正向为 i+1 → i。"""
        delta = self.z[i + 1] - self.z[i]
        if self.overlaps is not None:
            d = float(self.overlaps[i])
        else:
            scale = 0.5 * (np.linalg.norm(self.z[i]) + np.linalg.norm(self.z[i + 1])) + 1e-9
            r = np.linalg.norm(delta) / scale
            d = float(np.exp(-(r ** 2) / (2.0 * self.d_sigma ** 2)))

        c1 = d >= self.d_threshold_c1
        c2 = bool(self.leg_c2[i])
        if not c1:
            return np.zeros(self.dim), 0.0, False

        c = self.gain * self.pot_grad * d
        J = c * delta
        ghost = 0.0
        if not c2:
            J_sec = J[1:] * self.ghost_gain_ratio - self.ghost_base_ratio * abs(c)
            J = np.concatenate([[J[0]], J_sec])
            ghost = float(np.linalg.norm(J[1:] + c * delta[1:]))
        return J, ghost, True

    def step(self, t: float, dt: float) -> dict:
        fluxes, ghosts, conducts = [], [], []
        for i in range(self.n_legs):
            J, g, ok = self._leg_flux(i)
            fluxes.append(J)
            ghosts.append(g)
            conducts.append(ok)

        # 节点更新：端点收/发，内部节点净通量
        rates = [np.zeros(self.dim) for _ in self.z]
        rates[0] += fluxes[0]
        for i in range(1, self.n_legs):
            rates[i] += fluxes[i] - fluxes[i - 1]
        rates[-1] -= fluxes[-1]

        for k, node in enumerate(self.z):
            if self.clamped_ends and k in (0, len(self.z) - 1):
                continue
            self.z[k] = node + dt * rates[k]
            if self.clip_negative:
                self.z[k] = np.maximum(self.z[k], 0.0)

        # 端点自演化（可选）
        if self.with_self_dynamics and not self.clamped_ends:
            self._kernel.L.state = self.z[0]
            self._kernel.E.state = self.z[-1]
            self.z[0] = np.maximum(self.z[0] + dt * self._kernel._life_rhs(self.z[0]), 0.0)
            self.z[-1] = np.maximum(self.z[-1] + dt * self._kernel._env_rhs(self.z[-1]), 0.0)

        d_primary = float(fluxes[0][0])
        self.history["t"].append(float(t))
        self.history["legs_d"].append([self._leg_d(i) for i in range(self.n_legs)])
        self.history["legs_c2"].append(list(self.leg_c2))
        self.history["flux_primary"].append(d_primary)
        self.history["ghost_magnitude"].append(float(np.sum(ghosts)))
        self.history["zL"].append(self.z[0].copy())
        self.history["zE"].append(self.z[-1].copy())
        return {"flux_primary": d_primary, "ghost": float(np.sum(ghosts)),
                "conducts": conducts}

    def _leg_d(self, i: int) -> float:
        if self.overlaps is not None:
            return float(self.overlaps[i])
        delta = self.z[i + 1] - self.z[i]
        scale = 0.5 * (np.linalg.norm(self.z[i]) + np.linalg.norm(self.z[i + 1])) + 1e-9
        r = np.linalg.norm(delta) / scale
        return float(np.exp(-(r ** 2) / (2.0 * self.d_sigma ** 2)))

    def d_end_to_end(self) -> float:
        return self._leg_d(0) if self.n_legs == 1 else min(
            self._leg_d(i) for i in range(self.n_legs))


def simulate_transport(overlaps: list[float], leg_c2: list[bool] | None = None,
                       gain: float = 6.0, pot_grad: float = 0.40,
                       dt: float = 0.01, max_steps: int = 100_000,
                       tol: float = 1e-11) -> dict:
    """夹持端点的纯输运实验：量测稳态端到端传输率 T。

    T = J_ss / (gain·Φ·1.0)，理想单段（d=1，间距 1）的通量为分母。
    返回 dict：T_measured, T_product(=∏dᵢ), T_harmonic(=1/Σ(1/dᵢ)),
              contamination(远端副模态偏移), ghost。
    """
    n = len(overlaps)
    dim = 3
    z_start = np.array([1.0, 0.0, 0.0])
    z_end = np.array([0.0, 0.0, 0.0])
    chain = RelayChain(n_relay=n - 1, dim=dim, gain=gain, pot_grad=pot_grad,
                       overlaps=overlaps, leg_c2=leg_c2, clamped_ends=True,
                       clip_negative=False,
                       z_start=z_start, z_end=z_end, seed=0)

    J_prev = np.inf
    J_ss = 0.0
    for step in range(max_steps):
        before = [nz.copy() for nz in chain.z]
        res = chain.step(step * dt, dt)
        J_ss = res["flux_primary"]
        dz_max = max(float(np.max(np.abs(b - a))) for b, a in zip(before, chain.z))
        if dz_max < tol:
            break
        J_prev = J_ss

    c_max = gain * pot_grad * 1.0
    T = abs(J_ss) / c_max
    T_prod = float(np.prod(overlaps))
    # C1 失败段（overlap 0）→ 整条链阻断
    if any(o < 0.25 for o in overlaps):
        T_harm = 0.0
    else:
        T_harm = 1.0 / float(np.sum([1.0 / o for o in overlaps]))

    # 污染：内部节点副模态偏离初始 0 的程度（ghost 底噪项只作用在副模态）
    interior = chain.z[1:-1]
    contamination = float(max(np.mean(np.abs(nz[1:])) for nz in interior)) if interior else 0.0
    return {"T_measured": T, "T_product": T_prod, "T_harmonic": T_harm,
            "contamination": contamination,
            "ghost": float(chain.history["ghost_magnitude"][-1]),
            "steps": step + 1}


def leg_overlap(a: np.ndarray, b: np.ndarray, d_sigma: float = 0.45) -> float:
    """两状态在公共基底上的重合度（与 LECouplingKernel.calc_d_overlap 同式）。"""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = min(a.size, b.size)
    delta = a[:n] - b[:n]
    scale = 0.5 * (np.linalg.norm(a[:n]) + np.linalg.norm(b[:n])) + 1e-9
    r = float(np.linalg.norm(delta) / scale)
    return float(np.exp(-(r ** 2) / (2.0 * d_sigma ** 2)))


def min_relay_count(z_start: np.ndarray, z_end: np.ndarray,
                    d_threshold_c1: float = 0.25, d_sigma: float = 0.45,
                    max_n: int = 16) -> int | None:
    """桥接给定端距所需的最小中继数（线性插值布点）。

    这是命题 2′/I1 的定量预测：中继数下限由**弱端尺度**决定——
    d 是尺度相对量，端距不变时，端点范数越小，所需中继越多。
    返回 None 表示 max_n 内不可桥接。
    """
    z0 = np.asarray(z_start, dtype=float)
    z1 = np.asarray(z_end, dtype=float)
    for n in range(0, max_n + 1):
        ok = True
        for i in range(n + 1):
            a = z0 * (1 - i / (n + 1)) + z1 * (i / (n + 1))
            b = z0 * (1 - (i + 1) / (n + 1)) + z1 * ((i + 1) / (n + 1))
            if leg_overlap(a, b, d_sigma) < d_threshold_c1:
                ok = False
                break
        if ok:
            return n
    return None


def relay_bridge_demo(total_t: float = 40.0, dt: float = 0.02,
                      use_relay: bool = True, n_relay: int = 3) -> dict:
    """F1 吸收态的恢复对比：直接恢复 vs 经中继链恢复。

    两端为已坍塌的 L/E 状态（直接 d≈0.03，直接通量为零）。
    中继节点取线性插值位置。注意：**中继数有下限**——d 是尺度相对量，
    弱端（E，范数小）附近每段的相对距离大，端距 1.09 时 n_min=3
    （n=1: min d=0.094；n=2: 0.232；n=3: 0.374，阈值 0.25）。
    这本身就是一个可检验预测：中继数下限由弱端尺度决定。
    """
    z_l = np.array([1.30, 0.20, 0.00])
    z_e = np.array([0.21, 0.25, 0.00])
    chain = RelayChain(n_relay=n_relay if use_relay else 0, dim=3, gain=6.0,
                       pot_grad=0.40, with_self_dynamics=True,
                       clamped_ends=False, z_start=z_l, z_end=z_e, seed=3)

    ts, d_end, gap = [], [], []
    for i in range(int(total_t / dt)):
        chain.step(i * dt, dt)
        ts.append(i * dt)
        d_end.append(min(chain._leg_d(j) for j in range(chain.n_legs)))
        gap.append(float(np.linalg.norm(chain.z[0] - chain.z[-1])))
    return {"t": np.asarray(ts), "d_end": np.asarray(d_end), "gap": np.asarray(gap)}
