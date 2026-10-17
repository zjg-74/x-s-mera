"""
X-S-MERA V3.4 · LE (Life-Environment) coupling kernel
=====================================================

生命-环境同尺度耦合动力学内核。本文件是 experimentation 的**唯一真值源**，
experiments/ 与 tests/ 均只依赖此处的公开接口。

理论锚点（与《修订登记 v0》对齐）
----------------------------------
1. **d 是边界重合度标量**（无量纲量，非时空尺度差）：
   d = exp(-r^2 / (2 sigma^2))，r 为两子系统在公共基底上的相对距离。
   尺度差只是张力的"发生器"之一，d 是它在系统边界的外显测度。
2. **D 是 d 的统计分布集合**（不是一个数）。每步在边界采一个 d 样本，
   D 的封装性（A10 边界封装条款）作为 M-完备性指标：
   可封装 = 单峰、低方差、近对称 → 封装为"边界宏"，不计入 M-不完备；
   D 偏离可封装形态（多峰/重尾/方差发散）= M-不完备的真信号。
3. **势-流分离（A11a/A11b/A11-R）**：`semantic_potential_grad`（语义势）只调制
   耦合增益，**不直接驱动状态更新**；状态演化完全由载体通量 flux 驱动。
4. **交互双条件（命题 2'）**：
   C1 = d >= d_threshold_c1（存在满足粒度约束的公共细化构造 I1）
   C2 = 接口算子 S-兼容（c2_compatibility）
   三区：high(C1&C2) / middle(C1&!C2) / low(!C1)
   失效：F1 = 无重合、直接交互被禁止，通量置零；
        F2 = 有重合但接口不兼容，产生 ghost 伪力（礁石伪收敛）。
"""

from __future__ import annotations

import numpy as np

__all__ = ["LEProjectionSystem", "LECouplingKernel", "KERNEL_VERSION"]

KERNEL_VERSION = "3.4-le.1"


class LEProjectionSystem:
    """单个投影子系统（L 或 E），持有状态向量与边界测度样本集 D。"""

    def __init__(self, dim_state: int = 3, name: str = "LE", noise_std: float = 0.01,
                 var_max: float = 0.10, min_samples: int = 16):
        self.dim_state = int(dim_state)
        self.name = name
        self.noise_std = float(noise_std)
        self.var_max = float(var_max)
        self.min_samples = int(min_samples)

        self.state = np.zeros(self.dim_state, dtype=float)
        self.d_samples: list[float] = []

    # ---------------- 样本采集 ----------------
    def record_d_sample(self, sample: float) -> None:
        self.d_samples.append(float(sample))

    def clear_samples(self) -> None:
        self.d_samples = []

    @property
    def D(self) -> np.ndarray:
        return np.asarray(self.d_samples, dtype=float)

    # ---------------- 分布诊断 ----------------
    def _stats(self, window: int | None = None) -> dict | None:
        a = self.D
        if window is not None and a.size > window:
            a = a[-int(window):]
        n = a.size
        if n < 3:
            return None
        mean = float(a.mean())
        var = float(a.var(ddof=0))
        sd = float(np.sqrt(var)) + 1e-12
        z = (a - mean) / sd
        skew = float(np.mean(z ** 3))
        kurt = float(np.mean(z ** 4))                      # Pearson kurtosis，正态 = 3
        # Sarle 双峰系数：BC = (skew^2 + 1) / kurtosis，> 5/9 提示非单峰
        bc = float((skew ** 2 + 1.0) / max(kurt, 1e-12))
        return {"n": n, "mean": mean, "var": var, "sd": sd,
                "skew": skew, "kurtosis": kurt, "bc": bc}

    def test_D_encapsulation(self, bc_max: float = 0.60, skew_max: float = 1.0,
                             window: int | None = None) -> dict:
        """A10 边界封装条款的检验入口。

        可封装 = 样本量足够 且 非多峰(BC < bc_max) 且 近对称(|skew| < skew_max)
                 且 方差有界(var < var_max)。

        window：仅用最近 window 个样本（滚动窗口），用于跟踪**当前尺度的**
        D 形态；window=None 时用全部样本。
        """
        st = self._stats(window=window)
        if st is None or st["n"] < self.min_samples:
            # 样本不足：不判失败，也不判通过（避免启动阶段的伪影被计入 M-不完备）
            return {"is_encapsulated": False, "status": "warmup",
                    "reason": "insufficient_samples",
                    "n": st["n"] if st else len(self.d_samples),
                    "mean": float("nan"), "var": float("nan"),
                    "sd": float("nan"), "skew": float("nan"),
                    "kurtosis": float("nan"), "bc": float("nan")}

        reasons = []
        if st["bc"] >= bc_max:
            reasons.append("multimodal")
        if abs(st["skew"]) > skew_max:
            reasons.append("skewed")
        if st["var"] > self.var_max:
            reasons.append("variance_unbounded")

        out = {"is_encapsulated": len(reasons) == 0,
               "status": "fail" if reasons else "ok",
               "reason": "+".join(reasons) if reasons else "ok"}
        out.update(st)
        return out


class LECouplingKernel:
    """L-E 双投影子系统耦合内核。

    参数
    ----
    dim_life, dim_env : 两侧状态维度（公共基底取 min 维）
    d_threshold_c1    : C1 判定阈值
    c2_compatibility  : 接口算子是否 S-兼容
    semantic_potential_grad : 语义势梯度（只调制耦合增益，不直接驱动状态）
    coupling_gain     : 载体通量的基础增益
    ghost_gain_ratio  : F2 伪力幅值失配比（>1 即非守恒，注入净能量）
    seed              : 随机种子（可复现性）
    """

    def __init__(self, dim_life: int = 3, dim_env: int = 3,
                 d_threshold_c1: float = 0.25, c2_compatibility: bool = True,
                 semantic_potential_grad: float = 0.40,
                 d_sigma: float = 0.45, coupling_gain: float = 1.20,
                 ghost_gain_ratio: float = 1.60, ghost_base_ratio: float = 0.25,
                 ghost_sigma_gain: float = 0.60,
                 ghost_bimodal_offset: float = 0.05,
                 d_sample_noise: float = 0.02, d_var_max: float = 0.01,
                 d_window: int = 200,
                 seed: int | None = None):

        self.dim_life = int(dim_life)
        self.dim_env = int(dim_env)
        self.d_threshold_c1 = float(d_threshold_c1)
        self.c2_compatibility = bool(c2_compatibility)
        self.semantic_potential_grad = float(semantic_potential_grad)
        self.d_sigma = float(d_sigma)
        self.coupling_gain = float(coupling_gain)
        self.ghost_gain_ratio = float(ghost_gain_ratio)
        self.ghost_base_ratio = float(ghost_base_ratio)
        self.ghost_sigma_gain = float(ghost_sigma_gain)
        self.ghost_bimodal_offset = float(ghost_bimodal_offset)
        self.d_sample_noise = float(d_sample_noise)
        self.d_var_max = float(d_var_max)
        self.d_window = int(d_window)

        self._rng = np.random.default_rng(seed)

        # 两个投影子系统
        self.L = LEProjectionSystem(self.dim_life, "L", noise_std=0.01, var_max=self.d_var_max)
        self.E = LEProjectionSystem(self.dim_env, "E", noise_std=0.01, var_max=self.d_var_max)
        self._reset_states()

        # 载体侧生命/环境参数（LE 体系 F/E 公理的数值化原型）
        self.r_life = 0.90          # 增殖率
        self.k_life = 1.60          # 承载容量
        self.uptake = 0.55          # 摄食强度
        self.half_sat = 0.50        # 半饱和常数
        self.maint = 0.35           # 维持阈值
        self.damage_rate = 0.80     # 损伤累积率
        self.repair = 0.60          # 修复率
        self.env_supply = 0.35      # 资源供给率
        self.env_cap = 0.45         # 资源再生上限
        self.env_consum = 0.30      # 资源被消耗系数

        self.t = 0.0
        self.history = {k: [] for k in (
            "t", "d", "D_L_var", "D_E_var", "D_L_mean", "D_E_mean",
            "region_label", "failure_mode", "ghost_magnitude",
            "L_state", "E_state", "semantic_potential_grad",
            "L_encapsulated", "E_encapsulated")}

    # ---------------- 属性别名 ----------------
    @property
    def semantic_pot_grad(self) -> float:
        """兼容 experiment 脚本中的简写名。"""
        return self.semantic_potential_grad

    @semantic_pot_grad.setter
    def semantic_pot_grad(self, value: float) -> None:
        self.semantic_potential_grad = float(value)

    def _reset_states(self) -> None:
        self.L.state = np.array([0.90, 0.20, 0.60][: self.dim_life], dtype=float)
        self.E.state = np.array([1.00, 0.25, 0.55][: self.dim_env], dtype=float)

    # ---------------- 核心：重合度 ----------------
    def calc_d_overlap(self, z_life, z_env) -> float:
        """边界重合度标量 d ∈ [0,1]。

        在公共基底上取两侧状态差，用双侧尺度均值无量纲化后取高斯核。
        相对距离 r 的含义：0 = 完全重合，1 = 差异与自身尺度同量级。
        """
        a = np.asarray(z_life, dtype=float).ravel()
        b = np.asarray(z_env, dtype=float).ravel()
        n = min(a.size, b.size)
        if n == 0:
            return 0.0
        a, b = a[:n], b[:n]
        diff = a - b
        scale = 0.5 * (np.linalg.norm(a) + np.linalg.norm(b)) + 1e-9
        r = float(np.linalg.norm(diff) / scale)
        return float(np.exp(-(r ** 2) / (2.0 * self.d_sigma ** 2)))

    # ---------------- 核心：C1/C2 三区判定 ----------------
    def judge_interaction_conditions(self, d: float):
        """返回 (C1, C2, region, failure_mode)。"""
        C1 = bool(d >= self.d_threshold_c1)
        C2 = bool(self.c2_compatibility)
        if not C1:
            region, fail = "low", "F1"
        elif not C2:
            region, fail = "middle", "F2"
        else:
            region, fail = "high", None
        return C1, C2, region, fail

    # ---------------- 核心：载体通量 ----------------
    def compute_coupling_flux(self, d: float, C1: bool, C2: bool):
        """返回 (flux_L, flux_E)。

        F1（C1 不成立）：通量恒为零——无公共细化构造，直接交互被禁止。
        正常耦合（C1&C2）：作用-反作用，幅值相反，通量守恒。
        F2（C1 成立、C2 失效）：主模态勉强对齐（第 0 分量守恒），
            但副模态（数值通量/边力量）幅值失配 ghost_gain_ratio 倍，
            净能量注入 → ghost 伪力 → 礁石伪收敛。
        """
        n = min(self.L.state.size, self.E.state.size)
        zeros = np.zeros(n, dtype=float)

        if not (C1 and d >= self.d_threshold_c1):
            return zeros, zeros.copy()

        delta = self.E.state[:n] - self.L.state[:n]
        gain = self.coupling_gain * float(self.semantic_potential_grad) * float(d)
        J = gain * delta

        if C2:
            return J, -J

        # F2：幅值失配的 ghost 伪力。
        # 副模态额外叠加一项与状态差无关的底噪项（对应 MMS 中 patch test 失败：
        # 即使主模态对齐，界面算子不一致仍产生伪力），这正是"伪收敛"的动力学来源。
        flux_L = J.copy()
        flux_E = -J.copy()
        if n > 1:
            flux_E[1:] = -self.ghost_gain_ratio * J[1:] - self.ghost_base_ratio * abs(gain)
        return flux_L, flux_E

    # ---------------- 载体动力学 ----------------
    def _life_rhs(self, s: np.ndarray) -> np.ndarray:
        s = np.asarray(s, dtype=float)
        bio = max(s[0], 0.0)
        growth = self.r_life * bio * (1.0 - bio / self.k_life)
        uptake = self.uptake * bio * max(self.E.state[0], 0.0) / (self.half_sat + max(self.E.state[0], 0.0))
        d0 = growth - uptake
        if s.size >= 3:
            dmg = self.damage_rate * max(0.0, uptake - self.maint) - self.repair * s[1]
            div = 0.25 * (1.0 - s[2] / 0.60) - 0.10 * bio
            return np.array([d0, dmg, div])
        return np.array([d0] + [0.0] * (s.size - 1))

    def _env_rhs(self, s: np.ndarray) -> np.ndarray:
        s = np.asarray(s, dtype=float)
        res = max(s[0], 0.0)
        supply = self.env_supply * (self.env_cap - res)
        consum = self.env_consum * max(self.L.state[0], 0.0) * res
        d0 = supply - consum
        if s.size >= 3:
            return np.array([d0, 0.30 * (0.25 - s[1]), 0.30 * (0.55 - s[2])])
        return np.array([d0] + [0.0] * (s.size - 1))

    # ---------------- 单步推进 ----------------
    def step(self, t: float, dt: float) -> dict:
        # 1) 边界处测量 d，判定条件
        d = self.calc_d_overlap(self.L.state, self.E.state)
        C1, C2, region, fail = self.judge_interaction_conditions(d)

        # 2) 载体通量（势只进增益，不进状态导数）
        flux_L, flux_E = self.compute_coupling_flux(d, C1, C2)
        ghost_mag = float(np.linalg.norm(flux_E + flux_L))          # 守恒时恒为 0
        gain_ref = abs(self.coupling_gain * self.semantic_potential_grad * d) + 1e-12
        ghost_w = float(np.clip(ghost_mag / gain_ref, 0.0, 1.0))    # 无量纲伪力占比

        # 3) 状态更新：载体动力学 + 通量
        f_L = self._life_rhs(self.L.state) + flux_L
        f_E = self._env_rhs(self.E.state) + flux_E
        self.L.state = np.maximum(self.L.state + dt * f_L, 0.0)
        self.E.state = np.maximum(self.E.state + dt * f_E, 0.0)

        # 4) D 分布采样（根张力 d0 + 边界测度噪声；F2 时方差随伪力抬升）
        if fail == "F2":
            sigma = self.d_sample_noise + self.ghost_sigma_gain * ghost_w
            shift = self.ghost_bimodal_offset * (1.0 if self._rng.random() < 0.5 else -1.0)
            s_L = d + self._rng.normal(0.0, sigma) + shift
            s_E = d + self._rng.normal(0.0, sigma) + shift
        else:
            sigma = self.d_sample_noise
            s_L = d + self._rng.normal(0.0, sigma)
            s_E = d + self._rng.normal(0.0, sigma)

        self.L.record_d_sample(float(np.clip(s_L, 0.0, 1.0)))
        self.E.record_d_sample(float(np.clip(s_E, 0.0, 1.0)))

        # 5) A10 边界封装检验（滚动窗口，跟踪当前尺度的 D 形态）
        enc_L = self.L.test_D_encapsulation(window=self.d_window)
        enc_E = self.E.test_D_encapsulation(window=self.d_window)

        # 6) 历史登记
        h = self.history
        h["t"].append(float(t))
        h["d"].append(float(d))
        h["D_L_var"].append(enc_L["var"] if np.isfinite(enc_L["var"]) else 0.0)
        h["D_E_var"].append(enc_E["var"] if np.isfinite(enc_E["var"]) else 0.0)
        h["D_L_mean"].append(enc_L["mean"] if np.isfinite(enc_L["mean"]) else d)
        h["D_E_mean"].append(enc_E["mean"] if np.isfinite(enc_E["mean"]) else d)
        h["region_label"].append(region)
        h["failure_mode"].append(fail)
        h["ghost_magnitude"].append(ghost_mag)
        h["L_state"].append(self.L.state.copy())
        h["E_state"].append(self.E.state.copy())
        h["semantic_potential_grad"].append(float(self.semantic_potential_grad))
        h["L_encapsulated"].append(enc_L.get("status", "fail") != "fail")
        h["E_encapsulated"].append(enc_E.get("status", "fail") != "fail")

        self.t = float(t)
        return {
            "t": float(t),
            "d": float(d),
            "C1": C1,
            "C2": C2,
            "region": region,
            "fail_mode": fail,
            "ghost_magnitude": ghost_mag,
            "L_state": self.L.state.copy(),
            "E_state": self.E.state.copy(),
            "D_L_encapsulation": enc_L,
            "D_E_encapsulation": enc_E,
        }
