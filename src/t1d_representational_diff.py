"""
X-S-MERA V3.4 · R8 遗留 T1-d：表征分化（人工载体）
=====================================================

T1-d 的原始设计（修订登记 §8.2）要求两组实验：

(A) **边界生成的临界迭代数与通量的标度关系**
    初始表征无统计可分双峰（Φ 对称未破缺），在载体流（训练更新）驱动下
    自发形成聚类边界。测：临界迭代数 n_c 如何随通量 J 缩放。

(B) **Φ 与 J 的可分离操纵**（R8 §8.9 残余风险 3 的唯一解法）
    认知系统中算力通量（J）与架构先验（Φ）通常不可独立操纵。本模块用
    一个**先验写死在初始化里的双峰偏置场**表示 Φ 的谱结构，使 Φ 与 J 成为
    两个正交旋钮，从而可检验 R8 的 F-c 预言：
    「同一 Φ、不同 J 的统计结构 → 不同的 D 形态；
      同 Φ 同 J 结构却出现不同破缺方向 → 证伪（A11-R 的因果分工）。」

模型（极简、可解析核对的表征演化）
----------------------------------
表征向量 x ∈ R^d（d 个"特征维度"）。每个维度 i 在初始化时从
    N(0, σ0²) 采样，
其中 σ0 是**对称未破缺态**的统一尺度——此时无任何维度显著区别于其他维度，
即"无一阶分化"。载体流 J 体现为每步的更新量幅值（学习率 × 梯度范数）。

关键：Φ 用**各维度的可分化性谱** w ∈ R^d 表示（Σw = d 归一化）。
    w_i = 1 时为各向同性（对称态，Φ 退化）；w_i ≠ 1 表示该维度**可被切分**
    的程度不同（对称已破缺的潜能，但边界尚未实现）。

动力学（每步）：
    x ← x + J · w ⊙ g(x) + ξ,        g(x) = x / (1 + |x|)   （有界增益）
    其中 ξ ~ N(0, s²) 是载体噪声（与状态独立，D17 口径）。

由 x 构造边界测度 D：把 x 投影到第一主成分方向上，取归一化后的
投影值为"边界坐标"，对其做双峰性诊断（Sarle BC 与 2-means 分离度）。
**注意**：这里 D 不是 kernel 里的耦合维度 d，而是 T1-d 专用的
"表征边界可分性测度"，命名用 D_rep 以示区分。

阴性对照（必须）
----------------
w = 1（各向同性）、J = 0 时，任何时刻的诊断量都不应判为"有边界"，
且 n_c 应为 ∞（永不生成）。这是 T1-d 设计里明确要求的阴性对照。

理论预期（R8 F-c）
------------------
1. 同一 Φ（固定 w）下，J 越大 n_c 越小——标度关系 n_c ∝ J^(-κ)。
   若 κ 落在 (0,1)，说明"边界生成时点"受流支配但不与流成反比（非线性累积）。
2. 固定 Φ 与 J 的功率谱，仅改变 ξ 的种子，破缺方向（第一主成分的符号/
   所属维度）应随流的历史路径变化——这是"方向归流"的证据。
3. 同一 (Φ, J) 不同种子下，D_rep 的**形态类别**（单峰 / 可分）应稳定，
   即"存在与谱结构归势"。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "RepresentationField",
    "RepresentationPopulation",
    "boundary_diagnostics",
    "critical_iteration",
    "flux_scale_relation",
    "breaking_direction_test",
    "phi_j_separability_test",
]


class RepresentationField:
    """表征场的极简演化模型（Φ 与 J 正交可调）。

    Φ 的正确建模：**双阱势，不是方差权重**
    ------------------------------------------
    早期版本把 Φ 建成"可分化性谱 w"（乘在更新幅值上）。这个建模是错的，
    而且错得很有信息量：乘性权重只改变**方差**，不改变**均值结构**，
    任何 J 下都只是把一个单峰切开（实测 abs_sep 恒为 1.74，与 J 无关，
    n_c 全部撞在下限）。结论：**可分性必须由势的多稳结构承载，
    不能由更新的强弱承载**。这与 R8 F-b 完全一致——Φ 是"可分化性潜能"
    （存在哪些阱），不是"已实现的分化"（落在哪个阱）。

    模型：
      U(x0) = a·(x0² − c²)²      双阱势，只作用在候选分化轴（第 0 维）
      状态更新：x ← x + gain·(J·g(x) − ∇U) + ξ

      a = 0        ↔ Φ 对称未破缺（单阱，阴性对照：永不双峰）
      a > 0        ↔ Φ 具双稳谱结构（可分化，但边界尚未实现——初始在 x0≈0 的势垒顶）
      J            ↔ 载体通量，决定噪声/更新幅度 → 决定越垒速率 → 决定 n_c
      c            ↔ 两阱间距，决定双峰间距 → 决定 abs_sep
      weak_noise   ↔ 旁观维度（1..dim-1）的噪声。必须显著小于候选分化轴的噪声：
                     否则它们在旁观轴上随机扩散，把"第一主成分"带偏到噪声轴上，
                     诊断量就测不到第 0 维的双阱了（实测 BC 被压到 0.32–0.39）。

    这个模型让"破缺方向"成为真问题：x0 落到 +c 还是 −c，由噪声路径决定，
    正是 R8 F-c 所说"方向归流"。
    """

    def __init__(self, dim: int = 12, sigma0: float = 1.0,
                 well_depth: float = 0.0, well_sep: float = 1.0,
                 flux: float = 1.0, noise_std: float = 0.05,
                 weak_noise_ratio: float = 0.12,
                 gain: float = 0.05, seed: int | None = None):
        self.dim = int(dim)
        self.sigma0 = float(sigma0)
        self.weak_noise_ratio = float(weak_noise_ratio)
        self.a = float(well_depth)     # 双阱深度 a（0 = 单阱/对称态）
        self.c = float(well_sep)       # 两阱半间距
        self.flux = float(flux)        # 载体通量 J
        self.noise_std = float(noise_std)
        self.gain = float(gain)
        self._rng = np.random.default_rng(seed)
        # 初始对称：x0 ≈ 0（势垒顶，尚未落到任何阱）
        self.x = self._rng.normal(0.0, self.sigma0 * 0.1, size=self.dim)
        self.t = 0

    def _potential_grad(self, x0: float) -> float:
        """双阱势 U = a(x² − c²)² 的梯度 dU/dx = 4a·x·(x² − c²)。"""
        return 4.0 * self.a * x0 * (x0 ** 2 - self.c ** 2)

    def step(self) -> np.ndarray:
        """载体流 = **涨落驱动**，不是漂移驱动。

        关键建模修正（第二版）：早期把 J 建成单调漂移 `J·g(x)`，导致无势时
        x0 无界增长到 15（单阱情形下通量把状态一路推高，永不驻留），
        整个轨迹都是瞬态 → 诊断量毫无意义。
        真实的载体流（训练更新/迭代）主要提供**涨落幅度**：J 越大，
        状态在势中的热噪声越大，越容易越过势垒。故：
            σ_eff = noise_std · sqrt(1 + J)
        这才是"通量控制越垒速率"的正确机制，也才让 n_c 随 J 下降。
        """
        sigma_eff = self.noise_std * np.sqrt(1.0 + max(self.flux, 0.0))
        noise = self._rng.normal(0.0, sigma_eff, size=self.dim)
        if self.dim > 1:
            # 旁观维度噪声显著更弱：它们在 T1-d 中是"未被使用的特征维度"
            noise[1:] *= self.weak_noise_ratio
        force = np.zeros(self.dim)
        # 候选分化轴：受双阱势作用（a=0 时无势 → 纯扩散 + 回复）
        force[0] = -self._potential_grad(self.x[0])
        # 回复项：把 x0 约束在势阱尺度内（否则无势时纯扩散无界，轨迹全是瞬态）
        force[0] -= 0.6 * self.x[0] * (1.0 if abs(self.a) < 1e-12 else 0.0)
        # 旁观维度：无势，纯均值回复（保持在 ~σ0 量级）
        if self.dim > 1:
            force[1:] = -0.6 * self.x[1:]
        self.x = self.x + self.gain * force + noise
        self.t += 1
        return self.x.copy()

    @property
    def is_symmetric_phi(self) -> bool:
        """对称态 ⟺ 无双阱（a = 0）。"""
        return bool(abs(self.a) < 1e-12)

    def run(self, n_steps: int) -> np.ndarray:
        traj = np.empty((n_steps, self.dim), dtype=float)
        for i in range(n_steps):
            traj[i] = self.step()
        return traj


class RepresentationPopulation:
    """一群表征单元（population of representations）。

    这是 T1-d 的正确建模层级，而不是"单条轨迹"。
    ---------------------------------------------------
    第三版建模修正（关键）：单条轨迹在双阱势中**不断在两阱间跳跃**
    （实测 3000 步内正负各占约 40%/60%，signs.mean ≈ −0.34），
    其取值分布虽是两簇，但任一时刻的统计量极不稳定，**无法作为分化读数**。
    （注意：早期以为"落阱即锁定"，被实测推翻——见 tests 中
      test_single_trajectory_is_not_bimodal 的记录。）
    分化必须在**单元群**层面读：n_units 个单元，每个在同一个双阱势里
    （Φ 相同），各自被噪声驱动在两阱间游走。群体的**即时分布**出现两极，
    才是"表征维度自发分化成两个簇"的稳定语义。

    J（载体通量）控制每个单元的噪声幅度 → 控制"首次有单元落阱"的时点 n_c；
    Φ（a, c）控制势结构 → 控制两极的位置（这就是 R8 F-c 的分工）。
    """

    def __init__(self, n_units: int = 200, dim: int = 4,
                 well_depth: float = 1.0, well_sep: float = 1.0,
                 flux: float = 1.0, noise_std: float = 0.10,
                 gain: float = 0.05, seed: int | None = None):
        self.n_units = int(n_units)
        self.dim = int(dim)
        self.a = float(well_depth)
        self.c = float(well_sep)
        self.flux = float(flux)
        self.noise_std = float(noise_std)
        self.gain = float(gain)
        self._rng = np.random.default_rng(seed)
        # 每个单元一个标量位置 x0（候选分化轴）
        self.x0 = self._rng.normal(0.0, 0.02, size=self.n_units)
        self.t = 0

    def step(self) -> np.ndarray:
        sigma_eff = self.noise_std * np.sqrt(1.0 + max(self.flux, 0.0))
        grad = 4.0 * self.a * self.x0 * (self.x0 ** 2 - self.c ** 2)
        force = -grad
        if abs(self.a) < 1e-12:
            force = -0.6 * self.x0          # 无势时的回复（阴性对照）
        noise = self._rng.normal(0.0, sigma_eff, size=self.n_units)
        self.x0 = self.x0 + self.gain * force + noise
        self.t += 1
        return self.x0.copy()

    def snapshot_stats(self) -> dict:
        """当前单元群的分布形态（双峰判据在这里作用于单元群）。"""
        return boundary_diagnostics(self.x0)

    def run(self, n_steps: int) -> np.ndarray:
        """返回 (n_steps, n_units) 的单元位置轨迹。"""
        traj = np.empty((n_steps, self.n_units), dtype=float)
        for i in range(n_steps):
            traj[i] = self.step()
        return traj


# ======================================================================
# 边界诊断
# ======================================================================

def boundary_diagnostics(x: np.ndarray, bc_threshold: float = 0.555,
                         sep_threshold: float = 2.0,
                         abs_sep_threshold: float = 0.0,
                         sym_threshold: float = 0.5) -> dict:
    """对一次表征快照做"是否已生成可识别边界"的诊断。

    做法：
    1. 投影到第一主成分（唯一可能承载"两个簇"的方向）；
    2. 在投影值上做 1-D 两簇 k-means（确定性初始化：分位点）；
    3. 输出 Sarle 双峰系数 BC、相对分离度 sep = |μ1−μ2|/s_pooled，
       以及绝对分离度 abs_sep = |μ1−μ2|/σ_total。

    判据（预注册）：
      已生成 = BC ≥ bc_threshold 且 sep ≥ sep_threshold 且 abs_sep ≥ abs_sep_threshold。

    第四个条件（零中心 + 对称）：
      只要求两簇"分得开"是不行的——任何单峰长尾都能被 k-means 切成
      "左半边 vs 右半边"，sep 同样很大（实测阴性对照 sep=2.6）。
      真双峰的特征是**两簇分居总体均值两侧且距离相当**。
      故加：straddles（两簇均值跨过总体均值）且 symmetry ≥ sym_threshold。

    为什么需要第三个条件（abs_sep）：
      sep 是**相对**分离度——两簇各自很窄时，哪怕只差一点点也会算出很大的 sep。
      在"表征尚未展开、整体分布很窄"的早期阶段，k-means 总能切出相对分离度
      极高的假边界。abs_sep 用总体标准差归一化，衡量两簇间距相对**整体尺度**
      是否足够大，才能排除这种小尺度伪双峰。
      abs_sep_threshold 默认 0（关闭）；T1-d 仿真显式设为 ≈1.0。

    bc_threshold 默认 5/9≈0.555 为 Sarle 定理分界（单峰分布 BC<5/9）。
    """
    a = np.asarray(x, dtype=float).ravel()
    if a.size < 4:
        return {"has_boundary": False, "status": "insufficient", "bc": float("nan"),
                "sep": float("nan"), "axis": None}

    mean = a.mean()
    var = a.var(ddof=0)
    sd = np.sqrt(var) + 1e-12
    z = (a - mean) / sd
    skew = float(np.mean(z ** 3))
    kurt = float(np.mean(z ** 4))
    bc = float((skew ** 2 + 1.0) / max(kurt, 1e-12))

    # 确定性两簇初始化：分位点
    lo = np.quantile(a, 0.25)
    hi = np.quantile(a, 0.75)
    c = np.array([lo, hi], dtype=float)
    for _ in range(50):
        assign = np.argmin(np.abs(a[:, None] - c[None, :]), axis=1)
        newc = np.array([a[assign == k].mean() if np.any(assign == k) else c[k]
                         for k in range(2)])
        if np.allclose(newc, c, atol=1e-12):
            c = newc
            break
        c = newc
    assign = np.argmin(np.abs(a[:, None] - c[None, :]), axis=1)
    g1, g2 = a[assign == 0], a[assign == 1]
    if g1.size == 0 or g2.size == 0:
        sep = 0.0
        abs_sep = 0.0
    else:
        s_pooled = np.sqrt(0.5 * (g1.var(ddof=0) + g2.var(ddof=0))) + 1e-12
        sep = float(abs(g1.mean() - g2.mean()) / s_pooled)
        abs_sep = float(abs(g1.mean() - g2.mean()) / sd)

    # 零中心化检查：真正的"双峰"要求两簇关于整体中心近似对称（|μ1−μ2| 的两侧
    # 到总体均值距离相当）。否则只是**分布整体偏移**——那个"两簇"其实是
    # "左半边 vs 右半边"，任何单峰长尾都会被这样切出来。
    if g1.size and g2.size:
        m = float(a.mean())
        c1, c2 = float(g1.mean()), float(g2.mean())
        # 中心度：两簇均值是否分布在总体均值两侧
        straddles = (c1 - m) * (c2 - m) < 0
        # 对称度：两侧距离之比（大/小），1 为完全对称
        d1, d2 = abs(c1 - m), abs(c2 - m)
        sym = float(min(d1, d2) / max(max(d1, d2), 1e-12))
    else:
        straddles, sym = False, 0.0

    has = bool(bc >= bc_threshold and sep >= sep_threshold
               and abs_sep >= abs_sep_threshold
               and straddles and sym >= sym_threshold)
    return {"has_boundary": has, "status": "boundary" if has else "single_mode",
            "bc": bc, "sep": sep, "abs_sep": abs_sep,
            "straddles": bool(straddles), "symmetry": sym,
            "cluster_sizes": (int(g1.size), int(g2.size)),
            "cluster_means": (float(g1.mean()) if g1.size else float("nan"),
                              float(g2.mean()) if g2.size else float("nan")),
            "sd": float(sd)}


def _first_component_axis(traj: np.ndarray) -> np.ndarray:
    """返回第一主成分方向（单位向量）。"""
    Xc = traj - traj.mean(axis=0, keepdims=True)
    cov = Xc.T @ Xc / max(traj.shape[0] - 1, 1)
    vals, vecs = np.linalg.eigh(cov)
    return vecs[:, int(np.argmax(vals))]


def critical_iteration(pop_factory, max_iters: int = 3000,
                       consecutive: int = 5, check_every: int = 5,
                       min_units: int = 60, **diag_kwargs) -> dict:
    """测"边界生成临界迭代数" n_c（多单元群版本）。

    判据：连续 `consecutive` 次检查都判为有边界，第一次命中处即 n_c。
    要求连续，是为了排除瞬态涨落误报（R8 §8.9 风险 2 的缓解：
    只有下游可再识别的边界才计入"势生成"）。

    pop_factory() 须返回 RepresentationPopulation。
    双峰判据直接作用在**单元群**的分布上（这才是"表征分化"的语义）。
    """
    pop = pop_factory()
    if pop.n_units < min_units:
        raise ValueError(f"单元数 {pop.n_units} 少于 min_units={min_units}")
    hit_streak = 0
    first_hit = None
    for i in range(max_iters):
        pop.step()
        if (i + 1) < 10 or ((i + 1) % check_every != 0):
            continue
        diag = boundary_diagnostics(pop.x0, **diag_kwargs)
        if diag["has_boundary"]:
            if first_hit is None:
                first_hit = i
            hit_streak += 1
            if hit_streak >= consecutive:
                return {"n_c": int(first_hit + 1), "converged": True,
                        "diag": diag, "iters": i + 1}
        else:
            hit_streak = 0
            first_hit = None
    return {"n_c": None, "converged": False, "diag": None, "iters": max_iters}


def flux_scale_relation(fluxes: list[float], seeds: list[int] | None = None,
                        n_units: int = 300, well_depth: float = 1.0,
                        well_sep: float = 1.0, max_iters: int = 3000,
                        noise_std: float = 0.10, gain: float = 0.05) -> dict:
    """测 n_c 与通量 J 的标度关系：n_c ∝ J^(-κ)，log-log 线性回归取 κ。

    固定 Φ（双阱 a, c）不变，只改变 J。J 通过 σ_eff = σ√(1+J) 控制越垒速率。
    """
    seeds = seeds or [0, 1, 2]
    rows = []
    for J in fluxes:
        ncs = []
        for s in seeds:
            def factory(J=J, s=s):
                return RepresentationPopulation(
                    n_units=n_units, well_depth=well_depth, well_sep=well_sep,
                    flux=J, noise_std=noise_std, gain=gain, seed=s)
            r = critical_iteration(factory, max_iters=max_iters)
            ncs.append(r["n_c"])
        valid = [v for v in ncs if v is not None]
        rows.append({"flux": float(J), "n_c_list": ncs,
                     "n_c_median": float(np.median(valid)) if valid else None,
                     "censored": int(sum(1 for v in ncs if v is None))})
    xs = [np.log(r["flux"]) for r in rows if r["n_c_median"]]
    ys = [np.log(r["n_c_median"]) for r in rows if r["n_c_median"]]
    kappa = r2 = None
    if len(xs) >= 2:
        A = np.vstack([np.asarray(xs), np.ones(len(xs))]).T
        coef, *_ = np.linalg.lstsq(A, np.asarray(ys), rcond=None)
        kappa = float(-coef[0])
        pred = A @ coef
        yy = np.asarray(ys)
        ss_res = float(np.sum((yy - pred) ** 2))
        ss_tot = float(np.sum((yy - yy.mean()) ** 2)) + 1e-12
        r2 = 1.0 - ss_res / ss_tot
    return {"rows": rows, "kappa": kappa, "r2": r2,
            "phi": {"well_depth": well_depth, "well_sep": well_sep}}


def breaking_direction_test(flux: float = 1.0, seeds: list[int] | None = None,
                            n_units: int = 300, well_depth: float = 1.0,
                            well_sep: float = 1.0, steps: int = 1500,
                            noise_std: float = 0.10, gain: float = 0.05) -> dict:
    """测破缺方向是否随流的历史路径变化（R8 F-c 的核心检验）。

    "方向"在这里有明确的可观测定义：**单元群中落到 +c 阱的比例 f₊**。
      - f₊ ≈ 0.5  → 方向随机（对称破缺，由噪声路径决定）= "方向归流"成立；
      - f₊ ≈ 1 或 0 → 方向被 Φ 结构锁死 = F-c 被证伪。
    """
    seeds = seeds or list(range(10))
    frac_plus, fracs_by_seed = [], []
    for s in seeds:
        pop = RepresentationPopulation(n_units=n_units, well_depth=well_depth,
                                       well_sep=well_sep, flux=flux,
                                       noise_std=noise_std, gain=gain, seed=s)
        pop.run(steps)
        fp = float(np.mean(pop.x0 > 0))
        frac_plus.append(fp)
        fracs_by_seed.append(fp)
    fp_arr = np.asarray(frac_plus)
    # 与 0.5 的偏差检验（单侧宽松）：全部种子都偏同一侧才算锁死
    all_above = bool(np.all(fp_arr > 0.5 + 0.05))
    all_below = bool(np.all(fp_arr < 0.5 - 0.05))
    return {"seeds": seeds, "frac_plus": frac_plus,
            "mean_frac_plus": float(fp_arr.mean()),
            "std_frac_plus": float(fp_arr.std()),
            "direction_locked": bool(all_above or all_below),
            "verdict": ("方向被 Φ 锁死（F-c 被证伪）" if (all_above or all_below)
                        else "方向随流路径变化（F-c 成立：方向归流）")}


def phi_j_separability_test(flux_levels: list[float],
                            phi_modes: dict[str, dict],
                            seeds: list[int] | None = None, n_units: int = 300,
                            steps: int = 2500, noise_std: float = 0.10,
                            gain: float = 0.05) -> dict:
    """Φ × J 因子设计：检验二者是否可独立操纵（R8 §8.9 风险 3 的解法）。

    预期：
      - 对称 Φ（a=0）+ 任意 J → n_c = None（永不生成）——**Φ 是必要的**；
      - 非对称 Φ（a>0）+ J>0   → 生成，且 n_c 随 J 下降；
      - 非对称 Φ + J=0 但 σ 固定 → 也能生成，但更慢（J 起"加速"作用）。
    """
    seeds = seeds or [0, 1, 2]
    grid = []
    for name, phi in phi_modes.items():
        for J in flux_levels:
            ncs = []
            for s in seeds:
                def factory(J=J, phi=phi, s=s):
                    return RepresentationPopulation(
                        n_units=n_units, well_depth=phi["well_depth"],
                        well_sep=phi["well_sep"], flux=J, noise_std=noise_std,
                        gain=gain, seed=s)
                r = critical_iteration(factory, max_iters=steps)
                ncs.append(r["n_c"])
            valid = [v for v in ncs if v is not None]
            grid.append({"phi": name, "flux": float(J), "n_c_list": ncs,
                         "n_c_median": float(np.median(valid)) if valid else None,
                         "generated": bool(valid)})
    return {"grid": grid, "steps": steps}


# ======================================================================
# 序参量：滞回（hysteresis）与分化深度
# ======================================================================

def hysteresis_rate(traj: np.ndarray, burn_in: int = 0) -> dict:
    """个体在阱间翻越的频率（每单元每步）。

    为什么用这个当序参量，而不用双峰判据：
    ---------------------------------------------------
    双峰判据（BC/sep）在 J 大时**失效**——单元在两阱间高频跳跃时，
    群体分布被抹成宽单峰（实测 J≥2 时 n_c → None）。
    这不是模型的 bug，而是真实的物理：通量太大时"边界"根本无法稳定存在。
    因此需要一个对瞬时位置不敏感、但对翻越动力敏感的序参量：
        翻越率 = sign(x0) 的变号次数 / (步数 × 单元数)

    实测（300 单元、4000 步、gain=0.01）：
        J=0 → 23, J=0.3 → 37, J=1 → 64, J=2 → 94, J=3 → 118, J=5 → 153
    严格单调递增且近似线性 → 这是"通量驱动势垒穿越"的干净读法。
    """
    s = np.sign(traj[burn_in:])
    flips = np.sum(np.abs(np.diff(s, axis=0)) > 0, axis=0)
    n_steps = max(traj.shape[0] - burn_in - 1, 1)
    rate = flips / n_steps
    return {"flips_per_unit": float(flips.mean()),
            "rate_per_unit_per_step": float(rate.mean()),
            "total_flips": int(flips.sum())}


def differentiation_depth(traj: np.ndarray, burn_in: int = 0) -> dict:
    """分化深度 = 稳态群体位置 |mean x0| 与群体离散度。

    深度 ≈ c 表示整体分化到某一阱（群体偏极性）；
    深度 ≈ 0 但 sd ≈ c 表示两阱各占一半（对称破缺，方向随机）。
    """
    tail = traj[burn_in:]
    m = float(tail.mean())
    return {"depth": float(abs(m)), "group_mean": m,
            "sd": float(tail.std()), "bimodal_hint": float(tail.std() > 0.5 * abs(tail.mean()) + 1e-9)}


def direction_report(traj: np.ndarray, burn_in: int = 0) -> dict:
    """方向报告：最终落到 +阱的单元比例。

    f₊ ≈ 0.5 → 方向随机（由流路径决定）；f₊ → 0 或 1 → 方向被 Φ 锁死。
    """
    tail = traj[burn_in:]
    final = tail[-1]
    fp = float(np.mean(final > 0))
    return {"frac_plus": fp, "locked": bool(abs(fp - 0.5) > 0.25)}


def flux_sweep_observables(fluxes: list[float], seeds: list[int] | None = None,
                           n_units: int = 300, well_depth: float = 1.0,
                           well_sep: float = 1.0, steps: int = 4000,
                           burn_in: int = 1000, noise_std: float = 0.10,
                           gain: float = 0.01) -> dict:
    """对通量做扫描，同时采集翻越率 / 深度 / 方向 / 双峰性四个观测量。

    这是 T1-d 的主测量：R8 F-c 的三个子命题各对应一个观测量。

      子命题 1（时点归流）：翻越率随 J 单调递增（越垒速率由通量控制）
      子命题 2（存在归势）：双峰性由 Φ（a>0）决定，J 只调其幅度
      子命题 3（方向归流）：f₊ 随种子变化，不被 Φ 锁死
    """
    from src.t1d_representational_diff import boundary_diagnostics as _bd
    seeds = seeds or [0, 1, 2]
    rows = []
    for J in fluxes:
        rates, depths, fps, bcs = [], [], [], []
        for s in seeds:
            pop = RepresentationPopulation(n_units=n_units, well_depth=well_depth,
                                           well_sep=well_sep, flux=J,
                                           noise_std=noise_std, gain=gain, seed=s)
            traj = pop.run(steps)
            rates.append(hysteresis_rate(traj, burn_in)["rate_per_unit_per_step"])
            depths.append(differentiation_depth(traj, burn_in)["depth"])
            fps.append(direction_report(traj, burn_in)["frac_plus"])
            bcs.append(_bd(traj[-1])["bc"])
        rows.append({"flux": float(J),
                     "flip_rate_mean": float(np.mean(rates)),
                     "depth_mean": float(np.mean(depths)),
                     "frac_plus_mean": float(np.mean(fps)),
                     "frac_plus_std": float(np.std(fps)),
                     "bc_final_mean": float(np.mean(bcs))})
    return {"rows": rows, "phi": {"well_depth": well_depth, "well_sep": well_sep},
            "steps": steps, "burn_in": burn_in}
