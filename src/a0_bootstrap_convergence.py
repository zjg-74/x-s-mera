"""
X-S-MERA V3.4 · R7 遗留：A0 封闭方案甲（回归收敛）仿真
=========================================================

背景（修订登记 §五 第 4 条）
---------------------------
A0（设定的语义空间）与 A0*（认知者实际把握的语义空间）之间的差是根张力 d₀。
三种封闭方案：
  - 方案甲（回归收敛）：A0 的认知张力逐层衰减（粗粒化塔式），回归无害但需证明收敛；
  - 方案乙（自举封闭）：A0 的有效性由其下游全部成功投影回溯授证；
  - 方案丙（施为封闭）：A0 不可内证、只可外效。
R7 判定：**方案甲在人工系统四类中均无实例，零实证，降为待仿真项。**

本模块把"回归收敛"写成可证伪的动力学命题，并给出它的真实失效条件。

命题（形式化）
--------------
设认知者 k 的语义把握 A0*_k 由一次投影 P_k 从上一层 A0*_{k-1} 得到：
    A0*_k = Π(A0*_{k-1} ; θ_k) + η_k
其中 Π 是投影算子（方差压缩 g_k、均值偏移 b_k），η_k 是每层的投影噪声。

两种"回归"读法必须分开——这是方案甲此前被混淆的地方：

  **读法 R1（层级回归，塔式往上）**：k → ∞ 时 A0*_k → A0？
      在**有噪声**时：投影算子若收缩（g<1），噪声被压缩，d_k 收敛但收敛到
      **常驻偏置** b/(1−g)，不是 0。若 g→1 且 Σσ² < ∞，则 d_k → 0。
      结论：R1 收敛**依赖噪声可和**，而非"粗粒化本身"。

  **读法 R2（自我回归，自我指涉）**：A0*_{k+1} = Π(A0*_k)，问迭代是否收敛。
      这是压缩映射定理的直接应用：
        - ‖Π‖ < 1（严格收缩）→ 唯一不动点，且距离按 g^k 几何衰减 → **收敛**；
        - ‖Π‖ = 1（保距，如纯平移）→ 不收敛到 0，只保留初始偏置（等价方案丙）；
        - ‖Π‖ > 1（扩张）→ **发散**，张力逐层放大 → 方案甲在此类系统上被证伪。

因此方案甲的**真实适用域**是：投影算子在相关的张力度量下严格收缩。
R7 说"人工系统四类均无实例"，本仿真要回答的是更强的问题：
**这四类系统的投影算子，其收缩性如何？** —— 若存在一类可被视为收缩的，
方案甲在该类上就有实例；若全都不收缩，则方案甲应正式退役（降为数学趣例）。

四类系统的算子模型（与 R7 的类划分一一对应）
--------------------------------------------
  类1 纯语法机       ：Π = 恒等（语法变换不改变语义距离）→ g=1（保距）
  类2 被解释系统     ：Π = 规范-实现映射，误差有界但**跨层累积** → g≈1, 噪声可加
  类3 LLM 型         ：Π = 自回归自举，**承诺主体悬空** → 熵/张力随层数增（g>1）
  类4 演化功能语义   ：Π = 选择算子，**收敛到适应峰但非 A0**（不动点 ≠ A0）

这个模型给出的判定与 R7 的定性结论一致，但把"零实证"换成了可复算的判据：
**判据（预注册）**：
  C-1 收缩性：g < 1 且噪声方差有界 → 方案甲成立（张力几何衰减）；
  C-2 临界：g = 1 → 不收敛到零，但也不发散（等价于方案丙"登记 d₀"）；
  C-3 发散：g > 1 → 方案甲被证伪（张力随层级放大）。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "ProjectionChain",
    "run_chain",
    "fixed_point_iteration",
    "classify_system",
    "SYSTEM_MODELS",
]


class ProjectionChain:
    """层级回归（R1）：ε_k = g·ε_{k-1} + e_k，张力 d_k = |ε_k|。

    ε_k 是第 k 层认知者相对 A0 的偏差。关键语义：
      **每一层投影各注入一个独立的投影误差 e_k**（塔式回归的本义：
      每往上走一层就要多投影一次），而已有的偏差只被**部分保留**
      （保留率 g = 投影算子把上一层偏差原样传下来的比例）。

    于是 d_∞ 的行为（无噪声、偏置 b 每层叠加）：
      g < 1 : d → b/(1−g)  —— 有界，是**常驻偏置**
      g = 1 : d → b·k      —— 线性发散
      g > 1 : d → 几何发散
    有噪声时（e_k ~ N(0,σ²) 独立）：
      g < 1 : Var(d) → σ²/(1−g²)  —— 有界
      g ≥ 1 : 方差无界（g=1 时随 k 线性增）

    这个模型的意义：**"粗粒化塔式"本身不保证收敛到 0**。要让 d→0，
    必须同时满足 g<1 **且** 每层投影误差趋于零（b→0, σ→0）。
    这正是方案甲此前被当成"显然收敛"而漏掉的条件。
    """

    def __init__(self, gain: float = 0.95, layer_error: float = 0.0,
                 noise_std: float = 0.10, seed: int | None = None):
        self.gain = float(gain)
        self.layer_error = float(layer_error)   # 每层系统性投影误差 b
        self.noise_std = float(noise_std)       # 每层随机误差 σ
        self._rng = np.random.default_rng(seed)

    def run(self, n_layers: int, a0: float = 1.0) -> dict:
        d = np.empty(n_layers + 1, dtype=float)
        eps = 0.0                # 初始时 A0*_0 精确把握 A0 → 偏差为 0
        d[0] = 0.0
        noise = self._rng.normal(0.0, self.noise_std, size=n_layers)
        for k in range(1, n_layers + 1):
            eps = self.gain * eps + self.layer_error + noise[k - 1]
            d[k] = abs(eps)
        g = self.gain
        if abs(1.0 - g) > 1e-9:
            theo_bias = self.layer_error / (1.0 - g)
        else:
            theo_bias = float("inf") if self.layer_error > 0 else 0.0
        theo_var = (self.noise_std ** 2 / (1.0 - g ** 2)) if abs(g) < 1 else float("inf")
        tail_seg = d[n_layers // 2:]
        tail = float(tail_seg.mean())
        return {"layers": np.arange(n_layers + 1), "d": d,
                "theo_residual_bias": theo_bias, "theo_residual_var": theo_var,
                "tail_mean": tail,
                "tail_sq_mean": float((tail_seg ** 2).mean()),
                "tail_std": float(tail_seg.std()),
                "converged": bool(np.isfinite(d[-1]) and abs(d[-1] - d[-5]) < 1e-3)}


def run_chain(gains: list[float], n_layers: int = 60, layer_error: float = 0.02,
              noise_std: float = 0.10, seeds: list[int] | None = None) -> dict:
    """对一组增益跑层级回归，返回尾段张力（判定 C-1/C-2/C-3）。"""
    seeds = seeds or [0, 1, 2, 3, 4]
    rows = []
    for g in gains:
        tails, curves = [], []
        for s in seeds:
            pc = ProjectionChain(gain=g, layer_error=layer_error, noise_std=noise_std, seed=s)
            r = pc.run(n_layers)
            tails.append(r["tail_mean"])
            curves.append(r["d"])
        tail = float(np.mean(tails))
        if abs(g) < 1:
            verdict = "C-1 收缩（张力有界）"
        elif abs(g) == 1:
            verdict = "C-2 临界（不衰减至零）"
        else:
            verdict = "C-3 发散（方案甲失效）"
        rows.append({"gain": float(g), "tail_mean": tail,
                     "tail_max": float(np.max(tails)),
                     "verdict": verdict,
                     "curves": np.asarray(curves)})
    return {"rows": rows, "n_layers": n_layers}


def fixed_point_iteration(operator, x0: float, n_iters: int, tol: float = 1e-10) -> dict:
    """读法 R2（自我回归）的通用迭代器。

    operator(x) 返回 (new_x, local_lipschitz_estimate)。
    """
    xs = [float(x0)]
    lips = []
    for _ in range(n_iters):
        nx, L = operator(xs[-1])
        lips.append(float(L))
        xs.append(float(nx))
        if abs(xs[-1] - xs[-2]) < tol:
            break
    arr = np.asarray(xs)
    diffs = np.abs(np.diff(arr))
    ratio = (diffs[1:] / np.maximum(diffs[:-1], 1e-15)) if diffs.size > 1 else np.array([])
    return {"xs": arr, "iters": len(arr) - 1,
            "final_drift": float(diffs[-1]) if diffs.size else 0.0,
            "observed_ratio_median": float(np.median(ratio)) if ratio.size else float("nan"),
            "lipschitz_median": float(np.median(lips)) if lips else float("nan"),
            "converged": bool(diffs.size and diffs[-1] < tol)}


# ----------------------------------------------------------------------
# R7 四类人工系统的投影算子模型
# ----------------------------------------------------------------------

SYSTEM_MODELS = {
    # 类1 纯语法机：语法变换保距（g=1），不衰减；R7 判定 A0 条款不适用
    "class1_syntax": {"gain": 1.00, "noise_std": 0.05,
                      "note": "保距：张力不衰减（C-2）"},
    # 类2 被解释系统：规范→实现映射有界误差，但每层加新误差（g<1 但仍累积偏置）
    "class2_interpreted": {"gain": 0.98, "noise_std": 0.12,
                           "note": "近保距：收敛到常驻偏置而非 0"},
    # 类3 LLM 型：自举 → 语义漂移放大（g>1）
    "class3_llm": {"gain": 1.03, "noise_std": 0.15,
                   "note": "扩张：张力逐层放大（C-3）"},
    # 类4 演化功能语义：选择算子把表征压向局部适应峰（强收缩，但峰≠A0）
    "class4_evol": {"gain": 0.70, "noise_std": 0.10,
                    "note": "强收缩：收敛到一个'适应峰'，非 A0"},
}


def classify_system(name: str, n_layers: int = 200,
                    seeds: list[int] | None = None) -> dict:
    """对 R7 四类系统给出收缩性判定。"""
    if name not in SYSTEM_MODELS:
        raise KeyError(f"未知系统类型 {name}；可选 {list(SYSTEM_MODELS)}")
    cfg = SYSTEM_MODELS[name]
    seeds = seeds or [0, 1, 2, 3, 4]
    tails = []
    for s in seeds:
        pc = ProjectionChain(gain=cfg["gain"], layer_error=0.02,
                             noise_std=cfg["noise_std"], seed=s)
        tails.append(pc.run(n_layers)["tail_mean"])
    tail = float(np.mean(tails))
    g = cfg["gain"]
    if g < 0.95:
        cls = "C-1 强收缩"
    elif g < 1.0:
        cls = "C-1' 弱收缩（收敛到常驻偏置）"
    elif abs(g - 1.0) < 1e-12:
        cls = "C-2 临界（不衰减）"
    else:
        cls = "C-3 发散"
    return {"system": name, "gain": g, "tail_mean": tail,
            "tail_max": float(np.max(tails)), "class": cls, "note": cfg["note"]}
