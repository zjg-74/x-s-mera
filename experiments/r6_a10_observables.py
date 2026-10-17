"""R6: A10 封装性失败的可观测现象清单 —— 六场景判别力测试.

本脚本只做一件事：**检验预注册判据表（src/a10_observables.py）有没有判别力**。

它不"证明" A10 成立。它回答三个更小的问题：
  1. 每条判据在良性系统上会不会误报（特异性）；
  2. 两类失效（礁石 / 分辨率不足）与真坍塌在判据下能不能分开；
  3. 这份清单有没有盲区，盲区是什么。

场景（除 BIAS 外全部由内核真实动力学产生）
------------------------------------------
OK       稳定耦合，C1+C2 全成立              → 期望：无信号
F2       接口不兼容（ghost 伪力）            → 期望：礁石（修接口即恢复）
F1       d 跌破阈值（d_threshold=0.99）      → 期望：区域门排除，A10 不适用
SUB_OU   亚网格隐变量：OU 过程（有限相关时间）→ 期望：分辨率伪影（变粗即平均掉）
SUB_1F   亚网格隐变量：1/f 过程（无特征尺度） → 期望：M-不完备候选（坍塌）
BIAS     亚网格隐变量：恒定偏置（整体平移）   → 期望：**假阴性**（判据盲区）

BIAS 场景是故意加进来的阴性对照：P2–P5 全是位置不变量（IQR/BC/std/斜率），
对整体平移天生失明。若它被判为"M-完备"，说明清单的盲区与预注册所述一致，
这比"全对"更有价值——一份声称无盲区的清单本身就不可信。
"""

from __future__ import annotations

import copy
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.a10_observables import (
    CRITERIA, DELTA_ABS, DELTA_QUOTED, BETA_MAX, VERDICTS,
    classify, criterion_table, evaluate, evaluate_persistent,
    granularity_scan, repair_encapsulated, window_stats,
)
from src.le_coupling_kernel import LECouplingKernel

FIG_DIR = "figs"
os.makedirs(FIG_DIR, exist_ok=True)

DT = 0.05
T_TOTAL = 800.0
STEADY_FRAC = 0.40          # 只取后 60% 作稳态评估
REPAIR_STEPS = 400
REPAIR_TAIL = 200

COLORS = {
    "OK": "#27ae60", "F2": "#f39c12", "F1": "#e74c3c",
    "SUB_OU": "#2980b9", "SUB_1F": "#8e44ad", "BIAS": "#7f8c8d",
}

SCENARIOS = {
    "OK": dict(c2=True, drv=None, d_thr=0.25),
    "F2": dict(c2=False, drv=None, d_thr=0.25),
    "F1": dict(c2=True, drv=None, d_thr=0.99),
    "SUB_OU": dict(c2=True, drv=("ou", 0.030), d_thr=0.25),
    "SUB_1F": dict(c2=True, drv=("1f", 0.030), d_thr=0.25),
    # 偏置取负：d 稳态在 0.98 附近，正向偏置会被观测层 clip 成满格常数，
    # 那样测到的是裁剪伪影而不是盲区。
    "BIAS": dict(c2=True, drv=("bias", -0.060), d_thr=0.25),
}


# ---------------------------------------------------------------- 亚网格驱动
class SubgridDriver:
    """未建模的亚网格自由度：投影看不到它，但它渗进边界测量。

    kind:
      "ou"   Ornstein–Uhlenbeck，相关时间 TAU_OU=0.10s = 2 步（有限特征尺度）
      "1f"   幂律谱 S(f) ∝ 1/f，带宽跨整个记录长度（无特征尺度）
      "bias" 恒定偏置（整体平移；用于暴露判据盲区）

    注入方式：**不进入状态方程**，只叠加到边界观测量上——
    ``D_observed = d(state) + eps(t)``。这是"亚网格"这个词的字面含义：
    Ω 算子把它截断在粗粒化之外，所以它不能改变宏观轨迹，但边界测度
    里留着它的残差。这正是 M-不完备的定义（未闭合的残差出现在观测量上），
    也是多尺度建模里被丢弃的高频自由度在握手区的表现。
    """

    TAU_OU = 0.10

    def __init__(self, kind: str, amplitude: float, dt: float, n: int, seed: int):
        self.kind = kind
        self.amp = float(amplitude)
        self.dt = float(dt)
        self.n = int(n)
        rng = np.random.default_rng(seed)
        self.series = self._make(rng)

    def _make(self, rng) -> np.ndarray:
        n, dt = self.n, self.dt
        if self.kind == "bias":
            return np.full(n, self.amp, float)
        if self.kind == "ou":
            phi = np.exp(-dt / self.TAU_OU)
            sd = np.sqrt(1.0 - phi ** 2)
            e = np.empty(n, float)
            e[0] = rng.normal()
            for i in range(1, n):
                e[i] = phi * e[i - 1] + sd * rng.normal()
            return e * self.amp
        if self.kind == "1f":
            # 1/f 幅度谱（功率 ∝ 1/f）。下截止取 1/n（不设特征尺度）——
            # 若像初版那样切在 1/512，块平均在 ~1.5 个数量级内仍会拿掉一大截
            # 对数功率，斜率被拉到 -0.22，把"无特征尺度"误判成"有限尺度"。
            f = np.fft.rfftfreq(n, d=1.0)
            f[0] = f[1]
            band = f <= 0.45
            amp = np.zeros_like(f)
            amp[band] = 1.0 / np.sqrt(f[band])
            phase = rng.uniform(0, 2 * np.pi, f.size)
            x = np.fft.irfft(amp * np.exp(1j * phase), n)
            return x / x.std() * self.amp
        raise ValueError(f"unknown driver kind: {self.kind}")


def _observe(d: np.ndarray, drv: np.ndarray | None, offset: int = 0) -> np.ndarray:
    """把亚网格残差叠加到观测量上（不进状态方程）。"""
    d = np.asarray(d, float)
    if drv is None:
        return d
    seg = drv[offset:offset + d.size]
    return np.clip(d + seg, 0.0, 1.0)


def sweep_scale_separation(seed: int = 5) -> list[dict]:
    """τ 扫描：把 OU 的相关时间从亚采样步扫到超出扫描窗,看 beta 连续过渡。

    这是 P5 的**判别力边界**：判据只能在"特征尺度落在扫描窗口之外"时
    给出干净的符号；特征尺度落在窗口内部（τ ≈ 1–128 步）时斜率连续变化，
    任何固定阈值都只是把这条连续谱切成两半,而不是识别出两种机制。
    """
    n = 8000
    out = []
    for tau in (0.02, 0.05, 0.10, 0.25, 0.50, 1.0, 2.0, 4.0):
        phi = np.exp(-DT / tau)
        rng = np.random.default_rng(seed)
        e = np.empty(n, float)
        e[0] = rng.normal()
        for i in range(1, n):
            e[i] = phi * e[i - 1] + np.sqrt(1 - phi ** 2) * rng.normal()
        x = 0.03 * e + rng.normal(0, DELTA_QUOTED, n)
        sc = granularity_scan(x, max_level=7)
        out.append({"tau": tau, "tau_steps": tau / DT, "beta": sc["beta"],
                    "fires_P5": bool(sc["beta"] >= BETA_MAX),
                    "b": sc["b"], "a_norm": [v / sc["a"][0] for v in sc["a"]]})
    return out


# ---------------------------------------------------------------- 场景运行
def run_scenario(name: str, seed: int = 11) -> dict:
    spec = SCENARIOS[name]
    n_steps = int(T_TOTAL / DT)
    kernel = LECouplingKernel(
        dim_life=3, dim_env=3,
        d_threshold_c1=spec["d_thr"],
        c2_compatibility=spec["c2"],
        semantic_potential_grad=0.40,
        coupling_gain=6.0,          # 与 sim_le_coupling 同一标定
        d_sample_noise=DELTA_QUOTED,   # 申报分辨率 = 观测噪声底
        seed=seed,
    )
    drv = (None if spec["drv"] is None
           else SubgridDriver(spec["drv"][0], spec["drv"][1], DT, n_steps + REPAIR_STEPS, seed))

    d_all, c1_all = [], []
    for i in range(n_steps):
        res = kernel.step(i * DT, DT)
        d_all.append(kernel.L.d_samples[-1])   # 边界测度样本（含 ghost 污染）
        c1_all.append(res["C1"])

    d_all = _observe(np.asarray(d_all, float),
                     None if drv is None else drv.series, 0)
    c1_all = np.asarray(c1_all, bool)
    start = int(n_steps * STEADY_FRAC)
    d_ss, c1_ss = d_all[start:], c1_all[start:]

    # ---- 接口修复检验：恢复 C2，继续跑，看封装性是否回来
    k2 = copy.deepcopy(kernel)
    k2.c2_compatibility = True
    d_rep, c1_rep = [], []
    for j in range(REPAIR_STEPS):
        res = k2.step((n_steps + j) * DT, DT)
        d_rep.append(k2.L.d_samples[-1])
        c1_rep.append(res["C1"])
    d_rep = _observe(np.asarray(d_rep, float),
                     None if drv is None else drv.series, n_steps)[-REPAIR_TAIL:]
    c1_rep = np.asarray(c1_rep, bool)[-REPAIR_TAIL:]

    persist = evaluate_persistent(d_ss, c1_ss, n_windows=3, min_fire=2)
    scan = granularity_scan(d_ss, max_level=7)
    repair = repair_encapsulated(d_rep, c1_rep)
    return {"name": name, "d_ss": d_ss, "c1_ss": c1_ss,
            "persist": persist, "scan": scan, "repair": repair,
            "stats": window_stats(d_ss), "kernel": kernel}


# ---------------------------------------------------------------- 主流程
def main() -> dict:
    print("=" * 104)
    print("R6 · A10 可观测现象清单：预注册判据（运行前写定，阈值不随结果调整）")
    print("=" * 104)
    print(f"  申报精度 Δ_quoted = {DELTA_QUOTED:g}   |   绝对阈值 Δ_abs = 3Δ_quoted = {DELTA_ABS:g}")
    print(f"  标度律判据：beta >= {BETA_MAX:g} 即「无特征尺度」（残差不随粒度变粗而衰减）")
    print()
    print(f"{'code':5s}{'判据':34s}{'规则':16s}{'门':4s}")
    print("-" * 104)
    for row in criterion_table():
        print(f"{row['code']:5s}{row['name']:34s}{row['rule']:16s}{row['gate']:4s}")
    print("-" * 104)
    print("  已知盲区（预注册）：P2–P5 均为位置不变量 → 系统性偏置不可见，须靠跨实现一致率（C2）。")
    print()

    results = {}
    for name in SCENARIOS:
        results[name] = run_scenario(name)

    # ---------------- 表 1：逐场景统计
    print("=" * 104)
    print("表 1  逐场景判据统计（稳态窗口，绝对标度）")
    print("=" * 104)
    print(f"{'scenario':10s}{'C1占比':>8s}{'sd':>8s}{'IQR':>8s}{'BC':>8s}"
          f"{'尾质量':>8s}{'兑现比':>8s}{'beta':>8s}{'修接口恢复':>11s}")
    print("-" * 104)
    for name, r in results.items():
        s, sc = r["stats"], r["scan"]
        print(f"{name:10s}{r['persist']['c1_fraction']:8.2f}{s['sd']:8.4f}{s['iqr']:8.4f}"
              f"{s['bc']:8.3f}{s['tail_abs']:8.3f}{s['commit']:8.2f}{sc['beta']:8.3f}"
              f"{str(r['repair']['encapsulated']):>11s}")
    print("-" * 104)

    # ---------------- 表 2：判据触发矩阵 + 三联排除
    print()
    print("=" * 104)
    print("表 2  判据触发矩阵（3 窗口中 ≥2 触发记为 ●，不可评估记为 –）与三联排除裁决")
    print("=" * 104)
    codes = [c.code for c in CRITERIA]
    print(f"{'scenario':10s}" + "".join(f"{c:>6s}" for c in codes)
          + f"{'  裁决':4s}  说明")
    print("-" * 104)
    verdicts = {}
    for name, r in results.items():
        p = r["persist"]
        cells = ""
        for c in codes:
            f = p["fires"][c]
            cells += f"{'●' if f else ('·' if f is False else '–'):>6s}"
        v = classify(p, r["repair"], r["scan"], order="repair_first")
        verdicts[name] = v
        print(f"{name:10s}{cells}  {v['verdict']:<8s}{VERDICTS[v['verdict']]}")
    print("-" * 104)

    # ---------------- 表 3：顺序置换检验
    print()
    print("=" * 104)
    print("表 3  三联排除的顺序是不是 load-bearing？（repair_first vs granularity_first）")
    print("=" * 104)
    print(f"{'scenario':10s}{'repair_first':>16s}{'granularity_first':>20s}{'一致':>8s}")
    print("-" * 104)
    flips = 0
    for name, r in results.items():
        a = classify(r["persist"], r["repair"], r["scan"], order="repair_first")["verdict"]
        b = classify(r["persist"], r["repair"], r["scan"], order="granularity_first")["verdict"]
        flips += int(a != b)
        print(f"{name:10s}{a:>16s}{b:>20s}{('是' if a == b else '否 ← 顺序敏感'):>8s}")
    print("-" * 104)
    print(f"  顺序敏感场景数 = {flips}。若 >0，说明「先修接口、再看标度律」不是排版方便，")
    print("  而是**必要条件**：系统性接口伪影在标度律上与真坍塌同形（beta 都 ≈ 0），")
    print("  只有接口修复能把二者分开。")

    # ---------------- 表 4：判别力边界
    print()
    print("=" * 104)
    print("表 4  P5 的判别力边界：有限相关时间 τ 从 0.4 步扫到 80 步")
    print("=" * 104)
    print(f"{'τ (s)':>8s}{'τ (步)':>9s}{'beta':>9s}{'P5 触发':>9s}   判读")
    print("-" * 104)
    sweep = sweep_scale_separation()
    for row in sweep:
        if row["tau_steps"] > 128.0:
            tag = "尺度超出扫描窗 → 与无尺度不可分（保守触发）"
        elif row["fires_P5"]:
            tag = "窗内偏上：有限尺度被**误判为坍塌**"
        elif row["tau_steps"] < 2.0:
            tag = "尺度低于分辨率 → 正确不触发"
        else:
            tag = "正确不触发，但余量随 τ 增大而收窄"
        print(f"{row['tau']:8.2f}{row['tau_steps']:9.1f}{row['beta']:9.3f}"
              f"{('是' if row['fires_P5'] else '否'):>9s}   {tag}")
    print("-" * 104)
    flip = next((r["tau_steps"] for r in sweep if r["fires_P5"]), float("nan"))
    print(f"  判读：β 在 τ≈{flip:.0f} 步附近穿过阈值 {BETA_MAX:g}，且是**连续过渡**，")
    print("        不是两类机制之间的跃变。因此 P5 的真实身份是")
    print("        「特征尺度是否落在扫描窗下部」的探测器，而非「有无特征尺度」的判据。")
    print("        后果：τ 落在窗内偏上（约 25–128 步）的**有限尺度过程会被误判为坍塌**；")
    print("        τ 超出窗上界则与无尺度不可分（二者都触发，属保守方向）。")
    print("        使用纪律：先把扫描范围对齐到疑似尺度区间，再读 β，否则判据无效。")

    # ---------------- 图
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    fig.suptitle("R6 · A10 pre-registered observables: can D-form criteria separate "
                 "reef / resolution / collapse?", fontsize=13.5)

    # (1) ridgeline: D centered at median (absolute scale)
    ax = axes[0, 0]
    edges = np.linspace(-0.30, 0.30, 61)
    for i, (name, r) in enumerate(results.items()):
        x = r["d_ss"] - np.median(r["d_ss"])
        h, _ = np.histogram(x, bins=edges, density=True)
        h = h / max(h.max(), 1e-12) * 0.85
        c = edges[:-1]
        ax.fill_between(c, i, i + h, color=COLORS[name], alpha=0.55, lw=0)
        ax.plot(c, i + h, color=COLORS[name], lw=0.9)
        ax.text(-0.295, i + 0.18, f"{name}  sd={r['stats']['sd']:.3f}  "
                                  f"BC={r['stats']['bc']:.2f}", fontsize=7.6, color="#333333")
    ax.axvline(DELTA_ABS, color="#888888", ls="--", lw=0.9)
    ax.axvline(-DELTA_ABS, color="#888888", ls="--", lw=0.9, label=f"±Δ_abs = ±{DELTA_ABS:g}")
    ax.set_yticks([]); ax.set_xlabel("d − median(D)   [absolute scale]")
    ax.set_title("(1) D shape per scenario (median-centered, absolute scale)", fontsize=10.5, loc="left")
    ax.legend(fontsize=7.5, loc="upper right", frameon=False)
    ax.grid(alpha=0.2, axis="x")

    # (2) granularity scan（归一化：斜率就是判据本身）
    ax = axes[0, 1]
    for row in sweep:
        ax.plot(row["b"], row["a_norm"], "-", color="#b3b3b3", lw=0.7, alpha=0.7)
    ax.plot([], [], "-", color="#b3b3b3", lw=0.7, label="OU sweep τ=0.02…4 s")
    for name, r in results.items():
        sc = r["scan"]
        a0 = sc["a"][0]
        ax.plot(sc["b"], [v / a0 for v in sc["a"]], "o-", ms=3.5, lw=1.3,
                color=COLORS[name], label=f"{name}  β={sc['beta']:+.2f}")
    bb = np.asarray(results["OK"]["scan"]["b"], float)
    ax.plot(bb, bb ** BETA_MAX, "--", color="#c0392b", lw=1.2,
            label=f"critical slope β={BETA_MAX:g}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("description granularity B  [block length, samples]")
    ax.set_ylabel("a(B) / a(1)   (normalized residual amplitude)")
    ax.set_title("(2) granularity scan: above the red line → declared scale-free → collapse",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=7.0, frameon=False, ncol=2, loc="lower left")
    ax.grid(alpha=0.25, which="both")

    # (3) criteria matrix
    ax = axes[1, 0]
    for j, name in enumerate(results):
        p = results[name]["persist"]
        for i, c in enumerate(codes):
            f = p["fires"][c]
            val = 1.0 if f else (0.0 if f is False else 0.35)
            ax.add_patch(plt.Rectangle((i, j), 0.94, 0.94,
                                       color=("#c0392b" if f else
                                              ("#ecf0f1" if f is False else "#f7dc6f")),
                                       alpha=0.85, lw=0))
            ax.text(i + 0.47, j + 0.47, "●" if f else ("·" if f is False else "–"),
                    ha="center", va="center", fontsize=12,
                    color="white" if f else "#555555")
    ax.set_xlim(0, len(codes)); ax.set_ylim(0, len(results))
    ax.set_xticks([i + 0.47 for i in range(len(codes))]); ax.set_xticklabels(codes)
    ax.set_yticks([j + 0.47 for j in range(len(results))]); ax.set_yticklabels(list(results))
    ax.invert_yaxis()
    ax.set_title("(3) criterion firing matrix (● fires in ≥2/3 windows, – not evaluable)",
                 fontsize=10.5, loc="left")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    # (4) verdict table
    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("(4) elimination ladder verdict", fontsize=10.5, loc="left")
    col = {"ok": "#d5f5e3", "gate": "#d6eaf8", "reef": "#fdebd0",
           "resolution": "#e8daef", "collapse": "#f5b7b1"}
    rows = []
    for name, r in results.items():
        v = verdicts[name]
        rows.append((name, " → ".join(v["steps"]), v["verdict"]))
    y = 0.92
    ax.text(0.02, y + 0.06, "scenario", fontsize=8.6, weight="bold")
    ax.text(0.20, y + 0.06, "ladder steps", fontsize=8.6, weight="bold")
    ax.text(0.72, y + 0.06, "verdict", fontsize=8.6, weight="bold")
    for name, steps, vd in rows:
        ax.add_patch(plt.Rectangle((0.70, y - 0.035), 0.29, 0.085,
                                   color=col.get(vd, "#eeeeee"), transform=ax.transAxes, lw=0))
        ax.text(0.02, y, name, fontsize=8.4)
        ax.text(0.20, y, steps, fontsize=7.4, color="#333333")
        ax.text(0.72, y, vd, fontsize=8.4)
        y -= 0.135
    ax.text(0.02, y - 0.02,
            "BIAS is the pre-registered blind spot: all four criteria are\n"
            "location-invariant, so a pure bias is declared M-complete\n"
            "(false negative). Needs cross-route agreement, not D shape.",
            fontsize=7.8, color="#7b241c", va="top")

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "r6_a10_observables.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print()
    print(f"[ok] figure saved: {out}")

    return {"results": results, "verdicts": verdicts, "flips": flips}


if __name__ == "__main__":
    main()
