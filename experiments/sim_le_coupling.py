"""
X-S-MERA V3.4 · LE Life-Environment coupling simulation
=======================================================

Depend: src/le_coupling_kernel.py

四张子图：
1. L / E 状态时序（生命、环境子系统状态分量）
2. d（边界重合度标量）时序 + C1 判定阈值
3. Var(D)：两个子系统 D 分布方差（A10 边界封装 / M-完备性预警指标）+ ghost 伪力幅值
4. 三区判定与 F1/F2 失效标记

阶段设计（注意顺序，见 docs/experiment_manual.md §5.3）：
    S1 0–15s   Φ=0.40, C2=兼容   → high 稳定耦合
    S2 15–30s  Φ=0.45, C2=不兼容 → middle，F2 ghost 伪力（d 仍高，Var(D) 抬升）
    S3 30–45s  Φ=0.08, C2=兼容   → 跌出 C1，F1（无直接交互）
    S4 45–60s  Φ=0.40, C2=兼容   → 恢复尝试：**失败**（F1 是吸收态，滞后不可逆）
先 F2 后 F1 是必须的：F1 一旦成立，d 已低于阈值，后续再关 C2 也只会被判定为 F1，
F2 数据在动力学上不可能产出（原方案的阶段顺序有误）。

图像自动保存至 figs/le_coupling_timeseries.png

运行：
    python experiments/sim_le_coupling.py
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from src.le_coupling_kernel import LECouplingKernel, KERNEL_VERSION  # noqa: E402

FIG_DIR = os.path.join(_REPO_ROOT, "figs")

# (t_start, t_end, semantic_potential_grad, c2_compatibility)
STAGES = [
    (0.0, 15.0, 0.40, True),    # S1 稳定耦合
    (15.0, 30.0, 0.45, False),  # S2 F2 ghost 伪力（礁石伪收敛）
    (30.0, 45.0, 0.08, True),   # S3 F1 跌出重合区
    (45.0, 60.0, 0.40, True),   # S4 恢复尝试（滞后：失败）
]

STAGE_LABEL = ["S1 stable (C1+C2)", "S2 F2: ghost force",
               "S3 F1: overlap collapse", "S4 recovery attempt"]
STAGE_BG = ["#eaf5ee", "#fdf3e3", "#fdecea", "#f2f3f4"]

REGION_COLOR = {"high": "#1f7a45", "middle": "#d68910", "low": "#c0392b"}


def _apply_stage(kernel: LECouplingKernel, t: float) -> None:
    for t0, t1, grad, c2 in STAGES:
        if t0 <= t < t1:
            kernel.semantic_pot_grad = grad
            kernel.c2_compatibility = c2
            return


def run_le_simulation(t_total: float = 60.0, dt: float = 0.05, seed: int = 7) -> LECouplingKernel:
    kernel = LECouplingKernel(
        dim_life=3,
        dim_env=3,
        d_threshold_c1=0.25,
        c2_compatibility=True,
        semantic_potential_grad=0.40,
        coupling_gain=6.0,
        seed=seed,
    )

    for i in range(int(round(t_total / dt))):
        t = i * dt
        _apply_stage(kernel, t)
        kernel.step(t, dt)
    return kernel


def plot_timeseries(kernel: LECouplingKernel) -> str:
    os.makedirs(FIG_DIR, exist_ok=True)
    h = kernel.history
    t_arr = np.asarray(h["t"])
    d_arr = np.asarray(h["d"])
    var_L = np.asarray(h["D_L_var"])
    var_E = np.asarray(h["D_E_var"])
    ghost = np.asarray(h["ghost_magnitude"])
    L0 = np.asarray([s[0] for s in h["L_state"]])
    E0 = np.asarray([s[0] for s in h["E_state"]])
    L1 = np.asarray([s[1] for s in h["L_state"]])
    E1 = np.asarray([s[1] for s in h["E_state"]])
    regions = h["region_label"]
    fails = h["failure_mode"]
    enc_L = np.asarray(h["L_encapsulated"])
    enc_E = np.asarray(h["E_encapsulated"])

    fig, axes = plt.subplots(4, 1, figsize=(13.5, 12.5), sharex=True)
    fig.suptitle(f"X-S-MERA {KERNEL_VERSION} · LE Life-Environment Coupling Dynamics",
                 fontsize=15, y=0.985)

    for ax in axes:
        for i, (t0, t1, _, _) in enumerate(STAGES):
            ax.axvspan(t0, t1, color=STAGE_BG[i], alpha=0.7, lw=0, zorder=0)
            if i > 0:
                ax.axvline(t0, color="#b0b0b0", lw=0.8, ls=":", zorder=1)
        ax.grid(alpha=0.25, zorder=2)

    # ---- (1) 状态 ----
    axes[0].plot(t_arr, L0, color="#1f5f99", lw=1.8, label="Life L · state[0] (biomass)", zorder=3)
    axes[0].plot(t_arr, E0, color="#b03a2e", lw=1.8, label="Env E · state[0] (resource)", zorder=3)
    axes[0].plot(t_arr, L1, color="#5d8fc0", lw=1.1, ls="--", label="Life L · state[1] (damage)", zorder=3)
    axes[0].plot(t_arr, E1, color="#c98a83", lw=1.1, ls="--", label="Env E · state[1] (gradient)", zorder=3)
    axes[0].set_ylabel("State value")
    axes[0].set_title("(1) L / E subsystem state  —  E·state[1] collapse in S2 is a ghost-force artifact",
                      fontsize=11, loc="left")
    axes[0].legend(loc="upper right", fontsize=8.5, frameon=False, ncol=2)

    # ---- (2) d ----
    axes[1].plot(t_arr, d_arr, color="#1e8449", lw=1.9, label="d (boundary overlap)", zorder=3)
    axes[1].axhline(kernel.d_threshold_c1, color="#ca6f1e", ls="--", lw=1.3,
                    label=f"C1 threshold = {kernel.d_threshold_c1}", zorder=3)
    axes[1].fill_between(t_arr, 0, d_arr, where=(d_arr >= kernel.d_threshold_c1),
                         color="#1e8449", alpha=0.10, zorder=2)
    axes[1].set_ylabel("d")
    axes[1].set_ylim(0, 1.04)
    axes[1].set_title("(2) Boundary overlap d — scale-free overlap, NOT a length-scale ratio",
                      fontsize=11, loc="left")
    axes[1].legend(loc="upper right", fontsize=9, frameon=False)

    # ---- (3) Var(D) + ghost ----
    axes[2].plot(t_arr, var_L, color="#6c3483", lw=1.7, label="Var(D_L)", zorder=3)
    axes[2].plot(t_arr, var_E, color="#b9770e", lw=1.7, label="Var(D_E)", zorder=3)
    axes[2].axhline(kernel.d_var_max, color="#7b241c", ls=":", lw=1.4,
                    label=f"encapsulation bound = {kernel.d_var_max}", zorder=3)
    bad = ~(enc_L & enc_E)
    if bad.any():
        axes[2].fill_between(t_arr, 0, max(var_L.max(), var_E.max()) * 1.1, where=bad,
                             color="#c0392b", alpha=0.07, zorder=1,
                             label="A10 encapsulation failure (M-incomplete)")
    axes[2].set_ylabel("Var(D)")
    axes[2].set_title("(3) Var(D) — A10 boundary-encapsulation / M-completeness indicator",
                      fontsize=11, loc="left")

    ax3 = axes[2].twinx()
    ax3.plot(t_arr, ghost, color="#7f8c8d", lw=1.2, ls="-.",
             label="ghost force magnitude", zorder=3)
    ax3.set_ylabel("|ghost flux|", color="#5f6a6a")
    ax3.tick_params(axis="y", colors="#5f6a6a")
    ax3.set_ylim(0, max(ghost.max() * 1.6, 0.1))

    h1, l1 = axes[2].get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    axes[2].legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8.5, frameon=False)

    # ---- (4) 三区与失效 ----
    for tt, r, f in zip(t_arr, regions, fails):
        y = 1 if f == "F1" else (2 if f == "F2" else 3)
        axes[3].scatter(tt, y, s=10, color=REGION_COLOR[r], marker="s", lw=0, zorder=3)
    axes[3].set_yticks([1, 2, 3])
    axes[3].set_yticklabels(["F1  no overlap\n(flux = 0)",
                             "F2  ghost force\n(reef pseudo-convergence)",
                             "C1+C2  stable\ncoupling"])
    axes[3].set_ylim(0.4, 3.6)
    axes[3].set_ylabel("Regime")
    axes[3].set_xlabel("t")
    axes[3].set_xlim(t_arr[0], t_arr[-1])
    axes[3].set_title("(4) Region and failure mode — S4 cannot return to high: F1 is an absorbing state",
                      fontsize=11, loc="left")

    for i, (t0, t1, _, _) in enumerate(STAGES):
        axes[3].annotate(STAGE_LABEL[i], xy=((t0 + t1) / 2, 3.42), ha="center", va="center",
                         fontsize=8.5, color="#4a4a4a")

    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = os.path.join(FIG_DIR, "le_coupling_timeseries.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def summarize(kernel: LECouplingKernel) -> None:
    h = kernel.history
    t_arr = np.asarray(h["t"])
    d_arr = np.asarray(h["d"])
    var_L = np.asarray(h["D_L_var"])
    var_E = np.asarray(h["D_E_var"])
    ghost = np.asarray(h["ghost_magnitude"])
    n = t_arr.size

    print("\n==== Simulation summary ====")
    print(f"steps = {n}, dt = {t_arr[1] - t_arr[0]:.3f}, total t = {t_arr[-1]:.2f}")
    print(f"d: min={d_arr.min():.4f}  max={d_arr.max():.4f}  final={d_arr[-1]:.4f}")
    print(f"Var(D_L): max={var_L.max():.5f}  final={var_L[-1]:.5f}")
    print(f"Var(D_E): max={var_E.max():.5f}  final={var_E[-1]:.5f}")
    print(f"ghost flux magnitude: max={ghost.max():.4f}")

    for region in ("high", "middle", "low"):
        mask = np.array(h["region_label"]) == region
        print(f"time in {region:6s} region: {mask.mean() * 100:5.1f}%  ({mask.sum()} steps)")

    print(f"A10 encapsulation failures: L={int((~np.asarray(h['L_encapsulated'])).sum())} steps, "
          f"E={int((~np.asarray(h['E_encapsulated'])).sum())} steps")

    print("\n-- per stage --")
    for i, (t0, t1, grad, c2) in enumerate(STAGES):
        mask = (t_arr >= t0) & (t_arr < t1)
        if not mask.any():
            continue
        reg = np.array(h["region_label"])[mask]
        dom = max(set(reg), key=list(reg).count)
        fm = np.array([str(x) for x in np.array(h["failure_mode"])[mask]])
        modes = dict(zip(*np.unique(fm, return_counts=True)))
        print(f"  {STAGE_LABEL[i]:28s} grad={grad:.2f} C2={str(c2):5s} "
              f"mean d={d_arr[mask].mean():.3f} Var(D_L)={var_L[mask].mean():.5f} "
              f"ghost={ghost[mask].mean():.3f} dominant={dom:6s} {modes}")

    print("\n-- 判读 --")
    print("  S1: d 高、Var(D) 低、ghost=0     → C1+C2 成立，边界张力被封裝为边界宏（A10 不计入 M-不完备）")
    print("  S2: d 仍高（C1 成立）、ghost>0、Var(D) 越界 → F2 礁石伪收敛：伪力把副模态压到 0，"
          "但 D 已偏离可封装形态 = M-不完备真信号")
    print("  S3: d 跌破阈值 → F1，通量置零（无直接交互）")
    print("  S4: 势梯度恢复至 0.40 仍无法重建耦合 → F1 是吸收态：d 出现在通量里，d≈0 使通路自锁")


if __name__ == "__main__":
    kern = run_le_simulation()
    path = plot_timeseries(kern)
    print(f"[ok] figure saved: {path}")
    summarize(kern)
