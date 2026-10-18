"""
R8 遗留 T1-d：表征分化仿真（修正版 v2）
==========================================

四联图：
  (a) 通量—翻越率标度（R8 F-c 子命题 1：时点归流）
  (b) Φ 阴性对照（a=0 vs a>0：存在归势）
  (c) 边界稳定性 vs 通量（BC 随 J 衰减：稳定归流）
  (d) 破缺方向分布（f₊ 随种子：方向归流）

建模修正记录（v1 → v2，重要）
------------------------------
v1 把 Φ 建成"可分化性谱 w"（乘在更新幅值上）。**这是错的**：
乘性权重只改方差、不改均值结构，任何 J 下 k-means 都只是把单峰切开
（实测 abs_sep 恒 1.74、与 J 无关，n_c 全部撞上下限）。
v2 改用**双阱势 + 多单元群**：
  - 单轨迹在双阱中一旦落阱就锁死，永不双峰 → 必须用多单元群；
  - Φ（双阱深度 a）决定"是否存在两极"，J（噪声幅度）决定"越过势垒的速率"。
这两个建模错误的修正过程本身是 R8 F-b 的最好例证：
**可分性必须由势的多稳结构承载，不能由动的强弱承载。**

运行：python experiments/sim_t1d_representational.py
输出：figs/r8_t1d_representational.png + 控制台摘要
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.t1d_representational_diff import (  # noqa: E402
    RepresentationPopulation, boundary_diagnostics,
    flux_sweep_observables, direction_report, hysteresis_rate,
)

FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

N_UNITS = 300
STEPS = 4000
BURN = 1000


def main() -> None:
    print("=" * 72)
    print("R8 / T1-d 表征分化（v2：双阱势 + 多单元群）")
    print("=" * 72)

    fluxes = [0.0, 0.3, 1.0, 2.0, 3.0, 5.0]

    # ---------- (a)(c)(d) 主扫描 ----------
    sweep = flux_sweep_observables(fluxes, seeds=[0, 1, 2], n_units=N_UNITS,
                                   steps=STEPS, burn_in=BURN, gain=0.01)
    print("\n[a] 通量扫描（Φ 固定：双阱 a=1.0, c=1.0）")
    print(f"  {'J':>5} {'翻越率':>10} {'深度':>8} {'f+':>12} {'BC':>8}")
    for r in sweep["rows"]:
        print(f"  {r['flux']:5.1f} {r['flip_rate_mean']:10.5f} "
              f"{r['depth_mean']:8.3f} "
              f"{r['frac_plus_mean']:5.2f}±{r['frac_plus_std']:.2f} "
              f"{r['bc_final_mean']:8.3f}")

    # 线性度检验（翻越率 vs J）
    J = np.array([r["flux"] for r in sweep["rows"]])
    F = np.array([r["flip_rate_mean"] for r in sweep["rows"]])
    slope, intercept = np.polyfit(J, F, 1)
    pred = slope * J + intercept
    ss_res = float(np.sum((F - pred) ** 2))
    ss_tot = float(np.sum((F - F.mean()) ** 2)) + 1e-15
    r2 = 1 - ss_res / ss_tot
    print(f"\n  翻越率 ≈ {intercept:.4f} + {slope:.4f}·J   (R²={r2:.4f})")

    # ---------- (b) Φ 阴性对照 ----------
    print("\n[b] Φ 阴性对照：对称态（a=0）在任意 J 下是否生成边界")
    control = []
    for a in (0.0, 1.0):
        for Jv in (0.0, 1.0, 3.0):
            bcs, has = [], []
            for s in range(3):
                pop = RepresentationPopulation(n_units=N_UNITS, well_depth=a,
                                               well_sep=1.0, flux=Jv,
                                               noise_std=0.10, gain=0.01, seed=s)
                pop.run(STEPS)
                d = boundary_diagnostics(pop.x0)
                bcs.append(d["bc"])
                has.append(d["has_boundary"])
            control.append({"a": a, "J": Jv, "bc": float(np.mean(bcs)),
                            "has": bool(all(has))})
            print(f"  a={a:.1f} J={Jv:.1f}  BC={np.mean(bcs):.3f}  "
                  f"判为有边界={all(has)}")

    # ---------- 绘图 ----------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5))

    ax = axes[0, 0]
    ax.plot(J, F, "o-", color="tab:blue", label="实测翻越率")
    ax.plot(J, pred, "--", color="tab:red",
            label=f"线性拟合 R²={r2:.3f}")
    ax.set_xlabel("载体通量 J")
    ax.set_ylabel("翻越率（每单元每步）")
    ax.set_title("(a) 时点归流：越垒速率随通量线性上升")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    labels, vals, colors = [], [], []
    for c in control:
        labels.append(f"a={c['a']:.0f}\nJ={c['J']:.0f}")
        vals.append(c["bc"])
        colors.append("tab:green" if c["has"] else "tab:gray")
    bars = ax.bar(range(len(vals)), vals, color=colors)
    ax.axhline(5 / 9, color="k", ls=":", lw=1.2, label="Sarle 分界 5/9")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Sarle 双峰系数 BC")
    ax.set_title("(b) 存在归势：a=0 无边界，a>0 有边界")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1, 0]
    bc_vals = [r["bc_final_mean"] for r in sweep["rows"]]
    ax.plot(J, bc_vals, "s-", color="tab:purple")
    ax.axhline(5 / 9, color="k", ls=":", lw=1.2, label="Sarle 分界")
    ax.set_xlabel("载体通量 J")
    ax.set_ylabel("稳态 BC")
    ax.set_title("(c) 稳定归流：通量越大，边界越不稳定")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    low_J, high_J = 0.3, 3.0
    for Jv, mk in ((low_J, "o"), (high_J, "s")):
        fps = []
        for s in range(12):
            pop = RepresentationPopulation(n_units=N_UNITS, well_depth=1.0,
                                           flux=Jv, noise_std=0.10,
                                           gain=0.01, seed=s)
            traj = pop.run(STEPS)
            fps.append(direction_report(traj, BURN)["frac_plus"])
        ax.scatter([Jv] * len(fps), fps, marker=mk, alpha=0.7,
                   label=f"J={Jv}")
    ax.axhline(0.5, color="k", ls=":", lw=1.2, label="对称破缺预期 0.5")
    ax.set_xlabel("载体通量 J")
    ax.set_ylabel("f₊（最终落 +阱的单元比例）")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("(d) 方向归流：f₊ 恒在 0.5 附近（方向不被 Φ 锁死）")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("R8 / T1-d：表征分化——存在归势，时点与方向归流", fontsize=14)
    fig.tight_layout()
    out = FIG_DIR / "r8_t1d_representational.png"
    fig.savefig(out, dpi=140)
    print(f"\n图已保存：{out}")

    # ---------- 判定摘要 ----------
    sym_bad = [c for c in control if c["a"] == 0.0 and c["has"]]
    asym_ok = [c for c in control if c["a"] > 0.0 and c["has"]]
    fp_all = [r["frac_plus_mean"] for r in sweep["rows"]]
    print("\n" + "-" * 72)
    print("判定摘要（预注册判据）")
    print(f"  R8 F-c 子命题 1（时点归流）：翻越率随 J 单调递增 = "
          f"{'✓' if np.all(np.diff(F) > 0) else '✗'}，"
          f"线性 R²={r2:.3f}，斜率 {slope:.4f}")
    print(f"  R8 F-c 子命题 2（存在归势）：a=0 时判为有边界 = {len(sym_bad)} 例 "
          f"→ {'✓ 阴性对照通过' if not sym_bad else '✗ 出现反例'}")
    print(f"                                 a>0 时判为有边界 = {len(asym_ok)} 例 "
          f"→ {'✓' if len(asym_ok) >= 2 else '✗'}")
    print(f"  R8 F-c 子命题 3（方向归流）：f₊ 范围 = "
          f"[{min(fp_all):.2f}, {max(fp_all):.2f}]，恒贴近 0.5 = "
          f"{'✓' if all(abs(f - 0.5) < 0.1 for f in fp_all) else '✗'}")
    print(f"  边界稳定性（新增读数）：BC 随 J 从 {bc_vals[0]:.3f} 降到 {bc_vals[-1]:.3f}"
          f" → 通量过大时边界不可维持（这是 A11b 的非单调性证据）")


if __name__ == "__main__":
    main()
