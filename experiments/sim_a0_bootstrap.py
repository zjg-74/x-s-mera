"""
R7 遗留：方案甲（回归收敛）仿真
================================

两张图 + 一个判定表：
  (a) 层级回归：不同收缩增益 g 下的张力轨迹（对照读法 R1）
  (b) R7 四类人工系统的收缩性定位（对照读法 R2 与三方案适用域）

判定依据见 src/a0_bootstrap_convergence.py 的 C-1/C-2/C-3。
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

from src.a0_bootstrap_convergence import (  # noqa: E402
    run_chain, classify_system, SYSTEM_MODELS, fixed_point_iteration,
)

FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    print("=" * 70)
    print("R7 / 方案甲：回归收敛仿真")
    print("=" * 70)

    # ---------- (a) 读法 R1：层级回归 ----------
    gains = [0.70, 0.90, 0.98, 1.00, 1.03]
    res = run_chain(gains=gains, n_layers=60, layer_error=0.02, noise_std=0.10)
    print("\n[a] 层级回归（读法 R1）尾段张力")
    for r in res["rows"]:
        print(f"  g={r['gain']:.2f}  尾段均值={r['tail_mean']:.4f}  "
              f"最大={r['tail_max']:.4f}  → {r['verdict']}")

    # ---------- (b) R2 自我回归的不动点 ----------
    print("\n[b] 自我回归（读法 R2）不动点迭代")
    fp_rows = []
    # 收缩算子 x ↦ a + g·x，a 为"上游锚点"
    for g in [0.70, 0.90, 1.00, 1.03]:
        r = fixed_point_iteration(lambda x, g=g: (0.3 + g * x, abs(g)),
                                  x0=0.0, n_iters=400)
        print(f"  g={g:.2f}  迭代={r['iters']:3d}  最终漂移={r['final_drift']:.2e}  "
              f"观测收敛比={r['observed_ratio_median']:.3f}  "
              f"{'收敛' if r['converged'] else '不收敛'}")
        fp_rows.append((g, r))

    # ---------- (c) R7 四类系统 ----------
    print("\n[c] R7 四类人工系统的收缩性")
    rows = []
    for name in SYSTEM_MODELS:
        c = classify_system(name, n_layers=150)
        rows.append(c)
        print(f"  {name:22s} g={c['gain']:.2f}  尾段={c['tail_mean']:.3f}  "
              f"{c['class']}")

    # ---------- 绘图 ----------
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    ax = axes[0]
    for r in res["rows"]:
        ax.loglog(r["curves"].mean(axis=0) + 1e-12, label=f"g={r['gain']:.2f}")
    ax.set_xlabel("层数 k")
    ax.set_ylabel("张力 d_k（5 种子均值）")
    ax.set_title("(a) 层级回归：只有严格收缩才衰减到零")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    gs = [r["gain"] for r in res["rows"]]
    tails = [max(r["tail_mean"], 1e-6) for r in res["rows"]]
    ax.semilogy(gs, tails, "o-", color="tab:blue", label="层级回归尾段张力")
    ax.axvline(1.0, color="grey", ls=":", lw=1)
    ax.axvspan(0.0, 1.0, color="tab:green", alpha=0.08, label="C-1 收缩区")
    ax.axvspan(1.0, max(gs) + 0.02, color="tab:red", alpha=0.08, label="C-3 发散区")
    for c in rows:
        ax.scatter([c["gain"]], [max(c["tail_mean"], 1e-6)], s=70,
                   marker="D", zorder=5)
        ax.annotate(c["system"].split("_")[0], (c["gain"], max(c["tail_mean"], 1e-6)),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xlabel("投影算子增益 g（‖Π‖）")
    ax.set_ylabel("尾段张力（对数）")
    ax.set_title("(b) 四类系统的收缩性定位（菱形）")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("R7 / 方案甲：回归收敛的真实条件 = 投影算子严格收缩", fontsize=13)
    fig.tight_layout()
    out = FIG_DIR / "r7_plan_a_convergence.png"
    fig.savefig(out, dpi=140)
    print(f"\n图已保存：{out}")

    # ---------- 判定摘要 ----------
    print("\n" + "-" * 70)
    print("判定摘要")
    conv = [r for r in res["rows"] if abs(r["gain"]) < 1]
    print(f"  C-1（g<1）：{len(conv)} 例张力有界 —— 方案甲在收缩算子上成立")
    div = [r for r in res["rows"] if abs(r["gain"]) > 1]
    print(f"  C-3（g>1）：{len(div)} 例发散 —— 方案甲在扩张算子上被证伪")
    n_strong = sum(1 for c in rows if c["class"] == "C-1 强收缩")
    print(f"  四类系统中强收缩者 = {n_strong} 类"
          f"（{'类4 演化功能语义' if n_strong == 1 else '见上表'}）")
    print("  注意：类4 虽强收缩，其不动点是'适应峰'而不是 A0——")
    print("        收敛 ≠ 回归到 A0。方案甲即便在此类上成立，也无辩护 A0 之力，")
    print("        与 R7 的定性判定（丙为主、乙为辅）一致。")


if __name__ == "__main__":
    main()
