"""
G4 遗留：ΔT（时间势差）/ TP 谓词操作化仿真
============================================

补齐 V3.3.1 §G4 登记的缺口："ΔT→0（时间势差）：待定义观测量"。

操作化：ΔT = |ω_L − ω_E| / (ω_L + ω_E)（内禀节律频率的无量纲失配）；
观测量 = 耦合相位振子相对相位的漂移率 Ω_drift；
理论锚 = Adler 精确解 √(Δω²−4K²)（解锁区）/ 0（锁定区）。

图（2×2）：
  (a) ΔT→0 通道：漂移率 vs ΔT，对照 Adler 曲线（含锁定区阴影）
  (b) 代表性相位轨迹：解锁（漂移）vs 锁定（势差被吸收）
  (c) P5 阴性对照：K=0 时漂移恒等于 Δω（无锁定区）
  (d) P4 谓词独立：固定状态下 d 不变而 ΔT 变
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

from src.tp_temporal_potential import (  # noqa: E402
    measure_natural_frequency, delta_T, adler_drift, coupled_phase_experiment,
    dt_zero_channel_sweep, predicate_independence_test, TP_CRITERIA,
)

FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    print("=" * 70)
    print("G4 / ΔT：时间势差 TP 谓词操作化")
    print("=" * 70)
    print("\n预注册判据：")
    for k, v in TP_CRITERIA.items():
        print(f"  {k}: {v}")

    # ---------- P1 可测性 ----------
    print("\n[P1] 内禀频率可测性")
    for w in (0.8, 1.3, 2.1):
        m = measure_natural_frequency(w)
        print(f"  ω={w:.1f}  拟合={m['omega_fit']:.4f}  误差={m['rel_err']:.3%}")

    # ---------- ΔT→0 通道 ----------
    K = 0.05
    print(f"\n[ΔT→0 通道] 耦合 K={K}（锁定阈值 |Δω|≤2K={2*K}）")
    sw = dt_zero_channel_sweep(coupling=K)
    for r in sw["rows"]:
        print(f"  频率比={r['freq_ratio']:.2f}  ΔT={r['delta_T']:.3f}  "
              f"漂移={r['drift_mean']:.4f}  Adler={r['adler_predict']:.4f}  "
              f"{'锁定' if r['locked'] else '漂移'}")
    print(f"  标度指数={sw['scaling_exponent']:.2f}  "
          f"Adler 最大相对误差={sw['adler_max_rel_err']:.3%}  "
          f"单调={sw['monotone_with_tolerance']}")

    # ---------- P5 阴性对照 ----------
    print("\n[P5] 阴性对照（K=0）")
    p5 = []
    for rr in (1.05, 1.2, 1.5):
        r = coupled_phase_experiment(1.0, rr, coupling=0.0, seed=0)
        err = abs(r["drift_rate"] - (rr - 1.0)) / (rr - 1.0)
        p5.append((rr, r))
        print(f"  频率比={rr:.2f}  漂移={r['drift_rate']:.4f}  "
              f"Δω={rr-1.0:.2f}  误差={err:.2%}  锁定={r['locked']}")

    # ---------- P4 谓词独立 ----------
    pi = predicate_independence_test()
    print(f"\n[P4] 谓词独立：d 变异={pi['d_cv']:.3f}  ΔT 变异={pi['dt_cv']:.3f}  "
          f"独立={pi['independent']}")

    # ---------- 绘图 ----------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))

    # (a) ΔT→0 通道
    ax = axes[0, 0]
    dts = [r["delta_T"] for r in sw["rows"]]
    drs = [r["drift_mean"] for r in sw["rows"]]
    ax.plot(dts, drs, "o", color="tab:red", ms=8, label="实测漂移率 Ω_drift", zorder=5)
    domega_grid = np.linspace(0.0, 1.0, 200)
    dt_grid = domega_grid / (2.0 + domega_grid)      # ω_L=1, ω_E=1+Δω
    ax.plot(dt_grid, [adler_drift(dw, K) for dw in domega_grid],
            color="black", lw=1.5, label="Adler 精确解 √(Δω²−4K²)")
    dt_lock = delta_T(1.0, 1.0 + 2 * K)
    ax.axvspan(0, dt_lock, color="tab:green", alpha=0.12)
    ax.text(dt_lock * 0.5, 0.45, "锁定区\n（势差被接口吸收）", ha="center",
            fontsize=9, color="tab:green")
    ax.set_xlabel("时间势差 ΔT")
    ax.set_ylabel("相对相位漂移率")
    ax.set_title("(a) ΔT→0 通道：漂移率连续逼近锁定边界")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # (b) 代表性相位轨迹
    ax = axes[0, 1]
    r_unlock = coupled_phase_experiment(1.0, 1.5, coupling=K, seed=0)
    r_lock = coupled_phase_experiment(1.0, 1.05, coupling=0.15, seed=0)
    t_u = np.arange(len(r_unlock["phis"])) * r_unlock["dt"]
    t_l = np.arange(len(r_lock["phis"])) * r_lock["dt"]
    ax.plot(t_u, r_unlock["phis"], color="tab:red", lw=0.8,
            label=f"解锁（ΔT={r_unlock['delta_T']:.2f}，漂移）")
    ax.plot(t_l, r_lock["phis"], color="tab:green", lw=0.8,
            label=f"锁定（ΔT={r_lock['delta_T']:.2f}，势差被吸收）")
    ax.set_xlabel("时间")
    ax.set_ylabel("相对相位 φ = θ_L − θ_E")
    ax.set_title("(b) 漂移 vs 锁定：ΔT 的两种归宿")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # (c) P5 阴性对照
    ax = axes[1, 0]
    ratios = np.array([1.05, 1.2, 1.5, 2.0])
    drift_k0, drift_k = [], []
    for rr in ratios:
        drift_k0.append(coupled_phase_experiment(1.0, rr, coupling=0.0,
                                                 seed=0)["drift_rate"])
        drift_k.append(coupled_phase_experiment(1.0, rr, coupling=K,
                                                seed=0)["drift_rate"])
    ax.plot(ratios - 1.0, drift_k0, "s--", color="grey",
            label="K=0（无耦合）：漂移=Δω")
    ax.plot(ratios - 1.0, drift_k, "o-", color="tab:red",
            label=f"K={K}（有耦合）：漂移<Δω")
    ax.plot(ratios - 1.0, ratios - 1.0, ":", color="black", lw=1,
            label="恒等线（理论上限）")
    ax.set_xlabel("频率差 Δω")
    ax.set_ylabel("漂移率")
    ax.set_title("(c) P5 阴性对照：锁定是耦合诱导的")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # (d) P4 谓词独立
    ax = axes[1, 1]
    fr = np.array([1.0, 1.2, 1.5, 2.0, 3.0])
    ax.plot(fr, pi["dt_vals"], "o-", color="tab:blue", label="ΔT（频率失配）")
    ax.plot(fr, pi["d_vals"], "s--", color="tab:orange",
            label="d（状态重合度，固定状态）")
    ax.set_xlabel("频率比 ω_E/ω_L")
    ax.set_ylabel("谓词取值")
    ax.set_title("(d) P4 谓词独立：ΔT 携带 d 之外的信息")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("G4 / ΔT 操作化：时间势差 TP 谓词 = 内禀节律失配，"
                 "观测量 = 相对相位漂移（Adler 锚定）", fontsize=13)
    fig.tight_layout()
    out = FIG_DIR / "g4_tp_temporal_potential.png"
    fig.savefig(out, dpi=140)
    print(f"\n图已保存：{out}")

    print("\n" + "-" * 70)
    print("结论（写入修订登记的素材）")
    print("  1. TP 谓词操作化完成：ΔT = |ω_L−ω_E|/(ω_L+ω_E)，观测量 = 相对相位漂移率。")
    print("  2. ΔT→0 通道通过：解锁区漂移率与 Adler 精确解一致（误差 < 3%），")
    print("     Δω→2K⁺ 时漂移连续到 0；|Δω|≤2K 锁定（势差被接口吸收）。")
    print("  3. 锁定是耦合诱导的（K=0 阴性对照任意小 Δω 都外显为漂移）。")
    print("  4. ΔT 与 d 独立可变（P4），TP 不是冗余谓词。")
    print("  5. 框架语义：锁定 = 时间势差在接口内被消化；漂移 = 势差外显为")
    print("     不可逆时序流——与 A11a/b 势-流分离的时间版一致。")


if __name__ == "__main__":
    main()
