"""
X-S-MERA V3.4 · Relay-chain experiment（命题 2′ 中继部分）
=========================================================

三个问题：
1. F1 吸收态能否经中继链恢复？（直接恢复 vs n=1 不足 vs n=3 成功）
2. 串联中继的复合律：乘积律 ∏dᵢ 还是调和/短板律 1/Σ(1/dᵢ)？
3. F2 段在链中的角色：阻断还是污染？

输出：figs/relay_recovery_and_law.png

运行：python experiments/sim_relay_recovery.py
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

from src.relay_chain import (  # noqa: E402
    min_relay_count, relay_bridge_demo, simulate_transport)

FIG_DIR = os.path.join(_REPO_ROOT, "figs")

Z_L = np.array([1.30, 0.20, 0.00])
Z_E = np.array([0.21, 0.25, 0.00])


def run() -> dict:
    out = {}

    # ---- (1) 恢复对比 ----
    out["recovery"] = {
        "no relay": relay_bridge_demo(use_relay=False),
        "n=1 (insufficient)": relay_bridge_demo(use_relay=True, n_relay=1),
        "n=3 (sufficient)": relay_bridge_demo(use_relay=True, n_relay=3),
    }
    out["n_min"] = min_relay_count(Z_L, Z_E)

    # ---- (2) 复合律 ----
    d_leg = 0.9
    ns = [1, 2, 3, 4]
    out["law"] = {"ns": ns,
                  "measured": [simulate_transport([d_leg] * n)["T_measured"] for n in ns],
                  "harmonic": [d_leg / n for n in ns],
                  "product": [d_leg ** n for n in ns]}
    out["unequal"] = simulate_transport([0.95, 0.7])

    # ---- (3) 失效模式 ----
    out["failure"] = {
        "all compatible": simulate_transport([0.9, 0.9]),
        "one F2 leg": simulate_transport([0.9, 0.9], leg_c2=[True, False]),
        "one F1 leg": simulate_transport([0.9, 0.1]),
    }

    # 中继数下限随端距/弱端尺度的变化
    out["n_min_scan"] = {
        "collapsed pair (gap 1.09)": min_relay_count(Z_L, Z_E),
        "closer pair (gap 0.6)": min_relay_count(np.array([0.80, 0.22, 0.0]),
                                                 np.array([0.21, 0.25, 0.0])),
    }
    return out


def plot(out: dict) -> str:
    os.makedirs(FIG_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
    fig.suptitle("X-S-MERA · Relay chain: recovery, composition law, failure modes",
                 fontsize=14, y=0.99)

    # ---- (1) 恢复对比 ----
    ax = axes[0]
    colors = {"no relay": "#c0392b", "n=1 (insufficient)": "#ca6f1e",
              "n=3 (sufficient)": "#1f7a45"}
    for label, res in out["recovery"].items():
        ax.plot(res["t"], res["gap"], lw=1.8, label=f"{label}  (gap→{res['gap'][-1]:.2f})",
                color=colors[label])
    ax.axhline(1.091, color="#888", lw=0.8, ls=":", label="initial gap = 1.09")
    ax.set_xlabel("t")
    ax.set_ylabel("‖L − E‖  (end-to-end distance)")
    ax.set_title(f"(1) F1 recovery — min relays needed: n = {out['n_min']}",
                 fontsize=11, loc="left")
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.3)

    # ---- (2) 复合律 ----
    ax = axes[1]
    law = out["law"]
    ax.plot(law["ns"], law["product"], "o--", color="#c0392b", lw=1.4, ms=5,
            label=r"product law  $\prod d_i$")
    ax.plot(law["ns"], law["harmonic"], "s-", color="#1f5f99", lw=1.8, ms=5,
            label=r"harmonic law  $1/\Sigma(1/d_i)$")
    ax.plot(law["ns"], law["measured"], "D", color="#111", ms=7, zorder=5,
            label="measured T")
    ue = out["unequal"]
    ax.scatter([2], [ue["T_measured"]], marker="*", s=130, color="#6c3483", zorder=6,
               label=f"unequal [0.95, 0.7]: T={ue['T_measured']:.3f} "
                     f"(harm {ue['T_harmonic']:.3f} / prod {ue['T_product']:.3f})")
    ax.set_yscale("log")
    ax.set_xticks(law["ns"])
    ax.set_xlabel("number of legs n (d = 0.9 each)")
    ax.set_ylabel("end-to-end transmission T")
    ax.set_title("(2) Composition law — harmonic (weak-leg) law,\nNOT the product law",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.3, which="both")

    # ---- (3) 失效模式 ----
    ax = axes[2]
    labels = list(out["failure"].keys())
    T = [out["failure"][k]["T_measured"] for k in labels]
    C = [out["failure"][k]["contamination"] for k in labels]
    x = np.arange(len(labels))
    bars = ax.bar(x, T, width=0.55, color=["#1f7a45", "#d68910", "#c0392b"], alpha=0.85)
    for xi, t_val in zip(x, T):
        ax.text(xi, t_val + 0.015, f"T={t_val:.2f}", ha="center", fontsize=9)
    ax2 = ax.twinx()
    ax2.plot(x, C, "D", color="#111", ms=7, label="contamination (secondary-mode offset)")
    for xi, c_val in zip(x, C):
        ax2.text(xi + 0.12, c_val + 0.004, f"{c_val:.3f}", fontsize=8.5)
    ax2.set_ylabel("ghost contamination", color="#444")
    ax2.tick_params(axis="y", colors="#444")
    ax2.set_ylim(0, max(C) * 1.6 + 0.02)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("transmission T")
    ax.set_ylim(0, 0.6)
    ax.set_title("(3) Failure modes — F1 blocks; F2 transmits\nbut contaminates",
                 fontsize=10.5, loc="left")
    ax.grid(alpha=0.3, axis="y")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(FIG_DIR, "relay_recovery_and_law.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def summarize(out: dict) -> None:
    print("\n==== Relay chain summary ====")
    print(f"min relays to bridge collapsed pair (gap 1.09): n = {out['n_min']}")
    for k, v in out["n_min_scan"].items():
        print(f"  n_min[{k}] = {v}")
    print("\ncomposition law (d = 0.9 per leg):")
    print(f"  n      : {out['law']['ns']}")
    print(f"  meas   : {[round(v, 4) for v in out['law']['measured']]}")
    print(f"  harmonic: {[round(v, 4) for v in out['law']['harmonic']]}")
    print(f"  product : {[round(v, 4) for v in out['law']['product']]}")
    print("\nfailure modes (n = 2, d = 0.9):")
    for k, v in out["failure"].items():
        print(f"  {k:16s} T={v['T_measured']:.4f}  contamination={v['contamination']:.4f}")
    print("\n-- 判读 --")
    print("  1) 串联中继服从调和/短板律 T=1/Σ(1/dᵢ)，乘积律 ∏dᵢ 被否定（命题 2′ 需修正：")
    print("     ∏ 只适用于各段独立必要的并联/事件语义；串联传输由最弱段支配）")
    print("  2) C2 不是通行必要条件：F2 段主模态照常传输，但注入 ghost 污染")
    print("     → C2 是'可信度'条件而非'连通性'条件（F2=伪收敛的又一证据）")
    print("  3) 中继数下限由弱端尺度决定（n_min 随端距/尺度差增大）")
    print("     → I1 的可构造性有定量代价，可当作跨尺度桥接的预测")


if __name__ == "__main__":
    results = run()
    path = plot(results)
    print(f"[ok] figure saved: {path}")
    summarize(results)
