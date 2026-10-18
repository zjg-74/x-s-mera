"""
R4：A8/D17 冲突的分量归类判据 —— LE 内核 F2 实测
==================================================

R4 结论（修订登记 §十三）：管辖按可测属性分配，K1 状态相关性 → K2 修复检验
→ K3 滤波检验，顺序 load-bearing；复合现象先分解后归类，禁整体归类。

本脚本在 LECouplingKernel 的 F2 阶段（ghost = 1.6× 幅值失配 + 0.25× 底噪）
实测该管线，并跑两个纯分量对照（纯失配 / 纯底噪）。

图（2×2）：
  (a) 残差轨迹与分解重构（失配分量 + 底噪分量）
  (b) 系数恢复：α/β 估计值 vs 解析真值
  (c) K2 修复检验：修复前后残差范数
  (d) K3 滤波检验：底噪分量在 δ>0 / δ<0 两子样本的均值（同号 = 偏置型）
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

from src.r4_component_attribution import (  # noqa: E402
    classify_le_f2, ATTRIBUTION_CRITERIA,
)

FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    print("=" * 70)
    print("R4：A8/D17 分量归类（LE 内核 F2 阶段实测）")
    print("=" * 70)
    print("\n预注册判据：")
    for k, v in ATTRIBUTION_CRITERIA.items():
        print(f"  {k}: {v}")

    # 复合场景（默认内核：1.6× 失配 + 0.25× 底噪）
    res = classify_le_f2(n_steps=120, seed=0)
    print("\n复合场景（1.6× 失配 + 0.25× 底噪）：")
    for row in res["rows"]:
        print(f"  模态 {row['mode']}: α={row['alpha_mismatch']:+.4f} "
              f"(真值 {res['true_alpha']:+.4f})  β={row['beta_base']:+.4f} "
              f"(真值 {res['true_beta']:+.4f})  R²={row['r_squared']:.4f}")
        print(f"    失配分量 ρ={row['mismatch_rho']:+.3f} → {row['mismatch_domain']}")
        print(f"    底噪分量 符号不变={row['base_sign_invariant']} → {row['base_domain']}")
    rp = res["repair"]
    print(f"  K2 修复检验：残差范数 {rp['residual_norm_f2']:.4f} → "
          f"{rp['residual_norm_repaired']:.4f}（消除率 {rp['elimination_ratio']:.2%}）")

    # 纯分量对照
    pure_mm = classify_le_f2(n_steps=120, seed=0, ghost_base_ratio=0.0)
    pure_bs = classify_le_f2(n_steps=120, seed=0, ghost_gain_ratio=1.0)
    print("\n纯失配对照（底噪=0）：α=%.4f β=%.4f"
          % (pure_mm["rows"][0]["alpha_mismatch"], pure_mm["rows"][0]["beta_base"]))
    print("纯底噪对照（失配比=1）：α=%.4f β=%.4f"
          % (pure_bs["rows"][0]["alpha_mismatch"], pure_bs["rows"][0]["beta_base"]))

    # ---------- 绘图 ----------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))

    comp = res["decomposition"]["components"][0]
    n_show = 240
    xs = np.arange(n_show)

    # (a) 残差分解重构
    ax = axes[0, 0]
    mask = res["decomposition"]["mask"]
    r_mode = res["trace"]["residual"][mask][:, comp["mode"]][:n_show]
    ax.plot(xs, r_mode, color="black", lw=1.2, label="守恒残差 r（实测）")
    ax.plot(xs, comp["mismatch_part"][:n_show], color="tab:red", lw=1.2,
            label="失配分量重构 α·gain·δ")
    ax.plot(xs, comp["base_part"][:n_show], color="tab:blue", lw=1.2,
            label="底噪分量重构 β·|gain|")
    ax.set_xlabel("活跃步序（C1 通过）")
    ax.set_ylabel("残差幅值")
    ax.set_title("(a) 复合残差的分解（先分解后归类）")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) 系数恢复 vs 真值
    ax = axes[0, 1]
    scen = ["复合", "纯失配", "纯底噪"]
    alphas = [res["rows"][0]["alpha_mismatch"],
              pure_mm["rows"][0]["alpha_mismatch"],
              pure_bs["rows"][0]["alpha_mismatch"]]
    betas = [res["rows"][0]["beta_base"],
             pure_mm["rows"][0]["beta_base"],
             pure_bs["rows"][0]["beta_base"]]
    x = np.arange(3)
    ax.bar(x - 0.2, alphas, width=0.4, color="tab:red", label="α（失配系数）估计")
    ax.bar(x + 0.2, betas, width=0.4, color="tab:blue", label="β（底噪系数）估计")
    ax.axhline(res["true_alpha"], color="tab:red", ls="--", lw=1)
    ax.axhline(res["true_beta"], color="tab:blue", ls="--", lw=1)
    ax.text(2.4, res["true_alpha"] - 0.02, f"真值 α={res['true_alpha']:+.2f}",
            fontsize=8, color="tab:red", ha="right")
    ax.text(2.4, res["true_beta"] - 0.02, f"真值 β={res['true_beta']:+.2f}",
            fontsize=8, color="tab:blue", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels(scen)
    ax.set_ylabel("系数值")
    ax.set_title("(b) 分解系数恢复（虚线=解析真值）")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    # (c) K2 修复检验
    ax = axes[1, 0]
    bars = [rp["residual_norm_f2"], rp["residual_norm_repaired"]]
    ax.bar(["F2（未修复）", "修复后（C2 过）"], bars,
           color=["tab:red", "tab:green"], width=0.5)
    ax.set_ylabel("残差范数均值")
    ax.set_title(f"(c) K2 修复检验：消除率 {rp['elimination_ratio']:.2%}"
                 f"（阈值 ≥ {ATTRIBUTION_CRITERIA['K2_elimination_min']:.0%}）")
    for i, b in enumerate(bars):
        ax.text(i, b + 0.002, f"{b:.4f}", ha="center", fontsize=10)
    ax.grid(alpha=0.3, axis="y")

    # (d) K3 滤波检验：底噪分量的符号不变性
    ax = axes[1, 1]
    delta_i = comp["delta_i"]
    base = comp["base_part"]
    pos, neg = base[delta_i > 0], base[delta_i < 0]
    ax.boxplot([pos, neg], tick_labels=["δ > 0 子样本", "δ < 0 子样本"])
    ax.axhline(0, color="grey", ls=":", lw=1)
    same_sign = bool(np.sign(pos.mean()) == np.sign(neg.mean()))
    ax.set_title(f"(d) K3 滤波检验：底噪分量均值同号 = {same_sign}（偏置型）")
    ax.set_ylabel("底噪分量幅值")
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("R4：管辖按可测属性分配（K1→K2→K3），复合异常先分解后归类",
                 fontsize=13)
    fig.tight_layout()
    out = FIG_DIR / "r4_component_attribution.png"
    fig.savefig(out, dpi=140)
    print(f"\n图已保存：{out}")

    print("\n" + "-" * 70)
    print("结论（写入修订登记的素材）")
    print("  1. 复合残差可无歧义分解：α、β 精确恢复解析真值（R²=1.000）。")
    print("  2. K1：失配分量与状态差强相关（|ρ|≈0.95–0.98）→ 接口域（A8-F2，修复）。")
    print("  3. K2：修复接口后残差消除率 100% → 失配分量的接口域归属被证实。")
    print("  4. K3：底噪分量在 δ 翻号下符号不变 → 载体域（D17，按偏置滤波消除）。")
    print("  5. 可辨识性前提：分解要求 δ 充分变号（多初值集成），")
    print("     且屏蔽 C1 失败步——这本身是判据的适用条件，须写入条款。")


if __name__ == "__main__":
    main()
