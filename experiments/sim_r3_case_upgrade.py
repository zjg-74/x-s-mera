"""
R3-v2 案例库定量升级：S3（流感）与 S5（鳕鱼）
==============================================

数据源与判据见 src/r3_case_upgrade.py 模块文档串（全部公开来源 + 预注册判据）。

图（1×2）：
  (a) CDC 逐季疫苗有效性（2004-05→2025-26），标注 WHO 换株季与 2014-15 漂移错配季；
      附 2014-15 窗口时间线（决策时点 >90% 匹配 → 季中 <50% 抑制 → VE 19%）
  (b) 北方鳕 SSB 锚点序列（1960→2024），标注禁捕令、危急区出口、崩溃/恢复时标
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

from src.r3_case_upgrade import (  # noqa: E402
    FLU_VE, FLU_H3N2, DRIFT_SEASON, FLU_2014_TIMELINE,
    flu_verdict, cod_verdict, COD_SSB_ANCHORS, COD_SSB_RECENT, COD_FACTS,
    CASE_CRITERIA,
)

FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    print("=" * 70)
    print("R3-v2 案例定量升级：S3 流感抗原漂移 / S5 鳕鱼崩溃")
    print("=" * 70)
    print("\n预注册判据：")
    for k, v in CASE_CRITERIA.items():
        print(f"  {k}: {v}")

    fv, cv = flu_verdict(), cod_verdict()
    r1, r2, r3 = fv["C-S3-1"], fv["C-S3-2"], fv["C-S3-3"]
    print("\n[S3 流感]")
    print(f"  C-S3-1 漂移季 VE={r1['drift_ve']:.0f}% vs 其余紧 CI 季均值 "
          f"{r1['others_mean']:.1f}%（z={r1['z_score']:.2f}，比值 {r1['ratio']:.2f}，"
          f"系列最低={r1['is_series_min_tight']}）→ {'过' if r1['pass'] else '否'}")
    print(f"  C-S3-2 WHO 换株率 {r2['n_updates']}/{r2['n_transitions']} = "
          f"{r2['update_rate']:.0%} → {'过' if r2['pass'] else '否'}")
    print(f"  C-S3-3 窗口：决策时点匹配 {r3['match_at_decision']:.0%} → 季中抑制 "
          f"{r3['inhibition_mid_season']:.0%} → VE {r3['ve']:.0%} → "
          f"{'过' if r3['pass'] else '否'}")

    print("\n[S5 鳕鱼]")
    print(f"  C-S5-1 F>M（1950s-80s）+ 1990-95 SSB 净降 "
          f"{cv['C-S5-1']['ssb_drop_1990_95']:.0%} → {'过' if cv['C-S5-1']['pass'] else '否'}")
    print(f"  C-S5-2 禁捕后低于 LRP {cv['C-S5-2']['years_below_lrp_after_moratorium']} 年 "
          f"→ {'过' if cv['C-S5-2']['pass'] else '否'}")
    print(f"  C-S5-3 恢复/崩溃时标 {cv['C-S5-3']['t_recover_yr']}/"
          f"{cv['C-S5-3']['t_collapse_yr']} = {cv['C-S5-3']['ratio']:.1f} "
          f"（2024 年 SSB 仍仅 1990 年的 {cv['C-S5-3']['ssb_2024_vs_1990']:.0%}）"
          f"→ {'过' if cv['C-S5-3']['pass'] else '否'}")
    print(f"  C-S5-4 崩溃后 M 均值 {COD_FACTS['m_post_mean']} = 崩溃前假设的 "
          f"{cv['C-S5-4']['m_ratio_post_pre']:.1f}×；Capelin 仅剩 "
          f"{cv['C-S5-4']['capelin_pct']:.0f}% → {'过' if cv['C-S5-4']['pass'] else '否'}")
    print(f"\nS3 全过 = {fv['all_pass']}；S5 全过 = {cv['all_pass']}")

    # ---------- 绘图 ----------
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))

    # (a) 流感 VE 序列 + 换株标记
    ax = axes[0]
    seasons = [s for s in FLU_VE if FLU_VE[s] is not None]
    ves = [FLU_VE[s] for s in seasons]
    xs = np.arange(len(seasons))
    colors = ["tab:red" if s == DRIFT_SEASON else "tab:blue" for s in seasons]
    ax.bar(xs, ves, color=colors, alpha=0.85)
    # 换株季标记（瞬时读法报警）
    prev = None
    for i, s in enumerate(seasons):
        if prev is not None and FLU_H3N2.get(s) != FLU_H3N2.get(prev):
            ax.axvline(i - 0.5, color="grey", ls=":", lw=0.9)
        prev = s
    ax.axhline(r1["others_mean"], color="tab:green", ls="--", lw=1.2)
    ax.text(len(xs) - 1, r1["others_mean"] + 1.5,
            f"其余紧 CI 季均值 {r1['others_mean']:.1f}%", ha="right",
            fontsize=8, color="tab:green")
    i_drift = seasons.index(DRIFT_SEASON)
    ax.annotate("2014-15 漂移错配季\nVE=19%（z=−2.76）\n窗口：>90%→<50% 抑制",
                (i_drift, 19), textcoords="offset points", xytext=(-10, 45),
                fontsize=8, color="tab:red",
                arrowprops=dict(arrowstyle="->", color="tab:red"))
    ax.set_xticks(xs[::2])
    ax.set_xticklabels([seasons[i] for i in range(0, len(xs), 2)],
                       rotation=45, ha="right", fontsize=7.5)
    ax.set_ylabel("调整后总 VE（%）")
    ax.set_title("(a) S3：CDC 逐季疫苗有效性；灰虚线 = WHO 换株（瞬时读法报警，"
                 "换株率 71%）")
    ax.grid(alpha=0.3, axis="y")

    # (b) 鳕鱼 SSB 锚点 + 事件
    ax = axes[1]
    yrs = sorted(COD_SSB_ANCHORS) + sorted(COD_SSB_RECENT)
    ssb = [COD_SSB_ANCHORS[y] for y in sorted(COD_SSB_ANCHORS)] + \
          [COD_SSB_RECENT[y] for y in sorted(COD_SSB_RECENT)]
    ax.plot(yrs, ssb, "o-", color="tab:blue", lw=1.5, ms=5,
            label="SSB 锚点（DFO-FW 图 4 粗读数 + SAR 2024 表 1）")
    ax.axvline(1992, color="tab:red", ls="--", lw=1.2)
    ax.text(1992.5, 1750, "1992 禁捕令\n（交互停止）", fontsize=8, color="tab:red")
    ax.axvline(2016, color="tab:green", ls="--", lw=1.2)
    ax.text(2016.5, 1500, "2016 出危急区\n（低于 LRP 24 年）", fontsize=8,
            color="tab:green")
    ax.annotate("", xy=(1995, 700), xytext=(1990, 1000),
                arrowprops=dict(arrowstyle="->", color="black"))
    ax.text(1987.5, 720, "崩溃 5 年 −95%", fontsize=8)
    ax.annotate("", xy=(2024, 430), xytext=(1995, 60),
                arrowprops=dict(arrowstyle="->", color="tab:orange"))
    ax.text(2002, 250, "恢复 29 年仅至 ~1.2×LRP\n（2024 = 1990 年的 34%）",
            fontsize=8, color="tab:orange")
    ax.set_xlabel("年份")
    ax.set_ylabel("SSB（kt）")
    ax.set_title("(b) S5：北方鳕 SSB——耗竭期 F>M、吸收态 24 年、迟滞比 5.8、"
                 "M 均值升至 2.35×")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    fig.suptitle("R3-v2 案例库定量升级：S3（流感）与 S5（鳕鱼）预注册判据全过",
                 fontsize=13)
    fig.tight_layout()
    out = FIG_DIR / "r3_case_upgrade.png"
    fig.savefig(out, dpi=140)
    print(f"\n图已保存：{out}")

    print("\n" + "-" * 70)
    print("结论（回写修订登记 §17.3/§17.5）")
    print("  S3：增量合规 ∧ 绝对越界在真实序列上成立——2014-02 决策时点 >90% 匹配，")
    print("      同一生产窗口内漂移株把抑制率压到 <50%，当季 VE 19%（z=−2.76）；")
    print("      WHO 换株率 71% = 瞬时读法（绝对基准）是现役协议，且 2014-09 改荐")
    print("      Switzerland 株正是报警触发的即时纠正。")
    print("  S5：耗竭期 F>M（1950s-80s），崩溃 5 年 −95%；禁捕后低于 LRP 24 年；")
    print("      恢复/崩溃时标比 5.8；崩溃后 M 均值 0.47=2.35× 崩溃前假设、")
    print("      Capelin 仅 10%——再生机制本身受损（S5 的 phi_min 通道在数据中找到对应）。")


if __name__ == "__main__":
    main()
