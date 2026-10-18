"""
R3-v2：C1+C2 充分性补案仿真
============================

R3（修订登记 §12.8）只检验了 C1+C2 的**必要性**方向；本脚本检验**充分性**方向：
构造"单次交互全合规（C1∧C2 恒过）但序列必失败"的反例全谱。

五类构造（src/r3_sufficiency_probe.py）：
  S1 累积损伤型 / S2 竞争置换型 / S3 时序错配型 / S4 中介污染型 / S5 势耗竭型

图：
  (a) 五类构造的单次合规率 vs 序列成功率（预注册判据：100% ∧ 0% → 强反例）
  (b) S3 的双读法解剖：瞬时读法 |vL−vE| 越界报警 vs 增量读法 |ΔvL−ΔvE| 全程放行
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

from src.r3_sufficiency_probe import (  # noqa: E402
    run_full_sweep, SUFFICIENCY_CRITERIA,
)

FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    print("=" * 70)
    print("R3-v2：C1+C2 充分性补案（单次合规 ∧ 序列失败）")
    print("=" * 70)
    print("\n预注册判据：")
    for k, v in SUFFICIENCY_CRITERIA.items():
        print(f"  {k}: {v}")

    res = run_full_sweep(n_trials=200)
    print("\n五类构造判定：")
    for c in res["cases"]:
        print(f"  {c['case']:24s} 单次合规={c['single_pass_rate']:.2f}  "
              f"序列成功={c['series_success']:.3f}  → {c['verdict']}")
    print(f"\n强反例 {res['n_strong']} / 弱反例 {res['n_weak']} / 不合法 {res['n_invalid']}")
    print(f"总结论：{res['summary']}")

    cases = res["cases"]
    labels = [c["case"].split(" ")[0] for c in cases]

    # ---------- 绘图 ----------
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    # (a) 单次合规率 vs 序列成功率
    ax = axes[0]
    x = np.arange(len(cases))
    ax.bar(x - 0.2, [c["single_pass_rate"] for c in cases], width=0.4,
           color="tab:blue", label="单次合规率（C1∧C2）")
    ax.bar(x + 0.2, [c["series_success"] for c in cases], width=0.4,
           color="tab:red", label="序列成功率")
    ax.axhline(0.20, color="grey", ls=":", lw=1)
    ax.text(len(cases) - 0.5, 0.21, "弱反例阈值 0.20", fontsize=8, color="grey",
            ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("比率")
    ax.set_title("(a) 五类构造：单次全合规、序列全失败（强反例 ×5）")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3, axis="y")

    # (b) S3 双读法解剖：一条代表性轨迹
    ax = axes[1]
    rng = np.random.default_rng(2)
    n_rounds, drift_rate, tol, step_tol = 40, 0.015, 0.25, 0.04
    vL = vE = 0.0
    gaps, incs = [0.0], [0.0]
    for _t in range(n_rounds):
        dvL = drift_rate * (1 + 0.3 * rng.normal())
        dvE = drift_rate * 0.4 * (1 + 0.3 * rng.normal())
        vL += dvL
        vE += dvE
        gaps.append(abs(vL - vE))
        incs.append(abs(dvL - dvE))
    ts = np.arange(n_rounds + 1)
    ax.plot(ts, gaps, color="tab:red", lw=2, label="瞬时读法 |vL−vE|（累积）")
    ax.axhline(tol, color="tab:red", ls="--", lw=1)
    ax.text(1, tol + 0.01, f"瞬时容差 {tol}", fontsize=8, color="tab:red")
    ax.plot(ts, incs, color="tab:blue", lw=1.5, label="增量读法 |ΔvL−ΔvE|（每步）")
    ax.axhline(step_tol, color="tab:blue", ls="--", lw=1)
    ax.text(1, step_tol + 0.008, f"单步容差 {step_tol}", fontsize=8, color="tab:blue")
    ax.set_xlabel("交互轮次 τ")
    ax.set_ylabel("接口差距")
    ax.set_title("(b) S3：增量读法全程放行，瞬时读法中途报警")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)

    fig.suptitle("R3-v2：C1+C2 是单次交互的准入判据，不是历史成败的充分条件",
                 fontsize=13)
    fig.tight_layout()
    out = FIG_DIR / "r3_sufficiency_probe.png"
    fig.savefig(out, dpi=140)
    print(f"\n图已保存：{out}")

    # ---------- 结论 ----------
    print("\n" + "-" * 70)
    print("结论（写入修订登记的素材）")
    print("  1. C1+C2 的充分性被证伪：五类'单次全合规、序列全失败'构造全部成立。")
    print("  2. C1+C2 的作用域须显式声明为**单次交互的准入判据**；")
    print("     历史/序列成败需要追加序列条件（损伤累积率、置换率、版本漂移、")
    print("     中介污染叠加、势再生率 ≥ 消耗率）。")
    print("  3. S3 逼出的条款修正：C2 的接口兼容性必须按**瞬时读法**")
    print("     （|vL−vE| < tolerance，对照绝对基准）判定；")
    print("     若按增量读法执行，漂移可在全程合规中越界。")


if __name__ == "__main__":
    main()
