"""R5: 推论 3（演化无向择优）的 A9 资格审计.

核心对照：判据②的**确定性读法** vs **统计版读法**，在 Ns 梯度上各考一次。
原降级理由用的是确定性读法；A9 修订（判据②统计版 + 判据③ d 重述）之后
读法变了，降级是否还成立，是本实验要裁决的唯一问题。

场景： Wright–Fisher 选择-漂移，N ∈ {50, 200, 2000, 20000}（Ns = 1, 4, 40, 400），
      s = 0.02，干预幅度 Δs = 0.02，T = 200 代，p0 = 0.10。
"""

from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.a9_qualification import (  # noqa: E402
    DS, N_BINS, P0, REPS, T_GENS, TH_D, TH_EFF,
    closure_leak, coarse_predict, d_overlap, frozen_transfer,
    intervention_response, p_fix, qualify, wright_fisher,
)

FIG_DIR = "figs"
os.makedirs(FIG_DIR, exist_ok=True)

NS_LIST = [50, 200, 2000, 20000]
S = 0.02
DET_BAND = 0.05      # 确定性读法的容差带（预注册：|p_T − p_det| < 0.05 占比 ≥ 0.8）
DET_TH = 0.80


def deterministic_path(s: float = S, T: int = T_GENS, p0: float = P0) -> float:
    p = p0
    for _ in range(T):
        p = p * (1.0 + s) / (1.0 + s * p)
    return p


def main() -> dict:
    print("=" * 104)
    print("R5 · 推论 3（演化无向择优）A9 资格审计 —— 判据②统计版 / 判据③ d 重述版")
    print("=" * 104)
    print(f"  预注册阈值：θ_d = {TH_D}（预测/迁移重合度），θ_eff = {TH_EFF}（干预响应），"
          f"Δs = {DS}，T = {T_GENS} 代，p0 = {P0}，reps = {REPS}，bins = {N_BINS}")
    print(f"  确定性读法（原降级理由的读法）：|p_T − p_det| < {DET_BAND} 的实现占比 ≥ {DET_TH}")
    print()

    p_det = deterministic_path()

    # ---------------- 逐语境审计
    rows = []
    for N in NS_LIST:
        rng = np.random.default_rng(N)
        real = wright_fisher(P0, S, N, T_GENS, REPS, rng)
        pred = coarse_predict(P0, S, N, T_GENS, REPS, np.random.default_rng(9000 + N))
        d_stat = d_overlap(pred, real)
        det_frac = float(np.mean(np.abs(real - p_det) < DET_BAND))
        iv = intervention_response(N, seed=N)
        rows.append({"N": N, "Ns": N * S, "d_stat": d_stat, "det_frac": det_frac,
                     "iv": iv, "pfix": p_fix(real), "passes": iv["passes"]})

    print("表 1  判据②：确定性读法 vs 统计版读法（同一组现实实现）")
    print("=" * 104)
    print(f"{'N':>8s}{'Ns':>7s}{'确定性读法通过率':>16s}{'统计版 d(预测,现实)':>20s}"
          f"{'干预响应 eff':>13s}{'②判决':>10s}")
    print("-" * 104)
    for r in rows:
        det_ok = r["det_frac"] >= DET_TH
        stat_ok = r["d_stat"] >= TH_D
        print(f"{r['N']:8d}{r['Ns']:7.0f}{r['det_frac']:16.3f}{r['d_stat']:20.3f}"
              f"{r['iv']['effect']:13.3f}"
              f"{'确定✗统计✓' if (not det_ok and stat_ok) else ('双双通过' if det_ok else '双双失败'):>10s}")
    print("-" * 104)
    print("  判读：确定性读法（原降级理由）在 Ns 小的语境失败；统计版读法下预测重合度全境通过。")
    print("        原降级理由不成立——它淘汰的不是'不可预测'，而是'点预测不准'，后者不是判据②。")
    print()

    # ---------------- 判据③ 迁移
    t3 = frozen_transfer(seed=3)
    print("表 2  判据③（d 重述版）：参照语境（N=20000）标定后冻结参数逐语境迁移")
    print("=" * 104)
    print(f"{'N':>8s}{'Ns':>7s}{'d(冻结)':>10s}{'d(重标定)':>11s}{'③判决':>8s}   P(优等位固定)")
    print("-" * 104)
    for r in t3:
        print(f"{r['N']:8d}{r['Ns']:7.0f}{r['d_frozen']:10.3f}{r['d_refit']:11.3f}"
              f"{'✓' if r['passes'] else '✗':>8s}   {r['pfix']:.3f}")
    print("-" * 104)
    print("  判读：冻结参数跨语境迁移在 N 与参照差一个数量级以上时崩塌；")
    print("        按语境重标定后 d 回升 → 失败原因是'参数没跟着语境走'，不是机制本身错了。")
    print()

    # ---------------- 判据① 闭合
    cl = closure_leak(N=2000, seed=5)
    print("表 3  判据①（闭合版，N=2000）：钳制外部池频率，下游不得改变")
    print("=" * 104)
    print(f"{'m':>7s}{'d(低池,高池)':>13s}{'闭合？':>7s}{'d(粗模型,现实)':>15s}")
    print("-" * 104)
    for c in cl:
        print(f"{c['m']:7.2f}{c['d_ext']:13.3f}{'✓' if c['closed'] else '✗':>7s}"
              f"{c['d_model']:15.3f}")
    print("-" * 104)
    print("  判读：m=0 时外部状态进不来（闭合）；m 增大后 p_ext 直接写进下游（边界泄漏），")
    print("        且忽略边界条件的粗模型预测同步劣化——①失败不是形式问题，是可测的。")
    print()

    # ---------------- 序层方向
    print("表 4  推论 3 的\"序层无向\"：D_seq = P(优等位固定) − p0")
    print("=" * 104)
    print(f"{'N':>8s}{'Ns':>7s}{'P_fix':>9s}{'D_seq':>9s}{'无向？':>8s}")
    print("-" * 104)
    dir_rows = list(rows)
    # 补一个 Ns→0 极限点（s=0.002, N=50, Ns=0.1）：检验"无向"是不是极限性质
    r_extra = wright_fisher(P0, 0.002, 50, T_GENS, REPS, np.random.default_rng(77))
    dir_rows.append({"N": 50, "Ns": 0.1, "pfix": p_fix(r_extra)})
    for r in dir_rows:
        dseq = r["pfix"] - P0
        print(f"{r['N']:8d}{r['Ns']:7.1f}{r['pfix']:9.3f}{dseq:9.3f}"
              f"{'≈无向' if dseq < 0.05 else '有向':>8s}")
    print("-" * 104)
    print("  判读：s=0.02 时即使 Ns=1 方向指数也有 0.106——\"无向\"在 s>0 时**不成立**；")
    print("        只有 Ns→0（Ns=0.1 实测 0.009）才落入无向带。无向性是 Ns→0 的极限性质，")
    print("        不是\"漂移语境\"的性质。推论 3 的定性内容随 Ns 连续变化——")
    print("        这正是它不能作为无条件宏 M 的内因：它自己就带着资格条件 Ns→0。")
    print()

    # ---------------- 资格裁决
    cls = [{"m": c["m"], "closed": c["closed"], "d_ext": c["d_ext"]} for c in cl]
    verdicts = qualify(rows, t3, cls)
    print("表 5  三判据汇总与 M 资格裁决")
    print("=" * 104)
    print(f"{'N':>8s}{'Ns':>7s}{'①':>5s}{'②':>5s}{'③':>5s}   裁决与资格条件")
    print("-" * 104)
    for v in verdicts:
        n = f"{v['N']:8d}" if v["N"] else f"{'—':>8s}"
        ns = f"{v['Ns']:7.0f}" if v["Ns"] is not None else f"{'—':>7s}"
        c2 = "✓" if v["C2"] else ("✗" if v["C2"] is False else "—")
        c3 = "✓" if v["C3"] else ("✗" if v["C3"] is False else "—")
        print(f"{n}{ns}{v['C1']:>5s}{c2:>5s}{c3:>5s}   {v['verdict']}  {v['qualifier']}")
    print("-" * 104)
    print("  总判决：")
    print("  (a) 原降级【理由】不成立——它依赖判据②的确定性读法，而统计版读法下")
    print("      预测重合度全境 ≥0.93；被淘汰的是\"点预测\"，不是\"可预测\"。")
    print("  (b) 原【结论】部分保留、须改写：该机制在每个合格语境**内部**是宏 M")
    print("      「宏 M | Ns ≳ 5 ∧ m = 0 ∧ (s, N) 逐语境显式申报」，")
    print("      但 ③（d 重述版）表明它**不是跨数量级的单一宏**——冻结参数迁移即崩塌。")
    print("      即：不是\"推论 3 降级为机制\"，而是\"宏 M 资格是逐语境的\"。")
    print("  (c) A9 修订落点：M 资格 = 三判据 ∧ 显式资格条件 C_q（Ns 区间 / 边界闭合 / 参数申报）；")
    print("      无 C_q 字段的 A9 会让同一机制在不同论文里得到相反的资格判定。")
    print()

    # ---------------- 图
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10))
    fig.suptitle("R5 · A9 qualification audit of Inference-3 (undirected variation, "
                 "selection filters): deterministic vs statistical reading", fontsize=12.5)

    # (1) 判据②两种读法
    ax = axes[0, 0]
    nss = [r["Ns"] for r in rows]
    ax.plot(nss, [r["det_frac"] for r in rows], "o-", color="#c0392b", lw=1.6,
            label=f"deterministic reading (band ±{DET_BAND})")
    ax.plot(nss, [r["d_stat"] for r in rows], "s-", color="#27ae60", lw=1.6,
            label="statistical reading d(pred, real)")
    ax.axhline(TH_D, color="#888888", ls="--", lw=0.9, label=f"θ = {TH_D}")
    ax.set_xscale("log"); ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel("Ns   (selection vs drift)")
    ax.set_ylabel("passing fraction / overlap d")
    ax.set_title("(1) Criterion-2: the original downgrade hangs on the DETERMINISTIC reading\n"
                 "statistical reading passes wherever selection is not swamped",
                 fontsize=10, loc="left")
    ax.legend(fontsize=7.6, frameon=False, loc="center right")
    ax.grid(alpha=0.25, which="both")

    # (2) 判据③迁移
    ax = axes[0, 1]
    ns3 = [r["N"] for r in t3]
    ax.plot(ns3, [r["d_frozen"] for r in t3], "o-", color="#8e44ad", lw=1.6,
            label="frozen params (calibrated at N=20000)")
    ax.plot(ns3, [r["d_refit"] for r in t3], "s-", color="#2980b9", lw=1.4, ls="--",
            label="re-instantiated per context")
    ax.axhline(TH_D, color="#888888", ls=":", lw=0.9)
    ax.set_xscale("log"); ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel("target context N")
    ax.set_ylabel("d(frozen pred, realization)")
    ax.set_title("(2) Criterion-3 (d-restated): transfer fails across ≥1 decade of N;\n"
                 "re-instantiation restores d → failure is parameter transfer, not mechanism",
                 fontsize=10, loc="left")
    ax.legend(fontsize=7.6, frameon=False, loc="center left")
    ax.grid(alpha=0.25, which="both")

    # (3) 判据①闭合
    ax = axes[1, 0]
    ms = [c["m"] for c in cl]
    ax.plot(ms, [c["d_ext"] for c in cl], "o-", color="#e67e22", lw=1.6,
            label="closure: d(downstream | p_ext low, p_ext high)")
    ax.plot(ms, [c["d_model"] for c in cl], "s--", color="#16a085", lw=1.3,
            label="coarse model that ignores m: d(pred, real)")
    ax.axhline(TH_D, color="#888888", ls=":", lw=0.9)
    ax.set_xscale("symlog", linthresh=0.005)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel("migration rate m  (boundary leak)")
    ax.set_ylabel("overlap d")
    ax.set_title("(3) Criterion-1 (closure): external state enters downstream once m > 0",
                 fontsize=10, loc="left")
    ax.legend(fontsize=7.6, frameon=False, loc="lower left")
    ax.grid(alpha=0.25, which="both")

    # (4) 序层方向 + 裁决
    ax = axes[1, 1]
    dseq = [r["pfix"] - P0 for r in rows]
    ax.plot(nss, dseq, "D-", color="#2c3e50", lw=1.6,
            label=r"$D_{seq}=P(\mathrm{fit\ allele\ fixes})-p_0$")
    ax.axhline(0.0, color="#888888", lw=0.8, ls=":")
    ax.axhline(0.05, color="#bbbbbb", lw=0.8, ls="-.",
               label="≈undirected band (D<0.05)")
    ax.set_xscale("log")
    ax.set_xlabel("Ns"); ax.set_ylabel("direction index")
    ax.set_title("(4) Inference-3's own \"undirected\" claim is regime-dependent:\n"
                 "true only where drift swamps selection",
                 fontsize=10, loc="left")
    ax.legend(fontsize=7.6, frameon=False, loc="upper left")
    ax.grid(alpha=0.25, which="both")
    ax.text(0.98, 0.04,
            "verdict: natural selection = macro M with qualifier\n"
            "          [macro M | Ns >= 5 AND m = 0 AND (s, N) declared per context]\n"
            "          original reason (deterministic unpredictability) rejected;\n"
            "          A9 needs an explicit qualifier field C_q",
            transform=ax.transAxes, fontsize=7.6, ha="right", va="bottom",
            color="#7b241c",
            bbox=dict(facecolor="#fdf2e9", edgecolor="#e59866", boxstyle="round,pad=0.35"))

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "r5_a9_qualification.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[ok] figure saved: {out}")
    return {"rows": rows, "t3": t3, "cl": cl, "verdicts": verdicts}


if __name__ == "__main__":
    main()
