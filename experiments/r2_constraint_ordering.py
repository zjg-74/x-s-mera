"""R2: 约束强度排序 vs D 形态排序 —— 物理 / 生物 / 人文三域对照。

判据（先立判据，避免事后裁量）
-----------------------------
1. tau 必须**外生**：取投影自身申报的相对精度 1/sigma_quoted。
   申报精度不是从 D 反推的，而是投影自己写进论文的承诺值。
2. D 的形态统计量必须在**申报 σ 标度**（z = r/sigma_quoted）下计算，
   因为文献里所有"4.2σ 反常""6σ 预测误差"都用这个口径。
3. 判据要能区分三种可能：
   (a) 排序一致          -> A4 的原表述可用；
   (b) 排序不一致        -> A4 的"正态性"版本被否；
   (c) 统计量不可辨识    -> A4 的"唯一决定"根本无法检验（本轮的实际结论）。

三域实例（文献锚定，参数为标定近似，见 docs/experiment_manual.md §7）
  物理   QED/精密测量投影  申报精度 5.4e-7（μ子 g-2 2025 理论白皮书 540 ppb）
  生物   数量遗传/分子演化  申报精度 ~1e-1（GWAS 效应量、分子钟节点）
  人文   宏观预测投影      申报精度 ~8e-1（GDP/通胀预测误差）

输出：表 1 排序对照；表 2 峰度构造成分塌缩；图 r2_d_form_ordering.png
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.d_form_audit import synthesize, core_width, tail_mass, spearman

FIG_DIR = "figs"
os.makedirs(FIG_DIR, exist_ok=True)

DOMAINS = {
    "physics": dict(tau=1.0 / 5.4e-7, coverage=0.98, n_paths=2,
                    a_unmodeled=2.0e-6, a_route=2.0e-6, color="#1f4e79",
                    note="mu g-2: exp 127 ppb / theory 540 ppb"),
    "biology": dict(tau=1.0 / 0.10, coverage=0.85, n_paths=2,
                    a_unmodeled=0.50, a_route=0.30, color="#2e7d32",
                    note="GWAS effect sizes: tail exponent 1<alpha<2"),
    "humanities": dict(tau=1.0 / 0.80, coverage=0.55, n_paths=1,
                       a_unmodeled=1.00, a_route=0.30, color="#b45309",
                       note="2008 forecast error ~6 sigma under Gaussian"),
}


def measure(dom: dict, n: int = 40000, seed: int = 7) -> dict:
    """按申报 σ 标度计算全部形态统计量。"""
    out = synthesize(tau=dom["tau"], coverage=dom["coverage"],
                     n_paths=dom["n_paths"], a_unmodeled=dom["a_unmodeled"],
                     a_route=dom["a_route"], n=n, seed=seed)
    r = out["r"]
    sigma_q = 1.0 / dom["tau"]          # 申报精度（外生给定）
    z = r / sigma_q
    v = z.var()
    exk = float(np.mean((z - z.mean()) ** 4) / v ** 2 - 3.0)
    return dict(
        tau=dom["tau"], sigma_q=sigma_q, coverage=dom["coverage"],
        n_paths=dom["n_paths"],
        rho=dom["a_unmodeled"] / sigma_q,          # 承诺-实际失配比
        core=core_width(r),
        overclaim=core_width(r) / sigma_q,         # 实测核心 / 申报精度
        exceed3=float(np.mean(np.abs(z) > 3.0)),
        kurt=exk,
        tail_abs=tail_mass(r, tol=0.1),
        r=r, color=dom["color"], note=dom["note"],
    )


def main() -> None:
    rows = {k: measure(v, seed=11 + i) for i, (k, v) in enumerate(DOMAINS.items())}

    print("=" * 108)
    print("表 1  tau 排序 vs D 形态排序（n=40000 / 域；z 以各域申报 sigma 为单位）")
    print("=" * 108)
    print(f"{'domain':12s}{'tau':>10s}{'rho':>7s}{'core/claim':>11s}"
          f"{'exceed3':>9s}{'kurt':>10s}{'tail_abs':>10s}")
    print("-" * 108)
    for k, v in rows.items():
        print(f"{k:12s}{v['tau']:10.3g}{v['rho']:7.2f}{v['overclaim']:11.3f}"
              f"{v['exceed3']:9.4f}{v['kurt']:10.1f}{v['tail_abs']:10.4f}")

    def order(key):
        return [k for k, _ in sorted(rows.items(), key=lambda kv: kv[1][key])]

    tau_ord = order("tau")[::-1]           # 强 -> 弱
    core_ord = order("core")               # 窄 -> 宽
    kurt_ord = order("kurt")               # 低 -> 高
    ex3_ord = order("exceed3")

    print("-" * 108)
    print("tau 排序（强→弱）      :", " > ".join(tau_ord))
    print("核心宽度排序（窄→宽）  :", " > ".join(core_ord), "| 与 tau 一致:", core_ord == tau_ord)
    print("峰度排序（低→高）      :", " > ".join(kurt_ord), "| 与 tau 一致:", kurt_ord == tau_ord)
    print("超标频率排序（低→高）  :", " > ".join(ex3_ord), "| 与 tau 一致:", ex3_ord == tau_ord)
    print()
    print("说明：核心宽度与 tau 的一致是**构造性的**（核心宽度由申报精度生成），")
    print("      不构成对 A4 的独立支持；只有峰度/超标频率携带额外信息。")

    # --------------------------------------- 表 2：等峰度线的多对一（不可辨识）
    rho_grid = 10.0 ** np.linspace(-1.0, 3.5, 26)
    kap_grid = np.linspace(0.50, 0.99, 15)
    K = np.zeros((len(kap_grid), len(rho_grid)))
    for i, kap in enumerate(kap_grid):
        for j, rho in enumerate(rho_grid):
            r = synthesize(tau=1.0, coverage=kap, n_paths=2,
                           a_unmodeled=rho, a_route=0.4 * rho,
                           n=12000, seed=5)["r"]
            v = r.var()
            K[i, j] = float(np.mean((r - r.mean()) ** 4) / v ** 2 - 3.0)

    targets = [1.0, 3.0, 5.0, 10.0]
    print()
    print("=" * 108)
    print("表 2  等峰度线：达到同一峰度所需的 rho（行=峰度目标，列=覆盖度 kappa）")
    print("=" * 108)
    print(f"{'kurt':>6s}" + "".join(f"{'k=%.2f' % k:>10s}" for k in (0.50, 0.60, 0.70, 0.80, 0.90, 0.99)))
    print("-" * 108)
    for tgt in targets:
        cells = []
        for kq in (0.50, 0.60, 0.70, 0.80, 0.90, 0.99):
            i = int(np.abs(kap_grid - kq).argmin())
            col = K[i]
            if col.max() >= tgt >= col.min():
                cells.append(f"{np.interp(tgt, col, rho_grid):10.2g}")
            else:
                cells.append(f"{'—':>10s}")
        print(f"{tgt:6.0f}" + "".join(cells))
    print("-" * 108)
    print("判读：每一行都能在 **多个 kappa 列** 上取到（rho 相应变化）→ 映射 (rho,kappa)->kurt")
    print("      是多对一的：读到峰度无法反推 kappa，也无法反推 rho。")
    print("      例：kurt≈3 既可由 (rho≈2.5, kappa=0.90) 也可由 (rho≈5.3, kappa=0.60) 给出。")

    # ------------------------------------------------------- 表 3 可辨识性
    recs = []
    for t in 10.0 ** np.linspace(0, 6, 9):
        for c in np.linspace(0.55, 0.99, 8):
            for p in (1, 2):
                m = synthesize(tau=t, coverage=c, n_paths=p,
                               a_unmodeled=1.0, a_route=0.35, n=8000, seed=0)
                r = m["r"]
                z = r * t
                v = z.var()
                recs.append(dict(log_tau=np.log10(t), coverage=c, n_paths=p,
                                 kurt=float(np.mean((z - z.mean()) ** 4) / v ** 2 - 3.0),
                                 core=core_width(r)))

    print()
    print("=" * 108)
    print("表 3  指标可辨识性：Spearman 秩相关（网格 9x8x2）")
    print("=" * 108)
    print(f"{'statistic':12s}{'log_tau':>11s}{'coverage':>11s}{'n_paths':>10s}")
    print("-" * 108)
    for st in ("core", "kurt"):
        vals = [spearman([d[st] for d in recs], [d[c] for d in recs])
                for c in ("log_tau", "coverage", "n_paths")]
        print(f"{st:12s}" + "".join(f"{v:11.3f}" for v in vals))
    print("-" * 108)
    print("判读：kurt 对 tau 与 coverage 同时敏感 -> 两参数不可分离；")
    print("      core 只随 tau 变（对 coverage/paths 近 0）-> 可读 tau，但构造上近恒真。")

    # ------------------------------------------------------------------ 图
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.7))
    fig.suptitle("R2  |  Does D-form order with constraint strength across domains?",
                 fontsize=13, x=0.008, ha="left")

    ax = axes[0]
    for k, v in rows.items():
        e = np.log10(np.abs(v["r"]) + 1e-14)
        ax.hist(e, bins=160, density=True, histtype="step", lw=1.7,
                color=v["color"], label=f"{k} (tau={v['tau']:.2g})")
    ax.set_title("(1) Residual distribution, log10 |r|\n(absolute dimensionless scale)",
                 fontsize=10.5, loc="left")
    ax.set_xlabel("log10 | residual |"); ax.set_ylabel("density")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25)

    ax = axes[1]
    for k, v in rows.items():
        ax.scatter(v["tau"], v["core"], s=120, color=v["color"], zorder=3,
                   edgecolor="white", linewidth=1.2)
        ax.annotate(k, (v["tau"], v["core"]), textcoords="offset points",
                    xytext=(9, 6), fontsize=8.5, color=v["color"])
    tt = np.array([0.4, 4e6])
    ax.plot(tt, 1.35 / tt, "--", color="#999999", lw=1.2,
            label=r"slope $-1$: core $\propto 1/\tau$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"constraint strength $\tau=1/\sigma_{quoted}$")
    ax.set_ylabel("core width (IQR)")
    ax.set_title("(2) Core width orders with tau — but\nby construction (nearly tautological)",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25, which="both")

    ax = axes[2]
    cf = ax.contourf(np.log10(rho_grid), kap_grid, K,
                     levels=np.linspace(-0.5, 14, 22), cmap="YlOrBr")
    cs = ax.contour(np.log10(rho_grid), kap_grid, K,
                    levels=[1, 3, 5, 10], colors="#333333", linewidths=0.9)
    ax.clabel(cs, fmt="%g", fontsize=7.5)
    ax.set_xlabel(r"$\log_{10}\rho$  with $\rho=A_{unmodeled}/\sigma_{quoted}$")
    ax.set_ylabel(r"coverage $\kappa$")
    ax.set_title("(3) Iso-kurtosis lines are diagonal:\n($\\rho$,$\\kappa$) $\\to$ kurt is many-to-one",
                 fontsize=10.5, loc="left")
    fig.colorbar(cf, ax=ax, pad=0.02, label="excess kurtosis")

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(FIG_DIR, "r2_d_form_ordering.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print()
    print(f"[ok] figure saved: {out}")


if __name__ == "__main__":
    main()
