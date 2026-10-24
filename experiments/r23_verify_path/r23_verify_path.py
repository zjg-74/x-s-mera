"""
R23 验证路径 · ρ_res 与 ρ_exp 大小关系的可执行判定（GAPS X2 步骤 1）

来源：
  - 补篇五 §1.6 跨层引用规则 / R23：仅在已验证 ρ_res > ρ_exp 前提下许可两层互引。
  - 推论 5.3（震荡的生产性）：震荡期（ψ̄ = 0）U 走向 ⟺ ρ_res 与 ρ_exp 大小关系。
  - 补篇四 §3.2 程序 P̂：真实系统上 U 需从解释分布估计（本脚本用模型真值 U）。

R23 操作程序（本脚本落实为可复现实现）：
  ① 令 A_D 系统进入震荡：控制 r = θ_E（⟹ sign(ψ̄) = sign(r−θ_E) = 0，Ψ ≈ 0）。
  ② 用程序 P̂ 估计震荡期 U 的走向。
  ③ 判据：U 降 ⟹ ρ_res > ρ_exp；U 升 ⟹ ρ_res < ρ_exp；U 不变 ⟹ ρ_res = ρ_exp。

分层说明（诚实边界）：
  - 本脚本验证"判据本身在 A_D 模型内正确"（模型内自洽，可复现）。
  - 真实系统上 U 不可直接观测，须先经补篇四程序 P̂ 从解释分布估计（本脚本不覆盖该真实测量步骤）。
  - A_D 是建模假设（B-B22），本脚本不使 A_D 为真，只在其下复现推论 5.3 的可达模式。

用法：python -u r23_verify_path.py
输出：experiments/r23_verify_path/results.json
"""
import json
import random
from pathlib import Path

# A_D 双速率系统：每步 = 展开（投入 a，dE+=a·e₊, dU+=a·u₊）+ 消解（投入 b，dE−=b·e₋, dU−=b·u₋）
# 记 r = a/b；θ_E = e₋/e₊（Ψ=0 临界比），θ_U = u₋/u₊（dU=0 临界比）
# ρ_exp = u₊/e₊（展开效率比），ρ_res = u₋/e₋（消解效率比）
# sign(ψ̄) = sign(r−θ_E)，sign(dU/dk) = sign(r−θ_U)；θ_E<θ_U ⟺ ρ_res>ρ_exp


def simulate_oscillation(rho_exp, rho_res, U0, steps, noise_std, rng):
    """A_D 系统在震荡（r=θ_E）下的 U 轨迹。

    归一化：e_plus = 1.0，e_minus = 1.0 ⟹ θ_E = e_minus/e_plus = 1.0。
    u_plus = ρ_exp·e_plus；u_minus = ρ_res·e_minus。
    震荡：r = a/b = θ_E ⟹ 取 a = θ_E, b = 1。
    每步 dU = a·u_plus − b·u_minus + ξ（ξ ~ N(0, noise_std)）。
    推论 5.3：震荡期 dU/dk = θ_E·u₊ − u₋ = u₊(θ_E − θ_U)，符号由 ρ_res − ρ_exp 决定。
    """
    e_plus = 1.0
    e_minus = 1.0
    u_plus = rho_exp * e_plus
    u_minus = rho_res * e_minus
    theta_E = e_minus / e_plus
    a, b = theta_E, 1.0          # r = a/b = θ_E ⟹ 震荡（Ψ ≈ 0）
    U = U0
    trace = [U]
    for _ in range(steps):
        dU = a * u_plus - b * u_minus
        if noise_std > 0:
            dU += rng.gauss(0.0, noise_std)
        U += dU
        trace.append(U)
    return trace


def classify(trace, bottom=0.05, top_mult=1.5):
    """震荡期 U 走向判别：
      触底 = min ≤ bottom（≈0，收敛完成）；
      触顶 = max ≥ U0·top_mult（持续上升）；
      净变化 = 末 − 初（含噪声下的稳健方向）。
    """
    U0 = trace[0]
    mn, mx = min(trace), max(trace)
    net = trace[-1] - U0
    direction = "down" if net < -1e-9 else ("up" if net > 1e-9 else "flat")
    return {"net": net, "min": mn, "max": mx,
            "touched_bottom": mn <= bottom, "touched_top": mx >= U0 * top_mult,
            "direction": direction}


def main():
    out_dir = Path(__file__).resolve().parent
    configs = {
        "I_rho_res_gt_rho_exp":  {"rho_exp": 0.30, "rho_res": 0.90, "truth": "rho_res > rho_exp"},
        "II_rho_res_lt_rho_exp": {"rho_exp": 0.90, "rho_res": 0.30, "truth": "rho_res < rho_exp"},
        "III_rho_res_eq_rho_exp":{"rho_exp": 0.50, "rho_res": 0.50, "truth": "rho_res = rho_exp"},
    }
    N_REP = 400
    STEPS = 60
    NOISE = 0.05
    U0 = 4.0
    seed = 20261024

    summary = {}
    for cid, cfg in configs.items():
        rng = random.Random(seed)
        outs = [classify(simulate_oscillation(cfg["rho_exp"], cfg["rho_res"],
                                               U0, STEPS, NOISE, rng))
                for _ in range(N_REP)]
        avg_net = sum(o["net"] for o in outs) / len(outs)
        bottom_rate = sum(o["touched_bottom"] for o in outs) / len(outs)
        top_rate = sum(o["touched_top"] for o in outs) / len(outs)
        down_frac = sum(1 for o in outs if o["direction"] == "down") / len(outs)
        up_frac = sum(1 for o in outs if o["direction"] == "up") / len(outs)
        if down_frac > 0.9:
            verdict = "rho_res > rho_exp"
        elif up_frac > 0.9:
            verdict = "rho_res < rho_exp"
        else:
            verdict = "rho_res = rho_exp"
        match = (verdict == cfg["truth"])
        summary[cid] = {
            "config": cfg["truth"],
            "avg_U_net_change": round(avg_net, 3),
            "U_bottom_rate": round(bottom_rate, 3),
            "U_top_rate": round(top_rate, 3),
            "down_frac": round(down_frac, 3),
            "up_frac": round(up_frac, 3),
            "verdict": verdict,
            "truth": cfg["truth"],
            "match": match,
        }
        print(f"[{cid}] avg_net={avg_net:+.3f} bottom={bottom_rate:.3f} "
              f"top={top_rate:.3f} down={down_frac:.3f} up={up_frac:.3f} "
              f"→ verdict={verdict} match={match}")

    ok = all(s["match"] for s in summary.values())
    out = {
        "model": "A_D double-rate (expansion/resolution), oscillation r=theta_E",
        "steps": STEPS, "reps": N_REP, "noise_std": NOISE, "U0": U0, "seed": seed,
        "method": ("R23: oscillate (control Psi≈0 by r=theta_E), read U trend "
                   "(program P-hat); U down => rho_res>rho_exp, up => <, flat => ="),
        "results": summary,
        "all_match": ok,
    }
    with open(out_dir / "results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("all_match:", ok)
    print("saved:", out_dir / "results.json")


if __name__ == "__main__":
    main()
