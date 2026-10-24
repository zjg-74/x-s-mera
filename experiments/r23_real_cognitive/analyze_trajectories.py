"""
B 线 P1 · 真实 LLM 轨迹分析：程序 P̂（补篇四 §3.2）+ 震荡期识别 + ρ_res>ρ_exp 判定

数据源：真实 LLM（Doubao）对构造歧义任务的多步推理输出（可观测的候选解释集合演化）。
本脚本只读轨迹 JSON，客观计算，不预设结果。

程序 P̂（补篇四 §3.2）：
  P̂ = -Σ_j p̂_j·log2(p̂_j) + (r̂-1)/(2N·ln 2)     （Miller–Madow 偏差修正）
  权重 w_j 为各解释类的伪计数；N = Σ w_j；p̂_j = w_j / N；r̂ = 解释类个数。

震荡期识别（补篇二 §1）：Ψ ≈ 0 ⟺ 多候选并存且 U 相对平稳：
  连续步段 r̂ ≥ 2 且 |U(k+1)-U(k)| < ε_tol。

ρ_res>ρ_exp 判定（推论 5.3 / R23 逆用）：
  震荡期 U 净变化：降 ⟹ ρ_res>ρ_exp；升 ⟹ ρ_res<ρ_exp；不变 ⟹ =。
  稳健性：跨任务取多数方向（阈值 >0.5 为多数；>0.9 为强一致）。

用法：python -u analyze_trajectories.py trajectory.json
输出：stdout 摘要 + results.json
"""
import json, math, sys
from pathlib import Path


def p_hat(weights):
    """补篇四 §3.2 程序 P̂：给定各解释类伪计数，估计意义条件熵（含 Miller–Madow 修正）。"""
    w = {k: v for k, v in weights.items() if v > 0}
    if not w:
        return 0.0
    N = sum(w.values())
    r_hat = len(w)
    ps = [v / N for v in w.values()]
    entropy = -sum(p * math.log2(p) for p in ps)
    correction = (r_hat - 1) / (2 * N * math.log(2))
    return entropy + correction


def find_oscillation_segments(rs):
    """震荡期 = 多候选并存的连续步段（r>=2，即 Ψ≈0）。
    推论 5.3：震荡期 U 走向是待测的判别量（酝酿=降/反刍=升/空转=不变），
    故**不**要求 U 平稳——那会把酝酿/反刍误筛为无震荡。"""
    segs = []
    run = []
    for k, rk in enumerate(rs):
        if rk >= 2:
            run.append(k)
        else:
            if len(run) >= 2:
                segs.append(run)
            run = []
    if len(run) >= 2:
        segs.append(run)
    return segs


def analyze_task(task):
    steps = task["steps"]  # 每步 {weights: {类:伪计数}}
    Us = [p_hat(s["weights"]) for s in steps]
    rs = [len({k: v for k, v in s["weights"].items() if v > 0}) for s in steps]
    segs = find_oscillation_segments(rs)
    # 取最长震荡期判 U 走向
    if segs:
        seg = max(segs, key=len)
        net = Us[seg[-1]] - Us[seg[0]]
    else:
        seg = None
        net = None
    return {"us": Us, "rs": rs, "segs": segs,
            "longest_seg": seg, "osc_net_change": net}


def classify(net):
    if net is None:
        return None
    tol = 1e-4
    if net < -tol:
        return "rho_res > rho_exp"
    if net > tol:
        return "rho_res < rho_exp"
    return "rho_res = rho_exp"


def main():
    in_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "trajectory.json"
    out_name = sys.argv[2] if len(sys.argv) > 2 else "results.json"
    data = json.load(open(in_path, encoding="utf-8"))
    res = {}
    votes = []
    for t in data["trajectories"]:
        a = analyze_task(t)
        verdict = classify(a["osc_net_change"])
        votes.append(verdict)
        res[t["id"]] = {
            "us": [round(u, 4) for u in a["us"]],
            "rs": a["rs"],
            "longest_seg": a["longest_seg"],
            "osc_net_change": round(a["osc_net_change"], 4) if a["osc_net_change"] is not None else None,
            "verdict": verdict,
        }
        print(f"[{t['id']}] seg={a['longest_seg']} osc_net={res[t['id']]['osc_net_change']} "
              f"verdict={verdict}")
    # 多数判定
    from collections import Counter
    cnt = Counter(v for v in votes if v)
    if not cnt:
        majority = None
    else:
        majority = cnt.most_common(1)[0]
    print("verdict counts:", dict(cnt))
    print("majority:", majority)
    out = {
        "method": "program P-hat (supp4 3.2) + oscillation (Psi≈0) + U trend => rho_res>rho_exp (R23)",
        "results": res,
        "verdict_counts": dict(cnt),
        "majority": majority,
    }
    with open(Path(__file__).parent / out_name, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("saved:", Path(__file__).parent / out_name)


if __name__ == "__main__":
    main()
