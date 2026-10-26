"""
R23 判别器稳健性：Bootstrap 扰动下多数方向是否稳定（纯内部方法验证）。

问题：R23 判别（震荡期 U 走向 ⟹ ρ_res>ρ_exp）读的是 LLM 自报告伪计数权重。
这些权重带报告噪声；稳健性验证 = 在权重上做 Bootstrap 重采样扰动，重跑判别，
统计"多数方向 = ρ_res>ρ_exp"在扰动下保持的比例。

方法（机械可复核，固定种子）：
1. 对每个任务每步权重 {w_j}，N=Σw_j，按 p̂_j=w_j/N 做 Multinomial(N) 重采样 → 扰动权重 {w'_j}。
2. 重跑 p_hat + 震荡识别 + U 走向判定（复用 analyze_trajectories）。
3. B=200 次扰动，统计：
   - 每任务 200 副本下判定分布（酝酿/反刍/无震荡）
   - 全局"多数方向 = ρ_res>ρ_exp"的比例
4. 稳健性判据：多数方向保持比例（参考 R23 单次随机 53.5%/46.5% 的边缘性）。

用法：python -u r23_robustness.py
"""
import json, random, sys
from pathlib import Path
from collections import Counter

random.seed(0)  # 机械可复核

# 复用判别逻辑（同 analyze_trajectories.py，避免重复实现）
import analyze_trajectories as at

BASE = Path(__file__).parent
IN = BASE / "trajectory_all.json"
OUT = BASE / "r23_robustness.json"
B = 200  # 扰动副本数


def bootstrap_weights(weights, rng):
    """按 p̂ 做 N 次独立采样（标准库，Multinomial 等价），返回扰动权重 dict。"""
    items = [(k, v) for k, v in weights.items() if v > 0]
    if not items:
        return weights
    N = sum(v for _, v in items)
    ps = [v / N for _, v in items]
    keys = [k for k, _ in items]
    counts = {k: 0 for k in keys}
    for _ in range(N):
        r = rng.random()
        acc = 0.0
        for k, p in zip(keys, ps):
            acc += p
            if r < acc:
                counts[k] += 1
                break
    return counts


def main():
    rng = random.Random(0)
    data = json.load(open(IN, encoding="utf-8"))
    tasks = data["trajectories"]
    print(f"tasks: {len(tasks)}, bootstrap B={B}")

    # 基准（无扰动）
    base_votes = []
    for t in tasks:
        base_votes.append(at.classify(at.analyze_task(t)["osc_net_change"]))
    base_cnt = Counter(v for v in base_votes if v)
    base_majority = base_cnt.most_common(1)[0] if base_cnt else None
    print("base verdict counts:", dict(base_cnt), "majority:", base_majority)

    # 逐任务：200 副本的判定分布
    task_verdict_dist = {}
    per_task_majority_stable = {}
    for t in tasks:
        votes = []
        for _ in range(B):
            # 扰动该任务全部步骤权重
            t_pert = {"steps": [
                {"k": s["k"], "weights": bootstrap_weights(s["weights"], rng)}
                for s in t["steps"]
            ]}
            a = at.analyze_task(t_pert)
            votes.append(at.classify(a["osc_net_change"]))
        c = Counter(v for v in votes if v)
        task_verdict_dist[t["id"]] = dict(c)
        # 该任务多数方向（若有）
        if c:
            m = c.most_common(1)[0][0]
            per_task_majority_stable[t["id"]] = {
                "majority": m,
                "share": c[m] / B,
                "n_votes": sum(c.values()),
            }
        print(f"[{t['id']}] dist={dict(c)}")

    # 全局：每副本的多数方向分布
    global_majority_dist = Counter()
    for _ in range(B):
        votes = []
        for t in tasks:
            t_pert = {"steps": [
                {"k": s["k"], "weights": bootstrap_weights(s["weights"], rng)}
                for s in t["steps"]
            ]}
            a = at.analyze_task(t_pert)
            v = at.classify(a["osc_net_change"])
            if v:
                votes.append(v)
        if votes:
            m = Counter(votes).most_common(1)[0][0]
            global_majority_dist[m] += 1
    global_majority_share = {
        k: v / B for k, v in global_majority_dist.items()
    }
    print("global majority share over B:", global_majority_share)

    out = {
        "method": "R23 判别器稳健性：Bootstrap 权重扰动(seeded) + 重跑判别",
        "tasks": len(tasks),
        "B": B,
        "seed": 0,
        "base_verdict_counts": dict(base_cnt),
        "base_majority": base_majority,
        "task_verdict_distribution": task_verdict_dist,
        "per_task_majority_stability": per_task_majority_stable,
        "global_majority_share": global_majority_share,
    }
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("saved:", OUT)


if __name__ == "__main__":
    main()
