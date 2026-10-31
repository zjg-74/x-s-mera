"""
B 线 P4 · S0 真人全样本 Bootstrap 稳健性（纯内部方法验证）

复用 r23_robustness 的方法：对现有全部 S0 真人轨迹的每步权重做 Bootstrap 重采样
扰动（固定种子，B=200），重跑判别，统计：
1. 每任务 200 副本下的判定分布（酝酿/反刍/空转/无震荡）
2. 全局"多数方向 = ρ_res>ρ_exp"在扰动下保持的比例

对应 README"R23 判别器稳健性"先例（当时仅对 S1 18 任务）。本脚本对 S0 真人全样本跑。
数据来自 analyze_s0_all 的 METADATA（同一批 S0 轨迹）。
"""
import json, random
from pathlib import Path
from collections import Counter
import analyze_trajectories as at
from analyze_s0_all import METADATA, load_task

random.seed(0)
rng = random.Random(0)
BASE = Path(__file__).parent
OUT = BASE / "s0_robustness.json"
B = 200


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
    tasks = []
    for fname, tid, display, solved, typ, batch in METADATA:
        tasks.append((display, load_task(fname, tid)))
    print(f"S0 tasks: {len(tasks)}, bootstrap B={B}")

    # 基准（无扰动）
    base_votes = [at.classify(at.analyze_task(t)["osc_net_change"]) for _, t in tasks]
    base_cnt = Counter(v for v in base_votes if v)
    base_majority = base_cnt.most_common(1)[0] if base_cnt else None
    print("base verdict counts:", dict(base_cnt), "majority:", base_majority)

    # 逐任务：B 副本判定分布
    task_dist, task_stable = {}, {}
    for disp, t in tasks:
        votes = []
        for _ in range(B):
            t_pert = {"steps": [
                {"k": s["k"], "weights": bootstrap_weights(s["weights"], rng)}
                for s in t["steps"]
            ]}
            votes.append(at.classify(at.analyze_task(t_pert)["osc_net_change"]))
        c = Counter(v for v in votes if v)
        task_dist[disp] = dict(c)
        if c:
            m = c.most_common(1)[0][0]
            task_stable[disp] = {"majority": m, "share": round(c[m] / B, 3), "n_votes": sum(c.values())}
        print(f"[{disp}] dist={dict(c)}")

    # 全局：每副本多数方向
    global_majority = Counter()
    for _ in range(B):
        votes = []
        for _, t in tasks:
            t_pert = {"steps": [
                {"k": s["k"], "weights": bootstrap_weights(s["weights"], rng)}
                for s in t["steps"]
            ]}
            v = at.classify(at.analyze_task(t_pert)["osc_net_change"])
            if v:
                votes.append(v)
        if votes:
            global_majority[Counter(votes).most_common(1)[0][0]] += 1
    global_share = {k: round(v / B, 3) for k, v in global_majority.items()}
    print("global majority share:", global_share)

    out = {
        "method": "S0 真人全样本 Bootstrap 稳健性(seeded B=200) + 重跑判别",
        "tasks": len(tasks), "B": B, "seed": 0,
        "base_verdict_counts": dict(base_cnt),
        "base_majority": base_majority,
        "task_verdict_distribution": task_dist,
        "per_task_majority_stability": task_stable,
        "global_majority_share": global_share,
    }
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("saved:", OUT)


if __name__ == "__main__":
    main()
