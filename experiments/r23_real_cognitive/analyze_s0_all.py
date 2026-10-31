"""
B 线 P4 · S0 真人全样本统一判别 + 解出/未解出混淆检验（纯内部，不新增采集）

把现有全部 S0 真人 think-aloud 轨迹（10-26 6q / 10-27 单题 / M 系列 / 10-31 2/2(1) /
10-31 更早批次）用同一判别器（analyze_trajectories）统一跑一遍，做两件事：

1. **全样本汇总**：全部可判任务的酝酿/反刍/空转分布与多数方向。
2. **解出/未解出混淆检验**（M 系列已验证的"解题成败"主混淆变量）：
   在数学题内按"是否解出"分组，对照两组的酝酿/反刍方向——
   若方向几乎由"是否解出"决定，则现有样本不能作为干净的酝酿/反刍认知信号。

候选类=解题路径（沿用各批次既有归约）。solved 标注取自各批次逐字稿/self_trace_note
（解出=得到明确结论；False=卡住/放弃/讲错未收敛；None=非数学无"解出"语义）。
"""
import json
from pathlib import Path
from collections import Counter
import analyze_trajectories as at

BASE = Path(__file__).parent
OUT = BASE / "results_s0_all.json"

# (文件, 文件内tid, 展示名(唯一), solved, type, batch)
# solved: True=数学题解出, False=数学题未解出, None=非数学
# 10-31a(2/2(1)) 与 10-31b(更早批次) 有同名 id(T-CUBE*/CORNER/SHAPES)，用展示名区分，保留全量 29 任务
METADATA = [
    # 10-26 六题（思辨，非数学）
    ("trajectory_human_6q.json", "H6q-1-zexiao", "H6q-1", None, "thinking", "10-26"),
    ("trajectory_human_6q.json", "H6q-2-xiaofei", "H6q-2", None, "thinking", "10-26"),
    ("trajectory_human_6q.json", "H6q-3-shijian", "H6q-3", None, "thinking", "10-26"),
    ("trajectory_human_6q.json", "H6q-4-beilun", "H6q-4", None, "thinking", "10-26"),
    ("trajectory_human_6q.json", "H6q-5-shadui", "H6q-5", None, "thinking", "10-26"),
    ("trajectory_human_6q.json", "H6q-6-jihe", "H6q-6", None, "thinking", "10-26"),
    # 10-27 单题（数学，解出：整除+三位数收敛）
    ("trajectory_human_1027.json", "T-1027-1", "T-1027-1", True, "math", "10-27"),
    # M 系列（数学）
    ("trajectory_human_1027_m1m2.json", "T-M1-1027b", "T-M1", False, "math", "10-27M"),
    ("trajectory_human_1027_m1m2.json", "T-M2-1027b", "T-M2", False, "math", "10-27M"),
    ("trajectory_human_1027_m3457.json", "T-M3-1027c", "T-M3", True, "math", "10-27M"),
    ("trajectory_human_1027_m3457.json", "T-M4-1027c", "T-M4", False, "math", "10-27M"),
    ("trajectory_human_1027_m3457.json", "T-M5-1027c", "T-M5", True, "math", "10-27M"),  # 收敛但结论错
    ("trajectory_human_1027_m3457.json", "T-M7-1027c", "T-M7", True, "math", "10-27M"),
    # 10-31 2/2(1)（数学应用题）
    ("trajectory_human_1031_2.json", "T-NORM-1031", "T-NORM-1031", True, "math", "10-31a"),
    ("trajectory_human_1031_2.json", "T-RATE-1031", "T-RATE-1031", True, "math", "10-31a"),
    ("trajectory_human_1031_2.json", "T-BUS-1031", "T-BUS-1031", True, "math", "10-31a"),
    ("trajectory_human_1031_2_1.json", "T-CUBE1-1031", "T-CUBE1-1031A", True, "math", "10-31a"),
    ("trajectory_human_1031_2_1.json", "T-CUBE2-1031", "T-CUBE2-1031A", True, "math", "10-31a"),
    ("trajectory_human_1031_2_1.json", "T-CORNER-1031", "T-CORNER-1031A", True, "math", "10-31a"),
    ("trajectory_human_1031_2_1.json", "T-SHAPES-1031", "T-SHAPES-1031A", True, "math", "10-31a"),
    ("trajectory_human_1031_2_1.json", "T-FISHTANK-1031", "T-FISHTANK-1031", True, "math", "10-31a"),
    # 10-31 更早批次（数学）
    ("trajectory_human_1031.json", "T-ROPE1-1031", "T-ROPE1-1031", True, "math", "10-31b"),
    ("trajectory_human_1031.json", "T-ROPE2-1031", "T-ROPE2-1031", True, "math", "10-31b"),
    ("trajectory_human_1031.json", "T-PRICE-1031", "T-PRICE-1031", True, "math", "10-31b"),
    ("trajectory_human_1031.json", "T-CREATE-1031", "T-CREATE-1031", True, "math", "10-31b"),
    ("trajectory_human_1031_1.json", "T-CUBE1-1031", "T-CUBE1-1031B", True, "math", "10-31b"),
    ("trajectory_human_1031_1.json", "T-CUBE2-1031", "T-CUBE2-1031B", True, "math", "10-31b"),
    ("trajectory_human_1031_1.json", "T-CORNER-1031", "T-CORNER-1031B", True, "math", "10-31b"),
    ("trajectory_human_1031_1.json", "T-SHAPES-1031", "T-SHAPES-1031B", False, "math", "10-31b"),  # 讲错重讲未收敛
]

_cache = {}
def load_task(fname, tid):
    key = (fname, tid)
    if key in _cache:
        return _cache[key]
    data = json.load(open(BASE / fname, encoding="utf-8"))
    for t in data["trajectories"]:
        if t["id"] == tid:
            _cache[key] = t
            return t
    raise KeyError(f"{tid} not in {fname}")

def main():
    per_task = {}
    for fname, tid, display, solved, typ, batch in METADATA:
        t = load_task(fname, tid)
        a = at.analyze_task(t)
        verdict = at.classify(a["osc_net_change"])
        per_task[display] = {
            "batch": batch, "type": typ, "solved": solved,
            "longest_seg": a["longest_seg"],
            "osc_net_change": round(a["osc_net_change"], 4) if a["osc_net_change"] is not None else None,
            "verdict": verdict,
        }
        print(f"[{display}] type={typ} solved={solved} seg={a['longest_seg']} "
              f"net={per_task[display]['osc_net_change']} verdict={verdict}")

    judgeable = {k: v for k, v in per_task.items() if v["verdict"]}
    cnt_all = Counter(v["verdict"] for v in judgeable.values())
    print("ALL judgeable:", len(judgeable), dict(cnt_all))
    print("ALL majority:", cnt_all.most_common(1)[0] if cnt_all else None)

    # 解出/未解出混淆检验（限数学题）
    math = {k: v for k, v in per_task.items() if v["type"] == "math"}
    grp_solved = Counter(v["verdict"] for v in math.values() if v["solved"] is True and v["verdict"])
    grp_unsolved = Counter(v["verdict"] for v in math.values() if v["solved"] is False and v["verdict"])
    print("MATH solved  verdicts:", dict(grp_solved))
    print("MATH unsolved verdicts:", dict(grp_unsolved))

    # 思辨题方向（对照，非数学）
    think = Counter(v["verdict"] for v in per_task.values() if v["type"] == "thinking" and v["verdict"])
    print("THINKING verdicts:", dict(think))

    out = {
        "method": "S0 真人全样本统一判别(analyze_trajectories) + 解出/未解出混淆检验",
        "n_tasks": len(METADATA),
        "per_task": per_task,
        "all_judgeable_counts": dict(cnt_all),
        "all_majority": cnt_all.most_common(1)[0] if cnt_all else None,
        "confusion_math_solved": dict(grp_solved),
        "confusion_math_unsolved": dict(grp_unsolved),
        "thinking_verdicts": dict(think),
    }
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("saved:", OUT)

if __name__ == "__main__":
    main()
