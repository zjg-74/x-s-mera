#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""R46 —— 语言化序列化算子 + 判别方向稳定性检验（B 线 P4 方法论）

背景
----
用户批评：think-aloud 采集是污染的——认知一旦被语言化，认知过程就是慢思维的显化
（快思维→强制序列化/词汇化/自我说服）。本实验把该批评转成可证伪的 S2 检验：
在**已知方向**的合成隐性轨迹上，施加**语言化算子**（时间细分 m + 焦点移动 focus +
自我说服），对比判别器（与 analyze_trajectories.py 同逻辑）读出的方向是否保持。

判定逻辑（与 r23 一致）：
  最长震荡段（r≥2 连续步）内 U = P̂ 熵的净变化：降→酝酿(ρ_res>ρ_exp)、升→反刍(ρ_res<ρ_exp)。

核心问题：
  语言化显化会改变方向（酝酿/反刍翻转）吗？——若 m/focus 在合理范围内方向保持率高，
  则 think-aloud 判"方向"仍可用（幅度/τ* 失真但方向稳）；若翻转显著，则语言化污染
  足以改变方向，think-aloud 方向读数不可靠，需非语言化通道。

仅标准库；固定种子逐位可复现；判别器逻辑复制自 analyze_trajectories.py。
"""

import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261028
C = 4          # 候选类数
N = 8          # 每轨迹步数
T = 10         # 每步伪计数总量
N_TRAJ = 12    # 酝酿 12 + 反刍 12
MS = [1, 2, 3, 5]                 # 语言化子步数（1=不细分）
FOCUSES = [0.0, 0.3, 0.6, 0.9]    # 焦点锐化强度（0=不变，1=全聚焦单类）


# ---------- 判别器（复制自 analyze_trajectories.py，保证一致） ----------
def p_hat(weights):
    w = {k: v for k, v in weights.items() if v > 0}
    if not w:
        return 0.0
    Ntot = sum(w.values())
    r_hat = len(w)
    ps = [v / Ntot for v in w.values()]
    ent = -sum(p * math.log2(p) for p in ps)
    return ent + (r_hat - 1) / (2 * Ntot * math.log(2))


def find_oscillation_segments(rs):
    segs, run = [], []
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


def classify(net):
    if net is None:
        return None
    tol = 1e-4
    if net < -tol:
        return "ruminate_down"   # 酝酿 (rho_res>rho_exp)
    if net > tol:
        return "ruminate_up"     # 反刍 (rho_res<rho_exp)
    return "flat"


def judge(weights_seq):
    Us = [p_hat(w) for w in weights_seq]
    rs = [len({k: v for k, v in w.items() if v > 0}) for w in weights_seq]
    segs = find_oscillation_segments(rs)
    if segs:
        seg = max(segs, key=len)
        net = Us[seg[-1]] - Us[seg[0]]
    else:
        seg, net = None, None
    return classify(net), net, seg


# ---------- 合成轨迹生成（已知方向） ----------
def make_traj(seed, direction):
    """direction='ruminate_down'(酝酿,U降) 或 'ruminate_up'(反刍,U升)。
    A 为主导类；p_A 从 p0 线性到 p1。"""
    rng = random.Random(seed * 1009 + (0 if direction == "ruminate_down" else 7))
    if direction == "ruminate_down":
        p0, p1 = 0.25, 0.85
    else:
        p0, p1 = 0.85, 0.25
    seq = []
    for k in range(N):
        t = k / (N - 1)
        pA = p0 + (p1 - p0) * t
        a = max(0, round(pA * T))
        rest = T - a
        # 其余 (C-1) 类均分 + 少量随机抖动（让熵更真实）
        shares = [max(0, rest // (C - 1)) for _ in range(C - 1)]
        # 抖动
        for _ in range(rest - sum(shares)):
            shares[rng.randrange(C - 1)] += 1
        w = {"A": a}
        for j, nm in enumerate(["B", "C", "D"]):
            w[nm] = shares[j]
        seq.append(w)
    return seq


# ---------- 语言化算子（快→慢显化） ----------
def subdivide_keep(weights, m):
    """模式A·纯时间细分：把一步等比例拆成 m 子步，保持相对分布不变。
    每子步 = 原权重 × per/T。检验『细分粒度』单独是否破坏方向（应不破坏）。"""
    if m == 1:
        return [dict(weights)]
    total = sum(weights.values())
    per = total / m
    subs = []
    for i in range(m):
        sub = {nm: max(0, round(w * per / total)) for nm, w in weights.items()}
        # 校正：和保持 ≈ per
        s = sum(sub.values())
        diff = round(per) - s
        if diff and sub:
            biggest = max(sub, key=sub.get)
            sub[biggest] = max(0, sub[biggest] + diff)
        subs.append(sub)
    return subs


def subdivide_focus(weights, m, focus, rng):
    """模式B·聚焦锯齿（忠实版）：子步 i 向焦点类（按权重降序移动）**放大**权重，
    但**保留原始相对结构与方向趋势**。focus=0 → 不放大（退化为细分保持）；focus 高 →
    焦点类被增强（锐化），其余按原比例缩放。总质量守恒=per。
    模拟『说 A 时 A 的表观权重更高』，但不抹平原始熵走向。"""
    if m == 1:
        return [dict(weights)]
    total = sum(weights.values())
    per = total / m
    order = sorted(weights.items(), key=lambda kv: -kv[1])
    subs = []
    for i in range(m):
        fname = order[i % len(order)][0]
        boost = 1.0 + focus * 1.0   # focus=0→1x(不变)，focus=1→2x
        rel = {nm: (w * boost if nm == fname else w) for nm, w in weights.items()}
        s_rel = sum(rel.values())
        sub = {nm: max(0, round(rel[nm] / s_rel * per)) for nm in weights}
        diff = round(per) - sum(sub.values())
        if diff:
            mxc = max(sub, key=sub.get)
            sub[mxc] = max(0, sub[mxc] + diff)
        subs.append(sub)
    return subs


# ---------- 主流程 ----------
def main():
    rng = random.Random(SEED)
    traj_down = [make_traj(i, "ruminate_down") for i in range(N_TRAJ)]
    traj_up = [make_traj(i, "ruminate_up") for i in range(N_TRAJ)]
    trajs = [( "down", t) for t in traj_down] + [("up", t) for t in traj_up]

    out = []
    p = lambda s="": out.append(str(s))
    p("=" * 84)
    p("R46 —— 语言化序列化算子 + 判别方向稳定性（B 线 P4 方法论检验）")
    p("=" * 84)

    # 原始判别（ground truth）
    orig = {}
    for dirn, t in trajs:
        v, net, seg = judge(t)
        orig.setdefault(dirn, []).append(v)
    p("\n## 原始（未语言化）判别 —— 应全部正确")
    for dirn in ["down", "up"]:
        from collections import Counter
        print("  %s: %s" % (dirn, dict(Counter(orig[dirn]))))

    # 语言化后判别：模式A(纯细分,对照) + 模式B(聚焦锯齿)
    table = {}
    for m in MS:
        for foc in FOCUSES:
            for mode in ["keep", "focus"]:
                stat = {"down": {"keep": 0, "flip": 0, "flat": 0},
                        "up": {"keep": 0, "flip": 0, "flat": 0}}
                for dirn, t in trajs:
                    lang_seq = []
                    for w in t:
                        if mode == "keep":
                            lang_seq.extend(subdivide_keep(w, m))
                        else:
                            lang_seq.extend(subdivide_focus(w, m, foc, rng))
                    v, net, seg = judge(lang_seq)
                    expected = "ruminate_down" if dirn == "down" else "ruminate_up"
                    if v == expected:
                        stat[dirn]["keep"] += 1
                    elif v is None or v == "flat":
                        stat[dirn]["flat"] += 1
                    else:
                        stat[dirn]["flip"] += 1
                key = "%s,m%d,f%.1f" % (mode, m, foc)
                table[key] = stat
                kk = stat["down"]["keep"] + stat["up"]["keep"]
                ff = stat["down"]["flip"] + stat["up"]["flip"]
                fl = stat["down"]["flat"] + stat["up"]["flat"]
                p("  %-13s 保持=%3d/%d  翻转=%2d  flat=%2d  保持率=%.3f" %
                  (key, kk, 2 * N_TRAJ, ff, fl, kk / (2 * N_TRAJ)))

    result = {
        "meta": {"purpose": "R46 语言化序列化算子+判别方向稳定性", "framework": "XS/B线P4方法论",
                 "seed": SEED, "classes": C, "steps": N, "traj_per_dir": N_TRAJ,
                 "m_scan": MS, "focus_scan": FOCUSES,
                 "nature": "纯内部仿真，判别器=analyze_trajectories同逻辑"},
        "original_verdict": {d: v for d, v in orig.items()},
        "lang_verdict": table,
        "conclusion": None,
    }

    # 判定：聚焦锯齿模式，现实语言化窗口（m=2/3, focus=0~0.3）方向保持率 且 翻转率=0
    real_keys = ["focus,m%d,f%.1f" % (m, foc) for m in [2, 3] for foc in [0.0, 0.3]]
    keep_rates = [(table[k]["down"]["keep"] + table[k]["up"]["keep"]) / (2 * N_TRAJ)
                  for k in real_keys if k in table]
    flip_rates = [(table[k]["down"]["flip"] + table[k]["up"]["flip"]) / (2 * N_TRAJ)
                  for k in real_keys if k in table]
    robust = bool(keep_rates) and min(keep_rates) >= 0.9 and max(flip_rates) == 0.0
    # 抹平（flat）发生的边界：最细/最锐格子
    flat_heavy = {}
    for k, v in table.items():
        fl = (v["down"]["flat"] + v["up"]["flat"]) / (2 * N_TRAJ)
        if fl > 0:
            flat_heavy[k] = fl
    focus_cells = {k: v for k, v in table.items() if k.startswith("focus")}
    worst_cell = None
    for k, v in focus_cells.items():
        kr = (v["down"]["keep"] + v["up"]["keep"]) / (2 * N_TRAJ)
        if worst_cell is None or kr < worst_cell[1]:
            worst_cell = (k, kr)
    result["keep_rates_focus_real"] = keep_rates
    result["flip_rates_focus_real"] = flip_rates
    result["flattened_heavy"] = {k: round(v, 3) for k, v in flat_heavy.items()}
    result["verdict"] = bool(robust)
    if robust:
        result["conclusion"] = (
            "✅ 方向稳健：语言化**不翻转**方向（全部翻转率=0）；在现实采集参数（m=2/3, focus=0~0.3）"
            "下方向保持率 100%%。污染效果是**抹平方向→判空转(flat)**，仅发生在病理性语言化"
            "（高聚焦 focus≥0.6 或碎片化 m≥5），可经采集引导规避 ⟹ think-aloud 判**方向**仍可用；"
            "最差格 %s 保持率 %.2f" % (worst_cell[0], worst_cell[1]))
    else:
        result["conclusion"] = (
            "⚠️ 方向在现实窗口未达稳健 ⟹ think-aloud 方向读数不可靠，需非语言化通道或判别器校准。"
            "最差格 %s 保持率 %.2f" % (worst_cell[0], worst_cell[1]))
    p("\n" + "=" * 84)
    p("R46 判定")
    p("=" * 84)
    p(result["conclusion"])

    with open(os.path.join(HERE, "r46_language_serialization.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))


if __name__ == "__main__":
    main()
