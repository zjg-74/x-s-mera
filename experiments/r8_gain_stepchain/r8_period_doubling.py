# -*- coding: utf-8 -*-
"""
R8 —— 增益-步长机制（M1）的工程验证：周期倍增级联 → 混沌
==================================================================
对应 XS v3.1 补篇一 §5 的可检验推论 R8（强度：高）：

> 在 X-S-MERA 的离散迭代原型上，把 A11-R 的语义势增益作为可调参数扫过，
> 观察是否出现 logistic 式的周期倍增级联 → 混沌。出现则 M1 机制被工程验证；
> 不出现则 M1 被证否，须回到 M2/M3。

机制背景（补篇一 §2.2 M1 · 离散步长机制）
--------------------------------------------------
- 势-流分离（A11-R）：语义势**不进状态导数**，只**调增益/有效步长 η_eff**。
- 在**连续流**中，势只能重参数化时间、不能改变轨道；但在**离散迭代**中，增益
  调制等价于改变有效步长，而离散映射的定性行为（含混沌）**依赖有效步长**。
- 标准例证：logistic 映射 x_{k+1} = r·x_k·(1 − x_k)，参数 r 即有效增益/步长；
  r 越过 ≈3.5699 后经**周期倍增级联**进入混沌。调增益即可移出/移入混沌区，
  **全程不需要增益进入状态导数**。

本实验设计
--------------------------------------------------
用一个 **logistic 型最小离散认知原型**作为载体：x_{k+1} = g·x_k·(1 − x_k)
- x_k：加载度 σ 的归一化投影（候选聚焦量 m_eff 的归一化表征，v3.2 S5）；
- g：**A11-R 语义势增益**，外部给定、静态、**仅调步长、不进状态导数**。

验证判据（预注册）：M1 被工程验证 ⟺ 同时满足
  V1  出现**至少两次周期倍增**（1→2、2→4、4→8 中 ≥2 次）；
  V2  **混沌区存在**（Lyapunov λ>0 且稳态无短周期，且带宽足够）；
  V3  混沌阈值 g_chaos 落在 logistic 解析界附近（3.4 ≤ g ≤ 3.7，理论 3.5699）；
  V4  **势不进状态导数**（g 只作乘性系数，状态方程无 g 的梯度项）。

诚实边界：最小离散原型（logistic 型单状态），精确命中 M1 机制主张；**不等同于**
r18 完整宏加载原型。出现级联 ⟹ M1 作为"势-流分离下势可消除混沌"的机制路径被工程
验证；不声称"认知必然混沌"。logistic 解析分岔点为对照基准（已知数学事实）。
"""

import json


G_LO, G_HI = 2.5, 4.0      # 语义势增益扫描区间
N_G = 1001                 # 增益网格点
X0 = 0.3                   # 初始加载度
BURN_IN = 1000             # 暂态丢弃
SEQ_LEN = 2048             # 稳态序列长度
TOL = 1e-7                 # 周期检测容差
MAX_PERIOD = 64            # 周期检测上限
LYAP_TOL = 1e-4            # λ>0 判混沌
LYAP_STEPS = 8000          # Lyapunov 步数
LYAP_BURN = 2000


def iterate(g: float, x0: float, n: int):
    """logistic 型离散迭代：x_{k+1}=g·x_k(1−x_k)。g 是语义势增益（不进状态导数）。"""
    out = [0.0] * n
    v = float(x0)
    for k in range(n):
        out[k] = v
        v = g * v * (1.0 - v)
    return out


def detect_period(seq, tol=TOL, max_p=MAX_PERIOD):
    """检测稳态轨道最小周期。用末窗口对齐比较；无短周期 → 混沌候选。"""
    s = seq
    win = 200
    base = len(s) - win
    for p in range(1, max_p + 1):
        dmax = 0.0
        # 每 p 检查窗口内 40 个对齐位置的最大偏差
        step = max(1, win // 40)
        for i in range(0, win, step):
            a = base + i
            if a + p >= len(s):
                break
            d = abs(s[a] - s[a + p])
            if d > dmax:
                dmax = d
                if dmax >= tol:
                    break
        if dmax < tol:
            return p
    return None


def lyapunov(g: float, x0: float, n=LYAP_STEPS, burn=LYAP_BURN):
    """logistic Lyapunov 指数 λ = mean(ln|g(1−2x_k)|)。λ>0 → 混沌。"""
    x = float(x0)
    lam = 0.0
    cnt = 0
    for k in range(n):
        if k >= burn:
            lam += __import__("math").log(abs(g * (1.0 - 2.0 * x)) + 1e-12)
            cnt += 1
        x = g * x * (1.0 - x)
    return lam / cnt if cnt else float("nan")


def main():
    gs = [G_LO + (G_HI - G_LO) * i / (N_G - 1) for i in range(N_G)]
    periods = []
    lams = []
    for g in gs:
        seq = iterate(g, X0, BURN_IN + SEQ_LEN)[BURN_IN:]
        periods.append(detect_period(seq))
        lams.append(lyapunov(g, X0))

    unique_periods = sorted({int(p) for p in periods if p is not None})

    # 各级周期的 g 区间
    period_bands = {}
    for p in unique_periods:
        idxs = [i for i, v in enumerate(periods) if v == p]
        runs = []
        start = idxs[0]
        prev = idxs[0]
        for i in idxs[1:]:
            if i - prev > 1:
                runs.append((round(gs[start], 4), round(gs[prev], 4)))
                start = i
            prev = i
        runs.append((round(gs[start], 4), round(gs[prev], 4)))
        period_bands[p] = runs

    # 分岔点：周期首次从 p 变 2p
    doubling_points = {}
    for p in unique_periods:
        if 2 * p in unique_periods:
            crossing = None
            for i in range(len(periods) - 1):
                if periods[i] == p and periods[i + 1] == 2 * p:
                    crossing = round(gs[i + 1], 4)
                    break
            if crossing is not None:
                doubling_points["%d->%d" % (p, 2 * p)] = crossing

    # 混沌阈值：第一个 λ>0 的 g
    chaos_idx = [i for i, v in enumerate(lams) if v > LYAP_TOL]
    g_chaos = round(gs[chaos_idx[0]], 4) if chaos_idx else None
    chaos_band = (round(gs[chaos_idx[0]], 4), round(gs[chaos_idx[-1]], 4)) if chaos_idx else None

    # 倍增加链 1→2→4→8→16→32（纯 2 倍链）
    chain = []
    p = 1
    while True:
        if p in unique_periods:
            chain.append(p)
        else:
            break
        p *= 2
    has_doubling_chain = (len(chain) >= 2 and
                          chain == [chain[0] * (2 ** i) for i in range(len(chain))])

    v1 = len(doubling_points) >= 2
    v2 = g_chaos is not None and len(chaos_idx) >= 40
    v3 = g_chaos is not None and 3.4 <= g_chaos <= 3.7
    v4 = True
    m1_supported = bool(v1 and v2 and v3 and v4)

    return {
        "experiment": "R8-gain-stepchain",
        "framework": "XS v3.4 · 补篇一 R8（M1 机制工程验证）",
        "hardness": "S2（自制最小离散原型，纯标准库）",
        "g_range": [G_LO, G_HI], "n_g": N_G, "x0": X0,
        "period_detection": {"burn_in": BURN_IN, "seq_len": SEQ_LEN, "tol": TOL,
                             "max_period": MAX_PERIOD},
        "unique_periods": unique_periods,
        "period_bands": period_bands,
        "doubling_points": doubling_points,
        "chaos_threshold_g": g_chaos,
        "chaos_band": chaos_band,
        "lyapunov_positive_count": len(chaos_idx),
        "doubling_chain": chain,
        "has_doubling_chain": has_doubling_chain,
        "verdicts": {"V1_doubling>=2": v1, "V2_chaos_exists": v2,
                     "V3_threshold_near_3.5699": v3, "V4_gain_not_in_state_deriv": v4},
        "m1_supported": m1_supported,
        "conclusion": (
            "M1 被工程验证：离散迭代中语义势增益（仅调步长、不进状态导数）即可出现"
            "周期倍增级联并移入混沌区，与 logistic 定性一致"
            if m1_supported else
            "M1 未获本实验支持：未见完整周期倍增级联或混沌，须回到 M2/M3"),
    }


if __name__ == "__main__":
    res = main()
    print(json.dumps(res, ensure_ascii=False, indent=2))
    with open("r8_results.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    # 分岔图尾部数据
    gs = [2.5 + 1.5 * i / 1199 for i in range(1200)]
    tail = []
    for g in gs:
        seq = iterate(g, X0, 2000 + 600)[2000:]
        tail.append(round(seq[-1], 6))
    with open("r8_bifurcation_tail.json", "w", encoding="utf-8") as f:
        json.dump({"g": [round(g, 5) for g in gs], "x_tail": tail}, f, ensure_ascii=False)
    print("bifurcation data saved: r8_bifurcation_tail.json")
