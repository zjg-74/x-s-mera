#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R28：微观回落的本体性、峰的尺度依赖、多峰（重开）、局部/全局收敛。

来源主张（用户 2026-10-19）
--------------------------
1. 本文中的"增 / 减"如果是基于统计的，那么**都可能存在微观回落**。
2. 有效认知 = 收敛，但**不保证排除伪收敛**。
3. 收敛可以是**过程的、局部的**，而全局可能是不收敛的。
4. 峰在理论上一定存在：X 从无到一个值，再收敛到一个认知。
   强调"有效认知"是因为**有效认知保证收敛 ⟹ 保证减相**；发散不一定有减的过程。

本实验只做前三条的**机制检验**（第 4 条是结构性论证，见补篇七 §1）。

与 R27 的关键差别（本轮最重要的建模改动）
----------------------------------------
R27 把"微观回落"当成**观测噪声**：dU 是确定性的，噪声只加在观测量上
（`Uobs = U + N(0, obs_sd)`）。于是"单峰性在 w=1 时检出率 0.045"被记为
一个**测量问题**（B-B27）。

用户第 1 条指出这是错的：如果增/减**基于统计**（a、b 是统计速率），
那么 dU 本身就是随机变量，**微观回落是过程的一部分，不是测量误差**。
本实验因此把噪声从观测层移到**过程层**：

    r(k) = r_bar(k) * exp(s*Z - s^2/2)        （乘性对数正态，E[r] = r_bar）

这直接给出一个**精确可解析**的预测（不需要正态近似）：

    dU < 0  <=>  r < theta_U
    P(微观回落) = Phi( (ln theta_U - ln r_bar + s^2/2) / s )

注意两点非平凡之处：
  * 在峰处（r_bar = theta_U）回落率**不是 1/2，而是 Phi(s/2) > 1/2**
    （因为中位数 r = r_bar * exp(-s^2/2) < r_bar）。这是尖锐的、可证否的预测。
  * 回落率随 |ln r_bar - ln theta_U| 增大而单调下降。

尺度问题因此被重新表述：**"单峰"是 E 层（统计层）命题，不是路径层命题**。
路径层的峰数依赖声明的粗粒化尺度 (w, h)（h = 显著度 / prominence）。

四大臂
------
A1  微观回落率：实测 vs 解析式 Phi((ln th - ln r_bar + s^2/2)/s)，扫 s
A2  峰数与单峰率 vs 粗粒化尺度 (w, h)；验证 E 层峰位不受微观回落影响
B   多峰（r 调度非单调）：**峰数 = 1 + 由减转增的相变次数**
C   m 过程的生灭模型（**不用 A_D**）：局部收敛 / 重开 / 全局收敛 / 微观回落

诚实声明
--------
* 本实验是**机制自洽性检验**，不是实证验证。所有系数与调度都是人造的。
* Arm C 的生灭模型是**辅助模型**，用于检验不依赖 A_D 的结构性主张；
  它本身不是 XS 的组成部分，也未与 A_D 对齐。
* 本实验**不触碰停止规则 / 接受阈值 / 判定函数**（R18-b 的死因）。
  Arm C 刻意**不加**停止规则，以便让"重开"可见——这是声明式分析选择，
  不是拟合结果。加了停止规则后重开会被截断，这一点在报告里写明。
* 陷阱（沿用 R27 的 T1-T4，全部已规避）：
  T1 用内置 hash() 做种子（str hash 按进程随机化）——本脚本不用 hash()
  T2 峰位点值比对（并列取小索引伪像）——用平台比对
  T3 用平台下界估峰位（噪声下系统性下偏）——用平滑序列 argmax
  T4 平滑后仍用固定容差——容差取 3*sd/sqrt(w)
"""

import json
import math
import random
import statistics as st
import zlib

# ------------------------------------------------------------------ 参数

K = 241
DT = 0.01
R_HI = 9.0
R_LO = 1.0 / 9.0
LN_HI = math.log(R_HI)
LN_LO = math.log(R_LO)

U_START = 5.0
SEED = 20261019

N_REP_A1 = 4000      # 回落率需要大样本（每步一个伯努利）
N_REP_A2 = 400
N_REP_C = 400

S_GRID = [0.10, 0.20, 0.40, 0.80]
WINDOWS = [1, 5, 11, 21, 41]
PROM_GRID = [0.0, 0.05, 0.10, 0.20, 0.40]

TAU_C_GRID = [1, 3, 5, 10]   # 局部收敛的持续步数门槛（声明式，做敏感性分析）


def make_config(name, e_plus, e_minus, u_plus, u_minus, note=""):
    return {
        "name": name,
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": e_minus / e_plus,
        "theta_U": u_minus / u_plus,
        "rho_exp": u_plus / e_plus,
        "rho_res": u_minus / e_minus,
        "note": note,
    }


CONFIG = make_config("res_efficient", 1.0, 1.0, 1.0, 3.0,
                     "theta_E=1 < theta_U=3, rho_res=3 > rho_exp=1")


# ------------------------------------------------------------------ 工具


def zlib_off(s):
    return zlib.crc32(s.encode("utf-8")) % 9973


def Phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def r_ramp(k):
    t = k / (K - 1.0)
    return math.exp(LN_HI + (LN_LO - LN_HI) * t)


def r_multi(k):
    """三段调度：降 -> 升 -> 降。用于制造多峰（B-B29）。"""
    if k < 80:
        t = k / 80.0
        return math.exp(LN_HI + (LN_LO - LN_HI) * t)
    if k < 160:
        t = (k - 80) / 80.0
        return math.exp(LN_LO + (LN_HI - LN_LO) * t)
    t = (k - 160) / (K - 1.0 - 160)
    return math.exp(LN_HI + (LN_LO - LN_HI) * t)


def predicted_kstar(theta, r_of_k):
    """数值求解 r(k*) = theta（二分）。仅用于单调段。"""
    if r_of_k(0) < theta or r_of_k(K - 1) > theta:
        return None
    lo, hi = 0.0, float(K - 1)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if r_of_k(mid) > theta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def dU_det(cfg, r):
    return DT * (r * cfg["u_plus"] - cfg["u_minus"])


def sample_r(rng, rbar, s):
    if s <= 0:
        return rbar
    return rbar * math.exp(s * rng.gauss(0.0, 1.0) - 0.5 * s * s)


def run_traj(rng, cfg, r_of_k, s=0.0, u_start=U_START):
    """过程层噪声版轨迹（无观测噪声——观测噪声是本轮要分离出去的量）。"""
    U = u_start
    out = []
    for k in range(K):
        r = sample_r(rng, r_of_k(k), s)
        U = U + dU_det(cfg, r)
        out.append(U)
    return out


# ---- 平滑 / 峰 / 单调性（沿用 R27，含 T2/T3/T4 的修正） ----

def moving_average(seq, w):
    if w <= 1:
        return list(seq)
    n = len(seq)
    half = w // 2
    out = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out.append(sum(seq[lo:hi]) / (hi - lo))
    return out


def argmax_of(seq):
    best_i, best_v = 0, seq[0]
    for i, v in enumerate(seq):
        if v > best_v:
            best_v, best_i = v, i
    return best_i


def argmax_plateau(seq, tol=1e-9):
    best_v = max(seq)
    idx = [i for i, v in enumerate(seq) if v >= best_v - tol]
    return idx[0], idx[-1]


def is_unimodal(seq, tol=0.0):
    p = argmax_of(seq)
    if p <= 0 or p >= len(seq) - 1:
        return False, p
    for i in range(1, p + 1):
        if seq[i] < seq[i - 1] - tol:
            return False, p
    for i in range(p + 1, len(seq)):
        if seq[i] > seq[i - 1] + tol:
            return False, p
    return True, p


def local_maxima(seq):
    """返回平台合并后的局部极大索引列表（端点不算内点峰）。"""
    n = len(seq)
    out = []
    i = 1
    while i < n - 1:
        if seq[i] >= seq[i + 1]:
            j = i
            while j + 1 < n and seq[j + 1] == seq[i]:
                j += 1
            left_ok = seq[i] > seq[i - 1]
            right_ok = (j + 1 < n and seq[j + 1] < seq[i]) or (j + 1 >= n)
            if left_ok and right_ok:
                out.append((i + j) // 2)
            i = j + 1
        else:
            i += 1
    return out


def prominence(seq, i):
    """标准 prominence（拓扑持久性）：
    从 i 向左/右下降，直到遇到更高点或边界，取两侧最低谷中较高者。"""
    h = seq[i]
    n = len(seq)
    left_min = h
    j = i
    while j >= 0:
        if seq[j] > h:
            break
        left_min = min(left_min, seq[j])
        j -= 1
    right_min = h
    j = i
    while j < n:
        if seq[j] > h:
            break
        right_min = min(right_min, seq[j])
        j += 1
    return h - max(left_min, right_min)


def count_peaks(seq, min_prom):
    return sum(1 for i in local_maxima(seq) if prominence(seq, i) > min_prom)


def local_minima(seq):
    return local_maxima([-v for v in seq])


def count_valleys(seq, min_prom):
    """回跳（由减转增）的计数：对偶于峰，用谷的 prominence 过滤。"""
    neg = [-v for v in seq]
    return sum(1 for i in local_maxima(neg) if prominence(neg, i) > min_prom)


def last_extremum_is_max(seq):
    """最后一个内点极值是极大吗？

    为什么必须判：**峰数 = 回跳数 + 1 只在过程止于下降段时成立**。
    若过程止于上升段，末端极大落在序列边界上，而 local_maxima 不计端点，
    于是恒等式退化为 峰数 = 回跳数。这不是公式错，是端点约定。
    """
    mx = local_maxima(seq)
    mn = local_minima(seq)
    if not mx:
        return False
    if not mn:
        return True
    return max(mx) > max(mn)


def tol_for(sd, w):
    """仅供对照。注意：**不要在过程层噪声下用它做主判据**——
    容差是按"每步"标定的，而随机游走的累积偏差远大于每步量级，
    会导致"既增又减"从而 100% 误判为单峰（本轮实测到的假阳性）。
    主判据一律用 peak_count == 1（不依赖容差）。"""
    if sd <= 0:
        return 1e-9
    return 3.0 * sd / math.sqrt(max(1, w))


# ------------------------------------------------------------------ 臂 A1


def arm_A1_fallback(cfg, s, n_rep=N_REP_A1, seed=SEED):
    """实测微观回落率 vs 解析式 Phi((ln th - ln r_bar + s^2/2)/s)。"""
    th = cfg["theta_U"]
    rng = random.Random(seed + zlib_off("A1") + int(s * 1000))
    cnt = [0] * K
    for _ in range(n_rep):
        for k in range(K):
            r = sample_r(rng, r_ramp(k), s)
            if r < th:
                cnt[k] += 1
    rows = []
    for k in range(K):
        rb = r_ramp(k)
        ana = Phi((math.log(th) - math.log(rb) + 0.5 * s * s) / s)
        emp = cnt[k] / float(n_rep)
        rows.append({"k": k, "r_bar": round(rb, 5),
                     "emp": round(emp, 4), "ana": round(ana, 4),
                     "dev": round(emp - ana, 4)})
    devs = [abs(r["dev"]) for r in rows]
    kstar = predicted_kstar(th, r_ramp)
    k_at = int(round(kstar)) if kstar is not None else None
    # MC 标准误（伯努利，p 最坏 0.5）
    se = 0.5 / math.sqrt(n_rep)
    return {
        "config": cfg["name"], "s": s, "n_rep": n_rep,
        "theta_U": th,
        "kstar_pred": None if kstar is None else round(kstar, 2),
        "max_abs_dev": round(max(devs), 4),
        "mean_abs_dev": round(st.mean(devs), 4),
        "mc_se": round(se, 4),
        "within_3se": bool(max(devs) <= 3.0 * se),
        "fallback_at_kstar": None if k_at is None else rows[k_at]["emp"],
        "analytic_at_kstar": None if k_at is None else rows[k_at]["ana"],
        "predicted_Phi_s_over_2": round(Phi(0.5 * s), 4),
        "rows_sample": rows[::20],
    }


# ------------------------------------------------------------------ 臂 A2


def arm_A2_scale(cfg, s, n_rep=N_REP_A2, seed=SEED):
    """峰数 / 单峰率 / 峰位 vs 粗粒化尺度 (w, h)。噪声是过程层的。"""
    th = cfg["theta_U"]
    kstar = predicted_kstar(th, r_ramp)
    rng = random.Random(seed + zlib_off("A2") + int(s * 1000))

    # 残余抖动的量级：单步 dU 的 sd
    # Var(dU) = (DT*u_plus)^2 * Var(r) ; Var(r) = r_bar^2 (e^{s^2}-1)
    u_p = cfg["u_plus"]
    sd_unit = DT * u_p * math.exp(LN_HI) * math.sqrt(math.exp(s * s) - 1.0)

    paths = []
    for _ in range(n_rep):
        paths.append(run_traj(rng, cfg, r_ramp, s=s))

    out_rows = []
    for w in WINDOWS:
        row = {"w": w}
        for h in PROM_GRID:
            peaks, uni = [], []
            for p in paths:
                sm = moving_average(p, w)
                n_pk = count_peaks(sm, h)
                peaks.append(n_pk)
                # 主判据：单峰 <=> 显著峰数恰为 1。不依赖任何容差。
                uni.append(1 if n_pk == 1 else 0)
            row["peaks_h_" + str(h)] = round(st.mean(peaks), 3)
            row["uni_h_" + str(h)] = round(sum(uni) / len(uni), 3)
            row["valleys_h_" + str(h)] = round(st.mean(
                [count_valleys(moving_average(p, w), h) for p in paths]), 3)
        # 峰位：用平滑序列 argmax（陷阱 T3：不能用平台下界）
        ps = []
        for p in paths:
            sm = moving_average(p, w)
            ps.append(argmax_of(sm))
        row["peak_mean"] = round(st.mean(ps), 2)
        row["peak_sd"] = round(st.pstdev(ps), 2)
        out_rows.append(row)

    # 单调性检查：峰数随 w 单调不增、随 h 单调不增
    mono_w = all(out_rows[i]["peaks_h_0.0"] >= out_rows[i + 1]["peaks_h_0.0"] - 1e-9
                 for i in range(len(out_rows) - 1))
    hs = ["peaks_h_" + str(x) for x in PROM_GRID]
    mono_h = all(row[hs[i]] >= row[hs[i + 1]] - 1e-9
                 for row in out_rows for i in range(len(hs) - 1))

    return {
        "config": cfg["name"], "s": s, "n_rep": n_rep,
        "kstar_pred": None if kstar is None else round(kstar, 2),
        "sd_dU_per_step": round(sd_unit, 5),
        "monotone_in_w": mono_w,
        "monotone_in_h": mono_h,
        "rows": out_rows,
    }


# ------------------------------------------------------------------ 臂 B


def arm_B_multipeak(cfg, s=0.0, n_rep=N_REP_A2, seed=SEED):
    """r 调度非单调 -> 多峰。检验 峰数 = 1 + #(由减转增的相变)。"""
    th = cfg["theta_U"]
    # E 层的相变（deterministic）
    signs = []
    for k in range(K - 1):
        d = dU_det(cfg, r_multi(k))
        signs.append(0 if abs(d) < 1e-15 else (1 if d > 0 else -1))
    # 统计由 - 到 + 的次数（回跳 / 重开）
    reopen = 0
    prev = signs[0]
    for sg in signs[1:]:
        if prev == -1 and sg == 1:
            reopen += 1
        if sg != 0:
            prev = sg
    det = run_traj(random.Random(seed), cfg, r_multi, s=0.0)
    det_sm = moving_average(det, 1)
    pred_peaks = count_peaks(det_sm, 0.0)
    maxima = local_maxima(det_sm)
    minima = []
    for i in range(1, len(det_sm) - 1):
        if det_sm[i] < det_sm[i - 1] and det_sm[i] <= det_sm[i + 1]:
            minima.append(i)

    out = {
        "config": cfg["name"],
        "theta_U": th,
        "deterministic_peak_count": pred_peaks,
        "deterministic_peak_idx": [int(i) for i in maxima],
        "deterministic_trough_idx": [int(i) for i in minima[:5]],
        "reopen_count_sign_flips": reopen,
        "identity_peaks_eq_1_plus_reopen": bool(pred_peaks == 1 + reopen),
    }

    # 带过程噪声的复核
    if s > 0:
        rng = random.Random(seed + zlib_off("B") + int(s * 1000))
        counts = []
        for _ in range(n_rep):
            p = run_traj(rng, cfg, r_multi, s=s)
            sm = moving_average(p, 21)
            counts.append(count_peaks(sm, 0.05))
        out["noisy_s"] = s
        out["noisy_peak_count_mean_w21_h005"] = round(st.mean(counts), 3)
    return out


# ------------------------------------------------------------------ 臂 C


def poisson(rng, lam):
    if lam <= 0:
        return 0
    if lam > 30:
        return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))
    L = math.exp(-lam)
    kk = 0
    p = 1.0
    while True:
        p *= rng.random()
        if p <= L:
            return kk
        kk += 1


def m_schedule_multi(k):
    """m 过程的 r_bar 调度：9 -> 1 -> 4 -> 1（制造一次重开）。"""
    if k < 80:
        t = k / 80.0
        return math.exp(LN_HI + (0.0 - LN_HI) * t)
    if k < 140:
        t = (k - 80) / 60.0
        return math.exp(0.0 + (math.log(4.0) - 0.0) * t)
    return 4.0 * math.exp((0.0 - math.log(4.0)) * min(1.0, (k - 140) / 60.0))


Q_DEATH = 0.05    # 每个候选被消解的概率（质量作用）


def run_m_process(rng, r_of_k, q=Q_DEATH, K_=K):
    """m = |A_k| 的生灭模型：加 1 ~ Bernoulli(p)，减 ~ Binomial(m, q)。

    E[dm | m] = p - q*m，故 m* = p/q = r_bar：**m 以时间常数 1/q 追踪 r_bar**。
    这个模型的意义：微观回落是**内生的**（Binomial 与 Bernoulli 都是随机的），
    不需要任何外加噪声。它是"基于统计的增减"的最小实现。
    """
    m = 0
    path = []
    for k in range(K_):
        rb = r_of_k(k)
        p = min(1.0, q * rb)
        add = 1 if rng.random() < p else 0
        rem = 0
        for _ in range(m):
            if rng.random() < q:
                rem += 1
        m = m + add - rem
        if m < 0:
            m = 0
        path.append(m)
    return path


def analyse_m_runs(path, tau_c):
    """局部收敛 / 重开 / 全局收敛 的事件计数。

    局部收敛段：连续的 m == 1，长度 >= tau_c
    重开：从 m == 1 的段离开到 m != 1
    全局收敛：终端段是 m == 1 且长度 >= tau_c
    """
    n = len(path)
    runs = []   # (start, end, value)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and path[j + 1] == path[i]:
            j += 1
        runs.append((i, j, path[i]))
        i = j + 1

    conv_runs = [(a, b) for (a, b, v) in runs if v == 1 and (b - a + 1) >= tau_c]
    # 重开：一个 m==1 段（>= tau_c）之后又出现 m >= 2
    reopen = 0
    for (a, b) in conv_runs:
        if b + 1 < n and path[b + 1] >= 2:
            reopen += 1
    terminal_global = bool(conv_runs and conv_runs[-1][1] == n - 1)
    empty_runs = sum(1 for (a, b, v) in runs if v == 0)
    return {
        "n_conv_runs": len(conv_runs),
        "n_reopen": reopen,
        "global_converged": terminal_global,
        "n_empty_runs": empty_runs,
        "m_max": max(path),
        "m_final": path[-1],
    }


def arm_C_mprocess(n_rep=N_REP_C, seed=SEED):
    """不依赖 A_D 的结构检验：局部/全局收敛 + 峰-重开对偶 + 内生微观回落。"""
    out = {"q": Q_DEATH, "n_rep": n_rep, "tau_c_grid": TAU_C_GRID}
    rng = random.Random(seed + zlib_off("C"))
    paths = [run_m_process(rng, m_schedule_multi) for _ in range(n_rep)]

    # (i) 内生微观回落率：在 m 还在上升的阶段，dm < 0 的比例
    down = tot = 0
    for p in paths:
        for k in range(1, 60):
            if p[k] < p[k - 1]:
                down += 1
            tot += 1
    out["micro_reversal_rate_expansion_phase"] = round(down / float(tot), 4)

    # (ii) 局部 / 全局收敛（扫 tau_c）
    table = []
    for tau in TAU_C_GRID:
        conv = glob = reop = 0
        ncr = 0
        for p in paths:
            a = analyse_m_runs(p, tau)
            conv += 1 if a["n_conv_runs"] > 0 else 0
            glob += 1 if a["global_converged"] else 0
            reop += a["n_reopen"]
            ncr += a["n_conv_runs"]
        table.append({
            "tau_c": tau,
            "frac_has_local_convergence": round(conv / float(n_rep), 3),
            "frac_global_converged": round(glob / float(n_rep), 3),
            "mean_conv_runs": round(ncr / float(n_rep), 3),
            "mean_reopen": round(reop / float(n_rep), 3),
        })
    out["convergence_table"] = table

    # (iii) 峰 / 回跳 / 重开 三者分离
    #   * 峰数、回跳数是**拓扑量**，在声明尺度 (w,h) 上定义，恒有 峰数 = 回跳数 + 1
    #     （本模型的轨迹起于上升、止于下降，故 #max = #min + 1）
    #   * 重开数是**结构量**：回落必须真正到达 m == 1（已决态）才算。
    #     因此 重开 ⊆ 回跳，二者**不等价**：可以无重开却多峰（未达已决态的震荡）。
    # 拓扑恒等式只在**未过滤**极值上成立（min_prom = -1 保留全部）。
    # 显著度过滤会破坏它：两个可比峰被一个浅谷分隔时，谷被滤掉而两峰都留下，
    # 于是 n_+ 不再等于 n_- + 1。这一点必须写明，否则"峰数"这个量没有定义。
    dual = []
    for w in [5, 11, 21]:
        pk = vl = ident = ident_desc = n_desc = alt = 0
        pk_f = vl_f = 0
        for p in paths:
            sm = [float(x) for x in moving_average([float(v) for v in p], w)]
            a = count_peaks(sm, -1.0)
            b = count_valleys(sm, -1.0)
            pk += a
            vl += b
            alt += 1 if abs(a - b) <= 1 else 0
            ident += 1 if a == b + 1 else 0
            if last_extremum_is_max(sm):
                n_desc += 1
                ident_desc += 1 if a == b + 1 else 0
            pk_f += count_peaks(sm, 0.5)
            vl_f += count_valleys(sm, 0.5)
        n = float(n_rep)
        dual.append({"w": w,
                     "mean_peaks_unfiltered": round(pk / n, 3),
                     "mean_valleys_unfiltered": round(vl / n, 3),
                     "frac_alternation_abs_le_1": round(alt / n, 3),
                     "frac_identity_peaks_eq_valleys_plus_1": round(ident / n, 3),
                     "n_ends_descending": n_desc,
                     "frac_identity_given_ends_descending":
                         (None if n_desc == 0 else round(ident_desc / float(n_desc), 3)),
                     "mean_peaks_prom0.5": round(pk_f / n, 3),
                     "mean_valleys_prom0.5": round(vl_f / n, 3)})
    out["peak_count_table"] = dual
    out["reference_1_plus_reopen_tau5"] = 1.0 + table[2]["mean_reopen"]

    # (iv) 代表性轨迹（用于报告，取第一条）
    out["example_path_head"] = paths[0][:40]
    out["example_path_max"] = max(paths[0])
    out["example_path_final"] = paths[0][-1]
    return out


# ------------------------------------------------------------------ main


def main():
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    res = {"K": K, "seed": SEED, "config": CONFIG}

    res["A1_fallback"] = [arm_A1_fallback(CONFIG, s) for s in S_GRID]
    res["A2_scale"] = [arm_A2_scale(CONFIG, s) for s in S_GRID]
    res["B_multipeak"] = [arm_B_multipeak(CONFIG, s=s) for s in [0.0, 0.20, 0.40]]
    res["C_mprocess"] = arm_C_mprocess()

    with open("r28_results.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)

    print("=" * 78)
    print("R28  微观回落的本体性 / 峰的尺度依赖 / 多峰 / 局部-全局收敛")
    print("=" * 78)
    print("config: {}  theta_U={}  theta_E={}".format(
        CONFIG["name"], CONFIG["theta_U"], CONFIG["theta_E"]))

    print("\n--- A1 微观回落率：实测 vs 解析 Phi((ln th - ln r_bar + s^2/2)/s) ---")
    for a in res["A1_fallback"]:
        print("  s={:<5} k*={:<7} max|dev|={:<7} mean|dev|={:<7} MC se={:<7} "
              "within_3se={}".format(a["s"], a["kstar_pred"], a["max_abs_dev"],
                                     a["mean_abs_dev"], a["mc_se"], a["within_3se"]))
        print("        峰处回落率 实测={} 解析={} 预测 Phi(s/2)={}".format(
            a["fallback_at_kstar"], a["analytic_at_kstar"],
            a["predicted_Phi_s_over_2"]))

    print("\n--- A2 峰数 / 单峰率 vs (w, h)；噪声在过程层 ---")
    for a in res["A2_scale"]:
        print("  s={}  k*_pred={}  sd(dU)/step={}  单调(w)={} 单调(h)={}".format(
            a["s"], a["kstar_pred"], a["sd_dU_per_step"],
            a["monotone_in_w"], a["monotone_in_h"]))
        print("     {:>4} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}".format(
            "w", "uni_h0", "uni_h.1", "pk_h0", "pk_h.1", "valley.1", "peak_pos"))
        for row in a["rows"]:
            print("     {:>4} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}".format(
                row["w"], row["uni_h_0.0"], row["uni_h_0.1"],
                row["peaks_h_0.0"], row["peaks_h_0.1"],
                row["valleys_h_0.1"], row["peak_mean"]))

    print("\n--- B 多峰：峰数 = 1 + 由减转增的相变次数 ---")
    for b in res["B_multipeak"]:
        print("  s={:<5} det_peaks={} @ {}  回跳数={}  恒等式成立={}".format(
            b.get("noisy_s", 0.0), b["deterministic_peak_count"],
            b["deterministic_peak_idx"], b["reopen_count_sign_flips"],
            b["identity_peaks_eq_1_plus_reopen"]))
        if "noisy_peak_count_mean_w21_h005" in b:
            print("        噪声下峰数均值(w=21,h=0.05)={}".format(
                b["noisy_peak_count_mean_w21_h005"]))

    print("\n--- C m 过程（生灭，不用 A_D）---")
    c = res["C_mprocess"]
    print("  扩张相内生微观回落率 = {}".format(c["micro_reversal_rate_expansion_phase"]))
    print("  局部收敛 / 全局收敛（扫 tau_c）:")
    print("     {:>6} {:>14} {:>14} {:>12} {:>10}".format(
        "tau_c", "有局部收敛%", "全局收敛%", "平均收敛段", "平均重开"))
    for t in c["convergence_table"]:
        print("     {:>6} {:>14} {:>14} {:>12} {:>10}".format(
            t["tau_c"], t["frac_has_local_convergence"],
            t["frac_global_converged"], t["mean_conv_runs"], t["mean_reopen"]))
    print("  峰 / 回跳 / 重开三者分离（1+重开@tau_c=5 参考={}）:".format(
        c["reference_1_plus_reopen_tau5"]))
    for d in c["peak_count_table"]:
        print("     w={:<3} 未过滤 峰={:<7} 回跳={:<7} |差|<=1 比例={:<7} "
              "止于下降[n={:<4}] 条件恒等式={:<7} | 过滤后 峰={:<6} 回跳={}".format(
                  d["w"], d["mean_peaks_unfiltered"], d["mean_valleys_unfiltered"],
                  d["frac_alternation_abs_le_1"], d["n_ends_descending"],
                  d["frac_identity_given_ends_descending"],
                  d["mean_peaks_prom0.5"], d["mean_valleys_prom0.5"]))

    print("\n结果已写入 r28_results.json")


if __name__ == "__main__":
    main()
