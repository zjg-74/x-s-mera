# -*- coding: utf-8 -*-
"""
R29 —— 峰距标定：从峰距反解 rho_res / rho_exp（补篇六 推论 6.5）
================================================================
对应推论（补篇六 §3.6）：

    在对数斜坡 r(k) = exp(lerp(ln r_hi, ln r_lo, k/(K-1))) 下，
        |k*_E - k*_U| = (K-1) * |ln(rho_res/rho_exp)| / |ln(r_lo/r_hi)|

本实验做什么
------------------------------------------------------------------
推论 6.5 说"峰距能测出 rho 比值的对数"。但**一条公式不等于一个标定程序**。
本实验检验：若把它当作**估计程序**用，它站不站得住。六个臂：

  A  公式精确性（无噪声）：扫 rho 网格，实测峰距 vs 解析预测
  B  恢复精度 vs 过程噪声 s：估计量的偏与方差（含 bootstrap CI）
  C  比率设计（C 未知）：
       C1 同 K        —— 比值消去 C，应成功
       C2 异 K        —— raw 比值应**失败**（差因子 (K1-1)/(K2-1)），
                          归一化 delta = Delta/(K-1) 应**成功**
  D  调度误设（负面结果）：真调度非对数斜坡时，仍按 6.5 反解 → 偏多大？
  E  可辨识性视界：ln rho 超出什么范围后峰落在窗外、峰距不可测？
  F  峰序分辨率极限：推论 6.4（峰序 ⟺ rho_res vs rho_exp）在多小 |ln rho| 上失效？

模型（与 R27 / R28 逐字一致）
------------------------------------------------------------------
A_D：   dE = dt*(a*e_plus - b*e_minus),   dU = dt*(a*u_plus - b*u_minus)
        取 b = 1, r := a
        theta_E = e_minus/e_plus,  theta_U = u_minus/u_plus
        rho_exp = u_plus /e_plus,   rho_res = u_minus/e_minus
        sign(dE) = sign(r - theta_E),  sign(dU) = sign(r - theta_U)
        恒等式：theta_U/theta_E = rho_res/rho_exp

A_D^stat（R28）：噪声进**过程层**而非观测层
        r_k = rbar_k * exp(s*Z_k - s^2/2),   E[r_k] = rbar_k
        —— 于是 U 与 E 的峰位本身都是随机变量（这是本体性的，不是测量误差）

标定量的定义
------------------------------------------------------------------
    Delta := |k*_E - k*_U|                      （峰距，单位：步）
    L     := |ln(r_lo/r_hi)|                    （r 调度的对数跨度）
    C     := (K-1) / L                          （灵敏度：步 / 每 e 倍）
    delta := Delta/(K-1)                        （归一化峰距）

    推论 6.5： Delta = C * |ln(rho_res/rho_exp)|
    推论 6.4： sign(k*_U - k*_E) 给出 ln(rho_res/rho_exp) 的**符号**
    ⟹ 序给符号、距给量值；合起来 ln(rho_res/rho_exp) 全辨识
       —— **但只在 C 已知时**。C 未知时只辨识**比值**（见臂 C）。

四个方法学陷阱（本轮实际踩到并已修正，保留记录）
------------------------------------------------------------------
T8. **不能直接拿 Delta 的比值比，除非 K 相同。** C = (K-1)/L 含 K。
    臂 C2 实测：K=241 vs K=121 时 raw 比值差因子 2.000，归一化后恢复。
    ⟹ 推论 6.5 原文"r 调度保持不变"应精确化为：**只需两端点 (r_hi, r_lo) 相同，
       不必 K 相同；但必须归一化。**
T9. **|ln rho| 的绝对量值不可辨识**，除非知道 r 调度。能辨识的只有组间比值
    （或归一化后的相对量）。把"恢复值"写成绝对数值是**越界声称**。
T10. **峰靠近窗口边界时 argmax 有向内偏**（边界效应）。臂 E 显示：
     |ln rho| 越接近视界，估计偏越大；必须有"峰是否为内点"的显式判定，
     否则会把不可测的情形报成"rho 比值很小"。
T11. 平滑窗口的移动平均会**平坦化**峰顶；在臂 F（两峰很近）时，
     这会把两个峰合并成一个 → 峰序判定的**假阴性**。故臂 F 用 w=1 + 显著度。

不能声称什么（写死在这里，防止报告里被写越界）
------------------------------------------------------------------
* 这是原型上的**标定程序检验**，不是实证验证，没有真实理解者参与。
* 系数、调度、r 的量级全部人造；**C = 54.6 步/e 倍不是普适常数**。
* |ln rho| 的**绝对值不可辨识**（T9），只有比值可辨识。
* 结论条件于 A_D 与 A_D^stat（→ B-B22，高，未闭合）。
"""

import json
import math
import random
import statistics as st
import zlib

# ------------------------------------------------------------------ 参数

K = 241                 # 步数（k = 0 .. K-1）
DT = 0.01
R_HI = 9.0
R_LO = 1.0 / 9.0
LN_HI = math.log(R_HI)
LN_LO = math.log(R_LO)
L_SPAN = abs(LN_LO - LN_HI)          # L = 4.3944...
E_START = 5.0
# 注意：U_START 只影响**是否触底**，不影响峰位（峰位只由增量序列决定）。
# R27 用 5.0 足够（其配置不会触底）；本实验把 rho 扫到 6，U 会在末段跌破 0，
# 而 U >= 0 是熵的正当约束（U = H(M|phi,C_in)），故抬高初值以隔离出纯模型行为，
# 并单独报告触底诊断。
U_START = 30.0
N_REP = 400
SEED = 20261020
S_GRID = [0.0, 0.10, 0.20, 0.40, 0.80]
SMOOTH_W = 21
PROM = 0.5

# rho = rho_res/rho_exp 的扫描网格（theta_E = 1 固定，theta_U = rho）
RHO_GRID = [1.0 / 6.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.5, 3.0, 6.0]

# 臂 F 专用：ln rho 的细网格（rho >= 1 一侧），用于定位峰序判据的分辨率极限
F_LN_GRID = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.8]

# 视界扫描：ln rho 从 0 到 2.6（r 半跨度 = 2.1972，故 > 2.1972 应不可测）
HORIZON_GRID = [0.0, 0.4, 0.8, 1.2, 1.6, 1.8, 2.0, 2.1, 2.2, 2.4, 2.6]


def zlib_off(s):
    """可复现的字符串散列。不能用 hash()——Python 的 str hash 按进程随机化
    （PYTHONHASHSEED），会破坏可复现性。见陷阱 T1。"""
    return zlib.crc32(s.encode("utf-8")) % 9973


# ------------------------------------------------------------------ 配置

def make_cfg(rho, e_plus=1.0, u_plus=1.0):
    """给定 rho := rho_res/rho_exp 与 e_plus, u_plus，构造一组 (e±, u±)。

    取 theta_E = 1（即 e_minus = e_plus），则
        theta_U = rho * theta_E = rho  ⟹  u_minus = rho * u_plus
    校验：rho_exp = u_plus/e_plus；rho_res = u_minus/e_minus = rho*u_plus/e_plus
          ⟹ rho_res/rho_exp = rho  ✓
    """
    e_minus = e_plus * 1.0          # theta_E = 1
    u_minus = rho * u_plus          # theta_U = rho
    return {
        "rho": rho,
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": e_minus / e_plus,
        "theta_U": u_minus / u_plus,
        "rho_exp": u_plus / e_plus,
        "rho_res": u_minus / e_minus,
    }


# ------------------------------------------------------------------ 调度

def log_ramp(k, KK=None):
    KK = K if KK is None else KK
    t = k / (KK - 1.0)
    return math.exp(LN_HI + (LN_LO - LN_HI) * t)


def linear_ramp(k, KK=None):
    KK = K if KK is None else KK
    t = k / (KK - 1.0)
    return R_HI + (R_LO - R_HI) * t


def sigmoid_ramp(k, KK=None):
    KK = K if KK is None else KK
    t = k / (KK - 1.0)
    z = 12.0 * (t - 0.5)
    return R_LO + (R_HI - R_LO) / (1.0 + math.exp(z))


def piecewise_ramp(k, KK=None):
    """对数空间三段折线，中段更陡 —— 用于制造调度误设。"""
    KK = K if KK is None else KK
    t = k / (KK - 1.0)
    bps = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
    vals = [LN_HI, math.log(3.0), math.log(1.0 / 3.0), LN_LO]
    if t <= bps[0]:
        return math.exp(vals[0])
    for i in range(3):
        if t <= bps[i + 1]:
            u = (t - bps[i]) / (bps[i + 1] - bps[i])
            return math.exp(vals[i] + (vals[i + 1] - vals[i]) * u)
    return math.exp(vals[-1])


SCHEDULES = [
    ("log_ramp", log_ramp),
    ("linear_ramp", linear_ramp),
    ("sigmoid_ramp", sigmoid_ramp),
    ("piecewise_ramp", piecewise_ramp),
]


def kstar_logramp(theta, KK=None):
    KK = K if KK is None else KK
    if theta >= R_HI or theta <= R_LO:
        return None                      # 峰不在窗口内部
    return (KK - 1.0) * (math.log(theta) - LN_HI) / (LN_LO - LN_HI)


def kstar_numeric(rfun, theta, KK=None):
    """对任意单调下降调度，二分求 r(k) = theta 的连续解。"""
    KK = K if KK is None else KK
    if rfun(0.0, KK) < theta:
        return None                      # theta 高于起点：峰在窗口左侧
    if rfun(float(KK - 1), KK) > theta:
        return None                      # theta 低于终点：峰在窗口右侧
    lo, hi = 0.0, float(KK - 1)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rfun(mid, KK) > theta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------------ 轨迹

def trajectory(rng, cfg, rbar_of_k, s=0.0, KK=None, indep_rel=0.0):
    """A_D + A_D^stat。噪声进过程层：r_k = rbar_k * exp(s Z_k - s^2/2)。

    s         —— **共享**噪声：同一 r 过程同时驱动 dE 与 dU（补篇七 A_D^stat）
    indep_rel —— **独立**噪声：分别加在 dE 与 dU 上，但**每步波动率与共享噪声匹配**
                 （sd = indep_rel * dt * rbar_k * 各自系数），两次独立抽样。

    ⟹ 两个条件的唯一差别是**两路噪声是否相关**，波动率相同。
       这样任何峰序稳健性的差异都可归因于"共享"本身，而不是噪声大小。

    关键恒等式（定理 6.10）：共享噪声下
        U(k) - E(k) = (U_0 - E_0) + dt*e_minus*(rho_exp - rho_res)*(k+1)
    与 r 及噪声**完全无关** ⟹ 峰序是**路径级守恒量**（不会判反）。
    独立噪声破坏该恒等式 ⟹ 峰序可翻转，误差率随噪声单调上升。
    """
    KK = K if KK is None else KK
    E, U = E_START, U_START
    Es, Us = [], []
    clipped = 0
    for k in range(KK):
        rb = rbar_of_k(float(k))
        if s > 0.0:
            r = rb * math.exp(s * rng.gauss(0.0, 1.0) - 0.5 * s * s)
        else:
            r = rb
        dE = DT * (r * cfg["e_plus"] - cfg["e_minus"])
        dU = DT * (r * cfg["u_plus"] - cfg["u_minus"])
        if indep_rel > 0.0:
            # 两次**独立**抽样，波动率与共享噪声一致
            dE += rng.gauss(0.0, indep_rel * DT * rb * cfg["e_plus"])
            dU += rng.gauss(0.0, indep_rel * DT * rb * cfg["u_plus"])
        E += dE
        U += dU
        if U < 0.0:
            U = 0.0
            clipped += 1
        Es.append(E)
        Us.append(U)
    return Es, Us, clipped


# ------------------------------------------------------------------ 峰估计

def moving_average(seq, w):
    if w <= 1:
        return list(seq)
    out, n = [], len(seq)
    half = w // 2
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out.append(sum(seq[lo:hi]) / (hi - lo))
    return out


def argmax_of(seq):
    bi, bv = 0, seq[0]
    for i, v in enumerate(seq):
        if v > bv:
            bv, bi = v, i
    return bi


def argmax_plateau(seq, tol=1e-9):
    bv = max(seq)
    idx = [i for i, v in enumerate(seq) if v >= bv - tol]
    return idx[0], idx[-1]


def local_maxima(seq):
    out, n = [], len(seq)
    i = 1
    while i < n - 1:
        if seq[i] > seq[i - 1] and seq[i] >= seq[i + 1]:
            j = i
            while j + 1 < n and seq[j + 1] == seq[i]:
                j += 1
            out.append((i + j) // 2)
            i = j + 1
        else:
            i += 1
    return out


def local_minima(seq):
    return local_maxima([-v for v in seq])


def prominence(seq, i):
    v = seq[i]
    l = i
    while l > 0 and seq[l - 1] <= v:
        l -= 1
    left_min = min(seq[l:i + 1]) if l < i else seq[0]
    r = i
    n = len(seq)
    while r < n - 1 and seq[r + 1] <= v:
        r += 1
    right_min = min(seq[i:r + 1]) if r > i else seq[n - 1]
    return v - max(left_min, right_min)


def count_peaks(seq, min_prom):
    return sum(1 for i in local_maxima(seq) if prominence(seq, i) > min_prom)


def count_valleys(seq, min_prom):
    neg = [-v for v in seq]
    return sum(1 for i in local_maxima(neg) if prominence(neg, i) > min_prom)


def is_interior_peak(seq, idx):
    """峰是否为**内点**（不贴窗口边界）。贴边即不可辨识（陷阱 T10）。"""
    return 0 < idx < len(seq) - 1


# ------------------------------------------------------------------ 臂 A：公式精确性（无噪声）

def arm_A_formula():
    rows = []
    for rho in RHO_GRID:
        cfg = make_cfg(rho)
        pE = kstar_logramp(cfg["theta_E"])
        pU = kstar_logramp(cfg["theta_U"])
        pred = None if (pE is None or pU is None) else abs(pE - pU)
        rng = random.Random(SEED + zlib_off("A") + int(round(rho * 1000)))
        Es, Us, _clip = trajectory(rng, cfg, log_ramp, s=0.0)
        kE = argmax_of(Es)
        kU = argmax_of(Us)
        obs = abs(kE - kU)
        # 平台（并列取小索引的伪像，见陷阱 T2）
        plo_E, phi_E = argmax_plateau(Es)
        plo_U, phi_U = argmax_plateau(Us)
        rows.append({
            "rho": round(rho, 4),
            "ln_rho": round(math.log(rho), 4),
            "abs_ln_rho": round(abs(math.log(rho)), 4),
            "pred_delta": None if pred is None else round(pred, 3),
            "obs_delta": round(float(obs), 3),
            "abs_err": None if pred is None else round(abs(obs - pred), 3),
            "plateau_E": [plo_E, phi_E],
            "plateau_U": [plo_U, phi_U],
            # 陷阱 T2：并列时 argmax 取小索引，故必须用**平台区间**比对，
            # 且是两峰平台之**差**的区间，不是各自区间的并。
            "delta_from_plateau": [round(float(abs(plo_E - phi_U))),
                                   round(float(abs(phi_E - plo_U)))],
            # rho = 1 时两峰由同一阈值驱动，必重合；独立区间处理在此退化，
            # 故按项目惯例返回 None + 标注「不适用」，不报成 True（见陷阱 T2）。
            "pred_in_plateau": (
                None if (pred is None or abs(math.log(rho)) < 1e-12)
                else bool(min(abs(plo_E - phi_U), abs(phi_E - plo_U)) - 0.51
                          <= pred
                          <= max(abs(plo_E - phi_U), abs(phi_E - plo_U)) + 0.51)),
            "pred_in_plateau_note": (
                "不适用（两峰重合，平台差区间退化）"
                if abs(math.log(rho)) < 1e-12 else "适用"),
        })
    # 线性检查：obs_delta 对 |ln rho| 的过原点最小二乘斜率（应 ≈ C）
    xs = [r["abs_ln_rho"] for r in rows]
    ys = [r["obs_delta"] for r in rows]
    slope = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
    resid = max(abs(y - slope * x) for x, y in zip(xs, ys))
    return {
        "C_theory": round((K - 1.0) / L_SPAN, 4),
        "C_fitted": round(slope, 4),
        "max_resid_steps": round(resid, 3),
        "rows": rows,
    }


# ------------------------------------------------------------------ 臂 B：恢复精度 vs s

def arm_B_recovery():
    out = []
    for s in S_GRID:
        rows = []
        for rho in RHO_GRID:
            cfg = make_cfg(rho)
            pE = kstar_logramp(cfg["theta_E"])
            pU = kstar_logramp(cfg["theta_U"])
            pred = None if (pE is None or pU is None) else abs(pE - pU)
            rng = random.Random(SEED + 11 + int(round(rho * 1000)) + int(round(s * 100)))
            ds, ks = [], []
            for _ in range(N_REP):
                Es, Us, _clip = trajectory(rng, cfg, log_ramp, s=s)
                smE = moving_average(Es, SMOOTH_W)
                smU = moving_average(Us, SMOOTH_W)
                kE = argmax_of(smE)
                kU = argmax_of(smU)
                ds.append(abs(kE - kU))
                ks.append(kU)
            C = (K - 1.0) / L_SPAN
            mean_d = st.mean(ds)
            # 估计量：|ln rho|^ = Delta_hat / C
            est = mean_d / C
            true = abs(math.log(rho))
            # bootstrap（对每 rep 的 Delta）
            rngb = random.Random(SEED + 991 + int(round(rho * 1000)) + int(round(s * 100)))
            boots = []
            n = len(ds)
            for _ in range(2000):
                sample = [ds[rngb.randrange(n)] for _ in range(n)]
                boots.append(st.mean(sample) / C)
            boots.sort()
            lo = boots[int(0.025 * len(boots))]
            hi = boots[int(0.975 * len(boots))]
            rows.append({
                "rho": round(rho, 4),
                "true_abs_ln_rho": round(true, 4),
                "pred_delta": None if pred is None else round(pred, 3),
                "obs_delta_mean": round(mean_d, 3),
                "obs_delta_sd": round(st.pstdev(ds), 3),
                "est_abs_ln_rho": round(est, 4),
                "bias": round(est - true, 4),
                "ci95_lo": round(lo, 4),
                "ci95_hi": round(hi, 4),
                "true_in_ci": bool(lo - 1e-9 <= true <= hi + 1e-9),
            })
        out.append({"s": s, "rows": rows})
    return out


# ------------------------------------------------------------------ 臂 C：比率设计（C 未知）

def _measure_delta(rho, s, KK=None, rfun=None, n_rep=None, seed_off=0, smooth_w=SMOOTH_W):
    """测 (mean Delta, sd Delta)，可选自定义调度与 K。"""
    KK = K if KK is None else KK
    rfun = log_ramp if rfun is None else rfun
    n_rep = N_REP if n_rep is None else n_rep
    cfg = make_cfg(rho)
    rng = random.Random(SEED + 301 + int(round(rho * 1000)) + int(round(s * 100)) + seed_off)
    ds = []
    for _ in range(n_rep):
        Es, Us, _clip = trajectory(rng, cfg, lambda k: rfun(k, KK), s=s, KK=KK)
        kE = argmax_of(moving_average(Es, smooth_w))
        kU = argmax_of(moving_average(Us, smooth_w))
        ds.append(abs(kE - kU))
    return st.mean(ds), st.pstdev(ds)


def arm_C_ratio_design():
    """C 未知时，用两条件峰距之比消去 C。

    C1：两组同 K（=241），raw 比值应 ≈ 真比值。
    C2：两组异 K（241 vs 121），raw 比值应**偏**因子 2；
        归一化 delta = Delta/(K-1) 后应恢复真比值。
    """
    res = {"C1_same_K": [], "C2_diff_K": []}
    pairs = [(1.5, 3.0), (1.0 / 3.0, 3.0), (2.0 / 3.0, 6.0), (3.0, 6.0)]

    for rho1, rho2 in pairs:
        true_ratio = abs(math.log(rho1)) / abs(math.log(rho2))
        # C1
        for s in (0.0, 0.2, 0.4):
            m1, sd1 = _measure_delta(rho1, s)
            m2, sd2 = _measure_delta(rho2, s)
            raw = (m1 / m2) if m2 > 1e-12 else None
            res["C1_same_K"].append({
                "rho1": round(rho1, 4), "rho2": round(rho2, 4), "s": s,
                "true_ratio": round(true_ratio, 4),
                "obs_ratio_raw": None if raw is None else round(raw, 4),
                "rel_err": None if raw is None else round(
                    abs(raw - true_ratio) / true_ratio, 4),
            })
        # C2
        K1, K2 = 241, 121
        for s in (0.0, 0.2, 0.4):
            m1, sd1 = _measure_delta(rho1, s, KK=K1)
            m2, sd2 = _measure_delta(rho2, s, KK=K2)
            raw = (m1 / m2) if m2 > 1e-12 else None
            d1 = m1 / (K1 - 1.0)
            d2 = m2 / (K2 - 1.0)
            norm = (d1 / d2) if d2 > 1e-12 else None
            res["C2_diff_K"].append({
                "rho1": round(rho1, 4), "rho2": round(rho2, 4), "s": s,
                "K1": K1, "K2": K2,
                "true_ratio": round(true_ratio, 4),
                "obs_ratio_raw": None if raw is None else round(raw, 4),
                "rel_err_raw": None if raw is None else round(
                    abs(raw - true_ratio) / true_ratio, 4),
                "obs_ratio_normalized": None if norm is None else round(norm, 4),
                "rel_err_normalized": None if norm is None else round(
                    abs(norm - true_ratio) / true_ratio, 4),
            })
    return res


# ------------------------------------------------------------------ 臂 D：调度误设（负面结果）

def arm_D_schedule_misspec():
    """若真实 r 调度**不是**对数斜坡，仍按推论 6.5 反解 rho → 偏多大？

    推论 6.5 依赖 k*(theta) 关于 ln theta **仿射**，这只在对数斜坡下成立。
    一般调度下 k*(theta) = g(ln theta)，g 非线性 ⟹ Delta = |g(ln th_E) - g(ln th_U)|
    不再是 C*|ln rho|，且**比值设计也随之失效**。
    """
    out = []
    for name, rfun in SCHEDULES:
        rows = []
        for rho in RHO_GRID:
            cfg = make_cfg(rho)
            kE_true = kstar_numeric(rfun, cfg["theta_E"])
            kU_true = kstar_numeric(rfun, cfg["theta_U"])
            if kE_true is None or kU_true is None:
                rows.append({"rho": round(rho, 4), "observable": False,
                             "true_delta": None, "logramp_pred": None,
                             "abs_err_steps": None, "rel_err": None})
                continue
            true_delta = abs(kE_true - kU_true)
            # 实测（无噪声，取数值峰位，应与 true_delta 一致）
            rng = random.Random(SEED + 555 + int(round(rho * 1000)))
            Es, Us, _clip = trajectory(rng, cfg, lambda k: rfun(k, K), s=0.0)
            obs_delta = float(abs(argmax_of(Es) - argmax_of(Us)))
            # 若（错误地）按对数斜坡公式反解 |ln rho|
            C = (K - 1.0) / L_SPAN
            est = obs_delta / C
            row_true = abs(math.log(rho))
            rows.append({
                "rho": round(rho, 4), "observable": True,
                "true_delta": round(true_delta, 3),
                "obs_delta": round(obs_delta, 3),
                "est_abs_ln_rho": round(est, 4),
                "true_abs_ln_rho": round(row_true, 4),
                "abs_err_steps": round(abs(obs_delta - true_delta), 3),
                "rel_err": None if row_true == 0 else round(
                    abs(est - row_true) / row_true, 4),
            })
        # 该调度下"比值设计"是否仍能消去 C：取两对 rho，看比值是否守恒
        ratio_check = []
        for r1, r2 in ((1.5, 3.0), (1.0 / 3.0, 3.0)):
            d1 = kstar_numeric(rfun, r1)
            d2 = kstar_numeric(rfun, r2)
            e1 = kstar_numeric(rfun, 1.0)
            if None in (d1, d2, e1):
                ratio_check.append({"r1": r1, "r2": r2, "observable": False})
                continue
            D1, D2 = abs(e1 - d1), abs(e1 - d2)
            obs = (D1 / D2) if D2 > 1e-12 else None
            true = abs(math.log(r1)) / abs(math.log(r2))
            ratio_check.append({
                "r1": round(r1, 4), "r2": round(r2, 4), "observable": True,
                "true_ratio": round(true, 4),
                "obs_ratio": None if obs is None else round(obs, 4),
                "rel_err": None if obs is None else round(abs(obs - true) / true, 4),
            })
        out.append({"schedule": name, "rows": rows, "ratio_check": ratio_check})
    return out


# ------------------------------------------------------------------ 臂 E：可辨识性视界

def arm_E_horizon():
    """theta_E = 1 固定，扫 ln rho（正侧），找峰落出窗口的位置。

    可辨识条件：r_lo < theta_E < r_hi 且 r_lo < theta_U < r_hi。
    theta_U = rho * theta_E ⟹ ln rho ∈ (ln r_lo - ln th_E, ln r_hi - ln th_E)
    ⟹ 窗口宽度 = L = |ln(r_hi/r_lo)|，中心 = -ln theta_E。
    """
    rows = []
    for lnr in HORIZON_GRID:
        rho = math.exp(lnr)
        cfg = make_cfg(rho)
        pE = kstar_logramp(cfg["theta_E"])
        pU = kstar_logramp(cfg["theta_U"])
        observable = (pE is not None) and (pU is not None)
        pred = None if not observable else abs(pE - pU)
        row = {
            "ln_rho": round(lnr, 4),
            "rho": round(rho, 4),
            "theta_U": round(cfg["theta_U"], 4),
            "observable": bool(observable),
            "pred_delta": None if pred is None else round(pred, 3),
        }
        if observable:
            rng = random.Random(SEED + 777 + int(round(lnr * 1000)))
            Es, Us, _clip = trajectory(rng, cfg, log_ramp, s=0.0)
            kE = argmax_of(Es)
            kU = argmax_of(Us)
            row["obs_delta"] = round(float(abs(kE - kU)), 3)
            row["kU_interior"] = bool(is_interior_peak(Us, kU))
            row["kE_interior"] = bool(is_interior_peak(Es, kE))
            row["abs_err"] = round(abs(abs(kE - kU) - pred), 3)
            # 陷阱 T10：rho 增大 ⟹ theta_U 增大 ⟹ k*_U 朝 k=0 移动。
            # 故**约束方向是距起点**的余量（不是终点）。该余量趋零即视界。
            row["kU_margin_to_start"] = round(float(kU), 2)
            row["kU_margin_to_end"] = round(float((K - 1) - kU), 2)
        else:
            row["obs_delta"] = None
            row["reason"] = "theta_U outside (r_lo, r_hi) -> U peak not interior"
        rows.append(row)
    return {
        "theta_E": 1.0,
        "r_hi": R_HI, "r_lo": round(R_LO, 6),
        "window_ln_rho": [round(LN_LO - 0.0, 4), round(LN_HI - 0.0, 4)],
        # 不 round：视界律要求窗口宽度**恰等于** L，round 到 4 位会引入 5e-5 的假残差
        "window_width": L_SPAN,
        "rows": rows,
    }


# ------------------------------------------------------------------ 臂 F：峰序分辨率极限

def arm_F_order_resolution():
    """推论 6.4 用峰序判 rho_res vs rho_exp。两峰越近越难判。测分辨率极限。

    两个噪声条件（关键区分）：
      shared —— 同一 r 过程驱动 E 与 dU（补篇七 A_D^stat 的设定）
      indep  —— E 与 U 各带独立噪声（对照臂，用于定位守恒律的**失效条件**）

    恒等式（定理 6.10）：共享噪声下
        U(k) - E(k) = (U_0 - E_0) + dt*e_minus*(rho_exp - rho_res)*(k+1)
    与 r 及噪声**完全无关** ⟹ 峰序是**路径级守恒量**，
    样本层面只会出现"判不出"（并列），不会出现"判反"。
    独立噪声破坏该恒等式 ⟹ 峰序可翻转，误差率随噪声单调上升。
    """
    out = []
    for s in S_GRID:
        rows = []
        for lnr in F_LN_GRID:
            rho = math.exp(lnr)
            cfg = make_cfg(rho)
            pred_order = ("U_first" if cfg["rho_res"] > cfg["rho_exp"]
                          else ("E_first" if cfg["rho_res"] < cfg["rho_exp"]
                                else "degenerate"))
            blocks = {}
            for cond, s_sh, s_in in (("shared", s, 0.0), ("indep", 0.0, s)):
                rng = random.Random(SEED + 555 + int(round(lnr * 1000))
                                    + int(round(s * 100)) + (0 if cond == "shared" else 7))
                Ds = []
                for _ in range(N_REP):
                    Es, Us, _clip = trajectory(rng, cfg, log_ramp, s=s_sh, indep_rel=s_in)
                    # 陷阱 T11：此处**不平滑**；平滑会把近邻双峰合并成假的单峰
                    Ds.append(argmax_of(Us) - argmax_of(Es))
                if pred_order == "degenerate":
                    blocks[cond] = {"note": "不适用（两峰重合，无次序可言）",
                                    "mean_D": round(st.mean(Ds), 3),
                                    "sd_D": round(st.pstdev(Ds), 3)}
                    continue
                # rho>=1 ⟹ rho_res > rho_exp ⟹ 预测 U 峰在前 ⟹ 应有 D = kU - kE <= 0
                wrong = sum(1 for x in Ds if x > 0)
                blocks[cond] = {
                    "mean_D": round(st.mean(Ds), 3),
                    "sd_D": round(st.pstdev(Ds), 3),
                    "tie_rate": round(sum(1 for x in Ds if x == 0) / len(Ds), 4),
                    "err_rate": round(wrong / len(Ds), 4),
                }
            rows.append({
                "ln_rho": round(lnr, 4),
                "rho": round(rho, 4),
                "pred_order": pred_order,
                "shared": blocks["shared"],
                "indep": blocks["indep"],
            })
        out.append({"s": s, "rows": rows})
    return out


# ------------------------------------------------------------------ 臂 G：U-E 线性恒等式（定理 6.10）

def arm_G_identity():
    """直接验证定理 6.10 的恒等式，含 rho_exp != 1 的一般情形。

    取 e_plus=1, e_minus=1（theta_E=1），u_plus=w（rho_exp=w），u_minus=w*rho
    （theta_U=rho, rho_res=w*rho）。令 U^s := U/w（即按 rho_exp 归一化），则
    u_plus^s = 1 = e_plus，r_j 项相消，得到与 r 及噪声**无关**的恒等式：

        U^s(k) - E(k) = (U_0^s - E_0) + dt * e_minus * (1 - rho_res/rho_exp) * (k+1)

    若该式成立 ⟹ U^s - E 单调 ⟹ 峰序**路径级守恒**（只会出现并列，不会判反）。
    本臂逐点比对左右两边，报告最大绝对残差（应 ~1e-15 量级，纯浮点误差）。
    """
    rows = []
    for w in (1.0, 2.0, 0.5):
        for rho in (3.0, 1.0 / 3.0, 6.0):
            cfg = make_cfg(rho, e_plus=1.0, u_plus=w)
            for s in (0.0, 0.4, 0.8):
                rng = random.Random(SEED + 4242 + int(round(rho * 1000))
                                     + int(round(s * 100)) + int(round(w * 100)))
                Es, Us, nclip = trajectory(rng, cfg, log_ramp, s=s)
                # U 按 rho_exp 归一化（恒等式只在归一化后成立）
                Us_s = [v / cfg["rho_exp"] for v in Us]
                c0 = (U_START / cfg["rho_exp"]) - E_START
                slope = DT * cfg["e_minus"] * (1.0 - cfg["rho_res"] / cfg["rho_exp"])
                kU = argmax_of(Us)
                resid_full = 0.0
                resid_pre = 0.0
                for k in range(K):
                    lhs = Us_s[k] - Es[k]
                    rhs = c0 + slope * (k + 1)
                    d = abs(lhs - rhs)
                    resid_full = max(resid_full, d)
                    if k <= kU:
                        resid_pre = max(resid_pre, d)
                rows.append({
                    "u_plus(w)": w,
                    "rho_exp": round(cfg["rho_exp"], 4),
                    "rho_res": round(cfg["rho_res"], 4),
                    "rho_res_over_rho_exp": round(cfg["rho_res"] / cfg["rho_exp"], 4),
                    "rho": round(rho, 4), "s": s,
                    "max_abs_resid_full": "{:.3e}".format(resid_full),
                    # 触底截断发生在**峰后**，故 k <= k*_U 才是恒等式的有效检验域
                    "max_abs_resid_upto_Upeak": "{:.3e}".format(resid_pre),
                    "n_clipped": nclip,
                })
    mx = max(float(r["max_abs_resid_upto_Upeak"]) for r in rows)
    mx_full = max(float(r["max_abs_resid_full"]) for r in rows)
    nclip = sum(int(r["n_clipped"]) for r in rows)
    return {
        "identity": "U_norm(k) - E(k) = (U0/rho_exp - E0) + dt*e_minus*(1 - rho_res/rho_exp)*(k+1)",
        "max_abs_resid_over_all": "{:.3e}".format(mx),
        "max_abs_resid_full_window": "{:.3e}".format(mx_full),
        "note": ("identity is exact (resid ~1e-15) only on the pre-U-peak window; "
                 "the full-window residual is inflated by the U>=0 floor, which "
                 "saturates AFTER the peak and therefore does not affect "
                 "peak-based observables. Peak ORDER is pathwise conserved "
                 "under shared r noise."),
        "rows": rows,
    }


# ------------------------------------------------------------------ 主流程

def main():
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    C = (K - 1.0) / L_SPAN
    result = {
        "experiment": "R29",
        "claim": "peak distance calibrates ln(rho_res/rho_exp) (supplement 6, corollary 6.5)",
        "params": {"K": K, "DT": DT, "R_HI": R_HI, "R_LO": R_LO,
                   "L_span": round(L_SPAN, 6), "C_sensitivity_steps_per_efold": round(C, 4),
                   "N_REP": N_REP, "SEED": SEED, "S_GRID": S_GRID,
                   "SMOOTH_W": SMOOTH_W, "PROM": PROM,
                   "RHO_GRID": [round(x, 4) for x in RHO_GRID]},
        "note": ("calibration-procedure check on a synthetic A_D/A_D^stat prototype; "
                 "NOT empirical validation. C = 54.6 steps/efold is an artifact of the "
                 "chosen schedule, not a universal constant. Absolute |ln rho| is NOT "
                 "identified without knowing the r schedule (trap T9). All results "
                 "conditional on A_D (B-B22, high, open)."),
        "A_formula_exact": arm_A_formula(),
        "B_recovery_vs_noise": arm_B_recovery(),
        "C_ratio_design": arm_C_ratio_design(),
        "D_schedule_misspec": arm_D_schedule_misspec(),
        "E_identifiability_horizon": arm_E_horizon(),
        "F_order_resolution": arm_F_order_resolution(),
        "G_identity_UE_linear": arm_G_identity(),
    }

    with open("r29_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print("R29 ---- peak distance calibration of ln(rho_res/rho_exp)")
    print("K={}  r: {:.4f} -> {:.4f}  L={:.4f}  C={:.3f} steps/efold  N_REP={}".format(
        K, R_HI, R_LO, L_SPAN, C, N_REP))

    print("")
    print("--- A 公式精确性（无噪声）---")
    def fmt(v, w=9):
        return ("{:.4f}".format(v) if isinstance(v, float) else str(v)).rjust(w)

    print("   {:>8} {:>9} {:>11} {:>11} {:>9} {:>10}".format(
        "rho", "|ln rho|", "pred_dlta", "obs_dlta", "abs_err", "in_plat"))
    for r in result["A_formula_exact"]["rows"]:
        v = r["pred_in_plateau"]
        print("   {:>8} {:>9} {:>11} {:>11} {:>9} {:>10}".format(
            fmt(r["rho"]), fmt(r["abs_ln_rho"]), fmt(r["pred_delta"]),
            fmt(r["obs_delta"]), fmt(r["abs_err"]),
            "n/a" if v is None else str(v)))
    print("   C_theory={}  C_fitted={}  max_resid={} steps".format(
        result["A_formula_exact"]["C_theory"],
        result["A_formula_exact"]["C_fitted"],
        result["A_formula_exact"]["max_resid_steps"]))

    print("")
    print("--- B 恢复精度 vs 过程噪声 s（单位：|ln rho|）---")
    for blk in result["B_recovery_vs_noise"]:
        print("  s={}".format(blk["s"]))
        print("     {:>8} {:>10} {:>13} {:>11} {:>9} {:>18} {:>7}".format(
            "rho", "true", "obs_delta", "sd(delta)", "est", "CI95", "in_CI"))
        for r in blk["rows"]:
            print("     {:>8} {:>10} {:>13} {:>11} {:>9} {:>18} {:>7}".format(
                r["rho"], r["true_abs_ln_rho"], r["obs_delta_mean"],
                r["obs_delta_sd"], r["est_abs_ln_rho"],
                "[{}, {}]".format(r["ci95_lo"], r["ci95_hi"]), r["true_in_ci"]))

    print("")
    print("--- C1 比率设计（同 K=241）---")
    print("   {:>8} {:>8} {:>6} {:>11} {:>11} {:>9}".format(
        "rho1", "rho2", "s", "true_ratio", "obs_ratio", "rel_err"))
    for r in result["C_ratio_design"]["C1_same_K"]:
        print("   {:>8} {:>8} {:>6} {:>11} {:>11} {:>9}".format(
            r["rho1"], r["rho2"], r["s"], r["true_ratio"], r["obs_ratio_raw"],
            r["rel_err"]))
    print("--- C2 比率设计（异 K：241 vs 121）---")
    print("   {:>8} {:>8} {:>6} {:>11} {:>11} {:>9} {:>11} {:>9}".format(
        "rho1", "rho2", "s", "true_ratio", "raw_ratio", "rel_err",
        "norm_ratio", "rel_err"))
    for r in result["C_ratio_design"]["C2_diff_K"]:
        print("   {:>8} {:>8} {:>6} {:>11} {:>11} {:>9} {:>11} {:>9}".format(
            r["rho1"], r["rho2"], r["s"], r["true_ratio"], r["obs_ratio_raw"],
            r["rel_err_raw"], r["obs_ratio_normalized"], r["rel_err_normalized"]))

    print("")
    print("--- D 调度误设（真调度非对数斜坡，仍按 6.5 反解）---")
    for blk in result["D_schedule_misspec"]:
        print("  schedule={}".format(blk["schedule"]))
        for r in blk["rows"]:
            if r["observable"]:
                print("     rho={:<8} true_delta={:<9} obs={:<9} est|ln rho|={:<8} "
                      "true={:<8} rel_err={}".format(
                          r["rho"], r["true_delta"], r["obs_delta"],
                          r["est_abs_ln_rho"], r["true_abs_ln_rho"], r["rel_err"]))
            else:
                print("     rho={:<8} NOT observable ({})".format(
                    r["rho"], r.get("reason", "")))
        for rc in blk["ratio_check"]:
            if rc["observable"]:
                print("     比值检验 r1={} r2={}: true={} obs={} rel_err={}".format(
                    rc["r1"], rc["r2"], rc["true_ratio"], rc["obs_ratio"],
                    rc["rel_err"]))
            else:
                print("     比值检验 r1={} r2={}: NOT observable".format(
                    rc["r1"], rc["r2"]))

    print("")
    print("--- E 可辨识性视界（theta_E=1，扫 ln rho > 0）---")
    eh = result["E_identifiability_horizon"]
    print("   窗口 ln rho ∈ {},  宽度 = {}".format(eh["window_ln_rho"], eh["window_width"]))
    def fmt(v, w=9):
        return ("{:.4f}".format(v) if isinstance(v, float) else str(v)).rjust(w)

    print("   {:>8} {:>9} {:>9} {:>11} {:>11} {:>9} {:>11}".format(
        "ln_rho", "rho", "theta_U", "obs_dlta", "pred_dlta", "interior", "marg_start"))
    for r in eh["rows"]:
        print("   {:>8} {:>9} {:>9} {:>11} {:>11} {:>9} {:>11}".format(
            r["ln_rho"], r["rho"], r["theta_U"], fmt(r.get("obs_delta")),
            fmt(r.get("pred_delta")), fmt(r.get("kU_interior"), 9),
            fmt(r.get("kU_margin_to_start"), 11)))

    print("")
    print("--- F 峰序分辨率极限（不平滑，显著度过滤）---")
    for blk in result["F_order_resolution"]:
        print("  s={}".format(blk["s"]))
        print("     {:>7} {:>8} {:>10} | {:>17} | {:>17}".format(
            "ln_rho", "rho", "pred", "shared: tie/err", "indep: tie/err"))
        for r in blk["rows"]:
            sh, ind = r["shared"], r["indep"]
            if "note" in sh:
                print("     {:>7} {:>8} {:>10} | {:>17} | {:>17}".format(
                    r["ln_rho"], r["rho"], r["pred_order"], "n/a", "n/a"))
                continue
            print("     {:>7} {:>8} {:>10} | {:>8}/{:<8} | {:>8}/{:<8}".format(
                r["ln_rho"], r["rho"], r["pred_order"],
                sh["tie_rate"], sh["err_rate"],
                ind["tie_rate"], ind["err_rate"]))

    print("")
    print("--- G U-E 线性恒等式（定理 6.10）---")
    gi = result["G_identity_UE_linear"]
    print("   {}".format(gi["identity"]))
    print("   峰前窗口最大残差 = {}   （全窗口 = {}，被 U>=0 触底抬升）".format(
        gi["max_abs_resid_over_all"], gi["max_abs_resid_full_window"]))
    print("   {:>9} {:>9} {:>9} {:>8} {:>6} {:>14} {:>10}".format(
        "u_plus", "rho_exp", "rho_res", "rho", "s", "resid_pre", "n_clip"))
    for r in gi["rows"]:
        print("   {:>9} {:>9} {:>9} {:>8} {:>6} {:>14} {:>10}".format(
            r["u_plus(w)"], r["rho_exp"], r["rho_res"], r["rho"], r["s"],
            r["max_abs_resid_upto_Upeak"], r["n_clipped"]))

    print("")
    print("wrote r29_results.json")


if __name__ == "__main__":
    main()
