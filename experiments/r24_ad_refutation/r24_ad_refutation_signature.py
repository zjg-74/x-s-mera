# -*- coding: utf-8 -*-
"""
R24 —— A_D 的证否签名：判别力分解与操作规程
================================================================
对应推论（补篇五 §6，B-B22 的证否签名）：

    R24：同一理解者在不同任务上不应同时出现 (+,-) 与 (-,+)。
         若出现，则 A_D 被证否。

本实验做什么
------------------------------------------------------------------
R24 是 A_D **唯一**的证否途径，但登记时只写了这一句。
本实验检验这句话的**判别力**：它究竟证否了什么？

结论先行（详见补篇九）：

    R24 原登记的跨任务并存签名，证否的是**合取**
    (A_D ∧ rho-任务不变性)，**不是 A_D 单独**。
    因为「rho 依赖任务」会产生**完全同样的**并存现象。
    分离二者需要**第三可观测量**：任务内线性失配残差
    ——它在 A_D 下**逐样本恒为 0**（定理 6.11）。

九个臂（★= 核心）

  A  预言复核：A_D + rho 任务不变 ⟹ 跨任务签名一致、归一化峰距不变
  B  混淆臂（★核心）：A_D **逐任务成立** + rho 任务依赖
        B1  跨任务签名判据（R24 原登记） ⟹ 报"证否"（**假阳**）
        B2  任务内线性残差判据（定理 6.11） ⟹ 报"A_D 成立"（**正确**）
        B3  任务内签名唯一性 ⟹ 也报"A_D 成立"（**正确**）
  C  分离力：注入 -eps*r^2（真违反 A_D），比较三个判据的灵敏度
  D  共享 vs 独立噪声：定理 6.11 的守恒依赖共享驱动
  E  可探测视界：eps 细网格，三个判据的探测率
  E2 阈值标度（★）：探测翻转点 eps_c = e_+/(4*theta_E) ∝ 1/theta_E
  F  操作规程端到端：三分类混淆矩阵
  G  四系数辨识（闭合 B-B23 的条件版）：回归而非峰距
  H  观测噪声地板与平均律：残差判据的 T16 修正

模型（与 R27 / R28 / R29 逐字一致）
------------------------------------------------------------------
A_D：    dE = dt*(r*e_plus - e_minus),   dU = dt*(r*u_plus - u_minus)
        theta_E = e_minus/e_plus,  theta_U = u_minus/u_plus
        rho_exp = u_plus/e_plus,    rho_res = u_minus/e_minus
        sign(dE) = sign(r - theta_E),  sign(dU) = sign(r - theta_U)
        恒等式：theta_U/theta_E = rho_res/rho_exp =: rho

A_D^stat（R28）：噪声进**过程层**
        r_k = rbar_k * exp(s Z_k - s^2/2),   E[r_k] = rbar_k

违反 A_D 的注入（与 R19 配置 IV 同构）
------------------------------------------------------------------
        dE <- dE - dt*eps*r^2
    Psi(r) := dE/dt = e_plus*(r - theta_E) - eps*r^2
    —— 由**直线**变为**开口向下的抛物线**：
       * eps = 0                    —— 恰 1 个零点（A_D）
       * 0 < eps < 1/(4 theta_E)，且两根都落在 r 窗口内 —— 2 个零点（可探测）
       * eps >= 1/(4 theta_E)       —— **无实根**，E 全程下降、无峰（另类签名）

    Psi = 0 处 dU/dk = dt*u_plus*(r - theta_U)：两个零点给出**相反**符号
    ⟹ 与推论 5.3（"Psi=0 处 U 走向只由 theta_E - theta_U 决定"）自相矛盾。

本轮的方法学陷阱（T12–T16，实际踩到并已修正，保留记录）
------------------------------------------------------------------
T12. **跨任务并存 ≠ 违反 A_D。** rho 任务依赖同样产生并存；
     R24 原签名只能证否合取 (A_D ∧ rho-任务不变)，不能证否 A_D 单独。
T13. **"是否报证否"是逐 rep 的随机变量**，必须报率；只跑一次就下结论是过度解读。
T14. **判别"内部极小"必须用显著度过滤**，否则噪声下的假变号淹没信号（承 T5/T11）。
T15. **eps 过大时 Psi 无实根**，此时 E 无峰——探测率关于 eps **非单调**。
T16. **观测层独立噪声也会使残差非零。** 定理 6.11 的"恒为 0"是
     A_D ∧ **共享驱动** 的结论；若 E、U 各自带独立观测噪声，残差出现
     噪声地板 κ·sigma_obs，会被误读为 ¬A_D。
     修正：**重复测量后平均**再算残差 ⟹ 噪声底按 1/sqrt(N) 下降，
     而系统性弯曲（eps）不降 ⟹ 收敛到残差的系统值（臂 H 验证）。

运行：  python r24_ad_refutation_signature.py
产物：  r24_results.json（写入当前工作目录）
"""

import json
import math
import random
import zlib

# ------------------------------------------------------------------ 参数

K = 241                 # 步数（k = 0 .. K-1）
DT = 0.01
R_HI = 9.0
R_LO = 1.0 / 9.0
LN_HI = math.log(R_HI)
LN_LO = math.log(R_LO)
L_SPAN = abs(LN_LO - LN_HI)          # L = |ln(r_lo/r_hi)| = 4.394449...
E_START = 5.0
# U_START 只影响是否触底，不影响峰位（峰位只由增量序列决定）。
# 抬高初值以隔离出纯模型行为（承 R29 的处置）。
U_START = 30.0
N_REP = 400
SEED = 20261021
S_GRID = [0.0, 0.10, 0.20, 0.40, 0.80]
SMOOTH_W = 21
PROM = 0.5

# 跨任务臂用的 rho 网格（保证 theta_U = rho*theta_E 落在 r 窗口内部）
RHO_GRID = [1.0 / 3.0, 0.5, 1.0, 2.0, 3.0]

# 任务的 theta_E 取值（"不同任务"的操作化）
TASK_THETA_E = [0.5, 1.0, 2.0]

# 违反强度网格
EPS_GRID = [0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40]

# rho 任务依赖强度（跨任务 ln rho 的变异半幅）
ETA_GRID = [0.0, 0.05, 0.10, 0.20, 0.40, 0.80]

# 臂 E 视界扫描用的细网格
EPS_FINE = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15,
            0.18, 0.20, 0.22, 0.25, 0.28, 0.32, 0.40, 0.60]

# 线性失配残差的判定阈值（定理 6.11：A_D + 共享驱动 ⟹ 残差恒为 0）
RESID_TOL = 1e-9


def zlib_off(s):
    """可复现的字符串散列。不能用 hash()——Python 的 str hash 按进程随机化
    （PYTHONHASHSEED），会破坏可复现性。见陷阱 T1。"""
    return zlib.crc32(s.encode("utf-8")) % 9973


def fmt(v, n=4):
    if v is None:
        return "n/a"
    return "%.*f" % (n, v)


def sig_str(d, tol=None):
    """有符号归一化峰距 -> 签名字符串。"""
    tol = 1.5 / (K - 1.0) if tol is None else tol
    if abs(d) < tol:
        return "degenerate"
    return "(+,-)" if d > 0 else "(-,+)"


# ------------------------------------------------------------------ 配置

def make_cfg(theta_E, theta_U, e_plus=1.0, u_plus=1.0):
    """给定 theta_E / theta_U 与 e_plus / u_plus，构造一组 (e±, u±)。

    校验：rho_exp = u_plus/e_plus
          rho_res = u_minus/e_minus = (theta_U*u_plus)/(theta_E*e_plus)
          ⟹ rho_res/rho_exp = theta_U/theta_E  ✓
    """
    e_minus = theta_E * e_plus
    u_minus = theta_U * u_plus
    return {
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": e_minus / e_plus,
        "theta_U": u_minus / u_plus,
        "rho_exp": u_plus / e_plus,
        "rho_res": u_minus / e_minus,
        "rho": (u_minus / e_minus) / (u_plus / e_plus),
    }


# ------------------------------------------------------------------ 调度

def log_ramp(k, KK=None):
    KK = K if KK is None else KK
    t = k / (KK - 1.0)
    return math.exp(LN_HI + (LN_LO - LN_HI) * t)


def psi_of_r(r, e_plus, theta_E, eps):
    """Psi(r) := dE/dt 在无噪声下的值。eps=0 即 e_plus*(r-theta_E)（A_D）。"""
    return e_plus * (r - theta_E) - eps * r * r


def psi_roots(e_plus, theta_E, eps):
    """Psi(r)=0 的解析根（升序）。eps=0 时恰 1 个；eps>=1/(4 theta_E) 时 0 个。"""
    if eps <= 0.0:
        return [theta_E]
    disc = e_plus * e_plus - 4.0 * eps * e_plus * theta_E
    if disc <= 0.0:
        return []
    sq = math.sqrt(disc)
    return sorted([(e_plus - sq) / (2.0 * eps),
                   (e_plus + sq) / (2.0 * eps)])


def eps_critical(e_plus, theta_E):
    """临界违反强度 eps_c = e_+/(4*theta_E)。

    Psi(r) = e_+*(r-theta_E) - eps*r^2 是开口向下的抛物线，
    顶点值 Psi_max = e_+^2/(4 eps) - e_+*theta_E。
        eps < eps_c  ⟹ Psi_max > 0  ⟹ 两个实根，E 有内部上升相
        eps >= eps_c ⟹ Psi_max <= 0 ⟹ **无实根**，E 全程下降、峰退到 k=0

    ⟹ eps_c 是一个**解析阈值**，且 eps_c ∝ 1/theta_E。
    """
    return e_plus / (4.0 * theta_E)


# ------------------------------------------------------------------ 轨迹

def trajectory(rng, cfg, rbar_of_k, s=0.0, KK=None, indep_rel=0.0,
             eps=0.0, r_out=None):
    """A_D + A_D^stat（噪声进过程层）。

    s         —— **共享**噪声：同一 r 过程同时驱动 dE 与 dU
    indep_rel —— **独立**噪声：分别加在 dE / dU 上，但**每步波动率与共享噪声匹配**
                 （sd = indep_rel * dt * rbar_k * 各自系数），两次独立抽样
    eps       —— 违反 A_D 的强度：dE <- dE - dt*eps*r^2
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
        if r_out is not None:
            r_out.append(r)
        dE = DT * (r * cfg["e_plus"] - cfg["e_minus"] - eps * r * r)
        dU = DT * (r * cfg["u_plus"] - cfg["u_minus"])
        if indep_rel > 0.0:
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


def count_valleys(seq, min_prom):
    neg = [-v for v in seq]
    return sum(1 for i in local_maxima(neg)
               if prominence(neg, i) > min_prom)


def is_interior_peak(seq, idx):
    """峰是否为**内点**（不贴窗口边界）。贴边即不可辨识（陷阱 T10）。"""
    return 0 < idx < len(seq) - 1


def increments(seq, start):
    out = [seq[0] - start]
    for i in range(1, len(seq)):
        out.append(seq[i] - seq[i - 1])
    return out


def linear_mismatch(Es, Us, rho_exp):
    """定理 6.11 的**线性失配**：min_v ||Pi^perp(v U - E)||。

    在 A_D 下，取 v = 1/rho_exp 时 v*U(k) - E(k) 关于 k **仿射**，
    故投影到 span{1,k} 的正交补后残差恒为 0。
    返回 (最小残差 RMS, 最优 v)。
    """
    n = len(Es)
    xm = (n - 1) * 0.5
    # 先把 U、E 分别中心化并对 k 去趋势（即 Pi^perp）
    def resid(seq):
        ym = sum(seq) / n
        sxx = sum((i - xm) ** 2 for i in range(n))
        sxy = sum((i - xm) * (seq[i] - ym) for i in range(n))
        b = sxy / sxx if sxx > 0 else 0.0
        a = ym - b * xm
        return [seq[i] - (a + b * i) for i in range(n)]

    rU, rE = resid(Us), resid(Es)
    uu = sum(v * v for v in rU)
    if uu <= 0.0:
        return 0.0, 1.0
    vstar = sum(a * b for a, b in zip(rU, rE)) / uu
    ss = sum((vstar * a - b) ** 2 for a, b in zip(rU, rE))
    return math.sqrt(ss / n), vstar


# ------------------------------------------------------------------ 单次测量

def measure_once(rng, theta_E, rho, s=0.0, indep_rel=0.0, eps=0.0,
                KK=None, rfun=None, u_plus=1.0, e_plus=1.0, dead=0.02):
    """在一个任务上跑一条轨迹，返回签名与诊断量。

    返回
      d        —— 有符号归一化峰距 (k*_E - k*_U)/(K-1)
      sig      —— "(+,-)" / "(-,+)" / "degenerate" / "non_interior"
      nmin     —— E 轨迹显著内部极小值个数（A_D ⟹ 0）
      resid    —— 任务内线性失配 RMS（A_D + 共享驱动 ⟹ 0）
      vstar    —— 最优缩放；A_D ⟹ 1/rho_exp
      sig_set  —— 该任务内观测到的全部 (sign dE, sign dU) 组合
    """
    KK = K if KK is None else KK
    rf = rfun or log_ramp
    cfg = make_cfg(theta_E, rho * theta_E, e_plus=e_plus, u_plus=u_plus)
    Es, Us, clipped = trajectory(
        rng, cfg, lambda k: rf(k, KK), s=s, KK=KK,
        indep_rel=indep_rel, eps=eps)

    kE_lo, kE_hi = argmax_plateau(Es)
    kU_lo, kU_hi = argmax_plateau(Us)
    kE = 0.5 * (kE_lo + kE_hi)
    kU = 0.5 * (kU_lo + kU_hi)
    d = (kE - kU) / (KK - 1.0)

    interior = (is_interior_peak(Es, kE_lo) and is_interior_peak(Es, kE_hi)
                 and is_interior_peak(Us, kU_lo) and is_interior_peak(Us, kU_hi))
    sig = "non_interior" if not interior else sig_str(d)

    smE = moving_average(Es, SMOOTH_W)
    nmin = count_valleys(smE, PROM)

    resid, vstar = linear_mismatch(Es, Us, cfg["rho_exp"])

    # 任务内观测到的全部符号组合（用平滑增量 + 死区过滤）
    dE = increments(Es, E_START)
    dU = increments(Us, U_START)
    smdE = moving_average(dE, SMOOTH_W)
    smdU = moving_average(dU, SMOOTH_W)
    scE = max(abs(v) for v in smdE) or 1.0
    scU = max(abs(v) for v in smdU) or 1.0
    sig_set = set()
    for a, b in zip(smdE, smdU):
        if abs(a) > dead * scE and abs(b) > dead * scU:
            sig_set.add(("+" if a > 0 else "-", "+" if b > 0 else "-"))

    return dict(d=d, sig=sig, nmin=nmin, resid=resid, vstar=vstar,
                clipped=clipped, kE=kE, kU=kU,
                sig_set=sorted("(" + a + "," + b + ")" for a, b in sig_set))


def rng_for(seed_off=0):
    return random.Random(SEED + seed_off)


# ------------------------------------------------------------------ 臂 A：预言复核（无噪声）

def arm_A_prediction(u_plus_grid=(1.0, 2.0, 0.5)):
    """A_D 成立时：签名只由 rho 决定，跨任务（不同 theta_E）不变。

    预言：
      (P1) 对每个 rho，三个不同 theta_E 的任务给出**同一个**签名
      (P2) 有符号归一化峰距 d = ln(rho)/L，与 theta_E 无关
      (P3) 任务内线性失配残差恒为 0（定理 6.11）
      (P4) E 轨迹无内部极小（单峰）
    且 (P2)(P3)(P4) 应对 u_plus 的一般取值也成立（不只 u_plus = e_plus = 1）
    """
    rows = []
    for up in u_plus_grid:
        for rho in RHO_GRID:
            pred_d = math.log(rho) / L_SPAN
            ptasks = []
            for th in TASK_THETA_E:
                rng = rng_for(zlib_off("A|%.4f|%.4f|%.4f" % (up, rho, th)))
                m = measure_once(rng, th, rho, s=0.0, eps=0.0, u_plus=up)
                ptasks.append({
                    "theta_E": th,
                    "theta_U": rho * th,
                    "d_obs": m["d"],
                    "sig_obs": m["sig"],
                    "n_interior_min": m["nmin"],
                    "lin_resid_rms": m["resid"],
                    "vstar": m["vstar"],
                    "U_clipped_steps": m["clipped"],
                })
            sigs = set(t["sig_obs"] for t in ptasks)
            rows.append({
                "u_plus": up,
                "rho_exp_true": up / 1.0,
                "rho": rho,
                "pred_d": pred_d,
                "pred_sig": sig_str(pred_d),
                "tasks": ptasks,
                "cross_task_consistent": len(sigs) == 1,
                "max_abs_d_err": max(abs(t["d_obs"] - pred_d)
                                     for t in ptasks),
                "max_lin_resid_rms": max(t["lin_resid_rms"]
                                         for t in ptasks),
                "max_interior_min": max(t["n_interior_min"]
                                        for t in ptasks),
                "max_vstar_err": max(abs(t["vstar"] - 1.0 / (up / 1.0))
                                     for t in ptasks),
            })
    return rows


# ------------------------------------------------------------------ 臂 B：混淆（★核心）

def arm_B_confound(n_rep=None, s_grid=None):
    """A_D **逐任务成立**但 rho 依赖任务 ⟹ 跨任务签名并存（R24 的假阳来源）。

    构造：理解者 Omega 在第 i 个任务上的比值
              rho_i = rho_0 * exp(eta * z_i),   z_i ∈ {-1, 0, +1}
          每个任务内部 A_D **严格成立**（线性叠加未被破坏）。

    三个判据：
      B1  跨任务签名并存率（R24 原登记）   —— 预期：eta 大时为 1.0（假阳）
      B2  任务内线性失配残差（定理 6.11）   —— 预期：恒为 0（正确判定 A_D 成立）
      B3  任务内签名唯一性                  —— 预期：恒唯一（正确判定 A_D 成立）
    """
    n_rep = n_rep or N_REP
    s_grid = s_grid or S_GRID
    blocks = []
    for eta in ETA_GRID:
        rows = []
        for rho0 in (1.0, 2.0):
            rhos = [rho0 * math.exp(eta * z) for z in (-1.0, 0.0, 1.0)]
            for s in s_grid:
                both = 0
                within_multi = 0
                resids = []
                for rep in range(n_rep):
                    rng = rng_for(zlib_off("B|%.4f|%.4f|%.4f|%d"
                                           % (eta, rho0, s, rep)))
                    sigs = []
                    multi = False
                    mx = 0.0
                    for th, rh in zip(TASK_THETA_E, rhos):
                        m = measure_once(rng, th, rh, s=s, u_plus=1.0)
                        sigs.append(m["sig"])
                        sset = set(m["sig_set"])
                        # T13 修正：sig_set 的元素已是字符串 "(+,-)" 形式，
                        # 早先用元组 ("+","-") 比对导致本判据恒为 0（假阴性）。
                        if "(+,-)" in sset and "(-,+)" in sset:
                            multi = True
                        mx = max(mx, m["resid"])
                    if multi:
                        within_multi += 1
                    resids.append(mx)
                    if "(+,-)" in sigs and "(-,+)" in sigs:
                        both += 1
                resids.sort()
                rows.append({
                    "rho0": rho0, "s": s,
                    "task_rhos": [round(x, 6) for x in rhos],
                    "r24_cross_task_both_rate": both / n_rep,
                    "within_task_multi_sig_rate": within_multi / n_rep,
                    "within_task_max_resid_rms": resids[-1],
                })
        blocks.append({"eta": eta, "rows": rows})
    return blocks


# ------------------------------------------------------------------ 臂 C：分离判据的判别力

def arm_C_separation(n_rep=None, s_grid=None):
    """注入违反 A_D 的二次项 -eps*r^2，比较两个判据的灵敏度。

    判据 1（R24 原登记：跨任务签名并存）—— 需要跨任务 rho 有变异；
           在**单任务**下**结构性无法触发**。这是它的缺陷。
    判据 2（定理 6.11：任务内线性失配残差）—— 单任务即可触发。
    判据 0（内部极小：E 出现显著内部极小值）—— 本注入族下**结构性不可检**，
           原因是 Psi 的较大根 r_2 随 eps 增大冲出 r 窗口顶端，
           把"下降段"压缩进轨迹最左端不足 1 步的区间（陷阱 T15 的量化）。
    判据 1（峰退到左边界：k*_E = 0）—— eps >= eps_c 时**可检**（阈值型）。
    判据 2（定理 6.11 任务内线性失配残差）—— 对**任意** eps > 0 可检，
           且量级 ∝ eps ⟹ 灵敏度最高，是**首选判据**。
    """
    n_rep = n_rep or N_REP
    s_grid = s_grid or S_GRID
    blocks = []
    for eps in EPS_GRID:
        roots = psi_roots(1.0, 1.0, eps)
        in_win = [r for r in roots if R_LO < r < R_HI]
        rows = []
        for s in s_grid:
            rate = 0
            both = 0
            edge = 0
            resids = []
            for rep in range(n_rep):
                rng = rng_for(zlib_off("C|%.4f|%.4f|%d" % (eps, s, rep)))
                m = measure_once(rng, 1.0, 1.0, s=s, eps=eps)
                if m["nmin"] > 0:
                    rate += 1
                if m["kE"] <= 0.5:
                    edge += 1
                sset = set(m["sig_set"])
                if "(+,-)" in sset and "(-,+)" in sset:
                    both += 1
                resids.append(m["resid"])
            resids.sort()
            rows.append({
                "s": s,
                "within_task_interior_min_rate": rate / n_rep,
                "within_task_both_sig_rate": both / n_rep,
                "peak_at_left_edge_rate": edge / n_rep,
                "lin_resid_rms_median": resids[len(resids) // 2],
            })
        blocks.append({
            "eps": eps,
            "roots_analytic": [round(r, 6) for r in roots],
            "n_roots": len(roots),
            "n_roots_in_window": len(in_win),
            "rows": rows,
        })
    return blocks


# ------------------------------------------------------------------ 臂 D：共享 vs 独立噪声

def arm_D_noise_regime(n_rep=None):
    """定理 6.11 的守恒依赖「共享驱动」。

    共享噪声：线性失配残差**逐样本恒为 0**
    独立噪声：残差出现噪声底，判据退化为统计判据
    """
    n_rep = n_rep or N_REP
    rows = []
    for s in S_GRID:
        for rho in (1.0 / 3.0, 1.0, 3.0):
            for regime, kw in (("shared", {"s": s, "indep_rel": 0.0}),
                               ("indep", {"s": 0.0, "indep_rel": s})):
                resids = []
                for rep in range(n_rep):
                    rng = rng_for(zlib_off("D|%.4f|%.4f|%s|%d"
                                           % (s, rho, regime, rep)))
                    m = measure_once(rng, 1.0, rho, **kw)
                    resids.append(m["resid"])
                resids.sort()
                rows.append({
                    "s": s, "rho": rho, "regime": regime,
                    "lin_resid_rms_max": resids[-1],
                    "lin_resid_rms_median": resids[len(resids) // 2],
                })
    return rows


# ------------------------------------------------------------------ 臂 E：可探测视界

def arm_E_horizon(n_rep=None, s_grid=None):
    """eps 过小 → 两根之一落在 r 窗口外 → 内部极小不可见（漏检）
       eps 过大 → Psi 无实根 → E 无峰（另类可探测签名）

    ⟹ 探测率关于 eps **非单调**：只有中间窗口能靠"内部极小"探测。
    """
    n_rep = n_rep or N_REP
    s_grid = s_grid or S_GRID
    blocks = []
    for s in s_grid:
        rows = []
        for eps in EPS_FINE:
            roots = psi_roots(1.0, 1.0, eps)
            in_win = [r for r in roots if R_LO < r < R_HI]
            rate = 0
            edge = 0
            resids = []
            for rep in range(n_rep):
                rng = rng_for(zlib_off("E|%.4f|%.4f|%d" % (s, eps, rep)))
                m = measure_once(rng, 1.0, 1.0, s=s, eps=eps)
                if m["nmin"] > 0:
                    rate += 1
                if m["kE"] <= 0.5:
                    edge += 1
                resids.append(m["resid"])
            resids.sort()
            rows.append({
                "eps": eps,
                "n_roots": len(roots),
                "n_roots_in_window": len(in_win),
                "detect_rate_via_interior_min": rate / n_rep,
                "peak_at_left_edge_rate": edge / n_rep,
                "lin_resid_rms_median": resids[len(resids) // 2],
            })
        blocks.append({"s": s, "rows": rows})
    return blocks


# ------------------------------------------------------------------ 臂 E2：阈值标度（★）

W_GRID = [0.20, 0.40, 0.60, 0.80, 0.90, 1.00, 1.10, 1.30, 1.60, 2.00]
THETA_E_GRID = [0.5, 1.0, 2.0]


def arm_E2_threshold(n_rep=None, s_grid=(0.0, 0.20), theta_grid=None):
    """探测翻转点 eps_c 的**解析值 vs 观测值**，以及标度不变性。

    解析（推论 6.8）：eps_c = e_+ / (4*theta_E)。
    观测：在归一化违反强度 w = eps/eps_c 的网格上，找第一个使
          「E 峰退到左边界（k*_E = 0）」的探测率 > 0.5 的 w。

    预言
      P6a（阈值位置）观测翻转点 w* = 1，即 eps* = e_+/(4*theta_E)
      P6b（标度不变）探测率作为 w 的函数**与 theta_E 无关**
                     ——theta_E = 0.5/1.0/2.0 三条曲线应重合
      P6c（残差线性）resid/eps 为常数（与 eps、theta_E 均无关）

    三条曲线重合是**结构性**预言：eps 与 theta_E 各自变化两个数量级，
    但归一化后的翻转位置不动 ⟹ 锁的是结构，不是数字。
    """
    n_rep = n_rep or N_REP
    theta_grid = theta_grid or THETA_E_GRID
    blocks = []
    for th in theta_grid:
        ec = eps_critical(1.0, th)
        rows = []
        for s in s_grid:
            for w in W_GRID:
                eps = w * ec
                edge = 0
                resids = []
                for rep in range(n_rep):
                    rng = rng_for(zlib_off(
                        "E2|%.4f|%.4f|%.4f|%d" % (th, w, s, rep)))
                    m = measure_once(rng, th, 1.0, s=s, eps=eps)
                    if m["kE"] <= 0.5:
                        edge += 1
                    resids.append(m["resid"])
                resids.sort()
                med = resids[len(resids) // 2]
                rows.append({
                    "s": s, "w": w, "eps": eps,
                    "peak_at_left_edge_rate": edge / n_rep,
                    "lin_resid_rms_median": med,
                    "resid_over_eps": (med / eps) if eps > 0 else None,
                })
        blocks.append({"theta_E": th, "eps_c_analytic": ec, "rows": rows})
    return blocks


# ------------------------------------------------------------------ 臂 H：观测噪声地板与平均律（T16 修正）

SIG_OBS_GRID = [0.0, 0.10, 0.40]
N_AVG_GRID = [1, 2, 4, 8, 16, 32]


def arm_H_obs_floor(n_rep=None, eps_grid=(0.0, 0.05),
                    sig_grid=None, n_grid=None):
    """定理 6.11 的"残差恒为 0"只在 **A_D ∧ 共享驱动** 下成立（陷阱 T16）。

    若 E、U 各自带**独立观测噪声**（sigma_obs），残差出现噪声地板：

        resid(N)² ≈ resid_sys(eps)² + kappa * sigma_obs² / N

    其中 N 是**重复测量后平均**的次数。

    预言
      P7a  eps = 0 时：resid → 0，且 log resid 对 log N 的斜率 ≈ -1/2
      P7b  eps > 0 时：resid → resid_sys(eps) > 0（**地板不降到 0**）
      P7c  sigma_obs = 0 时：resid 与 N 无关，恒等于系统值

    ⟹ 操作规程（补篇九 §7）必须先重复测量平均、再算残差，
       并把残差与**由重复测量的组内散度估出的噪声底**比较，
       而不是与一个绝对常数比较。
    """
    n_rep = n_rep or 120
    sig_grid = sig_grid or SIG_OBS_GRID
    n_grid = n_grid or N_AVG_GRID
    rows = []
    for eps in eps_grid:
        for sg in sig_grid:
            for N in n_grid:
                resids = []
                for rep in range(n_rep):
                    rng = rng_for(zlib_off("H|%.4f|%.4f|%d|%d"
                                           % (eps, sg, N, rep)))
                    cfg = make_cfg(1.0, 1.0)
                    Eacc = [0.0] * K
                    Uacc = [0.0] * K
                    for _ in range(N):
                        Es, Us, _ = trajectory(
                            rng, cfg, lambda k: log_ramp(k),
                            s=0.0, indep_rel=sg, eps=eps)
                        for i in range(K):
                            Eacc[i] += Es[i]
                            Uacc[i] += Us[i]
                    Ebar = [v / N for v in Eacc]
                    Ubar = [v / N for v in Uacc]
                    r_, _ = linear_mismatch(Ebar, Ubar, cfg["rho_exp"])
                    resids.append(r_)
                resids.sort()
                rows.append({
                    "eps": eps, "sigma_obs": sg, "N": N,
                    "lin_resid_rms_median": resids[len(resids) // 2],
                })
    return rows


def fit_loglog_slope(rows, eps, sig):
    """对 eps/sigma_obs 固定的一组 (N, resid) 做 log-log 斜率拟合。

    返回 (斜率, 截距)；点数不足或含非正值返回 None。
    """
    pts = [(r["N"], r["lin_resid_rms_median"]) for r in rows
           if r["eps"] == eps and r["sigma_obs"] == sig
           and r["lin_resid_rms_median"] > 0.0]
    if len(pts) < 3:
        return None
    xs = [math.log(n) for n, _ in pts]
    ys = [math.log(v) for _, v in pts]
    a, b = _ols(xs, ys)
    return b, a


# ------------------------------------------------------------------ 臂 F：操作规程端到端

def verdict(tasks, n_rep=None, s=0.0, indep_rel=0.0, eps=0.0):
    """完整操作规程的三分类判决。

    tasks : [(theta_E, rho), ...]

    判决规则（补篇九 §7）
      1. 任务内线性失配 RMS > 阈值  ⟹ reject_A_D（与任务结构无关，逐样本）
      2. 否则，若跨任务签名并存      ⟹ rho_task_dependent（A_D 仍成立）
      3. 否则                        ⟹ accept（与合取一致）
    """
    n_rep = n_rep or N_REP
    within_bad = 0
    cross_bad = 0
    resid_med = []
    for rep in range(n_rep):
        rng = rng_for(zlib_off("F|%s|%.4f|%.4f|%.4f|%d"
                               % (str(tasks), s, indep_rel, eps, rep)))
        sigs = []
        rmax = 0.0
        for th, rh in tasks:
            m = measure_once(rng, th, rh, s=s, indep_rel=indep_rel, eps=eps)
            sigs.append(m["sig"])
            rmax = max(rmax, m["resid"])
        resid_med.append(rmax)
        if rmax > RESID_TOL:
            within_bad += 1
        if "(+,-)" in sigs and "(-,+)" in sigs:
            cross_bad += 1
    resid_med.sort()
    within_rate = within_bad / n_rep
    cross_rate = cross_bad / n_rep
    if within_rate > 0.5:
        v = "reject_A_D"
    elif cross_rate > 0.5:
        v = "rho_task_dependent"
    else:
        v = "accept"
    return {
        "within_task_resid_rate": within_rate,
        "within_task_resid_rms_median": resid_med[len(resid_med) // 2],
        "r24_cross_task_both_rate": cross_rate,
        "verdict": v,
    }


def arm_F_sop(n_rep=None, s_grid=None):
    """三元组数据集上的混淆矩阵：

      D1  A_D 成立 + rho 任务不变        ⟹ 应判 accept
      D2  A_D 被违反 (eps = 0.20)        ⟹ 应判 reject_A_D
      D3  A_D 逐任务成立 + rho 任务依赖   ⟹ 应判 rho_task_dependent
         （R24 原签名会把它误判成 reject_A_D —— 这就是 T12）
    """
    n_rep = n_rep or N_REP
    s_grid = s_grid or S_GRID
    rows = []
    for s in s_grid:
        t_invar = [(th, 2.0) for th in TASK_THETA_E]
        r1 = verdict(t_invar, n_rep=n_rep, s=s, eps=0.0)
        r2 = verdict(t_invar, n_rep=n_rep, s=s, eps=0.20)
        t_dep = [(0.5, 0.45), (1.0, 1.0), (2.0, 2.2)]
        r3 = verdict(t_dep, n_rep=n_rep, s=s, eps=0.0)
        rows.append({
            "s": s,
            "D1_A_D_invariant": r1,
            "D2_A_D_violated": r2,
            "D3_rho_task_dependent": r3,
            "r24_alone_would_misjudge_D3_as_reject": (
                r3["r24_cross_task_both_rate"] > 0.5),
        })
    return rows


# ------------------------------------------------------------------ 臂 G：四系数辨识（闭合 B-B23 的条件版）

def arm_G_identification(n_rep=None, s_grid=None, u_plus_grid=(0.5, 1.0, 2.0)):
    """给定 r 调度，(e+, e-, u+, u-) 是否可单独辨识？

    方法：**不是**峰距（推论 6.5′ 已证其只给比值），而是
    把增量回归到 r 上：
        dE_k/dt = e_+ * r_k - e_-        ⟹ 斜率 = e_+，截距 = -e_-
        dU_k/dt = u_+ * r_k - u_-        ⟹ 斜率 = u_+，截距 = -u_-

    两种观测体制（这是 B-B23 能否闭合的**分水岭**）：
      obs = "r"     —— r_k **可观测**（过程层可知）：回归到**实现值** r_k
      obs = "rbar"  —— r_k **不可观测**：只能回归到**名义调度** rbar_k

    预言
      P8a  obs="r"    ⟹ 四个系数**精确**复原（相对误差 ~1e-15），与 s 无关
      P8b  obs="rbar" ⟹ 乘性对数正态噪声满足 E[r_k | rbar_k] = rbar_k，
                        OLS 仍**一致**，但估计量带 O(s/sqrt(K)) 的标准误
                        ⟹ 误差随 s 增长而**不随 N 消失**（系统性偏倚为零，
                           但每次实现的抖动不因重复而平均掉到 0 以下）

    ⟹ **在 r 调度已知/可控的前提下，B-B23 闭合。**
       若 r 不可观测，则仍只辨识 rho（承 推论 6.5′）。——登记为 B-B41。
    """
    n_rep = n_rep or N_REP
    s_grid = s_grid or S_GRID
    rows = []
    for up in u_plus_grid:
        for rho in (1.0 / 3.0, 1.0, 3.0):
            for s in s_grid:
                for obs in ("r", "rbar"):
                    e_acc = [0.0, 0.0]
                    u_acc = [0.0, 0.0]
                    for rep in range(n_rep):
                        rng = rng_for(zlib_off(
                            "G|%.4f|%.4f|%.4f|%s|%d" % (up, rho, s, obs, rep)))
                        cfg = make_cfg(1.0, rho, e_plus=1.0, u_plus=up)
                        rseq = []
                        Es, Us, _ = trajectory(
                            rng, cfg, lambda k: log_ramp(k),
                            s=s, eps=0.0, r_out=rseq)
                        dE = increments(Es, E_START)
                        dU = increments(Us, U_START)
                        xs = (list(rseq) if obs == "r"
                              else [log_ramp(float(i)) for i in range(K)])
                        ys_E = [v / DT for v in dE]
                        ys_U = [v / DT for v in dU]
                        aE, bE = _ols(xs, ys_E)
                        aU, bU = _ols(xs, ys_U)
                        e_acc[0] += bE
                        e_acc[1] += -aE
                        u_acc[0] += bU
                        u_acc[1] += -aU
                    e_hat = [v / n_rep for v in e_acc]
                    u_hat = [v / n_rep for v in u_acc]
                    cfg = make_cfg(1.0, rho, e_plus=1.0, u_plus=up)
                    rho_exp_hat = u_hat[0] / e_hat[0]
                    rho_res_hat = u_hat[1] / e_hat[1]
                    rows.append({
                        "u_plus": up, "rho": rho, "s": s, "obs": obs,
                        "e_plus_hat": e_hat[0], "e_minus_hat": e_hat[1],
                        "u_plus_hat": u_hat[0], "u_minus_hat": u_hat[1],
                        "rho_exp_hat": rho_exp_hat,
                        "rho_res_hat": rho_res_hat,
                        "rel_err_rho_exp": abs(rho_exp_hat - cfg["rho_exp"]) / cfg["rho_exp"],
                        "rel_err_rho_res": abs(rho_res_hat - cfg["rho_res"]) / cfg["rho_res"],
                    })
    return rows


LAM_GRID = [0.25, 0.5, 1.0, 2.0, 4.0]


def arm_G2_gauge(lams=None, rho=2.0, s=0.0, theta_E0=1.0):
    """定理 6.12：A_D 的**规范自由度**（标定不可辨识性）。

    考虑变换（lambda > 0）

        G_lambda :  r_k  ->  lambda * r_k
                    e_+  ->  e_+ / lambda
                    u_+  ->  u_+ / lambda
                    e_-, u_-  不变

    则
        dE = dt*(r*e_+ - e_-)  ->  dt*(lambda*r * e_+/lambda - e_-) = dE
        dU = dt*(r*u_+ - u_-)  ->  dt*(lambda*r * u_+/lambda - u_-) = dU

    ⟹ (E(k), U(k)) 序列**逐点不变**，而

        theta_E -> lambda*theta_E ,  theta_U -> lambda*theta_U
        e_+     -> e_+/lambda     ,  u_+     -> u_+/lambda

    ⟹ **四个系数 (e_+, e_-, u_+, u_-) 从 (E, U) 单独不可辨识**；
       可辨识的只有规范不变量
            rho_exp = u_+/e_+ ,  rho_res = u_-/e_- ,  rho = theta_U/theta_E
       （三者满足 rho = rho_res/rho_exp，故**独立不变量恰为 2 个**）。

    这与推论 6.4 / 6.5 一致：峰序与峰距只给 rho，不给四个系数。

    ⟹ **B-B23 不闭合。** 其闭合条件是「r 的绝对尺度被外部标定」，
       而这正是补篇四讨论的**标定通约性**问题。
       未标定 ⟹ 登记为 **B-B41**；已标定 ⟹ 臂 G 的回归法闭合 B-B23。
    """
    lams = lams or LAM_GRID
    cfg0 = make_cfg(theta_E0, rho * theta_E0, e_plus=1.0, u_plus=1.0)
    rng0 = rng_for(zlib_off("G2|base"))
    E0, U0, _ = trajectory(rng0, cfg0, lambda k: log_ramp(k), s=s, eps=0.0)
    rows = []
    for lam in lams:
        cfgL = make_cfg(lam * theta_E0, lam * rho * theta_E0,
                        e_plus=1.0 / lam, u_plus=1.0 / lam)
        rngL = rng_for(zlib_off("G2|%.4f" % lam))
        EL, UL, _ = trajectory(rngL, cfgL, lambda k: lam * log_ramp(k),
                               s=s, eps=0.0)
        dE = max(abs(a - b) for a, b in zip(EL, E0))
        dU = max(abs(a - b) for a, b in zip(UL, U0))
        rows.append({
            "lambda": lam,
            "theta_E": cfgL["theta_E"],
            "theta_U": cfgL["theta_U"],
            "e_plus": cfgL["e_plus"],
            "u_plus": cfgL["u_plus"],
            "e_minus": cfgL["e_minus"],
            "u_minus": cfgL["u_minus"],
            "rho_exp": cfgL["rho_exp"],
            "rho_res": cfgL["rho_res"],
            "rho": cfgL["rho"],
            "max_abs_E_diff": dE,
            "max_abs_U_diff": dU,
        })
    return rows


def _ols(xs, ys):
    n = len(xs)
    xm = sum(xs) / n
    ym = sum(ys) / n
    sxx = sum((x - xm) ** 2 for x in xs)
    if sxx <= 0.0:
        return ym, 0.0
    sxy = sum((x - xm) * (y - ym) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = ym - b * xm
    return a, b


# ------------------------------------------------------------------ main

def main():
    random.seed(SEED)

    print("=" * 70)
    print("R24  A_D 的证否签名：判别力分解与操作规程")
    print("=" * 70)

    print("")
    print("--- A  预言复核（无噪声）：签名只由 rho 决定，跨任务不变 ---")
    A = arm_A_prediction()
    print("  {:>7} {:>7} {:>10} {:>10} | {:>9} | {:>9} {:>10} {:>5}".format(
        "u_plus", "rho", "pred_d", "pred_sig", "consistent",
        "d_err", "resid_rms", "nmin"))
    for r in A:
        print("  {:>7} {:>7} {:>10} {:>10} | {:>9} | {:>9} {:>10} {:>5}".format(
            r["u_plus"], r["rho"], fmt(r["pred_d"]), r["pred_sig"],
            "yes" if r["cross_task_consistent"] else "NO",
            fmt(r["max_abs_d_err"]), fmt(r["max_lin_resid_rms"], 2),
            r["max_interior_min"]))

    print("")
    print("--- B  混淆臂（★）：A_D 逐任务成立 + rho 任务依赖 ---")
    print("    判据 B1 = R24 原登记（跨任务并存）   B2 = 任务内线性失配   B3 = 任务内签名唯一性")
    B = arm_B_confound(n_rep=200)
    print("  {:>6} {:>6} {:>5} | {:>10} {:>10} {:>10} | {}".format(
        "eta", "rho0", "s", "B1(假阳)", "B2", "B3", "R24 单独会误判?"))
    for blk in B:
        for r in blk["rows"]:
            print("  {:>6} {:>6} {:>5} | {:>10} {:>10} {:>10} | {}".format(
                blk["eta"], r["rho0"], r["s"],
                fmt(r["r24_cross_task_both_rate"]),
                fmt(r["within_task_max_resid_rms"], 2),
                fmt(r["within_task_multi_sig_rate"]),
                "YES" if r["r24_cross_task_both_rate"] > 0.5 else "no"))

    print("")
    print("--- C  分离判据的判别力（违反强度 eps）---")
    print("    root/inwin = 解析根数; min_rate = 内部极小法; both_rate = 任务内双象限法;")
    print("    edge_rate = 峰退到左边界法; resid = 任务内线性失配（定理 6.11）")
    C = arm_C_separation(n_rep=200)
    print("  {:>6} {:>5} {:>5} {:>5} | {:>9} {:>9} {:>9} | {:>12}".format(
        "eps", "root", "inwin", "s", "min_rate", "both_rate", "edge_rate",
        "resid_rms"))
    for blk in C:
        for r in blk["rows"]:
            print("  {:>6} {:>5} {:>5} {:>5} | {:>9} {:>9} {:>9} | {:>12}".format(
                blk["eps"], blk["n_roots"], blk["n_roots_in_window"],
                r["s"], fmt(r["within_task_interior_min_rate"]),
                fmt(r["within_task_both_sig_rate"]),
                fmt(r["peak_at_left_edge_rate"]),
                fmt(r["lin_resid_rms_median"], 8)))

    print("")
    print("--- D  共享 vs 独立噪声（线性失配残差）---")
    D = arm_D_noise_regime(n_rep=200)
    print("  {:>5} {:>6} {:>8} | {:>14} {:>14}".format(
        "s", "rho", "regime", "resid_median", "resid_max"))
    for r in D:
        print("  {:>5} {:>6} {:>8} | {:>14} {:>14}".format(
            r["s"], r["rho"], r["regime"],
            fmt(r["lin_resid_rms_median"], 8),
            fmt(r["lin_resid_rms_max"], 8)))

    print("")
    print("--- E  可探测视界（eps 的非单调探测窗口）---")
    E = arm_E_horizon(n_rep=120)
    print("  {:>5} | {:>24} | {:>8} | {:>24} | {}".format(
        "s", "内部极小法可探测区间", "命中", "峰退边界法可探测区间", "命中"))
    for blk in E:
        d1 = [r["eps"] for r in blk["rows"]
              if r["detect_rate_via_interior_min"] > 0.5]
        d2 = [r["eps"] for r in blk["rows"]
              if r["peak_at_left_edge_rate"] > 0.5]
        s1 = ("[%s, %s]" % (fmt(min(d1), 3), fmt(max(d1), 3))) if d1 else "全网格漏检"
        s2 = ("[%s, %s]" % (fmt(min(d2), 3), fmt(max(d2), 3))) if d2 else "全网格漏检"
        print("  {:>5} | {:>24} | {:>8} | {:>24} | {}".format(
            blk["s"], s1, len(d1), s2, len(d2)))

    print("")
    print("--- E2 阈值标度（★）：eps_c = e_+/(4*theta_E) 的解析值 vs 观测翻转点 ---")
    E2 = arm_E2_threshold(n_rep=150)
    print("  {:>6} {:>10} | {}".format("theta_E", "eps_c(解析)", "w -> 峰退到左边界的探测率"))
    for blk in E2:
        line = []
        for w in W_GRID:
            got = [r["peak_at_left_edge_rate"] for r in blk["rows"]
                   if r["w"] == w]
            line.append("%.2f:%s" % (w, fmt(max(got)) if got else "n/a"))
        print("  {:>6} {:>10} | {}".format(
            blk["theta_E"], fmt(blk["eps_c_analytic"], 4), "  ".join(line)))
    print("  预言 P6a：翻转点在 w = 1；P6b：三条曲线重合（与 theta_E 无关）")

    print("")
    print("--- F  操作规程端到端（三分类混淆矩阵）---")
    F = arm_F_sop(n_rep=200)
    print("  {:>5} | {:>14} | {:>14} | {:>18} | {}".format(
        "s", "D1 A_D-不变", "D2 A_D-违反", "D3 rho-依赖", "R24 误判 D3?"))
    for r in F:
        print("  {:>5} | {:>14} | {:>14} | {:>18} | {}".format(
            r["s"], r["D1_A_D_invariant"]["verdict"],
            r["D2_A_D_violated"]["verdict"],
            r["D3_rho_task_dependent"]["verdict"],
            "YES" if r["r24_alone_would_misjudge_D3_as_reject"] else "no"))

    print("")
    print("--- G  四系数辨识（给定 r 调度，回归法）---")
    G = arm_G_identification(n_rep=200)
    print("  {:>7} {:>6} {:>5} | {:>10} {:>10} | {:>10} {:>10}".format(
        "u_plus", "rho", "s", "e+_hat", "e-_hat", "relerr_exp", "relerr_res"))
    for r in G:
        print("  {:>7} {:>6} {:>5} | {:>10} {:>10} | {:>10} {:>10}".format(
            r["u_plus"], r["rho"], r["s"],
            fmt(r["e_plus_hat"], 6), fmt(r["e_minus_hat"], 6),
            fmt(r["rel_err_rho_exp"], 6), fmt(r["rel_err_rho_res"], 6)))

    print("")
    print("--- G2 规范自由度（定理 6.12）：r -> lambda*r, e+ -> e+/lambda ---")
    G2 = arm_G2_gauge()
    print("  {:>7} {:>9} {:>9} {:>9} {:>9} | {:>13} {:>13}".format(
        "lambda", "theta_E", "theta_U", "e+", "u+", "max|dE|", "max|dU|"))
    for r in G2:
        print("  {:>7} {:>9} {:>9} {:>9} {:>9} | {:>13} {:>13}".format(
            r["lambda"], fmt(r["theta_E"], 4), fmt(r["theta_U"], 4),
            fmt(r["e_plus"], 4), fmt(r["u_plus"], 4),
            fmt(r["max_abs_E_diff"], 10), fmt(r["max_abs_U_diff"], 10)))
    print("  预言：轨迹逐点不变（~1e-15），而四个系数变了 lam 倍 => 绝对值不可辨识")

    print("")
    print("--- H  观测噪声地板与平均律（T16 修正）---")
    H = arm_H_obs_floor(n_rep=120)
    print("  {:>6} {:>10} | {}".format("eps", "sigma_obs", "N -> resid_rms_median"))
    for eps in (0.0, 0.05):
        for sg in SIG_OBS_GRID:
            line = []
            for N in N_AVG_GRID:
                got = [r["lin_resid_rms_median"] for r in H
                       if r["eps"] == eps and r["sigma_obs"] == sg
                       and r["N"] == N]
                line.append("%2d:%.2e" % (N, got[0]) if got else "%2d:n/a" % N)
            print("  {:>6} {:>10} | {}".format(fmt(eps, 2), fmt(sg, 2),
                                               " ".join(line)))
    slopes = {}
    for eps in (0.0, 0.05):
        for sg in (0.10, 0.40):
            sl = fit_loglog_slope(H, eps, sg)
            slopes["eps=%.2f,sig=%.2f" % (eps, sg)] = (
                None if sl is None else round(sl[0], 4))
    print("  log-log 斜率 (resid ~ N^b)：")
    for k in sorted(slopes):
        print("    {} -> {}".format(k, slopes[k]))
    print("  预言 P7a：eps=0 时 b ≈ -0.5；P7b：eps>0 时 b ≈ 0（收敛到系统地板）")

    result = {
        "meta": {
            "experiment": "R24",
            "title": "A_D 证否签名的判别力分解与操作规程",
            "corresponds_to": "补篇五 §6 R24 / B-B22",
            "K": K, "DT": DT, "R_HI": R_HI, "R_LO": R_LO,
            "L_SPAN": L_SPAN, "N_REP": N_REP, "SEED": SEED,
            "S_GRID": S_GRID, "SMOOTH_W": SMOOTH_W, "PROM": PROM,
            "RHO_GRID": RHO_GRID, "TASK_THETA_E": TASK_THETA_E,
            "EPS_GRID": EPS_GRID, "ETA_GRID": ETA_GRID,
            "W_GRID": W_GRID, "THETA_E_GRID": THETA_E_GRID,
            "SIG_OBS_GRID": SIG_OBS_GRID, "N_AVG_GRID": N_AVG_GRID,
            "resid_tol": RESID_TOL,
            "nature": "A_D 原型上的推导结构检验，不是实证验证",
        },
        "A_prediction": A,
        "B_confound": B,
        "C_separation": C,
        "D_noise_regime": D,
        "E_horizon": E,
        "E2_threshold": E2,
        "F_sop": F,
        "G_identification": G,
        "G2_gauge": G2,
        "H_obs_floor": H,
        "H_loglog_slopes": slopes,
    }

    with open("r24_results.json", "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    print("")
    print("wrote r24_results.json")


if __name__ == "__main__":
    main()
