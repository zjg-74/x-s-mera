# -*- coding: utf-8 -*-
"""
R27 —— 认知过程的双相形状：X 单调 / M 单峰（补篇六）
====================================================
对应主张（用户 2026-10-19）：
  "完整有效的认知过程，不仅仅在于前期认知事实（X）的递增，
   也在于认知（M）的迭代，在增加到峰值的时候，又有个收敛的减的过程。"

本实验做什么
------------------------------------------------------------------
不验证"人是否如此"（那需要真实理解者），而是检验这句话在 A_D
（补篇五 §1.2 的双速率分解）下的**推导结构**是否自洽、是否给出
可判别的定量签名。具体五件事：

  (1) 单峰性：r(k) 单调下降时，U(k) 与 E(k) 是否均为单峰
  (2) 峰位解析式：k*(theta) = (K-1) * (ln theta - ln r_hi)/(ln r_lo - ln r_hi)
      实测峰位是否落在该预测点上
  (3) 峰序定理（推论 6.2）：sign(k*_U - k*_E) 是否由 rho_res vs rho_exp 决定
      —— 即"钻研相"与"空洞相"的中期性质是否可由峰序判别
  (4) 中期相签名：r 下降过程中的符号组合序列是否为
      (+,+) -> (+,-) -> (-,-)   [rho_res > rho_exp]
      (+,+) -> (-,+) -> (-,-)   [rho_res < rho_exp]
  (5) 对照臂：r 恒定（无消解相 / 无扩张相）时应无内点峰 —— 这是
      "只有增"与"只有减"两种残缺过程的签名

A_D 与关键量
------------------------------------------------------------------
    dE = dt * (a*e_plus - b*e_minus),   dU = dt * (a*u_plus - b*u_minus)
    取 b = 1, r := a
    theta_E = e_minus/e_plus   (dE = 0 的临界比)
    theta_U = u_minus/u_plus   (dU = 0 的临界比)
    rho_exp = u_plus /e_plus   (扩张效率比)
    rho_res = u_minus/e_minus  (消解效率比)
    sign(dE) = sign(r - theta_E),  sign(dU) = sign(r - theta_U)

为什么用对数斜坡
------------------------------------------------------------------
线性斜坡会把两个 theta 的穿越点都推到末段，峰位估计不可靠。
对数斜坡 r(k) = exp(lerp(ln r_hi, ln r_lo, k/(K-1))) 使
    theta = 1     -> k* = (K-1)/2
    theta = 3     -> k* = (K-1)/4
    theta = 1/3   -> k* = 3(K-1)/4
即峰位解析值为整数（60 / 120 / 180），便于直接比对。

三个方法学陷阱（本轮实际踩到并已修正，保留记录）
------------------------------------------------------------------
T1. 不能用 Python 内置 hash() 做字符串散列——str 的 hash 按进程随机化
    （PYTHONHASHSEED），会破坏可复现性。改用 zlib.crc32。
T2. "预测 60 / 实测 59"是**并列取小索引**的伪像：r(k*)=theta 处 dU=0，
    于是 U(k*-1)=U(k*)。点值比对必须改为**平台比对**，且平台只在无噪声下有意义。
T3. 峰位估计不能用"平台下界"：噪声下 tol 放宽后下界会被早期点拉走
    （实测 42 vs 真值 60，是**系统性下偏**）。正确做法是对**平滑序列**取 argmax，
    它在噪声下均值稳定、方差随噪声增大。
T4. 单峰检出率的容差必须随平滑窗口缩小：平滑只把噪声降到 ~ sigma/sqrt(w)，
    用固定 1e-9 会导致"平滑后仍检出 0%"。容差取 3*sigma/sqrt(w)。

不能声称什么（写死在这里，防止报告里被写越界）
------------------------------------------------------------------
* 这是原型上的**机制自洽性检验**，不是实证验证。
* 系数与调度都是人造的，"峰位 = 60"不是普适常数，只有**峰序关系**可迁移。
* U 与 E 的上下界在本实验中**被刻意抬高**，以便无截断地测形状；
  收敛/截断问题是 R19 的辖域，不在本实验。
"""

import json
import math
import random
import statistics as st
import zlib

# ------------------------------------------------------------------ 参数

K = 241              # 步数（k = 0 .. K-1）
DT = 0.01            # 步长缩放
R_HI = 9.0
R_LO = 1.0 / 9.0
E_START = 5.0
U_START = 5.0
U_MAX = 1e9          # 刻意抬高：本实验测形状，不测截断
E_MAX = 1e9
E_MIN = -1e9
N_REP = 400
SEED = 20261019
OBS_SD_GRID = [0.0, 0.05, 0.15, 0.40, 0.80]
SMOOTH_W = 21        # 声明式分析选择（非拟合）：噪声条件下的峰位估计窗口
WINDOWS = [1, 5, 11, 21, 41]

LN_HI = math.log(R_HI)
LN_LO = math.log(R_LO)


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


CONFIGS = [
    make_config("A_resolution_efficient", 1.0, 1.0, 1.0, 3.0,
                "theta_E=1 < theta_U=3, rho_res=3 > rho_exp=1"),
    make_config("B_expansion_efficient", 1.0, 3.0, 1.0, 1.0 / 3.0,
                "theta_E=3 > theta_U=1/3, rho_res=1/9 < rho_exp=1"),
]


# ------------------------------------------------------------------ 工具


def zlib_off(s):
    return zlib.crc32(s.encode("utf-8")) % 9973


def r_ramp(k):
    t = k / (K - 1.0)
    return math.exp(LN_HI + (LN_LO - LN_HI) * t)


def r_const(value):
    return lambda k: value


def predicted_kstar(theta):
    if theta >= R_HI or theta <= R_LO:
        return None
    return (K - 1.0) * (math.log(theta) - LN_HI) / (LN_LO - LN_HI)


def dE_of(cfg, r):
    return DT * (r * cfg["e_plus"] - cfg["e_minus"])


def dU_of(cfg, r):
    return DT * (r * cfg["u_plus"] - cfg["u_minus"])


def moving_average(seq, w):
    if w <= 1:
        return list(seq)
    out = []
    n = len(seq)
    half = w // 2
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


def is_monotone_inc(seq, tol=0.0):
    return all(seq[i] >= seq[i - 1] - tol for i in range(1, len(seq)))


def is_monotone_dec(seq, tol=0.0):
    return all(seq[i] <= seq[i - 1] + tol for i in range(1, len(seq)))


def run_trajectory(rng, cfg, r_of_k, obs_sd=0.0, xi_sd=0.0):
    E, U = E_START, U_START
    Eobs, Uobs = [], []
    clipped = False
    for k in range(K):
        r = r_of_k(k)
        dE = dE_of(cfg, r) + (rng.gauss(0.0, xi_sd) if xi_sd > 0 else 0.0)
        dU = dU_of(cfg, r)
        E = E + dE
        U = U + dU
        if U < 0.0:
            U = 0.0
            clipped = True
        Eobs.append(E + (rng.gauss(0.0, obs_sd) if obs_sd > 0 else 0.0))
        Uobs.append(U + (rng.gauss(0.0, obs_sd) if obs_sd > 0 else 0.0))
    return {"E": Eobs, "U": Uobs, "clipped": clipped}


def tol_for(obs_sd, w):
    """平滑后残余噪声量级 ~ obs_sd/sqrt(w)；容差取 3 倍。见陷阱 T4。"""
    if obs_sd <= 0:
        return 1e-9
    return 3.0 * obs_sd / math.sqrt(max(1, w))


# ------------------------------------------------------------------ (1)(2) 峰位与单峰性


def peak_locus(cfg, obs_sd=0.0, smooth_w=SMOOTH_W, n_rep=N_REP, seed=SEED):
    out = {"config": cfg["name"], "obs_sd": obs_sd, "smooth_w": smooth_w,
           "pred_kstar_E": predicted_kstar(cfg["theta_E"]),
           "pred_kstar_U": predicted_kstar(cfg["theta_U"])}
    for tag, key in (("E", "E"), ("U", "U")):
        pred = out["pred_kstar_" + tag]
        rng = random.Random(seed + zlib_off(cfg["name"]) + int(round(obs_sd * 1000)))
        pts, unis = [], []
        plo, phi = [], []
        for _ in range(n_rep):
            tr = run_trajectory(rng, cfg, r_ramp, obs_sd=obs_sd)
            s = moving_average(tr[key], smooth_w)
            pts.append(argmax_of(s))
            ok, _ = is_unimodal(s, tol=tol_for(obs_sd, smooth_w))
            unis.append(ok)
            if obs_sd == 0.0:
                a, b = argmax_plateau(tr[key])
                plo.append(a)
                phi.append(b)
        out["peak_mean_" + tag] = round(st.mean(pts), 2)
        out["peak_sd_" + tag] = round(st.pstdev(pts), 2)
        out["unimodal_rate_" + tag] = round(sum(unis) / len(unis), 3)
        if obs_sd == 0.0:
            lo, hi = st.mean(plo), st.mean(phi)
            out["plateau_" + tag] = [round(lo, 2), round(hi, 2)]
            out["pred_in_plateau_" + tag] = bool(lo - 0.51 <= pred <= hi + 0.51)
        if pred is not None:
            out["abs_err_" + tag] = round(abs(st.mean(pts) - pred), 2)
    return out


# ------------------------------------------------------------------ (3) 峰序


def peak_order(cfg, obs_sd=0.15, smooth_w=SMOOTH_W, n_rep=N_REP, seed=SEED):
    pred_E = predicted_kstar(cfg["theta_E"])
    pred_U = predicted_kstar(cfg["theta_U"])
    pred_order = "U_first" if pred_U < pred_E else "E_first"
    rng = random.Random(seed + 31)
    compared = 0
    correct = 0
    for _ in range(n_rep):
        tr = run_trajectory(rng, cfg, r_ramp, obs_sd=obs_sd)
        kE = argmax_of(moving_average(tr["E"], smooth_w))
        kU = argmax_of(moving_average(tr["U"], smooth_w))
        obs_order = "U_first" if kU < kE else ("E_first" if kE < kU else "tie")
        if obs_order != "tie":
            compared += 1
            if obs_order == pred_order:
                correct += 1
    return {
        "config": cfg["name"], "obs_sd": obs_sd, "smooth_w": smooth_w,
        "theta_E": cfg["theta_E"], "theta_U": cfg["theta_U"],
        "rho_exp": cfg["rho_exp"], "rho_res": cfg["rho_res"],
        "rho_res_gt_rho_exp": cfg["rho_res"] > cfg["rho_exp"],
        "pred_order": pred_order,
        "obs_order_acc": round(correct / compared, 3) if compared else None,
        "n_compared": compared,
    }


# ------------------------------------------------------------------ (4) 中期相签名


def phase_sequence(cfg):
    seq = []
    for k in range(K):
        r = r_ramp(k)
        dE, dU = dE_of(cfg, r), dU_of(cfg, r)
        sE = 0 if abs(dE) < 1e-12 else (1 if dE > 0 else -1)
        sU = 0 if abs(dU) < 1e-12 else (1 if dU > 0 else -1)
        key = "({},{})".format({1: "+", 0: "0", -1: "-"}[sE],
                               {1: "+", 0: "0", -1: "-"}[sU])
        if not seq or seq[-1][0] != key:
            seq.append([key, k, k])
        else:
            seq[-1][2] = k
    return [{"combo": s[0], "k_from": s[1], "k_to": s[2]} for s in seq]


# ------------------------------------------------------------------ (5) 对照臂


def control_arms(cfg, obs_sd=0.0, smooth_w=SMOOTH_W, n_rep=N_REP, seed=SEED):
    """obs_sd=0 给出干净的机制结论；obs_sd>0 只用来展示测量噪声本身的影响。

    注意：原始（未平滑）轨迹上的单调性检出在 obs_sd>0 时会被噪声直接破坏——
    相邻差分的噪声 sd 约 sqrt(2)*obs_sd，可超过每步 dU 的量级。那是测量问题，
    不是机制问题，两种都跑以避免误读。"""
    tol = tol_for(obs_sd, smooth_w)
    out = {}
    for arm, rfn in (("no_resolution_r=R_HI", r_const(R_HI)),
                     ("no_expansion_r=R_LO", r_const(R_LO))):
        rng = random.Random(seed + 77)
        u_inc = u_dec = u_uni = 0
        slopes = []
        for _ in range(n_rep):
            tr = run_trajectory(rng, cfg, rfn, obs_sd=obs_sd)
            s = moving_average(tr["U"], smooth_w)
            if is_monotone_inc(s, tol=tol):
                u_inc += 1
            if is_monotone_dec(s, tol=tol):
                u_dec += 1
            ok, _ = is_unimodal(s, tol=tol)
            if ok:
                u_uni += 1
            # 不依赖容差的统计量：整条轨迹的平均每步增量。
            # 逐步单调判定在噪声下会退化（容差 > 每步漂移时"既增又减"），
            # 平均斜率的符号不受该问题影响，因此作为主要判据。
            slopes.append((s[-1] - s[0]) / (len(s) - 1))
        pos = sum(1 for x in slopes if x > 0)
        out[arm] = {"U_monotone_inc_rate": round(u_inc / n_rep, 3),
                    "U_monotone_dec_rate": round(u_dec / n_rep, 3),
                    "U_unimodal_rate": round(u_uni / n_rep, 3),
                    "U_mean_slope": round(st.mean(slopes), 6),
                    "U_mean_slope_sd": round(st.pstdev(slopes), 6),
                    "U_slope_positive_rate": round(pos / n_rep, 3)}
    return {"config": cfg["name"], "obs_sd": obs_sd, "smooth_w": smooth_w,
            "primary_statistic": "U_mean_slope (tolerance-free); "
                                 "stepwise rates are secondary and degrade under noise",
            "arms": out}


# ------------------------------------------------------------------ 平滑恢复（负面结果的量化）


def smoothing_recovery(cfg, obs_sd=0.15, n_rep=200, seed=SEED):
    rows = []
    for w in WINDOWS:
        rng = random.Random(seed + 913 + w)
        unis, pts = [], []
        for _ in range(n_rep):
            tr = run_trajectory(rng, cfg, r_ramp, obs_sd=obs_sd)
            s = moving_average(tr["U"], w)
            ok, _ = is_unimodal(s, tol=tol_for(obs_sd, w))
            unis.append(ok)
            pts.append(argmax_of(s))
        rows.append({"window": w,
                     "unimodal_rate": round(sum(unis) / len(unis), 3),
                     "peak_mean": round(st.mean(pts), 2),
                     "peak_sd": round(st.pstdev(pts), 2)})
    return {"config": cfg["name"], "obs_sd": obs_sd,
            "pred_kstar_U": predicted_kstar(cfg["theta_U"]), "rows": rows}


def noise_robustness(cfg, smooth_w=SMOOTH_W):
    return {"config": cfg["name"], "smooth_w": smooth_w,
            "rows": [peak_locus(cfg, obs_sd=sd, smooth_w=smooth_w)
                     for sd in OBS_SD_GRID]}


# ------------------------------------------------------------------ 主流程


def main():
    result = {
        "experiment": "R27",
        "claim": "X monotone increasing / M unimodal (peak then converge)",
        "params": {"K": K, "DT": DT, "R_HI": R_HI, "R_LO": R_LO,
                   "E_START": E_START, "U_START": U_START,
                   "N_REP": N_REP, "OBS_SD_GRID": OBS_SD_GRID,
                   "SMOOTH_W": SMOOTH_W, "WINDOWS": WINDOWS},
        "note": ("mechanism self-consistency check on a synthetic A_D prototype; "
                 "NOT empirical validation. Peak indices 60/120/180 follow from the "
                 "log-ramp choice and are not universal constants. "
                 "Smoothing window is a declared analysis choice, not a fitted one."),
        "configs": [],
    }
    for cfg in CONFIGS:
        entry = {
            "config": cfg["name"],
            "derived": {
                "theta_E": cfg["theta_E"], "theta_U": cfg["theta_U"],
                "rho_exp": cfg["rho_exp"], "rho_res": cfg["rho_res"],
                "rho_res_gt_rho_exp": cfg["rho_res"] > cfg["rho_exp"],
                "pred_kstar_E": predicted_kstar(cfg["theta_E"]),
                "pred_kstar_U": predicted_kstar(cfg["theta_U"]),
            },
            "peak_locus_noiseless": peak_locus(cfg, obs_sd=0.0, smooth_w=1),
            "peak_locus_noisy": peak_locus(cfg, obs_sd=0.15, smooth_w=SMOOTH_W),
            "peak_order": peak_order(cfg),
            "phase_sequence": phase_sequence(cfg),
            "control_arms_noiseless": control_arms(cfg, obs_sd=0.0),
            "control_arms_noisy": control_arms(cfg, obs_sd=0.15),
            "smoothing_recovery": smoothing_recovery(cfg),
            "noise_robustness": noise_robustness(cfg),
        }
        result["configs"].append(entry)

    with open("r27_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print("R27 ---- X monotone / M unimodal")
    print("K={} DT={} r: {:.4f} -> {:.4f} (log ramp)  N_REP={}".format(
        K, DT, R_HI, R_LO, N_REP))
    for e in result["configs"]:
        d = e["derived"]
        print("")
        print("== {} ==".format(e["config"]))
        print("   theta_E={:.3f} theta_U={:.3f} rho_exp={:.3f} rho_res={:.3f} "
              "rho_res>rho_exp={}".format(d["theta_E"], d["theta_U"], d["rho_exp"],
                                          d["rho_res"], d["rho_res_gt_rho_exp"]))
        pl = e["peak_locus_noiseless"]
        for tag in ("E", "U"):
            print("   {} : pred k*={:<7.1f} plateau={} in_plateau={} "
                  "unimodal_rate={}".format(
                      tag, pl["pred_kstar_" + tag], pl["plateau_" + tag],
                      pl["pred_in_plateau_" + tag], pl["unimodal_rate_" + tag]))
        po = e["peak_order"]
        print("   peak order  pred={:<8} obs_acc={} (n={})".format(
            po["pred_order"], po["obs_order_acc"], po["n_compared"]))
        print("   phase seq   " + " -> ".join(
            "{}[{}..{}]".format(p["combo"], p["k_from"], p["k_to"])
            for p in e["phase_sequence"]))
        print("   control arms (noiseless):")
        for arm, v in e["control_arms_noiseless"]["arms"].items():
            print("      {} : slope={:+.5f} (pos_rate={}) inc={} dec={} unimodal={}".format(
                arm, v["U_mean_slope"], v["U_slope_positive_rate"],
                v["U_monotone_inc_rate"], v["U_monotone_dec_rate"],
                v["U_unimodal_rate"]))
        print("   control arms (obs_sd=0.15, w={}; 主判据为平均斜率):".format(SMOOTH_W))
        for arm, v in e["control_arms_noisy"]["arms"].items():
            print("      {} : slope={:+.5f} (sd {:.5f}, pos_rate={})".format(
                arm, v["U_mean_slope"], v["U_mean_slope_sd"],
                v["U_slope_positive_rate"]))
        sr = e["smoothing_recovery"]
        print("   smoothing recovery (obs_sd={}, pred k*_U={:.0f}):".format(
            sr["obs_sd"], sr["pred_kstar_U"]))
        for row in sr["rows"]:
            print("      w={:<3} unimodal_rate={:<6} peak={:<7} sd={}".format(
                row["window"], row["unimodal_rate"], row["peak_mean"], row["peak_sd"]))
        print("   noise robustness (peak U, w={}):".format(SMOOTH_W))
        for row in e["noise_robustness"]["rows"]:
            print("      obs_sd={:<5} peak_U={:<7} sd={:<6} unimodal_U={}".format(
                row["obs_sd"], row["peak_mean_U"], row["peak_sd_U"],
                row["unimodal_rate_U"]))
    print("")
    print("wrote r27_results.json")


if __name__ == "__main__":
    main()
