# -*- coding: utf-8 -*-
"""
R19 —— U 与 Psi 的桥（B-B19）：双速率分解 A_D 的内部一致性检验
================================================================
对应 XS v3.1 补篇五（B-B19 的定位）。

背景
------------------------------------------------------------------
补篇二（动力学层）：E(k+1) - E(k) = Psi(k) = psi_bar + xi(k)
    E = 系统熵（展开度），psi_bar < 0 -> 回落/收敛，psi_bar > 0 -> 熵增/跃迁
补篇四（信息层）：  U = H(M | phi, C_in) >= 0
    U 降 -> 走向确定；U 不降 -> 未决

补篇四猜想 5.2：  sign(psi_bar) = -sign(dU/dk)          （未证）

本实验不"经验验证"该猜想（那需要真实理解者），而是检验它的**推导结构**：
在显式假设 A_D（展开/消解双速率线性分解）下，猜想 5.2 变成定理 5.2'，
其成立条件为 rho_res > rho_exp。本实验做的是：
  (1) 代数层：四种符号组合的可达性是否恰由 rho_res vs rho_exp 决定
  (2) 轨迹层：震荡（psi_bar = 0）期内 U 的行为是否由同一不等式决定
      —— 推论 5.3（酝酿 vs 反刍的判别）
  (3) 噪声鲁棒性：加入涨落 xi 后上述结论是否仍成立
  (4) 证否条件的判别力：故意违反 A_D 的线性性（加二次交互项），
      检查"两种反号组合并存"这一证否签名是否真的能出现。
      若不能出现，则证否条件是空的，必须重写。

A_D 的定义
------------------------------------------------------------------
每步的演化由两个子过程线性叠加：
    展开（expansion） 投入 a >= 0：  dE += a * e_plus,   dU += a * u_plus
    消解（resolution）投入 b >= 0：  dE -= b * e_minus,  dU -= b * u_minus
四个系数均为正。记 r := a / b，
    theta_E := e_minus / e_plus        （Psi = 0 的临界比）
    theta_U := u_minus / u_plus        （dU = 0 的临界比）
    rho_exp := u_plus  / e_plus        （展开的效率比）
    rho_res := u_minus / e_minus       （消解的效率比）

则
    sign(psi_bar) = sign(dE/dk) = sign(r - theta_E)
    sign(dU/dk)                 = sign(r - theta_U)

四种组合的可达性
------------------------------------------------------------------
    r > max(theta_E, theta_U)                -> (+, +)  越想越乱
    r < min(theta_E, theta_U)                -> (-, -)  正常收敛
    theta_E < r < theta_U  （仅当 rho_res > rho_exp） -> (+, -)  钻研型  ★猜想 5.2
    theta_U < r < theta_E  （仅当 rho_res < rho_exp） -> (-, +)  空洞型  ★猜想 5.2 的反面

关键：情形 (+,-) 与 (-,+) **互斥**，由 rho_res vs rho_exp 单一不等式决定。
    猜想 5.2 成立的区间存在  <=>  theta_E < theta_U  <=>  rho_res > rho_exp

证否条件
------------------------------------------------------------------
若同一组 (e_plus, e_minus, u_plus, u_minus) 下**同时**出现 (+,-) 与 (-,+)，
则 A_D 的线性叠加假设本身被证否（不只是猜想 5.2）。
本实验用 Config IV 验证该签名**确实可出现**（否则证否条件为空）。
"""

import json
import random
import statistics as st
import zlib

# ------------------------------------------------------------------ 配置


def make_config(name, e_plus, e_minus, u_plus, u_minus, quad=0.0, note=""):
    """quad != 0 时，dE 含二次项 -quad*a^2 —— 故意违反 A_D 的线性性。"""
    return {
        "name": name,
        "e_plus": e_plus,
        "e_minus": e_minus,
        "u_plus": u_plus,
        "u_minus": u_minus,
        "quad": quad,
        "note": note,
        "theta_E": e_minus / e_plus,
        "theta_U": u_minus / u_plus,
        "rho_exp": u_plus / e_plus,
        "rho_res": u_minus / e_minus,
    }


# Config I  : rho_res(0.90) > rho_exp(0.30)  -> 只应出现 (+,-)，不应出现 (-,+)
# Config II : rho_res(0.30) < rho_exp(0.90)  -> 只应出现 (-,+)，不应出现 (+,-)
# Config III: rho_res == rho_exp(0.50)       -> 两种反号组合都不可达（退化）
# Config IV : 同 Config I 的系数，但 dE 加二次项 -0.2*a^2
#             -> 违反 A_D 线性性，预期 (+,-) 与 (-,+) **并存**
CONFIGS = [
    make_config("I_rho_res>rho_exp", 1.0, 1.0, 0.30, 0.90,
                note="消解效率高于展开：应可出现钻研型 (+,-)，不可出现 (-,+)"),
    make_config("II_rho_res<rho_exp", 1.0, 1.0, 0.90, 0.30,
                note="展开效率高于消解：应可出现空洞型 (-,+)，不可出现 (+,-)"),
    make_config("III_degenerate", 1.0, 1.0, 0.50, 0.50,
                note="两效率相等：theta_E == theta_U，两种反号均不可达"),
    make_config("IV_nonlinear_A_D_violated", 1.0, 1.0, 0.30, 0.90, quad=0.20,
                note="违反 A_D：dE 含 -0.2*a^2，dE>0 的区间有两个边界，"
                     "预期 (+,-) 与 (-,+) 并存 —— 用于验证证否条件非空"),
]

# 扫描的 r = a/b 网格（对数均匀）
R_GRID = [round(10 ** (-1.0 + 2.0 * i / 40), 4) for i in range(41)]

E_MAX = 10.0      # E 的上界（临界跃迁边界）
E_MIN = 0.0       # E 的下界（回落/固化边界）
U_MAX = 4.0       # U 的上界（最大未决量）
U_START = 0.8 * U_MAX
B_STEPS = 400     # 迭代预算
N_REP = 400       # 每个 (config, r) 的重复次数
XI_SD = 0.15      # 涨落 xi 的标准差（补篇二的局部回落来源）


def zlib_off(s):
    return zlib.crc32(s.encode("utf-8")) % 9973


def dE_of(cfg, a, b):
    return a * cfg["e_plus"] - b * cfg["e_minus"] - cfg["quad"] * a * a


def dU_of(cfg, a, b):
    return a * cfg["u_plus"] - b * cfg["u_minus"]


def sign(x, tol=1e-9):
    if abs(x) < tol:
        return 0
    return 1 if x > 0 else -1


def combo_key(sE, sU):
    m = lambda s: "+" if s > 0 else ("-" if s < 0 else "0")
    return "({},{})".format(m(sE), m(sU))


# ------------------------------------------------------------------ 开环符号（代数层，无噪声）


def openloop_signs(cfg, r):
    return sign(dE_of(cfg, r, 1.0)), sign(dU_of(cfg, r, 1.0))


def predicted_reachable(cfg):
    """
    解析预言：本配置下哪些符号组合可达。
    quad != 0 时 dE>0 的区间可能有两个边界，退化为对网格直接判定。
    """
    if cfg["quad"] == 0.0:
        tE, tU = cfg["theta_E"], cfg["theta_U"]
        reach = set()
        if max(tE, tU) < R_GRID[-1]:
            reach.add("(+,+)")
        if min(tE, tU) > R_GRID[0]:
            reach.add("(-,-)")
        if tE < tU:
            reach.add("(+,-)")
        if tU < tE:
            reach.add("(-,+)")
        return sorted(reach)
    # 非线性：直接用开环式逐点判定
    seen = set()
    for r in R_GRID:
        sE, sU = openloop_signs(cfg, r)
        seen.add(combo_key(sE, sU))
    return sorted(seen)


# ------------------------------------------------------------------ 带边界的轨迹（含噪声）


def trajectory(rng, cfg, r, xi_sd=XI_SD, n_steps=B_STEPS):
    E = 0.5 * E_MAX
    U = U_START
    a, b = r, 1.0
    dE_log, dU_log = [], []
    hit_U_floor = False
    hit_U_ceil = False
    for _ in range(n_steps):
        dE = dE_of(cfg, a, b) + rng.gauss(0.0, xi_sd)
        dU = dU_of(cfg, a, b)
        E_new = min(E_MAX, max(E_MIN, E + dE))
        U_new = min(U_MAX, max(0.0, U + dU))
        dE_log.append(E_new - E)
        dU_log.append(U_new - U)
        E, U = E_new, U_new
        if U <= 1e-12:
            hit_U_floor = True
        if U >= U_MAX - 1e-12:
            hit_U_ceil = True
        if E >= E_MAX or E <= E_MIN:
            break
    return {
        "E_end": E, "U_end": U,
        "absorbed_up": E >= E_MAX, "absorbed_down": E <= E_MIN,
        "mean_dE": st.mean(dE_log) if dE_log else 0.0,
        "mean_dU": st.mean(dU_log) if dU_log else 0.0,
        "U_floor": hit_U_floor, "U_ceil": hit_U_ceil,
    }


def run_config(cfg, n_rep=N_REP, seed=20261019, xi_sd=XI_SD):
    rows = []
    for ri, r in enumerate(R_GRID):
        rng = random.Random(seed + zlib_off(cfg["name"]) + ri * 7919)
        trs = [trajectory(rng, cfg, r, xi_sd=xi_sd) for _ in range(n_rep)]
        mE = st.mean(t["mean_dE"] for t in trs)
        mU = st.mean(t["mean_dU"] for t in trs)
        rows.append({
            "r": r,
            "mean_dE": round(mE, 5),
            "mean_dU": round(mU, 5),
            "combo": combo_key(sign(mE), sign(mU)),
            "absorb_up%": round(100 * sum(t["absorbed_up"] for t in trs) / n_rep, 1),
            "absorb_down%": round(100 * sum(t["absorbed_down"] for t in trs) / n_rep, 1),
            "U_floor%": round(100 * sum(t["U_floor"] for t in trs) / n_rep, 1),
            "U_ceil%": round(100 * sum(t["U_ceil"] for t in trs) / n_rep, 1),
            "U_end均值": round(st.mean(t["U_end"] for t in trs), 4),
        })
    combos = sorted({r_["combo"] for r_ in rows})
    real = [c for c in combos if "0" not in c]
    return {
        "config": cfg["name"],
        "note": cfg["note"],
        "theta_E": round(cfg["theta_E"], 4),
        "theta_U": round(cfg["theta_U"], 4),
        "rho_exp": round(cfg["rho_exp"], 4),
        "rho_res": round(cfg["rho_res"], 4),
        "quad": cfg["quad"],
        "预测可达组合": predicted_reachable(cfg),
        "实测出现组合": real,
        "一致": predicted_reachable(cfg) == real,
        "扫描": rows,
    }


# ------------------------------------------------------------------ 推论 5.3：震荡期内的 U 行为


def psi_zero_points(cfg):
    """
    Psi = 0 的 r 值。线性（quad = 0）时唯一，为 theta_E = e_minus/e_plus；
    非线性（quad > 0）时 dE = e_plus*r - e_minus - quad*r^2 是二次的，
    可能有 0 / 1 / 2 个实根 —— 这正是 Config IV 能同时出现两种反号的原因。
    """
    quad, ep, em = cfg["quad"], cfg["e_plus"], cfg["e_minus"]
    if quad == 0.0:
        return [em / ep]
    disc = ep * ep - 4.0 * quad * em
    if disc < 0:
        return []
    s = disc ** 0.5
    return sorted([(ep - s) / (2.0 * quad), (ep + s) / (2.0 * quad)])


def oscillation_probe_at(cfg, r, n_rep=N_REP, seed=20261019, steps=60):
    """
    把 r 精确设为一个 Psi = 0 的点（震荡判据），观察 U 的走向。
    推论 5.3：
        rho_res > rho_exp  ->  U 下降  -> 震荡是「酝酿」（生产性）
        rho_res < rho_exp  ->  U 上升  -> 震荡是「反刍」（退化性）
        rho_res = rho_exp  ->  U 不变  -> 纯空转
    注意：U 有界（0 <= U <= U_MAX），撞界本身是结论的一部分——
    Config I 的 U 触底意味着「在震荡中完成了确定」。
    """
    rng = random.Random(seed + zlib_off(cfg["name"]) + 12345)
    dUs, netEs, floors, ceils = [], [], 0, 0
    for _ in range(n_rep):
        U = U_START
        E = 0.5 * E_MAX
        tot_dE = 0.0
        for _ in range(steps):
            dE = dE_of(cfg, r, 1.0) + rng.gauss(0.0, XI_SD)
            dU = dU_of(cfg, r, 1.0)
            E = min(E_MAX, max(E_MIN, E + dE))
            U = min(U_MAX, max(0.0, U + dU))
            tot_dE += dE
        dUs.append(U - U_START)
        netEs.append(tot_dE)
        if U <= 1e-12:
            floors += 1
        if U >= U_MAX - 1e-12:
            ceils += 1
    mU = st.mean(dUs)
    expect_down = cfg["rho_res"] > cfg["rho_exp"] + 1e-12
    expect_up = cfg["rho_res"] < cfg["rho_exp"] - 1e-12
    linear_ok = cfg["quad"] == 0.0
    if linear_ok:
        ok = (mU < -1e-6 and expect_down) or (mU > 1e-6 and expect_up) \
            or (abs(mU) <= 1e-6 and not expect_down and not expect_up)
    else:
        # 非线性配置下推论 5.3 **不适用**（其推导依赖 A_D 的线性性），
        # 此处只报告观测值，不做符合性判定——写 True 会是假的。
        ok = None
    return {
        "config": cfg["name"],
        "r_set_to_Psi0": round(r, 4),
        "rho_res/rho_exp": round(cfg["rho_res"] / cfg["rho_exp"], 4),
        "震荡期U净变化均值": round(mU, 5),
        "U方向": "下降(酝酿)" if mU < -1e-6 else ("上升(反刍)" if mU > 1e-6 else "不变(空转)"),
        "U触底%": round(100 * floors / n_rep, 1),
        "U触顶%": round(100 * ceils / n_rep, 1),
        "E净漂移均值": round(st.mean(netEs), 5),
        "符合推论5.3": ok if ok is None else bool(ok),
        "推论5.3适用性": "适用" if linear_ok else "不适用（A_D 线性性已被违反）",
    }


def oscillation_probe(cfg, **kw):
    """对配置的每个 Psi=0 点各跑一次探针。"""
    pts = psi_zero_points(cfg)
    return [oscillation_probe_at(cfg, r, **kw) for r in pts]


# ------------------------------------------------------------------ 噪声鲁棒性


def xi_robustness(seed=20261019):
    out = []
    for cfg in CONFIGS:
        row = {"config": cfg["name"]}
        for xi in (0.0, 0.05, 0.15, 0.40, 0.80):
            res = run_config(cfg, n_rep=200, seed=seed, xi_sd=xi)
            row["xi=" + str(xi)] = res["实测出现组合"]
        row["预测"] = predicted_reachable(cfg)
        row["各xi下均与预测一致"] = all(
            row["xi=" + str(xi)] == row["预测"]
            for xi in (0.0, 0.05, 0.15, 0.40, 0.80)
        )
        out.append(row)
    return out


# ------------------------------------------------------------------ 主流程


if __name__ == "__main__":
    scans = [run_config(c) for c in CONFIGS]
    out = {
        "说明": "A_D 双速率分解的内部一致性检验；非经验验证",
        "参数": {"E_MAX": E_MAX, "U_MAX": U_MAX, "U_START": U_START,
                 "B_STEPS": B_STEPS, "N_REP": N_REP, "XI_SD": XI_SD,
                 "R_GRID范围": [R_GRID[0], R_GRID[-1]]},
        "扫描结果": scans,
        "推论5.3震荡探针": [p_ for c in CONFIGS for p_ in oscillation_probe(c)],
        "噪声鲁棒性": xi_robustness(),
    }
    print(json.dumps(out, ensure_ascii=False, indent=2))
    with open("r19_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
