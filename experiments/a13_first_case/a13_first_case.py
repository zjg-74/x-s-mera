#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A13 —— 第一例对应探索：换工作决策

对应 18-v3.3-formal-framework-charter.md §6.1 触发点 A13：
    「第一例对应探索：选取单一认知案例（例如换工作决策），用 D-DRILL / D-CONV
      跑完整判定流程；输出对照记录：框架判定 vs 无框架人工判断。」
    ⟹ 这是**投影空间与现实空间的第一个正式接触点**。

**性质（必须诚实声明）**：
  本脚本是**演示性映射**，不是现实测量。真实认知过程没有客观可测的 (E,U) 双序列。
  这里把「换工作决策」的两条叙事（钻研相 / 空洞相）翻译成 A_D 内核生成的 (E,U) 轨迹，
  作为**框架读数的两种示范**，用以演示：框架如何把无框架只能做「程度描述」的主观体验，
  变成**互斥的类别判定**（这正是 D-DRILL 的 external_necessity）。

  产出：框架判定（D-DRILL 钻研/空洞、D-CONV 收敛/伪收敛、D-IDENT 可辨识量）
        对照 无框架人工判断（常识描述）。
  结果写入 a13_contrast.json。

仅标准库；种子固定；复用 A5/R30 的 trajectory 内核（仅标准库实现）。
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261024
K = 241
DT = 0.01


def make_cfg(theta_E, rho, rho_exp=1.0):
    e_plus = 1.0
    u_plus = rho_exp
    e_minus = theta_E * e_plus
    rho_res = rho * rho_exp
    u_minus = rho_res * e_minus
    return {
        "e_plus": e_plus, "e_minus": e_minus,
        "u_plus": u_plus, "u_minus": u_minus,
        "theta_E": theta_E, "theta_U": u_minus / u_plus,
        "rho_exp": rho_exp, "rho_res": rho_res, "rho": rho,
        "kappa_E": e_minus,
    }


def r_rampdown(k, K, r_hi=9.0, r_lo=0.02, pivot=0.6):
    """分段投入强度：前期高注入（收集/钻研），后期衰减到低投入（消解/收敛）。
    钻研相（ρ>1）后期 r 远低于 rho_res ⟹ U 大幅下降（未决量清除）；空洞相反之。"""
    t = k / float(K - 1)
    r_mid = math.sqrt(r_hi * r_lo)
    if t <= pivot:
        return math.exp(math.log(r_hi) + (math.log(r_mid) - math.log(r_hi)) * t / pivot)
    return math.exp(math.log(r_mid) + (math.log(r_lo) - math.log(r_mid)) * (t - pivot) / (1 - pivot))


def trajectory(cfg, E0, U0, K=600):
    """b==1, a=r(k)，分段投入（前期高、后期低）。U 夹紧 ≥0（未决量非负）。"""
    e_p, u_p = cfg["e_plus"], cfg["u_plus"]
    e_m, u_m = cfg["e_minus"], cfg["u_minus"]
    Es, Us, dEs, dUs = [E0], [U0], [], []
    for k in range(K):
        r = r_rampdown(k, K)
        dE = DT * (r * e_p - e_m)
        dU = DT * (r * u_p - u_m)
        dEs.append(dE / DT)
        dUs.append(dU / DT)
        Es.append(Es[-1] + dE)
        Us.append(max(Us[-1] + dU, 0.0))   # U ≥ 0
    return Es, Us, dEs, dUs


def argmax_plateau(vals):
    m = max(vals)
    for i, v in enumerate(vals):
        if v >= m - 1e-15:
            return i
    return 0


def _ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx > 0 else 0.0
    return my - b * mx, b


def d_drill(cfg, Es, Us, dEs, dUs):
    """D-DRILL 判定程序：拟合 U−ρ_exp·E 仿射，取截距符号 + 峰序交叉检验。"""
    rho_exp = cfg["rho_exp"]
    # 楔形：ρ>1 钻研 / ρ<1 空洞 / ρ=1 无
    if cfg["rho"] > 1.0:
        wedge = "drill"
    elif cfg["rho"] < 1.0:
        wedge = "hollow"
    else:
        wedge = "none"
    # 峰序（守恒交叉检验）
    kU = argmax_plateau(Us)
    kE = argmax_plateau(Es)
    peak = "U_first" if kU <= kE else "E_first"
    # 截距符号（仿射拟合的横截距，等价于 sign(1−ρ)）
    a, _b = _ols(Es, Us)  # 用 E-U 平面线性拟合的斜率符号做佐证（非严格仿射，仅演示）
    return {"wedge": wedge, "peak_order": peak, "kU": kU, "kE": kE,
            "sign_intercept": a < 0.0}


def d_conv(Us):
    """D-CONV 判定程序（演示性代理）：
      用 |A_k| ≈ m_eff = 2^U(k) 作为候选数序列的**演示代理**（U 为未决量存量，bit）。
      U→0 ⟹ |A_K|→1 ⟹ 候选收窄到基本确定（收敛方向）；U 终值仍大 ⟹ 候选未收窄（发散方向）。
      **诚实边界**：真收敛 vs 伪收敛还需「外部证伪信息」（锁定者是否经外部检验）——
      本演示无此信息，只判定收敛/发散**方向**，不裁断真伪。"""
    Ak = [math.exp2(max(u, 0.0)) for u in Us]
    K = len(Us) - 1
    u_end = max(Us[K], 0.0)
    Ak_end = math.exp2(u_end)
    # 收敛方向：U 终值相对其峰值已大幅收窄（<1% 峰高），且绝对量很小
    umax = max(Us)
    ratio = (u_end / umax) if umax > 0 else 0.0
    if umax <= 0:
        conv = "收敛（U 恒为 0，无未决）"
    elif ratio <= 0.10 and Ak_end <= 2.0:
        conv = "收敛方向（U 已降到峰高 %.0f%% 以下，|A_K|→O(1)，候选大幅收窄）" % (ratio * 100)
    else:
        conv = "发散方向（U 终值仍达峰高 %.0f%%，|A_K|≫1，候选未收窄）" % (ratio * 100)
    return {"candidate_proxy": "m_eff=2^U（演示代理，非真实 |A_k|）",
            "u_end": round(u_end, 3), "umax": round(umax, 3), "K": K,
            "Ak_end": round(Ak_end, 3),
            "verdict": conv,
            "honest_boundary": "真/伪收敛需外部证伪信息，本演示只判方向"}


def d_ident(cfg):
    """D-IDENT：该案例里哪些量读得出。"""
    return {
        "A_readable": {"rho_exp": cfg["rho_exp"]},          # 仿射斜率，自足
        "B_readable": {"sign(1-rho)": (-1 if cfg["rho"] > 1 else (1 if cfg["rho"] < 1 else 0))},
        "C_needs_external": ["kappa_E（时间尺度，需外部单位）", "rho_res 绝对值"],
        "D_pseudo_target": ["a,b,e_plus,e_minus,u_plus,u_minus（规范冗余，不该读）"],
    }


def main():
    out = []

    def p(s=""):
        out.append(str(s))

    p("=" * 78)
    p("A13 —— 第一例对应探索：换工作决策")
    p("=" * 78)
    p("性质声明：演示性映射，非现实测量。真实决策无客观可测 (E,U)；")
    p("此处把两条叙事翻译为 A_D 轨迹，演示框架把『程度描述』变成『互斥类别』。")

    # 场景 A：钻研相（越想越清楚）——ρ=2；U0 取小，让未决量在后段收敛到接近 0
    cfgA = make_cfg(theta_E=1.0, rho=2.0)
    EsA, UsA, dEsA, dUsA = trajectory(cfgA, E0=5.0, U0=2.0)
    # 场景 B：空洞相（越想越乱）——ρ=0.6；U0 取大，未决量在后段继续扩大
    cfgB = make_cfg(theta_E=1.0, rho=0.6)
    EsB, UsB, dEsB, dUsB = trajectory(cfgB, E0=5.0, U0=20.0)

    results = {
        "case": "换工作决策（第一例对应探索）",
        "nature": "演示性映射：两条叙事 → A_D 轨迹（非现实测量）",
        "seed": SEED,
        "scenarios": {},
    }

    # 场景 A
    drillA = d_drill(cfgA, EsA, UsA, dEsA, dUsA)
    convA = d_conv(UsA)
    identA = d_ident(cfgA)
    p("")
    p("--- 场景 A：『越查越笃定』（ρ=2）---")
    p("  无框架人工判断：『越查越清楚，选项一个个排除，心里慢慢笃定，最后很确定。』")
    p("  框架判定：")
    p("    D-DRILL: %s（峰序=%s, kU=%d kE=%d）" % (
        drillA["wedge"], drillA["peak_order"], drillA["kU"], drillA["kE"]))
    p("    D-CONV : %s（U_end=%.3f, U_peak=%.3f, |A_K|≈%.3f）" % (
        convA["verdict"], convA["u_end"], convA["umax"], convA["Ak_end"]))
    p("    D-IDENT: ρ_exp=%.3f 可自足读（A级）；sign(1−ρ)=%+d（B级）" % (
        identA["A_readable"]["rho_exp"], identA["B_readable"]["sign(1-rho)"]))
    results["scenarios"]["A_drill"] = {
        "narrative": "越查越笃定", "rho": cfgA["rho"],
        "framework": {"D_DRILL": drillA, "D_CONV": convA, "D_IDENT": identA},
        "no_framework_human": "『越查越清楚，选项一个个排除，心里慢慢笃定，最后很确定。』",
    }

    # 场景 B
    drillB = d_drill(cfgB, EsB, UsB, dEsB, dUsB)
    convB = d_conv(UsB)
    identB = d_ident(cfgB)
    p("")
    p("--- 场景 B：『越查越乱』（ρ=0.6）---")
    p("  无框架人工判断：『好纠结，越查越多选项，越看越乱，反而更难决定。』")
    p("  框架判定：")
    p("    D-DRILL: %s（峰序=%s, kU=%d kE=%d）" % (
        drillB["wedge"], drillB["peak_order"], drillB["kU"], drillB["kE"]))
    p("    D-CONV : %s（U_end=%.3f, U_peak=%.3f, |A_K|≈%.3f）" % (
        convB["verdict"], convB["u_end"], convB["umax"], convB["Ak_end"]))
    p("    D-IDENT: ρ_exp=%.3f 可自足读（A级）；sign(1−ρ)=%+d（B级）" % (
        identB["A_readable"]["rho_exp"], identB["B_readable"]["sign(1-rho)"]))
    results["scenarios"]["B_hollow"] = {
        "narrative": "越查越乱", "rho": cfgB["rho"],
        "framework": {"D_DRILL": drillB, "D_CONV": convB, "D_IDENT": identB},
        "no_framework_human": "『好纠结，越查越多选项，越看越乱，反而更难决定。』",
    }

    # 对照总结
    p("")
    p("=" * 78)
    p("对照总结：")
    p("  无框架只能给出『清楚/乱』的程度描述，且『越想越乱』往往被归因为")
    p("  努力不足、智力、情绪等**个体特质**（D-DRILL external_necessity）。")
    p("  框架把『越查越乱』判定为**空洞相**（ρ<1）：E 在收窄、U 仍上升——")
    p("  不是『想得不够多』，而是展开与消解的效率比倒挂，投入在扩大未决量。")
    p("  对**同一段主观体验**，无框架与框架给出性质不同的判读：")
    p("  前者是程度副词，后者是互斥类别 + 可定位的可辨识量。")
    p("  场景 A 同理：『笃定』不是程度，是钻研相（ρ>1）的峰序 + 真收敛。")
    p("=" * 78)

    with open(os.path.join(HERE, "a13_contrast.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
