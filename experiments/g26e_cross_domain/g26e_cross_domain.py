#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""G-26e —— 非认知域对照：把 A13 的 D 系列判定程序移植到非流体/非认知域（流行病防控）

对应 `docs/xs-framework/26-v3.4-coupling-channel-from-mms-engineering.md` §7.4 G-26e 行：
    「单域（MMS）。跨域对照未做。【待验证】| A13 排期后加非流体力学域」
与 `27-v3.4-boundary-map.md` §5 表 B6 行（G-26e：非认知域对照）。

A13 已在**认知域**（换工作决策）演示了 D-DRILL / D-CONV / D-IDENT 判定程序。本项把同一套
判定程序**机械移植**到一个既非流体力学、也非认知过程的领域——**流行病防控**（公共卫生 /
系统动力学域），检验框架判读的**跨域可移植性**：
- E = 累计已干预量（确诊/隔离/处置累计）
- U = 未决风险存量（潜伏/未识别/未处理的风险，熵的演示代理）

**性质（必须诚实声明）**：
- 与 A13 同构：**演示性映射**，非现实疫情测量。真实疫情没有客观可测的 (E,U) 双序列。
- 这里把两条叙事（有效管控 / 防控失守）翻译成 A_D 轨迹，演示**同一判定程序**在非认知域
  上仍产出**互斥类别 + 可定位可辨识量**（external_necessity 跨域成立的一个演示证据）。
- 结论只支持「判定程序跨域可移植（演示层面）」，不构成对疫情机制或防控效果的断言。

复用：`experiments/a13_first_case/a13_first_case.py` 的 make_cfg / r_rampdown / trajectory /
d_drill / d_conv / d_ident（仅标准库，种子固定）。
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "a13_first_case"))
from a13_first_case import (  # noqa: E402
    make_cfg, trajectory, d_drill, d_conv, d_ident,
)

SEED = 20261024


def main():
    out = []

    def p(s=""):
        out.append(str(s))

    p("=" * 78)
    p("G-26e —— 非认知域对照：流行病防控（D 系列判定程序跨域移植）")
    p("=" * 78)
    p("性质声明：演示性映射，非现实疫情测量。真实疫情无客观可测 (E,U) 双序列；")
    p("此处把两条叙事翻译为 A_D 轨迹，演示同一判定程序在非认知域上仍产出互斥类别。")

    results = {
        "case": "流行病防控（非认知域对照，G-26e）",
        "domain": "public_health_epidemic（非流体、非认知域）",
        "nature": "演示性映射：两条叙事 → A_D 轨迹（非现实测量）；跨域可移植性演示",
        "seed": SEED,
        "method": "复用 A13 的 d_drill/d_conv/d_ident（跨域机械移植）",
        "interpretation": {
            "E": "累计已干预量（确诊/隔离/处置累计）",
            "U": "未决风险存量（潜伏/未识别/未处理风险，熵的演示代理）",
        },
        "scenarios": {},
    }

    # 场景 A：有效管控（钻研相 ρ=2）—— U0 小，后段收敛到接近 0
    cfgA = make_cfg(theta_E=1.0, rho=2.0)
    EsA, UsA, dEsA, dUsA = trajectory(cfgA, E0=5.0, U0=2.0)
    # 场景 B：防控失守（空洞相 ρ=0.6）—— U0 大，后段继续扩大
    cfgB = make_cfg(theta_E=1.0, rho=0.6)
    EsB, UsB, dEsB, dUsB = trajectory(cfgB, E0=5.0, U0=20.0)

    drillA = d_drill(cfgA, EsA, UsA, dEsA, dUsA)
    convA = d_conv(UsA)
    identA = d_ident(cfgA)
    p("")
    p("--- 场景 A：『严格管控，存量清零』（ρ=2，钻研相）---")
    p("  无框架人工判断：『管控有效，新增和存量都下来了，快清零了。』")
    p("  框架判定：")
    p("    D-DRILL: %s（峰序=%s, kU=%d kE=%d）" % (
        drillA["wedge"], drillA["peak_order"], drillA["kU"], drillA["kE"]))
    p("    D-CONV : %s（U_end=%.3f, U_peak=%.3f, |A_K|≈%.3f）" % (
        convA["verdict"], convA["u_end"], convA["umax"], convA["Ak_end"]))
    p("    D-IDENT: ρ_exp=%.3f 可自足读（A级）；sign(1−ρ)=%+d（B级）" % (
        identA["A_readable"]["rho_exp"], identA["B_readable"]["sign(1-rho)"]))
    results["scenarios"]["A_effective_control"] = {
        "narrative": "严格管控，存量清零", "rho": cfgA["rho"],
        "framework": {"D_DRILL": drillA, "D_CONV": convA, "D_IDENT": identA},
        "no_framework_human": "『管控有效，新增和存量都下来了，快清零了。』",
    }

    drillB = d_drill(cfgB, EsB, UsB, dEsB, dUsB)
    convB = d_conv(UsB)
    identB = d_ident(cfgB)
    p("")
    p("--- 场景 B：『防控失守，越防越多』（ρ=0.6，空洞相）---")
    p("  无框架人工判断：『越防越多，存量下不去，感觉失控了。』")
    p("  框架判定：")
    p("    D-DRILL: %s（峰序=%s, kU=%d kE=%d）" % (
        drillB["wedge"], drillB["peak_order"], drillB["kU"], drillB["kE"]))
    p("    D-CONV : %s（U_end=%.3f, U_peak=%.3f, |A_K|≈%.3f）" % (
        convB["verdict"], convB["u_end"], convB["umax"], convB["Ak_end"]))
    p("    D-IDENT: ρ_exp=%.3f 可自足读（A级）；sign(1−ρ)=%+d（B级）" % (
        identB["A_readable"]["rho_exp"], identB["B_readable"]["sign(1-rho)"]))
    results["scenarios"]["B_containment_lost"] = {
        "narrative": "防控失守，越防越多", "rho": cfgB["rho"],
        "framework": {"D_DRILL": drillB, "D_CONV": convB, "D_IDENT": identB},
        "no_framework_human": "『越防越多，存量下不去，感觉失控了。』",
    }

    p("")
    p("=" * 78)
    p("对照总结（跨域可移植性）：")
    p("  同一套判定程序（A13 的 d_drill/d_conv/d_ident）在**非认知、非流体**域（流行病防控）")
    p("  上仍产出：")
    p("    - 互斥类别（D-DRILL）：钻研相（有效管控）vs 空洞相（防控失守）")
    p("    - 收敛方向（D-CONV）：存量收敛 vs 存量发散")
    p("    - 可定位可辨识量（D-IDENT）：ρ_exp 自足（A）、sign(1−ρ)（B）、")
    p("      κ_E/ρ_res 需外部标定（C）、a/b/e/u 系数为伪目标（D）")
    p("  无框架只能给出『快清零 / 失控』的程度描述；框架给出互斥类别 + 可辨识量，")
    p("  且判定结构（A/B/C/D 分级、钻研/空洞楔形）不依赖具体领域——这是跨域可移植性的")
    p("  一个演示证据（演示层面，非经验断言）。")
    p("=" * 78)

    with open(os.path.join(HERE, "g26e_contrast.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main())
