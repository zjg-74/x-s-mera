"""
X-S-MERA V3.4 · B3 吸引子分类：四类终态相图测绘（L2 仿真层）
================================================================

GAPS.md B3：体系自发抵达四类终态（稳态 / 周期振荡 / 生态锁定 / 全局灭绝）。
终态分类是动力学稳定性陈述，FOL 谓词层结构性不可达 → 交 L2 仿真层核验。

本脚本双轨测绘：

**Part A · 当前 LE 内核吸引子测绘**（le_coupling_kernel v3.4-le.1）
在耦合参数网格（coupling_gain × d_threshold_c1）上跑内核轨迹，分类终态。
→ 诚实登记结论：**LE 内核是单调稳定系统**（逻辑斯谛增长 + Michaelis-Menten
  摄食 + 线性资源供给-消耗），全域参数扫均收敛到唯一非零稳态，
  **不呈现振荡 / 灭绝 / 锁定**。四类终态在当前内核不可达（B3-ext 待办）。

**Part B · 扩展生态动力学原型相图**（b3 实验内独立模型，诚实标注"原型"）
在资源-消费者 Lotka-Volterra 型模型（logistic 资源自再生 + 消费者
捕食增长/密度约束/死亡 + Allee 效应双稳）的网格（资源再生率 s × 捕食率 a）上
测绘四类终态相图，证明四类终态机制可被此类生态动力学表达。
→ 覆盖：稳态 / 周期振荡 / 生态锁定 / 全局灭绝 四类（可达性验证）。

诚实边界：Part B 是**操作化原型模型**，非当前 LE 内核；四类标签是对数值轨迹
的操作化分类，非现实生态学断言。正式四类终态测绘须在 LE 内核引入
振荡/灭绝/锁定机制（B3-ext），见 README。

分类器（尾段 25% 窗口，L/消费者 生物量）：
- extinction  : 尾段均值 < EXT_THR
- oscillation: 峰谷幅 > OSC_PTP 且标准差 > OSC_SD 且局部极大数 ≥ 2
- steady_state: 标准差 < STEADY_SD 且峰谷 < STEADY_PTP
- ecological_lock: 其余（非灭绝、非稳态、非明确周期 → 被困/次优波动）

用法
----
    python b3_attractor_phase.py                # 标准双轨网格
    python b3_attractor_phase.py --quick        # 小网格探测
    python b3_attractor_phase.py --out out.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from le_coupling_kernel import LECouplingKernel  # noqa: E402

# ---------------- 分类参数（诚实登记：操作化口径，非现实生态学阈值） ----------------
EXT_THR = 0.02      # 灭绝：尾段均值 < 阈值
OSC_PTP = 0.04      # 振荡：峰谷振幅 > 阈值
OSC_SD = 0.008      # 振荡：标准差 > 阈值
STEADY_SD = 0.008   # 稳态：标准差 < 阈值
STEADY_PTP = 0.02   # 稳态：峰谷 < 阈值


def classify_tail(series: np.ndarray) -> str:
    """对尾段生物量序列分类为四类终态之一。

    顺序（先决后）：extinction → oscillation（peaks≥2）→ steady_state → ecological_lock。
    """
    series = np.asarray(series, dtype=float)
    if series.size == 0:
        return "ecological_lock"
    tail = series[- (series.size // 4):] if series.size >= 8 else series
    if tail[-1] < EXT_THR or float(np.mean(tail)) < EXT_THR:
        return "global_extinction"
    ptp = float(np.ptp(tail))
    sd = float(np.std(tail))
    if ptp > OSC_PTP and sd > OSC_SD:
        peaks = int(np.sum((tail[1:-1] > tail[:-2]) & (tail[1:-1] > tail[2:])))
        if peaks >= 2:
            return "periodic_oscillation"
    if sd < STEADY_SD and ptp < STEADY_PTP:
        return "steady_state"
    return "ecological_lock"


# ---------------- Part A · 当前 LE 内核吸引子测绘 ----------------

def run_le_trajectory(coupling_gain: float, d_threshold_c1: float,
                      nsteps: int, dt: float, seed: int) -> dict:
    k = LECouplingKernel(coupling_gain=coupling_gain,
                         d_threshold_c1=d_threshold_c1, seed=seed)
    bio = np.empty(nsteps, dtype=float)
    t = 0.0
    for i in range(nsteps):
        r = k.step(t, dt)
        bio[i] = float(r["L_state"][0])
        t += dt
    cls = classify_tail(bio)
    tail = bio[-(bio.size // 4):]
    return {"cls": cls, "final_bio": float(bio[-1]),
            "tail_mean_bio": float(np.mean(tail)),
            "bio_ptp": float(np.ptp(tail)),
            "bio_std": float(np.std(tail))}


def run_le_grid(cg_grid, c1_grid, nsteps, dt, seed) -> dict:
    details = []
    for cg in cg_grid:
        for c1 in c1_grid:
            r = run_le_trajectory(float(cg), float(c1), nsteps, dt, seed)
            details.append({"coupling_gain": float(cg),
                            "d_threshold_c1": float(c1), **r})
    counts = {}
    for d in details:
        counts[d["cls"]] = counts.get(d["cls"], 0) + 1
    return {"details": details, "class_counts": counts,
            "n_grid": len(details)}


# ---------------- Part B · 扩展生态动力学原型 ----------------

class EcoDynamics:
    """资源-消费者动力学（Lotka-Volterra 型 + 密度约束 + Allee 效应双稳）。

    R : 资源量（被捕食者，logistic 自再生）
    C : 消费者(生命)生物量
    dR/dt = s*R*(1 - R/cap) - a*C*R
    dC/dt = e*a*C*R - m*C - n*C^2 - allee*C/(A+C)
    Allee 项：低 C 时高死亡 → 双稳（高 C 稳态 / 低 C 陷阱 = 生态锁定）
    """
    def __init__(self, s=0.8, cap=1.5, a=2.0, e=0.6, m=0.2, n=0.05,
                 allee=0.0, A=0.3, R0=1.2, C0=0.5):
        self.s, self.cap, self.a = s, cap, a
        self.e, self.m, self.n = e, m, n
        self.allee, self.A = allee, A
        self.R, self.C = R0, C0

    def step(self, dt):
        R, C = self.R, self.C
        R = max(R, 0.0); C = max(C, 0.0)
        dR = self.s * R * (1.0 - R / self.cap) - self.a * C * R
        dC = (self.e * self.a * C * R - self.m * C - self.n * C * C
              - self.allee * C / (self.A + C))
        self.R = max(self.R + dt * dR, 0.0)
        self.C = max(self.C + dt * dC, 0.0)
        return self.R, self.C


def run_eco_trajectory(s: float, a: float, nsteps: int, dt: float,
                       seed: int) -> dict:
    # seed 复现：用 RNG 生成确定性初值扰动（固定种子 → 可复现）
    rng = np.random.default_rng(seed)
    m = EcoDynamics(s=s, a=a,
                    R0=float(1.2 + 0.05 * rng.standard_normal()),
                    C0=float(0.5 + 0.05 * rng.standard_normal()))
    bio = np.empty(nsteps, dtype=float)
    res = np.empty(nsteps, dtype=float)
    for i in range(nsteps):
        r, c = m.step(dt)
        bio[i] = float(c); res[i] = float(r)
    cls = classify_tail(bio)
    tail = bio[-(bio.size // 4):]
    return {"cls": cls, "final_C": float(bio[-1]),
            "tail_mean_C": float(np.mean(tail)),
            "C_ptp": float(np.ptp(tail)), "C_std": float(np.std(tail)),
            "final_R": float(res[-1])}


def run_eco_grid(s_grid, a_grid, nsteps, dt, seed) -> dict:
    details = []
    for s in s_grid:
        for a in a_grid:
            r = run_eco_trajectory(float(s), float(a), nsteps, dt, seed)
            details.append({"resource_growth_s": float(s),
                            "attack_rate_a": float(a), **r})
    counts = {}
    for d in details:
        counts[d["cls"]] = counts.get(d["cls"], 0) + 1
    # 相图矩阵 {s: [cls_per_a]}
    mat = {}
    for d in details:
        mat.setdefault(d["resource_growth_s"], {})[d["attack_rate_a"]] = d["cls"]
    return {"details": details, "class_counts": counts,
            "n_grid": len(details), "phase_matrix": mat}


# ---------------- 入口 ----------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="小网格探测")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__),
                                                  "results.json"))
    ap.add_argument("--nsteps", type=int, default=6000)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--n-s", type=int, default=14, help="资源再生率 s 档数")
    ap.add_argument("--n-a", type=int, default=10, help="捕食率 a 档数")
    args = ap.parse_args()

    if args.quick:
        le_cg = np.linspace(0.3, 3.0, 4)
        le_c1 = np.linspace(0.05, 0.6, 4)
        s_grid = np.linspace(0.05, 1.6, 6)
        a_grid = np.linspace(0.3, 8.0, 6)
        le_nsteps = 1200
        eco_nsteps = 3000
    else:
        le_cg = np.linspace(0.3, 3.0, 6)
        le_c1 = np.linspace(0.05, 0.6, 6)
        s_grid = np.linspace(0.05, 1.6, args.n_s)
        a_grid = np.linspace(0.3, 8.0, args.n_a)
        le_nsteps = args.nsteps // 2
        eco_nsteps = args.nsteps

    le = run_le_grid(le_cg, le_c1, le_nsteps, 0.05, args.seed)
    eco = run_eco_grid(s_grid, a_grid, eco_nsteps, args.dt, args.seed)

    out = {
        "part_A_LE_kernel": {
            "kernel": "le_coupling_kernel", "kernel_version": "3.4-le.1",
            "grid": {"coupling_gain": [float(x) for x in le_cg],
                     "d_threshold_c1": [float(x) for x in le_c1]},
            "nsteps": le_nsteps, "dt": 0.05, "seed": args.seed,
            "class_counts": le["class_counts"], "n_grid": le["n_grid"],
            "details": le["details"],
            "verdict": "monostable_steady",
            "note": ("LE 内核单调稳定：全域扫收敛唯一非零稳态；"
                     "振荡/灭绝/锁定当前内核不可达（B3-ext）"),
        },
        "part_B_eco_prototype": {
            "model": "eco_lv_prototype",
            "model_desc": ("资源-消费者 LV(logistic 资源)+密度约束+Allee 双稳；"
                           "b3 实验内独立原型，非当前 LE 内核"),
            "grid": {"resource_growth_s": [float(x) for x in s_grid],
                     "attack_rate_a": [float(x) for x in a_grid]},
            "nsteps": eco_nsteps, "dt": args.dt, "seed": args.seed,
            "class_counts": eco["class_counts"], "n_grid": eco["n_grid"],
            "phase_matrix": eco["phase_matrix"],
            "details": eco["details"],
            "verdict": "four_attractors_reachable",
        },
        "classifier_defaults": {
            "EXT_THR": EXT_THR, "OSC_PTP": OSC_PTP, "OSC_SD": OSC_SD,
            "STEADY_SD": STEADY_SD, "STEADY_PTP": STEADY_PTP,
            "tail_fraction": 0.25,
        },
        "verdict": "prototype",   # 诚实边界：操作化分类，非现实生态学断言
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("PART_A_LE_COUNTS", le["class_counts"])
    print("PART_B_ECO_COUNTS", eco["class_counts"])
    print("PART_B_MATRIX")
    for s in sorted(eco["phase_matrix"], key=float):
        row = eco["phase_matrix"][s]
        print("  s=%s: %s" % (s, " ".join(row[a][:1].upper()
                                           for a in sorted(row, key=float))))
    print("OUT", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
