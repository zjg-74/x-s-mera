"""
X-S-MERA V3.4 · B3-ext：LE 内核四类终态机制扩展测绘
====================================================

B3（v1.6.45）Part A 发现当前 LE 内核（v3.4-le.1）是单调稳定系统——四类终态
（稳态/周期振荡/生态锁定/全局灭绝）在内核本身**不可达**，仅扩展生态原型可及，
故登记 **B3-ext 待办**。

本轮为 `src/le_coupling_kernel.py` 引入**可选扩展机制**（默认 None 保持向后
兼容的单调稳定内核），使四类终态在 **LE 内核本身**可达：

| 机制 | 键 | 原理 |
|---|---|---|
| 全局灭绝 | `death_threshold`/`death_rate` | life 低于死亡阈值 → 不可逆加速死亡（吸收态） |
| 生态锁定 | `allee_threshold`/`allee_strength` | Allee 双稳：bio 低于阈值趋 0、高于趋 k → 低态/高态锁定 |
| 周期振荡 | `oscillation_gain` + `env_r`/`consumer_mortality`/`half_override`/`uptake_override`/`conversion` | **纯消费者 R-M 结构**（bio 完全依赖资源摄食 + 线性死亡；资源 logistic 再生）→ Hopf 极限环 |
| 稳态 | （默认，无机制） | 原单调稳定内核 |

产出 results.json（四类各一代表用例 + 机制参数扫描摘要）。

诚实边界：机制为**操作化原型扩展**，非对现实的生态学断言；四类标签是对数值
轨迹的操作化分类；默认内核（无机制）行为保持 v3.4-le.1 不变（向后兼容）。
"""

import json
import os

import numpy as np

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
from le_coupling_kernel import LECouplingKernel  # noqa: E402

SEED = 20261029
T = 2000
DT = 0.1
WARM = 600


def run_trajectory(mechs, t=T, dt=DT, seed=SEED):
    """运行扩展内核，返回 bio（life 第 0 分量）时间序列。"""
    k = LECouplingKernel(mechanisms=mechs, seed=seed)
    for _ in range(int(t / dt)):
        k.step(k.t, dt)
    return [s[0] for s in k.history['L_state']]


def classify(bio, warm=WARM, tol=0.05, osc_thr=0.10):
    """操作化终态分类（与 b3 一致口径）。"""
    a = np.asarray(bio[warm:], dtype=float)
    if max(a[-100:]) < tol:
        return "extinction"
    tail_wave = max(a[-150:]) - min(a[-150:])
    if tail_wave > osc_thr:
        return "oscillation"
    mean_tail = float(np.mean(a[-150:]))
    if mean_tail < 0.5:
        return "lock"
    return "steady"


# 四类代表参数（固定种子下稳健）
REPRESENTATIVES = {
    "steady": {},
    "extinction": {"death_threshold": 0.90, "death_rate": 1.5},
    "lock": {"allee_threshold": 1.00, "allee_strength": 1.2},
    "oscillation": {"oscillation_gain": 1.0, "half_override": 0.05,
                    "uptake_override": 2.5, "env_r": 2.5,
                    "consumer_mortality": 0.5, "conversion": 0.8},
}


def sweep(label, variants):
    """对机制变体扫描，返回各变体终态分类。"""
    out = {}
    for name, mech in variants.items():
        bio = run_trajectory(mech)
        out[name] = {"mechanism": mech, "class": classify(bio),
                     "final_bio": float(bio[-1]),
                     "tail_wave": float(max(bio[WARM:]) - min(bio[WARM:]))}
    return out


def main() -> int:
    # 1) 四类代表用例
    reps = {}
    for name, mech in REPRESENTATIVES.items():
        bio = run_trajectory(mech)
        reps[name] = {"mechanism": mech, "class": classify(bio),
                      "final_bio": float(bio[-1]),
                      "tail_wave": float(max(bio[WARM:]) - min(bio[WARM:]))}

    # 2) 机制参数扫描摘要（各一维）
    extinction_sweep = sweep("extinction", {
        f"th{int(th*100):02d}": {"death_threshold": th, "death_rate": 1.5}
        for th in (0.60, 0.75, 0.90)})
    lock_sweep = sweep("lock", {
        f"th{int(th*100):02d}": {"allee_threshold": th, "allee_strength": 1.2}
        for th in (0.70, 0.90, 1.00, 1.20)})
    osc_sweep = sweep("oscillation", {
        f"conv{int(conv*10):02d}": {"oscillation_gain": 1.0, "half_override": 0.05,
                                    "uptake_override": 2.5, "env_r": 2.5,
                                    "consumer_mortality": 0.5, "conversion": conv}
        for conv in (0.5, 0.8, 1.2)})

    out = {
        "experiment": "B3-ext: LE-kernel four-attractor mechanism extension",
        "version": "v1.6.47",
        "kernel_version": "3.4-le.1+ext",
        "method": ("optional mechanisms on le_coupling_kernel; "
                   "death_threshold(extinction) / allee(lock) / "
                   "pure-consumer R-M(oscillation); default None = monotone stable"),
        "representatives": reps,
        "sweeps": {"extinction": extinction_sweep, "lock": lock_sweep,
                   "oscillation": osc_sweep},
        "honest_boundary": ("mechanisms are operationalized prototype extensions, "
                             "not claims about reality; default kernel (no mechanism) "
                             "behavior unchanged v3.4-le.1; four-class labels are "
                             "operational classification of numerical trajectories"),
        "seed": SEED,
    }

    path = os.path.join(os.path.dirname(__file__), "results.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("=== representatives ===")
    for name, r in reps.items():
        print(f"  {name:12s} -> {r['class']:11s} final={r['final_bio']:.3f}")
    print("=== sweeps ===")
    for sname, sw in (("extinction", extinction_sweep),
                      ("lock", lock_sweep), ("oscillation", osc_sweep)):
        for k, v in sw.items():
            print(f"  {sname}[{k}] -> {v['class']:11s} final={v['final_bio']:.3f}")
    print("OUT", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
