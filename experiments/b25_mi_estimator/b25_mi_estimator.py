"""
X-S-MERA V3.4 · B-B25 估计器演示
==================================

两场景验证 I(M;Ω|φ) 估计器：

- 场景 A（有偏倚）：两理解者 Ω=0/1 对候选 M 的先验分布不同 →
  混合估计 H(M|φ) 高于分群体 H(M|φ,Ω) → I(M;Ω|φ) > 0（偏倚可检出）。
- 场景 B（无偏倚）：Ω 与 M 分布无关 → I(M;Ω|φ) ≈ 0（跨群体可比）。

产出 results.json。合成标签数据，非真实认知数据（诚实边界）。
"""

import json
import os

import numpy as np

from mi_estimator import MIEstimator

SEED = 20261029


def gen_biased(n_per=2000, seed=SEED):
    """场景 A：Ω 影响 M 分布（有偏倚）。"""
    rng = np.random.default_rng(seed)
    M = np.concatenate([
        rng.choice([0, 1, 2, 3], n_per, p=[0.55, 0.25, 0.15, 0.05]),  # Ω=0
        rng.choice([0, 1, 2, 3], n_per, p=[0.05, 0.15, 0.30, 0.50]),  # Ω=1
    ])
    Omega = np.concatenate([np.zeros(n_per, int), np.ones(n_per, int)])
    phi = np.zeros(M.size, int)   # 单一语境
    return M, Omega, phi


def gen_unbiased(n_per=2000, seed=SEED):
    """场景 B：Ω 与 M 无关（无偏倚）。"""
    rng = np.random.default_rng(seed)
    p = [0.30, 0.30, 0.20, 0.20]
    M = rng.choice([0, 1, 2, 3], n_per * 2, p=p)
    Omega = np.concatenate([np.zeros(n_per, int), np.ones(n_per, int)])
    phi = np.zeros(M.size, int)
    return M, Omega, phi


def gen_phi_effect(n_per=1500, seed=SEED):
    """场景 C：语境 φ 影响 M，但 Ω 不（偏倚按 φ 校正后为 0）。"""
    rng = np.random.default_rng(seed)
    M, Omega, phi = [], [], []
    for ph in (0, 1):
        for om in (0, 1):
            p = ([0.6, 0.2, 0.1, 0.1] if ph == 0 else
                 [0.1, 0.1, 0.3, 0.5])
            M.append(rng.choice([0, 1, 2, 3], n_per, p=p))
            Omega.append(np.full(n_per, om, int))
            phi.append(np.full(n_per, ph, int))
    return (np.concatenate(M), np.concatenate(Omega), np.concatenate(phi))


def main() -> int:
    scenarios = {}
    for name, gen in [("A_biased_by_Omega", gen_biased),
                      ("B_Omega_independent", gen_unbiased),
                      ("C_phi_effect_Omega_irrelevant", gen_phi_effect)]:
        M, Omega, phi = gen()
        est = MIEstimator(M, Omega, phi)
        scenarios[name] = {"n": int(M.size),
                           "estimate": est.report()}

    out = {
        "experiment": "B-B25 mutual-information estimator I(M;Omega|phi)",
        "version": "v1.6.46",
        "method": ("conditional mutual information = H(M|phi) - H(M|phi,Omega);"
                   "Miller-Madow corrected empirical entropy"),
        "scenarios": scenarios,
        "expected": {
            "A_biased_by_Omega": "I > 0 (cross-understander bias detectable)",
            "B_Omega_independent": "I ~ 0 (comparable across groups)",
            "C_phi_effect_Omega_irrelevant": "I ~ 0 (bias removed once phi conditioned)",
        },
        "verdict": "prototype",   # 诚实边界：合成标签数据，非真实认知
    }

    out_path = os.path.join(os.path.dirname(__file__), "results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    for name, s in scenarios.items():
        e = s["estimate"]
        print(f"{name}: n={s['n']} H_M_phi={e['H_M_phi']:.4f} "
              f"H_M_phi_Omega={e['H_M_phi_Omega']:.4f} "
              f"I(M;Omega|phi)={e['I_M_Omega_phi']:+.4f}")
    print("OUT", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
