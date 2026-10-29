# -*- coding: utf-8 -*-
"""R2：F1 ↔ 发散（¬C1）/ F2 ↔ 伪收敛（C1∧¬C2）的条件分离（主文 §9 / 推论 R2）。

规格来源：`01-v3.1-main-dichotomy-macro-loading-trichotomy.md` §7.2 + §9 R2——
  F1 ↔ 发散：¬C1 时应观察到"发散现象：A_k = ∅ 且不可新构"；
  F2 ↔ 伪收敛：C1 ∧ ¬C2 时应观察到 ghost 伪力（加载成功但事后验证为错）。
  此二者的条件分离在工程原型上应当是可分离的。

本文在"宏-槽-事实"离散加载原型上验证该分离（原型确认，非现实证据）：
  - C1（事实闭合可构）：存在候选宏其加载度 σ ≥ 可加载阈值 σ_load。
  - C2（正确性）：收敛到的宏 == 预设正确宏（ground truth）。
  - 三态判定：
      C1 ∧  C2  → 真收敛
      ¬C1       → F1 发散（A_k = ∅，无候选可构）
      C1 ∧ ¬C2  → F2 伪收敛（A_k ≠ ∅ 且收敛，但收敛到错误宏 = ghost 伪力）
  - 分离性判据：F1 输出 A_k = ∅（发散），F2 输出 A_k ≠ ∅ 但收敛错误——
    二者在"是否已加载"这一可观察量上**条件分离**。

仅标准库，固定种子。
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from typing import Dict, List

SEED = 20261029
SIGMA_LOAD = 0.5       # 可加载阈值（C1 判据）
SIGMA_STAR = 0.8       # 完成加载阈值


@dataclass
class Macro:
    mid: str
    slots: List[Dict]
    is_correct: bool = False

    def sigma(self, fact_avail: Dict[str, int], w_min: int = 1) -> float:
        total = sum(s["w"] for s in self.slots)
        if total == 0:
            return 0.0
        ok = 0.0
        for s in self.slots:
            if all(fact_avail.get(f, 0) >= w_min for f in s["feats"]):
                ok += s["w"]
        return ok / total


def make_candidates(rng: random.Random) -> List[Macro]:
    """正确宏 + 错误宏 + 不可加载宏（槽期望无关特征）。"""
    pool = [f"f{i}" for i in range(40)]
    correct = Macro("M_correct", [{"feats": ["a", "b"], "w": 2},
                                  {"feats": ["c"], "w": 1}], is_correct=True)
    wrong = Macro("M_wrong", [{"feats": ["x", "y"], "w": 2},
                              {"feats": ["z"], "w": 1}], is_correct=False)
    # 不可加载宏：期望一组极稀有特征，事实批不会足额供给
    rare = rng.sample(pool, 4)
    unloadable = Macro("M_rare", [{"feats": rare, "w": 3}], is_correct=False)
    return [correct, wrong, unloadable]


def facts_for_true_conv():
    """C1∧C2：足额供给正确宏所需事实（a,b,c）。"""
    return {"a": 3, "b": 3, "c": 3, "x": 0, "y": 0, "z": 0}


def facts_for_f1_divergence():
    """¬C1：所有候选 σ < σ_load（0.5）。

    正确宏槽1 需 a∧b（w=2）、槽2 需 c（w=1）：b 缺 ⟹ σ=0；
    错误宏槽1 需 x∧y（w=2）、槽2 需 z（w=1）：全缺 ⟹ σ=0；
    稀有宏需 rare（w=3）：缺 ⟹ σ=0。⟹ A_k = ∅。
    """
    return {"a": 3, "b": 0, "c": 0, "x": 0, "y": 0, "z": 0}


def facts_for_f2_pseudoconvergence():
    """C1∧¬C2：错误宏足额（x,y,z），正确宏缺 c ⟹ 收敛到错误宏。"""
    return {"a": 3, "b": 3, "c": 0, "x": 3, "y": 3, "z": 3}


def classify(candidates: List[Macro], fact_avail: Dict[str, int],
             sigma_load: float = SIGMA_LOAD,
             sigma_star: float = SIGMA_STAR) -> Dict:
    """三态判定 + 可观察量（A_k、加载成功与否）。"""
    A_k = [m for m in candidates if m.sigma(fact_avail) >= sigma_load]
    loaded = [m for m in A_k if m.sigma(fact_avail) >= sigma_star]

    if not A_k:
        state = "F1_divergence"          # ¬C1：无候选可构 ⟹ 发散
    else:
        winner = max(A_k, key=lambda m: m.sigma(fact_avail))
        if winner.is_correct and loaded:
            state = "true_convergence"   # C1∧C2
        else:
            state = "F2_pseudoconvergence"  # C1∧¬C2：加载成功但收敛到错对象

    return {
        "A_k": [m.mid for m in A_k],
        "loaded": [m.mid for m in loaded],
        "winner": winner.mid if A_k else None,
        "state": state,
        "separability": {   # F1/F2 在"是否已加载 + 是否收敛到正确"上可分离
            "A_k_empty": not A_k,
            "converged_wrong": bool(A_k) and not winner.is_correct,
        },
    }


def run(seed: int = SEED):
    rng = random.Random(seed)
    candidates = make_candidates(rng)

    scenes = {
        "true_convergence": facts_for_true_conv(),
        "F1_divergence": facts_for_f1_divergence(),
        "F2_pseudoconvergence": facts_for_f2_pseudoconvergence(),
    }

    out = {}
    for name, facts in scenes.items():
        out[name] = classify(candidates, facts)

    # 分离性断言
    f1 = out["F1_divergence"]
    f2 = out["F2_pseudoconvergence"]
    tc = out["true_convergence"]
    sep_ok = (
        f1["state"] == "F1_divergence" and f1["separability"]["A_k_empty"]
        and f2["state"] == "F2_pseudoconvergence" and f2["separability"]["converged_wrong"]
        and not f2["separability"]["A_k_empty"]
        and tc["state"] == "true_convergence"
    )

    return {
        "scenes": out,
        "conditions_separated": sep_ok,
        "verdict": "HOLDS (F1↔发散/¬C1、F2↔伪收敛/C1∧¬C2 条件分离)" if sep_ok
                   else "FAIL",
        "evidence": "原型确认（结构分离，非现实证据）",
    }


if __name__ == "__main__":
    res = run()
    print(json.dumps(res, ensure_ascii=False, indent=2))
