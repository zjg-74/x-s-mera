# -*- coding: utf-8 -*-
"""R1：加载度 σ(M, F_k) 与工程边界重合度 d 的同构对偶（主文 §9 / 推论 R1）。

规格来源：`docs/xs-framework/01-v3.1-main-dichotomy-macro-loading-trichotomy.md`
  §4.4 + §9 R1——加载度 σ(M, F_k) 给出单个宏与其可用事实的贴合程度；工程版边界
  重合度 d = exp(-r^2/(2σ^2)) 给出两系统边界的重合程度。二者**同构，不同尺度**。
  推论 R1：**当多个宏的槽被同一批事实填充时，这些宏之间的 d 应显著高于随机配对
  的宏**。

本文在"宏-槽-事实"离散原型上验证该对偶的结构成立（原型确认，非现实证据）：

  - 事实批 F：特征空间上的一批实例（"同一批事实"）。
  - 宏 M：一组槽，每槽声明其期望特征子集。
  - 加载度 σ(M, F) := 被事实批满足的槽的加权比例（槽满足 = 槽期望特征在 F 中
    均有足额实例）。【约定】槽间独立加权，不模拟嵌套依赖。
  - 边界重合度 d(M1, M2) := 两宏在"槽满足模式"（加权 0/1 向量）上的余弦相似度。
    对偶语义：同批事实填充 → 共享依赖宏的加载模式高度重合 → d 高；随机配对宏
    无共享事实依赖 → 加载模式不相交 → d 低。

比较在**同一共享事实批 F** 下进行：相关宏（依赖 F 的核心）d 高，随机宏（依赖
F 之外特征）d 低。仅标准库，固定种子。
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from typing import Dict, List

SEED = 20261029


@dataclass
class Macro:
    mid: str
    slots: List[Dict]

    def sigma(self, fact_avail: Dict[str, int], w_min: int = 1) -> float:
        total = sum(s["w"] for s in self.slots)
        if total == 0:
            return 0.0
        ok = 0.0
        for s in self.slots:
            if all(fact_avail.get(f, 0) >= w_min for f in s["feats"]):
                ok += s["w"]
        return ok / total

    def load_vector(self, fact_avail: Dict[str, int], w_min: int = 1) -> List[float]:
        return [float(s["w"]) if all(fact_avail.get(f, 0) >= w_min for f in s["feats"])
                else 0.0 for s in self.slots]


def boundary_overlap(m1: Macro, m2: Macro, fact_avail: Dict[str, int],
                     w_min: int = 1) -> float:
    """边界重合度 d(M1,M2)：槽满足模式（加权 0/1 向量）的余弦相似度。"""
    l1 = m1.load_vector(fact_avail, w_min)
    l2 = m2.load_vector(fact_avail, w_min)
    dot = sum(a * b for a, b in zip(l1, l2))
    n1 = math.sqrt(sum(a * a for a in l1))
    n2 = math.sqrt(sum(b * b for b in l2))
    if n1 == 0 or n2 == 0:
        return 0.0
    return dot / (n1 * n2)


# ---------------------------------------------------------------------------
# 构造
# ---------------------------------------------------------------------------

def make_related_macros(rng: random.Random, core: List[str], excl: List[str],
                        n: int = 5, n_core: int = 3, n_uniq: int = 2):
    """相关宏：核心槽（依赖共享 core，同一批事实填充 ⟹ 均满足）+ 独有槽（依赖
    excl 的随机子集，仅当其特征被共享批供给时满足）。⟹ 各宏核心槽全 1、独有槽
    满足情况不同 → 加载模式部分重合，d 落在中间区间。"""
    macros = []
    for i in range(n):
        slots = [{"feats": sorted(core), "w": 1} for _ in range(n_core)]
        slots += [{"feats": rng.sample(excl, 2), "w": 1} for _ in range(n_uniq)]
        macros.append(Macro(f"rel{i}", slots))
    return macros


def make_random_macros(rng: random.Random, excl: List[str], n: int = 5,
                       slots_per_macro: int = 3):
    """随机宏：每槽依赖核心 base 之外的随机特征（同一共享批下不满足）。"""
    macros = []
    for i in range(n):
        slots = [{"feats": rng.sample(excl, 3), "w": 1}
                 for _ in range(slots_per_macro)]
        macros.append(Macro(f"rnd{i}", slots))
    return macros


def fact_batch(feats: List[str], instances: int = 3) -> Dict[str, int]:
    """同一批事实：给指定特征足额实例。"""
    return {f: instances for f in feats}


def run(seed: int = SEED):
    rng = random.Random(seed)
    pool = [f"f{i}" for i in range(100)]
    base = rng.sample(pool, 9)                 # 共享核心特征池
    excl = [f for f in pool if f not in base]  # 核心之外
    core = base[:4]                            # 相关宏共享依赖的核心事实
    excl_supply = rng.sample(excl, 6)          # 共享批对 excl 的部分供给

    rel_macros = make_related_macros(rng, core, excl)
    rnd_macros = make_random_macros(rng, excl)

    # 同一批共享事实 F：core（全供）∪ excl 的部分供给 ⟹ 相关宏核心槽全满足、
    # 独有槽满足各异 → 加载模式部分重合；随机宏（依赖 excl 外）不重合。
    shared_facts = fact_batch(list(core) + excl_supply)

    d_rel = [boundary_overlap(rel_macros[i], rel_macros[j], shared_facts)
             for i in range(len(rel_macros)) for j in range(i + 1, len(rel_macros))]
    d_rnd = [boundary_overlap(rnd_macros[i], rnd_macros[j], shared_facts)
             for i in range(len(rnd_macros)) for j in range(i + 1, len(rnd_macros))]

    mean_rel = sum(d_rel) / len(d_rel)
    mean_rnd = sum(d_rnd) / len(d_rnd)
    sep = mean_rel - mean_rnd

    return {
        "n_related_pairs": len(d_rel),
        "n_random_pairs": len(d_rnd),
        "d_related": [round(x, 4) for x in d_rel],
        "d_random": [round(x, 4) for x in d_rnd],
        "mean_d_related": round(mean_rel, 4),
        "mean_d_random": round(mean_rnd, 4),
        "separation": round(sep, 4),
        "verdict": "HOLDS (同批事实填充 → 宏间 d 显著高于随机配对)" if sep > 0.3
                   else "WEAK/FAIL",
        "evidence": "原型确认（结构对偶，非现实证据）",
    }


if __name__ == "__main__":
    res = run()
    print(json.dumps(res, ensure_ascii=False, indent=2))
