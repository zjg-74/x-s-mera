# -*- coding: utf-8 -*-
"""R1 加载度对偶测试（推论 R1 / experiments/r1_loading_duality/）。"""
import os
import sys
import json

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), "..", "experiments", "r1_loading_duality"))

import r1_loading_duality as r1


def test_structure_counts():
    res = r1.run()
    assert res["n_related_pairs"] == 10          # C(5,2)
    assert res["n_random_pairs"] == 10
    assert len(res["d_related"]) == 10
    assert len(res["d_random"]) == 10


def test_related_over_random():
    res = r1.run()
    # 每个相关宏对 d > 每个随机宏对 d（原型上应整体分离）
    assert all(dr > drnd for dr in res["d_related"] for drnd in res["d_random"])


def test_separation_holds():
    res = r1.run()
    assert res["separation"] > 0.3
    assert res["verdict"].startswith("HOLDS")


def test_deterministic_seed():
    a = r1.run(r1.SEED)
    b = r1.run(r1.SEED)
    assert a == b
