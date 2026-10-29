# -*- coding: utf-8 -*-
"""R2 失效分离测试（推论 R2 / experiments/r2_failure_separation/）。"""
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), "..", "experiments", "r2_failure_separation"))

import r2_failure_separation as r2


def test_scenes_classified():
    res = r2.run()
    assert res["scenes"]["true_convergence"]["state"] == "true_convergence"
    assert res["scenes"]["F1_divergence"]["state"] == "F1_divergence"
    assert res["scenes"]["F2_pseudoconvergence"]["state"] == "F2_pseudoconvergence"


def test_conditions_separated():
    res = r2.run()
    assert res["conditions_separated"] is True
    assert res["verdict"].startswith("HOLDS")


def test_f1_empty_ak():
    assert r2.run()["scenes"]["F1_divergence"]["separability"]["A_k_empty"] is True


def test_f2_converged_wrong_not_empty():
    s = r2.run()["scenes"]["F2_pseudoconvergence"]
    assert s["separability"]["converged_wrong"] is True
    assert s["separability"]["A_k_empty"] is False


def test_true_conv_loaded_correct():
    s = r2.run()["scenes"]["true_convergence"]
    assert s["loaded"] == ["M_correct"]


def test_deterministic_seed():
    assert r2.run(r2.SEED) == r2.run(r2.SEED)
