#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B3 —— 残差反演 g（R32 / 推论 R43）的结构性测试（锁结构，不锁数字）

对应：`experiments/b3_residual_inversion/`、
`docs/xs-framework/15-v3.2-gauge-fixed-framework.md` §5.7（推论 R43）。

锁的结构预言：
  * 推论 R43 的可检验形式：共线残差 ≡0 ⟺ g 沿路径为常数；残差形状 = g 变化的反演
  * 在**已知 g 的合成轨迹**上，反演值与真值的（中心化）相关 ≈1、归一化 RMSE 在判阈内
  * 负对照（g≡0 与 g≡const）不产生虚假结构（反演 std ≈0）
  * 反演只恢复**形状/变异**，绝对水平不可辨识：常数平移被 ∫(g−ḡ) 吸收
    ⟹ 反演输出是"中心化的 g"，其均值必为 0（规范自由度 B-B23 的体现）
  * 工具函数有牙齿：pearson(x,x)=1、pearson(x,−x)=−1、正交序列≈0；
    中心差分对线性序列还原斜率

不锁：任何具体 corr/RMSE 数值（只锁"高于/低于判阈"）、U 的具体取值范围。
诚实边界：这是**方法验证**（已知 g 的合成轨迹），不是经验验证；
「真实认知残差确实是 g 信号」不由此声称。
"""

import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "b3_residual_inversion"))

import b3_residual_inversion as R  # noqa: E402


def test_module_imports():
    for fn in ("GFUNCS", "ad_plus_trajectory", "invert_g",
               "linfit", "pearson", "center_diff", "main"):
        assert hasattr(R, fn), fn


def test_pearson_has_teeth():
    xs = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert abs(R.pearson(xs, xs) - 1.0) < 1e-12
    assert abs(R.pearson(xs, [-v for v in xs]) + 1.0) < 1e-12
    # 两列正交（和为零的对称取值）⟹ 相关 ≈0
    ys = [1.0, -1.0, 1.0, -1.0, 0.0]
    zs = [1.0, 1.0, -1.0, -1.0, 0.0]
    assert abs(R.pearson(ys, zs)) < 1e-9


def test_linfit_recovers_line():
    xs = [0.0, 1.0, 2.0, 3.0]
    ys = [5.0, 7.0, 9.0, 11.0]
    a, b = R.linfit(xs, ys)
    assert abs(a - 5.0) < 1e-12
    assert abs(b - 2.0) < 1e-12


def test_center_diff_recovers_linear_slope():
    step = 0.01
    xs = [2.0 + 3.0 * i * step for i in range(50)]
    d, start = R.center_diff(xs, step)
    assert start == 1
    assert all(abs(v - 3.0) < 1e-9 for v in d)


def test_g_functions_satisfy_domain_conditions():
    """g 的辖域条件：g(0)=0 且单调不减（在采样区间内）。"""
    for name, g in R.GFUNCS.items():
        assert abs(g(0.0)) < 1e-12, name
        vals = [g(0.25 * i) for i in range(1, 20)]
        assert all(v2 >= v1 - 1e-12 for v1, v2 in zip(vals, vals[1:])), name


def test_invert_g_recovers_shape():
    """已知 g 的合成轨迹上，反演值与真值一致（R32 判据原文）。"""
    for name, g in R.GFUNCS.items():
        d = R.invert_g(g)
        assert d["corr"] > 0.99, (name, d["corr"])
        assert d["norm_rmse"] < 0.05, (name, d["norm_rmse"])
        assert d["g_true_std"] > 0.0, name


def test_negative_controls_produce_no_structure():
    """g≡0 与 g≡const（无变异）不得反演出虚假结构。"""
    for name, g in (("g_zero", R.g_zero), ("g_const", R.g_const)):
        d = R.invert_g(g)
        assert d["g_true_std"] < 1e-12, (name, d["g_true_std"])
        assert d["g_inv_std"] < 0.02, (name, d["g_inv_std"])


def test_affine_part_of_residual_is_absorbed_but_not_curvature():
    """规范自由度 B-B23：残差的常数项与线性项被吸收 ⟹ g 的绝对水平不可辨识；
    但**非线性部分不被吸收**（有牙齿：吸收不是"一律不变"）。"""
    ks, Es, Us = R.ad_plus_trajectory(R.g_linear)
    resid = [u - R.RHO_EXP * e for u, e in zip(Us, Es)]
    a, b = R.linfit(ks, resid)
    base = [r - (a + b * k) for r, k in zip(resid, ks)]
    d0, _ = R.center_diff(base, R.STEP)

    for shift, slope in ((5.0, 0.0), (0.0, 0.3), (-2.0, 0.7)):
        resid2 = [r + shift + slope * k for r, k in zip(resid, ks)]
        a2, b2 = R.linfit(ks, resid2)
        d2, _ = R.center_diff([r - (a2 + b2 * k) for r, k in zip(resid2, ks)], R.STEP)
        assert max(abs(p - q) for p, q in zip(d0, d2)) < 1e-9, (shift, slope)

    resid3 = [r + 0.001 * k * k for r, k in zip(resid, ks)]
    a3, b3 = R.linfit(ks, resid3)
    d3, _ = R.center_diff([r - (a3 + b3 * k) for r, k in zip(resid3, ks)], R.STEP)
    assert max(abs(p - q) for p, q in zip(d0, d3)) > 1e-9


def test_r43_equivalence_holds_both_directions():
    """推论 R43 的 ⟺ 双向：g 沿路径为常数 ⟺ 共线残差无结构。"""
    d_const = R.invert_g(R.g_const)
    assert d_const["g_true_std"] < 1e-12
    assert d_const["g_inv_std"] < 1e-9
    for name, g in R.GFUNCS.items():
        d = R.invert_g(g)
        assert d["g_true_std"] > 1e-9, name
        assert d["g_inv_std"] > 1e-9, name


def test_main_records_cases_and_controls(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "HERE", str(tmp_path))
    R.main()
    import json
    with open(tmp_path / "b3_results.json", "r", encoding="utf-8") as fh:
        res = json.load(fh)
    assert set(res["cases"]) == set(R.GFUNCS)
    for name, d in res["cases"].items():
        assert d["corr"] > 0.99 and d["norm_rmse"] < 0.05, (name, d)
    assert set(res["negative_controls"]) == {"g_zero", "g_const"}
    for name, d in res["negative_controls"].items():
        assert d["g_inv_std"] < 0.02, (name, d)
    assert "方法" in res["meta"]["task"] or "合成" in res["meta"]["engine"]
