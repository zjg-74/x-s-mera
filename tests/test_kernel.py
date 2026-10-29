# -*- coding: utf-8 -*-
"""src/kernel.py 的结构性测试（仅标准库）。

锁：内核的四条第一步能力——骨架/仿射积分器、T-affine/T-gauge 自证伪、
D-SCALE 三值判别接口、R23 输出适配器。锁结构预言与恒等式（纯 A_D 共线性
残差≈0、规范不变、楔形方向、ρ_res 走向分类），不锁具体数值；调底层函数、
不调 main()、零文件写入。
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from kernel import (  # noqa: E402
    GaugeFixedParams, KernelState, g_saturated,
    step, t_affine, t_gauge,
    fit_beta_with_se, dscale_verdict, dscale_classify,
    to_r23_json, rho_res_verdict,
)


def _params(rho=1.0, rho_exp=1.3, kappa_E=0.5, s=1.0):
    return GaugeFixedParams(s=s, rho_exp=rho_exp, rho=rho, kappa_E=kappa_E)


def _run(states_n=20, rho=1.0, dE=0.5):
    p = _params(rho=rho)
    st = KernelState()
    for _ in range(states_n):
        step(st, p, dE=dE)
    return st


# ---- 骨架与参数 ----
def test_params_reject_nonpositive():
    for kwargs in ({"s": 1.0, "rho_exp": 0.0, "rho": 1.0, "kappa_E": 1.0},
                   {"s": 1.0, "rho_exp": 1.0, "rho": 1.0, "kappa_E": 0.0}):
        try:
            GaugeFixedParams(**kwargs)
            raise AssertionError("应拒绝非正参数")
        except ValueError:
            pass


def test_step_advances_and_records():
    p = _params()
    st = KernelState()
    step(st, p, dE=0.5)
    assert st.k == 1
    assert st.E == 0.5
    assert len(st.history) == 1
    assert st.history[0][0] == 1  # (k, E, U)
    step(st, p, dE=0.5, cand_weights={"c1": 2.0, "c2": 1.0})
    assert len(st.cand_weights) == 2
    assert st.cand_weights[1] == {"c1": 2.0, "c2": 1.0}


# ---- T-affine：纯 A_D（rho=1，drift=0）共线性残差≈0 ----
def test_t_affine_pure_ad_is_collinear():
    st = _run(states_n=20, rho=1.0)   # rho=1 -> drift=0 -> U 严格随 E
    states = [(s[1], s[2]) for s in st.history]   # 转 (E, U)
    # 纯 A_D 下 U_k = U_0 + rho_exp*E_k，残差应为 0（浮点）
    residual = t_affine(states, _params(rho=1.0))
    assert residual < 1e-9, residual


def test_t_affine_single_point_ok():
    assert t_affine([(0.0, 1.0)], _params()) == 0.0


# ---- T-gauge：相同轨迹偏差 0 ----
def test_t_gauge_identical_is_zero():
    traj = [(0.5 * k, 0.2 * k) for k in range(10)]   # (E, U)
    assert t_gauge(traj, traj) == 0.0


def test_t_gauge_detects_deviation():
    a = [(0.5 * k, 0.2 * k) for k in range(10)]           # (E, U)
    b = [(0.5 * k, 0.2 * k + 1.0) for k in range(10)]
    dev = t_gauge(a, b)
    assert dev > 0.999, dev   # 浮点舍入下应 ≈1.0


# ---- D-SCALE 三值判读 ----
def test_dscale_verdict_three_way():
    assert dscale_verdict(0.5, 0.1) == "mid_hollow"      # beta>0 空洞方向
    assert dscale_verdict(-0.5, 0.1) == "mid_drill"      # beta<0 钻研方向
    assert dscale_verdict(0.1, 0.1) == "vertex_band"     # |beta|<=z*SE 禁区
    assert dscale_verdict(0.0, 0.0) == "vertex_band"     # 无噪声顶点


def test_dscale_verdict_band_widens_with_z():
    # z 越大，禁区越宽：同 (beta=0.15, SE=0.1)，z=1.96 判禁区，z=1.0 判符号
    assert dscale_verdict(0.15, 0.1, z=1.96) == "vertex_band"
    assert dscale_verdict(0.15, 0.1, z=1.0) == "mid_hollow"


def test_dscale_classify_wedge_direction():
    # 非共线观测：E=k^2（E 与 k 非线性相关，避免回归奇异），U 线性于 k。
    # rho>1 语义 = U 随 k 升 -> mid_hollow；rho<1 = U 随 k 降 -> mid_drill。
    def obs(beta_k):
        return [(float(k), float(k * k), beta_k * k) for k in range(1, 12)]
    hollow = dscale_classify(obs(+0.5))
    drill = dscale_classify(obs(-0.5))
    assert hollow["verdict"] == "mid_hollow", hollow
    assert drill["verdict"] == "mid_drill", drill
    assert "beta_hat" in hollow["fit"] and "se_beta" in hollow["fit"]


def test_fit_beta_needs_4_points():
    try:
        fit_beta_with_se([(0, 0.0, 0.0), (1, 1.0, 1.0), (2, 2.0, 2.0)])
        raise AssertionError("应拒绝 <4 点")
    except ValueError:
        pass


# ---- R23 输出适配器 ----
def test_to_r23_json_shape():
    states = [(0, 0.0, 0.0), (1, 1.0, 0.5)]
    out = to_r23_json(states, task_id="k1")
    traj = out["trajectories"][0]
    assert traj["id"] == "k1"
    assert len(traj["steps"]) == 2
    assert traj["steps"][0]["weights"] == {}          # 未提供权重 -> 空（不虚构）


def test_to_r23_json_passthrough_weights():
    states = [(0, 0.0, 0.0), (1, 1.0, 0.5)]
    cw = [{"a": 3.0, "b": 1.0}, {"a": 2.0, "b": 0.0}]
    out = to_r23_json(states, cand_weights=cw)
    assert out["trajectories"][0]["steps"][0]["weights"] == {"a": 3.0, "b": 1.0}


def test_rho_res_verdict_classification():
    # 震荡期 U 走向：降 -> >，升 -> <，不变 -> =
    assert rho_res_verdict([1.0, 0.9, 0.8, 0.7, 0.6], [2, 2, 2, 2, 2]) == "rho_res > rho_exp"
    assert rho_res_verdict([0.6, 0.7, 0.8, 0.9, 1.0], [2, 2, 2, 2, 2]) == "rho_res < rho_exp"
    assert rho_res_verdict([0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2, 2, 2]) == "rho_res = rho_exp"


def test_rho_res_verdict_no_oscillation():
    # 无 r>=2 段 -> None（不可判）
    assert rho_res_verdict([0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1]) is None
