"""
LE coupling kernel — unit tests
===============================

覆盖：
1. d 重合度计算：值域 [0,1]、单调性（距离越大 d 越小）
2. C1/C2 判定与三区划分（high / middle / low）
3. D 分布封装检验（A10 边界封装条款）：单峰可封装，多峰/大方差不可封装
4. F1 分支：C1 不成立 → 耦合通量恒为零
5. F2 分支：C1 成立、C2 失效 → 非零 ghost 伪力，作用-反作用幅值失配
6. step 的历史登记、状态返回完整性
7. 势-流分离（A11b/A11-R）：势只进增益、流驱动状态；通量为零时 D 不漂移
8. 随机种子可复现

运行：
    pytest tests/test_le_kernel.py -v
    python tests/test_le_kernel.py         # 无 pytest 时的简易入口
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from src.le_coupling_kernel import LECouplingKernel, LEProjectionSystem  # noqa: E402


# ---------------------------------------------------------------- 1. d 重合度
def test_d_overlap_calc_range_and_monotonicity():
    kernel = LECouplingKernel()
    near = kernel.calc_d_overlap(np.array([0.1, 0.1, 0.1]), np.array([0.12, 0.12, 0.12]))
    far = kernel.calc_d_overlap(np.array([0.1, 0.1, 0.1]), np.array([8.0, 8.0, 8.0]))
    mid = kernel.calc_d_overlap(np.array([0.1, 0.1, 0.1]), np.array([0.25, 0.28, 0.22]))

    assert 0.0 <= near <= 1.0
    assert 0.0 <= mid <= 1.0
    assert 0.0 <= far <= 1.0
    assert near > 0.8, f"near states should overlap strongly, got {near}"
    assert far < 0.01, f"far states should not overlap, got {far}"
    assert near > mid > far, "d must decrease monotonically with distance"


def test_d_overlap_identical_states_is_one():
    kernel = LECouplingKernel()
    z = np.array([0.3, 0.7, 0.1])
    assert kernel.calc_d_overlap(z, z) == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------- 2. C1/C2 三区
def test_C1_C2_judge_three_region():
    kernel = LECouplingKernel(d_threshold_c1=0.25)

    # Case 1: C1=True, C2=True → high
    C1, C2, region, fail = kernel.judge_interaction_conditions(d=0.40)
    assert C1 is True and C2 is True
    assert region == "high" and fail is None

    # Case 2: C1=True, C2=False → middle / F2
    kernel.c2_compatibility = False
    C1, C2, region, fail = kernel.judge_interaction_conditions(d=0.40)
    assert C1 is True and C2 is False
    assert region == "middle" and fail == "F2"

    # Case 3: C1=False（不论 C2）→ low / F1
    kernel.c2_compatibility = True
    C1, C2, region, fail = kernel.judge_interaction_conditions(d=0.10)
    assert C1 is False
    assert region == "low" and fail == "F1"


def test_C1_threshold_boundary():
    kernel = LECouplingKernel(d_threshold_c1=0.25)
    assert kernel.judge_interaction_conditions(0.25)[0] is True       # 含等号
    assert kernel.judge_interaction_conditions(0.2499)[0] is False


# ---------------------------------------------------------------- 3. D 封装检验
def test_D_encapsulation_unimodal_passes():
    rng = np.random.default_rng(0)
    proj = LEProjectionSystem(dim_state=3, name="TestProj", noise_std=0.01)
    for _ in range(200):
        proj.record_d_sample(rng.normal(loc=0.5, scale=0.05))

    res = proj.test_D_encapsulation()
    assert res["is_encapsulated"] is True, res
    assert abs(res["mean"] - 0.5) < 0.05


def test_D_encapsulation_bimodal_rejected():
    rng = np.random.default_rng(1)
    proj = LEProjectionSystem(dim_state=3, name="TestProj", noise_std=0.01)
    for _ in range(150):
        proj.record_d_sample(rng.normal(-2, 0.2))
        proj.record_d_sample(rng.normal(2, 0.2))

    res = proj.test_D_encapsulation()
    assert res["is_encapsulated"] is False, res
    assert res["bc"] > 0.6


def test_D_encapsulation_variance_blowup_rejected():
    rng = np.random.default_rng(2)
    proj = LEProjectionSystem(dim_state=3, name="TestProj", var_max=0.02)
    for _ in range(200):
        proj.record_d_sample(rng.normal(0.5, 0.4))       # 方差 0.16 > 0.02

    res = proj.test_D_encapsulation()
    assert res["is_encapsulated"] is False
    assert "variance_unbounded" in res["reason"]


def test_D_encapsulation_insufficient_samples():
    proj = LEProjectionSystem(dim_state=3, name="TestProj")
    proj.record_d_sample(0.5)
    res = proj.test_D_encapsulation()
    assert res["is_encapsulated"] is False
    assert res["reason"] == "insufficient_samples"


# ---------------------------------------------------------------- 4/5. F1 / F2 分支
def test_F1_branch_flux_zero():
    kernel = LECouplingKernel(d_threshold_c1=0.25, semantic_potential_grad=0.5)
    d_f1 = 0.10
    C1, C2, region, fail = kernel.judge_interaction_conditions(d_f1)
    assert fail == "F1"
    flux_L, flux_E = kernel.compute_coupling_flux(d_f1, C1, C2)
    assert np.allclose(flux_L, 0.0)
    assert np.allclose(flux_E, 0.0)


def test_F2_branch_ghost_flux():
    kernel = LECouplingKernel(semantic_potential_grad=0.4)
    flux_L, flux_E = kernel.compute_coupling_flux(0.4, True, False)

    assert not np.allclose(flux_L, 0.0), "F2 must produce non-zero ghost flux"
    assert not np.allclose(flux_E, 0.0)
    assert np.sign(flux_L[0]) == -np.sign(flux_E[0]), "primary mode keeps action-reaction sign"

    # 主模态守恒，副模态幅值失配 → 通量和不为零（净能量注入）
    mismatch = float(np.linalg.norm(flux_L + flux_E))
    assert mismatch > 1e-6


def test_C2_compatible_flux_is_conservative():
    kernel = LECouplingKernel(semantic_potential_grad=0.4)
    flux_L, flux_E = kernel.compute_coupling_flux(0.4, True, True)
    assert not np.allclose(flux_L, 0.0)
    assert np.allclose(flux_L + flux_E, 0.0), "C2 兼容时通量必须守恒"


# ---------------------------------------------------------------- 6. step 历史
def test_full_step_records_history():
    kernel = LECouplingKernel(dim_life=2, dim_env=2, seed=3)
    res = kernel.step(t=0.0, dt=0.01)

    for key in ("d", "fail_mode", "L_state", "E_state",
                "D_L_encapsulation", "D_E_encapsulation"):
        assert key in res

    assert len(kernel.history["t"]) == 1
    assert len(kernel.history["d"]) == 1
    assert len(kernel.history["L_state"]) == 1
    assert kernel.history["L_state"][0].shape == (2,)
    assert set(kernel.history["region_label"]) <= {"high", "middle", "low"}


def test_history_length_matches_steps():
    kernel = LECouplingKernel(seed=4)
    dt, n = 0.05, 40
    for i in range(n):
        kernel.step(i * dt, dt)
    assert all(len(v) == n for v in kernel.history.values())


# ---------------------------------------------------------------- 7. 势-流分离
def test_potential_modulates_gain_not_state():
    """A11-R：势梯度为零时，直接耦合通量为零（流不产生），但载体自演化仍然进行。"""
    kernel = LECouplingKernel(semantic_potential_grad=0.0, seed=5)
    z0 = kernel.L.state.copy()
    flux_L, flux_E = kernel.compute_coupling_flux(0.9, True, True)
    assert np.allclose(flux_L, 0.0) and np.allclose(flux_E, 0.0)

    for i in range(50):
        kernel.step(i * 0.05, 0.05)
    assert not np.allclose(kernel.L.state, z0), "zero potential still allows carrier self-dynamics"


def test_zero_flux_limits_D_drift():
    """R8 的 F-a 检验：通量截断后，D 的漂移只剩载体自演化分量，且应随自演化收敛而衰减。

    注意：d 是状态函数，**不是通量函数**。J=0 并不意味着 D 立刻静止——
    只要两个子系统各自还在自演化，d 仍会漂移。R8 的"J→0 则 dD/dt→0"必须
    读作"耦合诱导分量→0"，本测试检验的是收敛后的残余漂移。
    """
    kernel = LECouplingKernel(d_threshold_c1=1.01, semantic_potential_grad=0.4, seed=6)
    for i in range(1200):
        kernel.step(i * 0.05, 0.05)

    tail = kernel.L.D[-200:]
    assert np.var(tail) < 1e-3, f"residual D drift too large: {np.var(tail)}"
    assert set(kernel.history["failure_mode"][-50:]) == {"F1"}


def test_potential_never_enters_state_equation():
    """A11-R：无耦合通道时，改变语义势梯度不得改变任何状态轨迹。"""
    def run(grad):
        k = LECouplingKernel(d_threshold_c1=1.01, semantic_potential_grad=grad, seed=8)
        for i in range(200):
            k.step(i * 0.05, 0.05)
        return np.array(k.history["d"]), np.array([s[0] for s in k.history["L_state"]])

    d_a, s_a = run(0.40)
    d_b, s_b = run(0.00)
    assert np.allclose(d_a, d_b)
    assert np.allclose(s_a, s_b)


def test_coupling_changes_trajectory():
    """对照：接通耦合（C1 成立）后轨迹必须与无耦合情形不同，证明通量确实驱动状态。"""
    def run(threshold):
        k = LECouplingKernel(d_threshold_c1=threshold, semantic_potential_grad=0.4,
                             coupling_gain=6.0, seed=9)
        for i in range(400):
            k.step(i * 0.05, 0.05)
        return np.array(k.history["d"])

    assert not np.allclose(run(0.25), run(1.01))


# ---------------------------------------------------------------- 8. 可复现性
def test_seed_reproducibility():
    def run():
        k = LECouplingKernel(seed=11)
        for i in range(100):
            k.step(i * 0.05, 0.05)
        return np.array(k.history["d"])

    assert np.allclose(run(), run())


def test_semantic_pot_grad_alias():
    kernel = LECouplingKernel(semantic_potential_grad=0.4)
    kernel.semantic_pot_grad = 0.12
    assert kernel.semantic_potential_grad == pytest.approx(0.12)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"[pass] {fn.__name__}")
    print(f"\nAll {len(fns)} LE kernel tests passed.")
