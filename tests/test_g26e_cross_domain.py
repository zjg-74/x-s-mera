#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""G-26e 结构性测试：D 系列判定程序跨域可移植性（锁结构，不锁数字）。"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, "..", "experiments", "g26e_cross_domain")
sys.path.insert(0, EXP)

from g26e_cross_domain import (  # noqa: E402
    make_cfg, trajectory, d_drill, d_conv, d_ident,
)


def _read(cfg, E0, U0):
    Es, Us, dEs, dUs = trajectory(cfg, E0, U0)
    return d_drill(cfg, Es, Us, dEs, dUs), d_conv(Us), d_ident(cfg)


def test_drill_effective_control():
    """流行病有效管控（ρ=2）：钻研相 + 收敛方向。"""
    cfg = make_cfg(1.0, 2.0)
    drill, conv, ident = _read(cfg, E0=5.0, U0=2.0)
    assert drill["wedge"] == "drill"
    assert drill["peak_order"] == "U_first"
    assert "收敛" in conv["verdict"]
    assert abs(ident["A_readable"]["rho_exp"] - 1.0) < 1e-9


def test_hollow_containment_lost():
    """流行病防控失守（ρ=0.6）：空洞相 + 发散方向。"""
    cfg = make_cfg(1.0, 0.6)
    drill, conv, ident = _read(cfg, E0=5.0, U0=20.0)
    assert drill["wedge"] == "hollow"
    assert drill["peak_order"] == "E_first"
    assert "发散" in conv["verdict"]
    assert ident["B_readable"]["sign(1-rho)"] == +1


def test_cross_domain_structure_invariant():
    """判定结构（A/B/C/D 分级 + 钻研/空洞楔形）跨域不变（不引用领域）。"""
    for rho in (2.0, 0.6):
        cfg = make_cfg(1.0, rho)
        ident = d_ident(cfg)
        assert set(ident) == {"A_readable", "B_readable", "C_needs_external", "D_pseudo_target"}
        assert "kappa_E" in ident["C_needs_external"][0]
    # 楔形与 sign(1-ρ) 一致
    assert make_cfg(1.0, 2.0)["rho"] > 1
    assert make_cfg(1.0, 0.6)["rho"] < 1
