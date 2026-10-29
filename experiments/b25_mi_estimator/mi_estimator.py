"""
X-S-MERA V3.4 · B-B25 互信息估计器 I(M; Ω | φ)
================================================================

B-B25 理论部分已闭合（v3.2 审计 S8）：跨理解者混合估计给的是 H(M|φ)，
而 H(M|φ) = H(M|φ,Ω) + I(M;Ω|φ)，故 **混合偏倚恰等于 I(M;Ω|φ)**，
是 P := H(M|φ,Ω) 定义的直接推论。剩余的是纯工程：**需要 I(M;Ω|φ) 的估计器**。

本文件实现该估计器（经验条件熵 + Miller–Madow 修正），并给出演示：

- `cond_mutual_information(M, Omega, phi)`：估计 I(M;Ω|φ) = H(M|φ) − H(M|φ,Ω)
- 分量为 H(M|φ)（混合）、H(M|φ,Ω)（分群体）、H(M|Ω,φ)（对称）
- Miller–Madow 修正：Ĥ_MM = Ĥ_ML + (K−1)/(2N)（K 非零箱数，N 样本数）

用途：P̂ 的数值**不可跨群体直接比较**，除非同时报告 I(M;Ω|φ)
（15-v3.2 §4.3 / 16-v3.2）。本估计器使该偏倚可计算、可报告。

诚实边界：估计器用**合成标签数据**演示与自测（正偏倚当 Ω 影响 M 分布、
≈0 当 Ω 无关），非真实认知数据；真实偏倚报告须待 B 线真实样本。
"""

from __future__ import annotations

import numpy as np

__all__ = ["entropy_mm", "conditional_entropy", "cond_mutual_information",
           "MIEstimator"]


def entropy_mm(series: np.ndarray) -> float:
    """经验熵 + Miller–Madow 修正（bit）。

    Ĥ_MM = −Σ p log₂p + (K−1)/(2N)，K 非零箱数，N 样本数。
    """
    series = np.asarray(series).ravel()
    n = int(series.size)
    if n == 0:
        return 0.0
    counts = np.bincount(series.astype(int))
    counts = counts[counts > 0]
    p = counts.astype(float) / n
    h = -float(np.sum(p * np.log2(p)))
    k = int(counts.size)
    correction = (k - 1.0) / (2.0 * n) if n > 0 else 0.0
    return h + correction


def conditional_entropy(M: np.ndarray, G: np.ndarray) -> float:
    """H(M | G) = Σ_g p(g) · H(M | g)，每组的 MM 修正熵按组样本量加权。"""
    M = np.asarray(M).ravel()
    G = np.asarray(G).ravel()
    assert M.size == G.size, "M 与 G 样本量须一致"
    n = int(M.size)
    if n == 0:
        return 0.0
    total = 0.0
    for g in np.unique(G):
        mask = G == g
        sub = M[mask]
        if sub.size == 0:
            continue
        total += (sub.size / n) * entropy_mm(sub)
    return total


def cond_mutual_information(M: np.ndarray, Omega: np.ndarray,
                            phi: np.ndarray) -> float:
    """I(M; Ω | φ) = H(M | φ) − H(M | φ, Ω)。"""
    M = np.asarray(M).ravel()
    Omega = np.asarray(Omega).ravel()
    phi = np.asarray(phi).ravel()
    assert M.size == Omega.size == phi.size
    # H(M|φ)：按 φ 分组的混合熵（跨 Ω 混合 → 混合偏倚的来源）
    h_m_phi = conditional_entropy(M, phi)
    # H(M|φ,Ω)：按 (φ,Ω) 联合分组的熵（无偏基准）
    phi_omega = phi.astype(int) * 1000 + Omega.astype(int)
    h_m_phi_omega = conditional_entropy(M, phi_omega)
    return h_m_phi - h_m_phi_omega


class MIEstimator:
    """I(M;Ω|φ) 估计器的面向对象接口（含分量报告）。"""

    def __init__(self, M=None, Omega=None, phi=None):
        if M is not None:
            self.fit(M, Omega, phi)

    def fit(self, M, Omega, phi):
        self.M = np.asarray(M).ravel()
        self.Omega = np.asarray(Omega).ravel()
        self.phi = np.asarray(phi).ravel()
        assert self.M.size == self.Omega.size == self.phi.size
        return self

    def report(self) -> dict:
        """返回偏倚估计的完整分量。"""
        phi_omega = (self.phi.astype(int) * 1000
                     + self.Omega.astype(int))
        h_m = entropy_mm(self.M)                     # 无条件 H(M)（整体）
        h_m_phi = conditional_entropy(self.M, self.phi)      # 混合 H(M|φ)
        h_m_phi_omega = conditional_entropy(self.M, phi_omega)  # 分群体 H(M|φ,Ω)
        imi = h_m_phi - h_m_phi_omega
        return {
            "H_M": h_m,
            "H_M_phi": h_m_phi,           # 混合估计（有偏）
            "H_M_phi_Omega": h_m_phi_omega,  # 分群体（无偏基准）
            "I_M_Omega_phi": imi,         # 偏倚 = 混合 − 分群体
            "bias_note": ("P̂ 不可跨群体比较，除非报告 I(M;Ω|φ)"
                          "（15-v3.2 §4.3）"),
        }
