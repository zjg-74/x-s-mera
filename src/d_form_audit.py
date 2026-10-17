"""R2: D-form audit utilities.

D 定义为「投影系统在其有效域边界处的残差分布」。本模块提供：
  - 两套标度下的形态统计量
      * 绝对标度（无量纲相对残差 r）  -> 核心宽度、尾部质量
      * 申报 σ 标度（z = r / sigma_core）-> 峰度
  - 合成 D 的生成器（tau / coverage / n_paths 三参数）
  - 秩相关（用于可辨识性检验，不依赖 scipy）

关键动机：物理、生物、人文三个域的 D 在 z 标度与绝对标度上
给出的排序**方向相反**。若不做标度分离，A4 的度量无法定义。
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------- statistics
def robust_sigma(x: np.ndarray) -> float:
    """1.4826 * MAD，对重尾/多峰稳健的核心尺度估计。"""
    x = np.asarray(x, float)
    med = np.median(x)
    return float(1.4826 * np.median(np.abs(x - med)))


def core_width(x: np.ndarray) -> float:
    """核心宽度：四分位距（绝对标度）。"""
    q1, q3 = np.percentile(np.asarray(x, float), [25, 75])
    return float(q3 - q1)


def excess_kurtosis(x: np.ndarray) -> float:
    """超额峰度（Pearson 峰度 - 3）。"""
    x = np.asarray(x, float)
    v = x.var()
    if v <= 0:
        return 0.0
    return float(np.mean((x - x.mean()) ** 4) / v ** 2 - 3.0)


def z_kurtosis(r: np.ndarray, sigma_core: float | None = None) -> float:
    """把残差折算成「以自身申报精度为单位」的 z，再算超额峰度。

    这是文献里最常用的口径（"4.2σ 反常"），但它对 tau 的依赖是四次的。
    """
    r = np.asarray(r, float)
    s = robust_sigma(r) if sigma_core is None else float(sigma_core)
    if s <= 0:
        return 0.0
    return excess_kurtosis(r / s)


def tail_mass(x: np.ndarray, tol: float) -> float:
    """绝对标度尾部质量：|r| > tol 的样本占比。

    注意 tol 是**绝对**阈值（无量纲相对误差），不是 σ 的倍数——
    这是让 A10 的失败信号不被 tau 放大的唯一办法。
    """
    return float(np.mean(np.abs(np.asarray(x, float)) > tol))


def bimodality_coefficient(x: np.ndarray) -> float:
    """BC = (skew^2 + 1) / pearson_kurtosis。BC > 0.555 提示双峰/多峰。"""
    x = np.asarray(x, float)
    v = x.var()
    if v <= 0:
        return np.nan
    s = np.sqrt(v)
    m3 = np.mean((x - x.mean()) ** 3)
    m4 = np.mean((x - x.mean()) ** 4)
    kurt_p = m4 / v ** 2
    if kurt_p <= 0:
        return np.nan
    skew = m3 / s ** 3
    return float((skew ** 2 + 1.0) / kurt_p)


# ---------------------------------------------------------------- generator
def synthesize(tau: float, coverage: float, n_paths: int = 1,
               a_unmodeled: float = 1.0, a_route: float | None = None,
               n: int = 20000, seed: int = 0) -> dict:
    """合成一个投影系统的边界残差 D。

    参数
    ----
    tau          : 约束强度 = 1/sigma_core（相对精度）。物理 ~1e6，人文 ~1。
    coverage     : 约束覆盖度 kappa = 残差落在已建模核心内的概率。
    n_paths      : 同一 X 下不可通约的 S 路径条数（模态来源）。
    a_unmodeled  : 未建模余项的**绝对**尺度（与 tau 无关）。
    a_route      : 不同 S 路径之间的绝对偏移尺度（默认 0.3*a_unmodeled）。

    返回 dict(r, sigma_core, overclaim, params...)

    overclaim = a_unmodeled / sigma_core 是「承诺精度与实际未建模尺度的
    失配比」——峰度基本由它决定，而不是由 tau 决定。
    """
    rng = np.random.default_rng(seed)
    sigma_core = 1.0 / float(tau)
    if a_route is None:
        a_route = 0.3 * a_unmodeled

    r = rng.normal(0.0, sigma_core, n)

    # 未建模余项：绝对尺度，与 tau 无关（这是关键建模假设）
    out = rng.random(n) > coverage
    if out.any():
        r[out] += rng.normal(0.0, a_unmodeled, int(out.sum()))

    # 多条 S 路径：离散模态偏移，同样是绝对尺度
    if n_paths > 1:
        k = rng.integers(0, n_paths, n)
        r += (k - (n_paths - 1) / 2.0) * a_route

    return {
        "r": r,
        "tau": float(tau),
        "coverage": float(coverage),
        "n_paths": int(n_paths),
        "sigma_core": sigma_core,
        "overclaim": float(a_unmodeled / sigma_core),
        "core_width": core_width(r),
        "z_kurtosis": z_kurtosis(r, sigma_core),
        "tail_mass": tail_mass(r, tol=0.1),
        "bc": bimodality_coefficient(r),
        "n_outlier": int(out.sum()),
    }


# ---------------------------------------------------------------- ranking
def _rank(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float)
    order = a.argsort()
    r = np.empty(len(a), float)
    r[order] = np.arange(len(a), dtype=float)
    return r


def spearman(a, b) -> float:
    """秩相关（不依赖 scipy）。"""
    ra, rb = _rank(a), _rank(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0
