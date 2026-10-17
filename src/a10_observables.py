"""R6: A10 封装性失败的**可观测现象清单**（预注册）与三联排除判定.

为什么需要本模块
----------------
A10 的原文是：M-完备性失败 ⇔ 存在"必须无限细分才能解释的稳定行为"。
这句话在纸面上无法证伪——"无限细分"是极限陈述，任何有限观测都用不上它。
A10-03 因此要求：**先写下"出现什么现象即判定失败"，再去观测**；
写不出来，就说明这条公理没有证伪能力。

本模块做三件事：

1. **预注册判据表** ``CRITERIA``：任何仿真运行之前写定，阈值不得在看到
   结果之后调整（若要调整，必须开新版本并登记原因）。
2. **粒度扫描** ``granularity_scan``：把"无限细分"换成有限可测的
   **残差–粒度标度律**。若残差幅度随描述粒度变粗而衰减（斜率 < 0），
   说明该异常有有限特征尺度，细分即可解释 → 不构成 M-不完备；
   若残差幅度对粒度**近似不变**（斜率 ≈ 0），说明没有特征尺度，
   任何有限细分都解释不了 → M-不完备候选。
3. **三联排除** ``classify``：把"失败信号"逐级排除成
   礁石（接口问题）/ 分辨率伪影（有限尺度）/ 真坍塌（无特征尺度）。

R2 的三条教训，本模块强制执行
------------------------------
* **禁用 σ 倍数口径**（跨域不可比）：所有阈值都是绝对标度量。
* **禁用峰度作主判据**：峰度对 (ρ, κ) 多对一，跨域不可辨识；
  这里只在"多峰"判据里用 Sarle 双峰系数 BC（形状量，不用峰度值）。
* **τ 必须外生**：申报精度 ``DELTA_QUOTED`` 来自投影系统自己声明的分辨率，
  不从 D 反推，否则判据循环。

判据的盲区（也是预注册的一部分）
--------------------------------
P2–P5 **全部是形态量或标度律量，对"整体平移"不变**（IQR、BC、std 都是
位置不变量）。因此系统性偏置（D 整体偏了，但形状完全正常）**不可能**
被本清单捕获，必须靠独立路径比对（C2 = 跨实现一致率，见 R2 结论）。
这一点在 R6 里被显式构造为阴性对照，登记为判据缺口 G2。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Criterion", "CRITERIA", "DELTA_QUOTED", "DELTA_ABS",
    "block_average", "granularity_scan", "window_stats", "evaluate",
    "evaluate_persistent", "repair_encapsulated", "classify", "VERDICTS",
]

# ---------------------------------------------------------------- 外生标度
# 申报精度：投影系统自己声明的分辨率（= 内核 d_sample_noise）。
# 这是外生量，不由 D 估计——否则"承诺兑现比"变成同义反复。
DELTA_QUOTED = 0.01
# 绝对标度阈值 = 3 倍申报精度（0.03）。注意是 3×Δ_quoted，不是 3σ：
# σ 从 D 估计出来会把判据自身的失败重新吸收进阈值。
DELTA_ABS = 3.0 * DELTA_QUOTED   # 0.06

# Sarle 双峰系数阈值（BC > 5/9 提示非单峰）
BC_MAX = 5.0 / 9.0
# 承诺兑现比上限：核心宽度（IQR）不得超过申报精度的这么多倍
COMMIT_MAX = 2.0
# 标度律斜率上限：比这更平（更接近 0）= 无特征尺度 = 不可细分解释
BETA_MAX = -0.15


@dataclass(frozen=True)
class Criterion:
    """一条可观测判据。字段 design 上是"给观测者看的作业指导"。"""
    code: str
    name: str
    statistic: str      # window_stats / granularity_scan 的键
    op: str             # "ge"（大于等于即触发）| "le"
    threshold: float
    gate: str           # "C1" = 仅在 C1 成立的窗口内评估；"none" = 无条件
    observable: str     # 现场/实验上怎么测
    limitation: str     # 已知盲区


# ============================================================ 预注册判据表
# 写定时间：2026-10-17，早于本模块所配的任何一次运行。
# 阈值依据：DELTA_QUOTED 取内核申报分辨率；BC 取 Sarle 标准值 5/9；
#          COMMIT_MAX 取"核心宽度不超过申报精度 2 倍"这一工程惯例；
#          BETA_MAX 取"一个数量级的粒度变化内残差至少衰减到 1/√10 ≈ 0.32"
#          对应的斜率 -0.5 的一半，留出余量。
CRITERIA: tuple[Criterion, ...] = (
    Criterion(
        code="P2", name="多峰（D 分裂为多个可分辨支）",
        statistic="bc", op="ge", threshold=BC_MAX, gate="C1",
        observable="D 的直方图/密度出现两个以上可分辨支，Sarle BC 越过 5/9",
        limitation="白噪声也能抬 BC；不能区分界面失配与真结构分歧",
    ),
    Criterion(
        code="P3", name="绝对标度尾部质量超界",
        statistic="tail_abs", op="ge", threshold=0.05, gate="C1",
        observable="|d - median(D)| > 3Δ_quoted 的样本占比超过 5%",
        limitation="对任何额外方差都敏感，特异性差",
    ),
    Criterion(
        code="P4", name="承诺兑现比失配",
        statistic="commit", op="ge", threshold=COMMIT_MAX, gate="C1",
        observable="核心宽度 IQR(D) / Δ_quoted 超过 2：申报精度没有兑现",
        limitation="位置不变量：整体平移（系统性偏置）看不见",
    ),
    Criterion(
        code="P5", name="残差无特征尺度（无法靠细分解释）",
        statistic="beta", op="ge", threshold=BETA_MAX, gate="C1",
        observable="粒度扫描：残差幅度 a(B) 对块长 B 的双对数斜率 ≥ -0.15",
        limitation="位置不变量；且需至少两个数量级的粒度范围才可估",
    ),
)

# 判定取值（避免字符串散落在各处）
VERDICTS = {
    "ok": "M-完备（无失败信号）",
    "gate": "A10 不适用（C1 不成立：交互缺失 F1，非尺度问题）",
    "reef": "礁石锁定（接口伪影：修接口即恢复，非尺度坍塌）",
    "resolution": "分辨率伪影（有限尺度过程：变粗即平均掉）",
    "collapse": "M-不完备候选（尺度坍塌：无特征尺度，需新公理/新尺度）",
    "blind": "假阴性（判据盲区：形态完好但存在系统性偏置）",
}


# ---------------------------------------------------------------- 基础算子
def block_average(x: np.ndarray, block: int) -> np.ndarray:
    """把序列按 block 长做块平均 = 用"粒度 block"的投影去看同一段历史。"""
    x = np.asarray(x, float).ravel()
    block = max(int(block), 1)
    n = (x.size // block) * block
    if n == 0:
        return np.asarray([x.mean()], float)
    return x[:n].reshape(-1, block).mean(axis=1)


def _slope(xs, ys) -> float:
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    ok = (xs > 0) & (ys > 0) & np.isfinite(ys)
    if ok.sum() < 2:
        return float("nan")
    lx, ly = np.log10(xs[ok]), np.log10(ys[ok])
    return float(np.polyfit(lx, ly, 1)[0])


def granularity_scan(x: np.ndarray, max_level: int = 7) -> dict:
    """粒度扫描：k=0..max_level，块长 B=2^k，测残差幅度 a(B)。

    返回
    ----
    b     : 块长序列
    a     : 各粒度下的幅度（块平均序列的标准差）
    bc    : 各粒度下的双峰系数
    beta  : log a ~ log b 的斜率。**这是 A10 的核心观测量**
            beta < 0  有特征尺度（异常会被更粗的描述平均掉）
            beta ≈ 0  无特征尺度（任何粒度都留着，细分永远解释不完）
    """
    levels = list(range(int(max_level) + 1))
    b, a, bc = [], [], []
    for k in levels:
        B = 2 ** k
        xb = block_average(x, B)
        b.append(float(B))
        a.append(float(xb.std(ddof=0)))
        v = float(xb.var(ddof=0))
        if v <= 0 or xb.size < 8:
            bc.append(float("nan"))
            continue
        s = np.sqrt(v)
        m3 = float(np.mean((xb - xb.mean()) ** 3))
        m4 = float(np.mean((xb - xb.mean()) ** 4))
        kp = m4 / v ** 2
        bc.append(float(((m3 / s ** 3) ** 2 + 1.0) / kp) if kp > 0 else float("nan"))
    return {"b": b, "a": a, "bc": bc, "beta": _slope(b, a), "levels": levels}


def window_stats(x: np.ndarray) -> dict:
    """单个评估窗口的全部判据统计量（绝对标度）。"""
    x = np.asarray(x, float).ravel()
    n = int(x.size)
    if n < 8:
        return {k: float("nan") for k in
                ("n", "mean", "median", "sd", "iqr", "bc", "tail_abs", "commit")}
    med = float(np.median(x))
    q1, q3 = np.percentile(x, [25, 75])
    iqr = float(q3 - q1)
    v = float(x.var(ddof=0))
    sd = float(np.sqrt(v))
    if v <= 0:
        bc = float("nan")
    else:
        m3 = float(np.mean((x - x.mean()) ** 3))
        m4 = float(np.mean((x - x.mean()) ** 4))
        kp = m4 / v ** 2
        bc = float(((m3 / sd ** 3) ** 2 + 1.0) / kp) if kp > 0 else float("nan")
    return {
        "n": float(n), "mean": float(x.mean()), "median": med, "sd": sd,
        "iqr": iqr, "bc": bc,
        "tail_abs": float(np.mean(np.abs(x - med) > DELTA_ABS)),
        "commit": float(iqr / DELTA_QUOTED),
    }


# ---------------------------------------------------------------- 判据评估
def evaluate(stats: dict, c1_fraction: float = 1.0) -> dict:
    """对单窗口统计量逐条判据裁决。

    gate="C1" 的判据只在 c1_fraction ≥ 0.5 的窗口内评估；
    否则记为 None（不可评估），而不是 False——避免"没测"被当成"通过"。
    """
    out = {"gate_open": bool(c1_fraction >= 0.5), "c1_fraction": float(c1_fraction),
           "fires": {}, "values": {}, "evaluable": {}}
    for c in CRITERIA:
        v = float(stats.get(c.statistic, float("nan")))
        out["values"][c.code] = v
        if c.gate == "C1" and not out["gate_open"]:
            out["evaluable"][c.code] = False
            out["fires"][c.code] = None
            continue
        if not np.isfinite(v):
            out["evaluable"][c.code] = False
            out["fires"][c.code] = None
            continue
        out["evaluable"][c.code] = True
        out["fires"][c.code] = bool(v >= c.threshold) if c.op == "ge" else bool(v <= c.threshold)
    return out


def evaluate_persistent(series: np.ndarray, c1_series: np.ndarray,
                        n_windows: int = 3, min_fire: int = 2,
                        include_beta: bool = True) -> dict:
    """持续性检验：判据必须在 3 个连续窗口里至少 2 个触发才算触发。

    单窗口触发可能只是瞬态——A10 说的是"稳定行为"，不是抖动。
    """
    x = np.asarray(series, float).ravel()
    c1 = np.asarray(c1_series, bool).ravel()
    n = x.size
    edges = np.linspace(0, n, n_windows + 1).astype(int)
    gates = [float(c1[edges[i]:edges[i + 1]].mean()) if edges[i + 1] > edges[i] else 0.0
             for i in range(n_windows)]

    per_code = {c.code: [] for c in CRITERIA}
    per_stats = []
    for i in range(n_windows):
        seg = x[edges[i]:edges[i + 1]]
        st = window_stats(seg)
        st["beta"] = float("nan")
        if include_beta:
            st["beta"] = granularity_scan(seg, max_level=6)["beta"]
        per_stats.append(st)
        ev = evaluate(st, gates[i])
        for code in per_code:
            per_code[code].append(ev["fires"][code])

    fires, counts = {}, {}
    for code, seq in per_code.items():
        truthy = [b for b in seq if b is not None]
        counts[code] = len(truthy)
        # 一个窗口都不可评估 → 返回 None（"没测"不等于"没触发"）
        fires[code] = (None if len(truthy) == 0
                       else bool(len(truthy) >= min_fire and sum(truthy) >= min_fire))
    gate_open = bool(np.mean(gates) >= 0.5)
    return {"fires": fires, "counts": counts, "gate_open": gate_open,
            "c1_fraction": float(np.mean(gates)), "windows": per_stats,
            "any_fire": bool(any(fires.values()))}


def repair_encapsulated(post_repair_series: np.ndarray,
                        c1_series: np.ndarray | None = None) -> dict:
    """接口修复检验：把接口条件（C2）恢复后，封装性是否回来。

    这是三联排除的第一步，也是唯一能把"礁石"从"坍塌"里分出来的步骤——
    因为系统性接口失配在形态上可以长得和真尺度坍塌一模一样。
    """
    c1f = 1.0 if c1_series is None else float(np.asarray(c1_series, bool).mean())
    # 只用形态判据（P2–P4）：标度律需要两个数量级的粒度范围，修复后的
    # 短窗口估不准 beta，硬算会制造假信号。
    ev = evaluate_persistent(post_repair_series, np.full_like(
        np.asarray(post_repair_series), c1f >= 0.5), n_windows=1, min_fire=1,
        include_beta=False)
    return {"encapsulated": not ev["any_fire"], "c1_fraction": c1f,
            "fires": ev["fires"]}


# ---------------------------------------------------------------- 三联排除
def classify(persist: dict, repair: dict, scan: dict,
             order: str = "repair_first") -> dict:
    """三联排除：把"失败信号"逐级归因为礁石 / 分辨率 / 坍塌。

    order="repair_first"（默认，推荐）：
        ① 区域门：C1 不成立 → A10 不适用（F1 是交互缺失，不是尺度问题）
        ② 无信号 → M-完备
        ③ 修接口后恢复 → 礁石锁定（F2 一级伪影）
        ④ 标度律斜率 < 0（异常会被粗描述平均掉）→ 分辨率伪影
        ⑤ 否则 → M-不完备候选（无特征尺度）

    order="granularity_first"：交换 ③④。保留此选项**只为了检验顺序本身
    是否 load-bearing**——若顺序无关，说明这两步是冗余的。
    """
    if not persist.get("gate_open", True):
        return {"verdict": "gate", "text": VERDICTS["gate"], "steps": ["gate_closed"]}
    if not persist.get("any_fire", False):
        return {"verdict": "ok", "text": VERDICTS["ok"], "steps": ["no_signal"]}

    beta = float(scan.get("beta", float("nan")))
    beta_decays = bool(np.isfinite(beta) and beta < BETA_MAX)
    repaired = bool(repair.get("encapsulated", False))

    if order == "repair_first":
        if repaired:
            return {"verdict": "reef", "text": VERDICTS["reef"],
                    "steps": ["signal", "repair_ok"]}
        if beta_decays:
            return {"verdict": "resolution", "text": VERDICTS["resolution"],
                    "steps": ["signal", "repair_fail", "beta_decays"]}
        return {"verdict": "collapse", "text": VERDICTS["collapse"],
                "steps": ["signal", "repair_fail", "beta_flat"]}

    # 反序：先看标度律，再看接口
    if beta_decays:
        return {"verdict": "resolution", "text": VERDICTS["resolution"],
                "steps": ["signal", "beta_decays"], "repaired": repaired}
    if repaired:
        return {"verdict": "reef", "text": VERDICTS["reef"],
                "steps": ["signal", "beta_flat", "repair_ok"]}
    return {"verdict": "collapse", "text": VERDICTS["collapse"],
            "steps": ["signal", "beta_flat", "repair_fail"]}


def criterion_table() -> list[dict]:
    """把预注册表导出为可打印的行，便于写进文档/附录。"""
    return [{"code": c.code, "name": c.name, "statistic": c.statistic,
             "rule": f"{c.statistic} {'≥' if c.op == 'ge' else '≤'} {c.threshold:g}",
             "gate": c.gate, "observable": c.observable,
             "limitation": c.limitation} for c in CRITERIA]
