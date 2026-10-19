"""R28 回归测试：微观回落的本体性 / 峰的尺度依赖 / 多峰 / 局部-全局收敛。

设计原则（沿用 R18 以来的一贯做法）：**锁结构，不锁数值**。
数值是原型参数（对数斜坡、生灭率）的产物，换成别的调度就会变；
可迁移的是**关系**——解析律、单调性、恒等式、以及"局部收敛 ⟹̸ 全局收敛"
这样的分离性断言。

每条测试都对应补篇七的一条命题或一条已登记的缺口。
"""

import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, os.pardir, "experiments", "r28_micro_reversal")
JSON = os.path.join(EXP, "r28_results.json")
SRC = os.path.join(EXP, "r28_micro_reversal.py")

with open(JSON, encoding="utf-8") as f:
    R = json.load(f)


# ------------------------------------------------------------------ 工具


def _code_lines(text):
    """剔除模块 docstring 与 # 注释，只留可执行代码行。

    必须这么做：脚本的 docstring 里**记录了**陷阱本身（含相关字样），
    若不做剔除，守卫会把自己的文档当成违规代码，产生假阳性
    （R27 的 test_source_does_not_use_builtin_hash_for_seeding 已踩过一次）。
    """
    out = []
    in_doc = False
    for ln in text.splitlines():
        if ln.count('"""') % 2 == 1:
            in_doc = not in_doc
            continue
        if in_doc:
            continue
        ln = "" if ln.lstrip().startswith("#") else ln.split("#", 1)[0]
        if ln.strip():
            out.append(ln)
    return out


# ------------------------------------------------------------------ 臂 A1


def test_A1_analytic_fallback_law_holds_within_mc_error():
    """命题 6.6：P(微观回落) = Phi((ln th - ln r_bar + s^2/2)/s)，精确成立。"""
    for a in R["A1_fallback"]:
        assert a["within_3se"], (
            "s={} 最大偏差 {} 超过 3 倍 MC 标准误 {}".format(
                a["s"], a["max_abs_dev"], 3 * a["mc_se"]))
        assert a["mean_abs_dev"] < a["mc_se"], a


def test_A1_fallback_at_peak_exceeds_one_half_and_grows_with_s():
    """推论：峰处回落率 = Phi(s/2) > 1/2，且随 s 单调增。

    这是本轮最尖锐的预测——直觉会说"峰顶增量为 0，回落率应是 1/2"，
    但对数噪声的中位数 < 均值，所以严格大于 1/2。
    """
    prev = 0.0
    for a in R["A1_fallback"]:
        assert a["fallback_at_kstar"] > 0.5, a
        assert abs(a["analytic_at_kstar"] - a["predicted_Phi_s_over_2"]) < 1e-6, a
        assert a["fallback_at_kstar"] > prev
        prev = a["fallback_at_kstar"]


# ------------------------------------------------------------------ 臂 A2


def test_A2_unimodality_is_false_on_path_level_and_recovers_with_scale():
    """B-B27 升级后的核心：单峰性在路径层**不成立**（这是过程事实，不是测不了）。"""
    for a in R["A2_scale"]:
        rows = {r["w"]: r for r in a["rows"]}
        if a["s"] >= 0.2:
            assert rows[1]["uni_h_0.0"] <= 0.05, a["s"]
        assert rows[41]["uni_h_0.0"] >= 0.8, a["s"]


def test_A2_peak_count_monotone_in_window_and_prominence():
    for a in R["A2_scale"]:
        assert a["monotone_in_w"], a["s"]
        assert a["monotone_in_h"], a["s"]


def test_A2_prominence_filtering_beats_smoothing():
    """显著度过滤在 w=1 就已恢复单峰，而平滑需要 w=41。

    含义：**定义"峰"应优先用显著度（拓扑持久性），而不是平滑**——
    平滑会移动峰位（见下一条测试），显著度不会。
    """
    for a in R["A2_scale"]:
        if a["s"] <= 0.4:
            rows = {r["w"]: r for r in a["rows"]}
            assert rows[1]["uni_h_0.1"] >= 0.99, a["s"]


def test_A2_peak_position_is_robust_but_smoothing_biases_it_forward():
    """峰位对微观回落鲁棒（E 层量）；但平滑窗口越大，峰位越向后偏。"""
    for a in R["A2_scale"]:
        pred = a["kstar_pred"]
        for row in a["rows"]:
            assert abs(row["peak_mean"] - pred) <= 5.0, (a["s"], row)
        rows = {r["w"]: r for r in a["rows"]}
        assert rows[41]["peak_mean"] >= rows[1]["peak_mean"] - 1e-9, a["s"]


# ------------------------------------------------------------------ 臂 B


def test_B_multipeak_identity_peaks_eq_reopen_plus_one():
    """峰数 = 由减转增的相变次数 + 1（补篇六 N2 的精确形式）。"""
    for b in R["B_multipeak"]:
        assert b["deterministic_peak_count"] == 2, b
        assert b["reopen_count_sign_flips"] == 1, b
        assert b["identity_peaks_eq_1_plus_reopen"], b


def test_B_multipeak_survives_process_noise():
    for b in R["B_multipeak"]:
        if "noisy_peak_count_mean_w21_h005" in b:
            assert abs(b["noisy_peak_count_mean_w21_h005"] - 2.0) < 0.2, b


# ------------------------------------------------------------------ 臂 C


def test_C_micro_reversal_is_endogenous():
    """生灭模型**没有任何外加噪声**，扩张相仍有 13% 的回落 —— 第 1 条成立。"""
    c = R["C_mprocess"]
    assert c["micro_reversal_rate_expansion_phase"] > 0.05, c


def test_C_local_convergence_does_not_imply_global():
    """命题 6.7：局部收敛 ⟹̸ 全局收敛。分离度必须在每个 tau_c 上都显著。"""
    c = R["C_mprocess"]
    for t in c["convergence_table"]:
        assert t["frac_has_local_convergence"] > 0.8, t
        assert t["frac_global_converged"] < 0.5, t
        assert (t["frac_has_local_convergence"]
                - t["frac_global_converged"]) > 0.3, t
        assert t["mean_reopen"] > 0, t


def test_C_tightening_tau_c_reduces_both_counts():
    """tau_c 是**声明式检测门槛**，不是内约束强度（R18-b 的教训）：
    收紧它会同时减少收敛段数与重开数，两者都趋近 0 —— 说明它是尺度参数。"""
    tab = R["C_mprocess"]["convergence_table"]
    for i in range(len(tab) - 1):
        assert tab[i]["mean_conv_runs"] > tab[i + 1]["mean_conv_runs"]
        assert tab[i]["mean_reopen"] > tab[i + 1]["mean_reopen"]


def test_C_topological_alternation_law_holds_exactly():
    """定理 6.8：|n_+ - n_-| <= 1 恒成立（极值交替的拓扑律）。"""
    for d in R["C_mprocess"]["peak_count_table"]:
        assert d["frac_alternation_abs_le_1"] == 1.0, d
        assert d["n_ends_descending"] > 0, d
        # 止于下降段时 n_+ = n_- + 1；残差来自端点并列平台（已知，见报告）
        assert d["frac_identity_given_ends_descending"] >= 0.85, d


def test_C_prominence_filtering_breaks_the_exact_identity():
    """显著度过滤**破坏**恒等式（两可比峰 + 一浅谷 -> 谷被滤掉、两峰保留）。
    这条测试锁住的是"不要假装过滤后计数仍满足恒等式"。
    """
    for d in R["C_mprocess"]["peak_count_table"]:
        assert d["mean_peaks_prom0.5"] < d["mean_peaks_unfiltered"], d
        assert d["mean_valleys_prom0.5"] < d["mean_valleys_unfiltered"], d


# ------------------------------------------------------------------ 源码守卫


def test_source_does_not_use_builtin_hash_for_seeding():
    src = open(SRC, encoding="utf-8").read()
    bad = [ln for ln in _code_lines(src)
           if re.search(r"\bhash\s*\(", ln) and "zlib" not in ln]
    assert not bad, "builtin hash() used for seeding: {}".format(bad)
    assert "zlib.crc32" in src


def test_source_keeps_noise_in_process_layer_not_observation_layer():
    """本轮的核心建模改动：噪声必须在**过程层**，不能退化回观测层。

    R27 把噪声加在观测量上（`Uobs = U + gauss`），于是"单峰性测不出来"
    被误记为测量问题。R28 必须把随机性放进 r 本身。
    """
    src = open(SRC, encoding="utf-8").read()
    code = _code_lines(src)
    assert any("def sample_r" in ln for ln in code), "sample_r missing"
    body = src.split("def sample_r", 1)[1].split("\ndef ", 1)[0]
    assert "exp(" in body and "gauss" in body, "r must be multiplicative-lognormal"
    # 轨迹函数不得再叠加观测噪声
    traj = src.split("def run_traj", 1)[1].split("\ndef ", 1)[0]
    assert "obs_sd" not in traj, "观测噪声不得回到 run_traj"
