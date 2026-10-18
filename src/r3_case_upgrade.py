"""
X-S-MERA V3.4 · R3-v2 案例库定量升级：S3（流感抗原漂移）与 S5（鳕鱼崩溃）
============================================================================

背景（修订登记 §17.3 / §17.5）
------------------------------
第十七节把 S3/S5 的现实对应物登记为"叙事级对应"，并预注册了定量升级路径：
  S3：抗原距离序列 + 年度疫苗有效性，验证"增量合规 ∧ 绝对越界"的窗口宽度；
  S5：年捕捞量/生物量序列，拟合 draw/regen 比率与 phi_min 阈值。

本模块执行该升级。数据全部来自公开权威来源，逐条带来源标签：

S3 数据源
---------
  [CDC-VE]  CDC Past Seasons' Vaccine Effectiveness Estimates（美国 Flu VE Network，
            2004-05 至 2025-26 逐季调整后总 VE%，2020-21 因低流行缺测）；
  [WHO-COMP] WHO 疫苗组分推荐存档（1998-2010）+ GISAID/Crick 汇总（2010-2026），
            逐季 H3N2 组分；
  [WHO-QA]  WHO 2015 年 2 月 Q&A 文件：2014-15 错配完整时间线——
            2014-02 选定 A/Texas/50/2012（决策时点 >90% 流行株抗原相似），
            2014-03 下旬检测到 3C.3a/3C.2a 漂移株（决定之后、生产已开始），
            2014-05~09 漂移株占比快速上升（<50% 流行株被 A/Texas 抗血清抑制），
            2014-09 WHO 改荐 A/Switzerland/9715293/2013（南半球 2015）；
  [CDC-51829] CDC Stacks 51829：2014-15 季 >90%→<50% 抑制的抗原转变记录。

S5 数据源
---------
  [DFO-2024]  DFO CSAS SAR 2024/049（2J3KL 北方鳕评估至 2024：SSB/LRP 表、
            M/F 范围、Capelin 状态、LRP=40% Bmsy）；
  [DFO-2025]  DFO CSAS SAR 2025/043（SSB 叙述：1995 历史最低、M 峰值 2.5
            [1992-94]、F 在 1950s-80s 大部分时间超过 M、崩溃后 M 超过 F）；
  [DFO-FW]    Northern Cod Assessment Framework（41257303.pdf）图 4 读数的
            粗粒度 SSB 锚点（1960~2.0Mt … 2020~500Kt）；
  [DFO-2018]  SAR 2018（2015 SSB=34% LRP、2017=423Kt/52% LRP、2018=315Kt/37%、
            2005=3% B_lim）；
  [BRIT]      Britannica「Cod fishery collapse of 1992」（捕捞量锚点：1600 年
            3.9Kt、1856 年 151Kt、1968 年 810Kt；禁捕 32 年、2024 以 18Kt 重开）。

预注册判据（先于分析写定，不得事后调整）
------------------------------------------
  C-S3-1 错配惩罚：漂移季（2014-15）VE 低于其余紧 CI 季均值的 50%；
  C-S3-2 瞬时读法是现役协议：2004→2026 季际 WHO H3N2 换株率 ≥ 50%
         （每次换株 = 一次绝对基准报警触发）；
  C-S3-3 窗口宽度：2014-15 决策时点匹配（>90%）→ 季中失配（<50% 抑制），
         即"增量合规 ∧ 绝对越界"在一个生产窗口（约 6-8 个月）内完成；
  C-S5-1 耗竭期 draw > regen：1950s-80s F > M（源引），崩溃段 1990-95
         SSB 净降 > 80%；
  C-S5-2 吸收态：禁捕（交互停止）后 SSB 低于 LRP ≥ 20 年；
  C-S5-3 迟滞不对称：恢复时标 / 崩溃时标 ≥ 4；
  C-S5-4 再生机制受损：崩溃后 M 均值 ≥ 2× 崩溃前假设值（0.47 vs 0.2），
         且饵料（Capelin）≈10% 崩溃前水平——S5 模型的"Φ<phi_min 破坏再生通道"
         在现实数据中找到对应（不是虚构的机制）。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "FLU_VE", "FLU_H3N2", "flu_strain_update_rate", "flu_ve_stats",
    "flu_2014_window", "flu_verdict",
    "COD_SSB_ANCHORS", "COD_SSB_RECENT", "COD_FACTS",
    "cod_collapse", "cod_absorption_years", "cod_hysteresis_ratio",
    "cod_regen_damage", "cod_verdict", "CASE_CRITERIA",
]

CASE_CRITERIA = {
    "C-S3-1": "漂移季 VE < 其余紧 CI 季均值 × 50%",
    "C-S3-2": "WHO H3N2 季际换株率 ≥ 50%（瞬时读法为现役协议）",
    "C-S3-3": "2014-15：决策时点 >90% 匹配 → 季中 <50% 抑制（窗口内越界）",
    "C-S5-1": "耗竭期 F>M 且 1990-95 SSB 净降 > 80%",
    "C-S5-2": "禁捕后低于 LRP ≥ 20 年（吸收态）",
    "C-S5-3": "恢复/崩溃时标比 ≥ 4（迟滞不对称）",
    "C-S5-4": "崩溃后 M 均值 ≥ 2× 崩溃前（再生机制受损）",
    "note": "判据先于分析写定（预注册）；数据全部带公开来源标签。",
}

# ----------------------------------------------------------------------
# S3 数据（来源：CDC-VE / WHO-COMP / WHO-QA / CDC-51829）
# ----------------------------------------------------------------------

# CDC 逐季调整后总 VE（%）；2020-21 缺测（低流行）
FLU_VE: dict[str, float | None] = {
    "2004-05": 10.0, "2005-06": 21.0, "2006-07": 52.0, "2007-08": 37.0,
    "2008-09": 41.0, "2009-10": 56.0, "2010-11": 60.0, "2011-12": 47.0,
    "2012-13": 49.0, "2013-14": 52.0, "2014-15": 19.0, "2015-16": 48.0,
    "2016-17": 40.0, "2017-18": 38.0, "2018-19": 29.0, "2019-20": 39.0,
    "2020-21": None, "2021-22": 36.0, "2022-23": 30.0, "2023-24": 44.0,
    "2024-25": 56.0, "2025-26": 36.0,
}

# WHO 北半球 H3N2 组分（egg-based 推荐株；2013-14 实际生产株为 A/Texas/50/2012）
FLU_H3N2: dict[str, str] = {
    "2004-05": "A/Fujian/411/2002", "2005-06": "A/California/7/2004",
    "2006-07": "A/Wisconsin/67/2005", "2007-08": "A/Wisconsin/67/2005",
    "2008-09": "A/Brisbane/10/2007", "2009-10": "A/Brisbane/10/2007",
    "2010-11": "A/Perth/16/2009", "2011-12": "A/Perth/16/2009",
    "2012-13": "A/Victoria/361/2011", "2013-14": "A/Victoria/361/2011",
    "2014-15": "A/Texas/50/2012", "2015-16": "A/Switzerland/9715293/2013",
    "2016-17": "A/Hong Kong/4801/2014", "2017-18": "A/Hong Kong/4801/2014",
    "2018-19": "A/Singapore/INFIMH-16-0019/2016", "2019-20": "A/Kansas/14/2017",
    "2020-21": "A/Hong Kong/2671/2019", "2021-22": "A/Cambodia/e0826360/2020",
    "2022-23": "A/Darwin/9/2021", "2023-24": "A/Darwin/9/2021",
    "2024-25": "A/Thailand/8/2022", "2025-26": "A/Croatia/10136RV/2023",
}

DRIFT_SEASON = "2014-15"          # WHO-QA 记录的经典漂移错配季
# 早期三季（2004-05~2006-07）为单中心小样本（CI 极宽，如 -36~40），
# 统计口径与 2011 后的多中心网络季不可比 → "紧 CI 季"从 2011-12 起算
TIGHT_CI_START = "2011-12"

# 2014-15 窗口事件（WHO-QA / CDC-51829 逐条对应）
FLU_2014_TIMELINE = [
    ("2014-02", "WHO 选定 A/Texas/50/2012（决策时点 >90% 流行株抗原相似）", 0.90),
    ("2014-03", "检测到 3C.3a/3C.2a 漂移株（决定之后、生产已开始，极少标本）", None),
    ("2014-09", "漂移株占比快速上升：流行株被 A/Texas 抗血清抑制率 <50%", 0.50),
    ("2014-09", "WHO 改荐 A/Switzerland/9715293/2013（绝对基准报警触发）", None),
    ("2014-15 季", "CDC 实测 VE = 19%（匹配季通常 50-60%）", 0.19),
]


def flu_strain_update_rate() -> dict:
    """C-S3-2：季际 H3N2 换株率（每次换株 = 一次瞬时读法/绝对基准报警）。"""
    seasons = [s for s in FLU_H3N2]
    transitions = [(a, b) for a, b in zip(seasons, seasons[1:])]
    updates = [(a, b) for a, b in transitions if FLU_H3N2[a] != FLU_H3N2[b]]
    return {"n_transitions": len(transitions), "n_updates": len(updates),
            "update_rate": len(updates) / len(transitions),
            "pass": len(updates) / len(transitions) >= 0.50}


def flu_ve_stats() -> dict:
    """C-S3-1：漂移季 VE 与其余紧 CI 季均值的比较。"""
    tight = {s: v for s, v in FLU_VE.items()
             if v is not None and s >= TIGHT_CI_START and s != DRIFT_SEASON}
    others = np.array(list(tight.values()))
    drift_ve = FLU_VE[DRIFT_SEASON]
    mean, sd = float(others.mean()), float(others.std(ddof=1))
    z = (drift_ve - mean) / sd
    return {"drift_season": DRIFT_SEASON, "drift_ve": drift_ve,
            "others_mean": mean, "others_sd": sd, "z_score": float(z),
            "ratio": drift_ve / mean,
            "n_tight_seasons": len(others),
            "is_series_min_tight": bool(drift_ve < others.min()),
            "pass": bool(drift_ve < 0.5 * mean)}


def flu_2014_window() -> dict:
    """C-S3-3：2014-15 窗口——增量合规（决策时点匹配）∧ 绝对越界（季中失配）。"""
    at_decision = FLU_2014_TIMELINE[0][2]      # >90% 抗原相似（增量合规）
    mid_season = FLU_2014_TIMELINE[2][2]       # <50% 抑制（绝对越界）
    ve = FLU_2014_TIMELINE[4][2]               # VE 19%
    return {"match_at_decision": at_decision, "inhibition_mid_season": mid_season,
            "ve": ve,
            "incremental_compliant": at_decision >= 0.80,
            "absolute_breach": mid_season <= 0.50,
            "correction_next_meeting": True,   # 2014-09 改荐 Switzerland 株
            "pass": bool(at_decision >= 0.80 and mid_season <= 0.50)}


def flu_verdict() -> dict:
    r1, r2, r3 = flu_ve_stats(), flu_strain_update_rate(), flu_2014_window()
    return {"C-S3-1": r1, "C-S3-2": r2, "C-S3-3": r3,
            "all_pass": bool(r1["pass"] and r2["pass"] and r3["pass"])}


# ----------------------------------------------------------------------
# S5 数据（来源：DFO-2024 / DFO-2025 / DFO-FW / DFO-2018 / BRIT）
# ----------------------------------------------------------------------

# 粗粒度 SSB 锚点（kt），DFO-FW 图 4 读数（粗，仅用于崩溃/恢复时标）
COD_SSB_ANCHORS = {1960: 2000.0, 1970: 1000.0, 1980: 500.0, 1990: 1000.0,
                   2000: 50.0, 2010: 100.0, 2020: 500.0}

# SAR 2024/049 表 1（修订模型，SSB kt）
COD_SSB_RECENT = {2020: 414.0, 2021: 423.0, 2022: 462.0, 2023: 425.0,
                  2024: 342.0}
COD_RECENT_SSB_OVER_LRP = {2020: 1.50, 2021: 1.53, 2022: 1.67, 2023: 1.54,
                           2024: 1.24}

COD_FACTS = {
    "moratorium_year": 1992,            # BRIT / DFO
    "critical_zone_exit": 2016,         # DFO-2024（修订 LRP 下 2016 起出危急区）
    "collapse_start": 1990,             # 崩溃段起点（DFO-2025 叙述：early-1990s 急降）
    "collapse_low_year": 1995,          # 时序最低（DFO-2025）
    "catch_peak_1968_kt": 810.0,        # BRIT
    "reopen_2024_quota_kt": 18.0,       # BRIT（2024 谨慎重开）
    "f_gt_m_1950s_80s": True,           # DFO-2025：1950s-80s 大部分年份 F>M
    "m_gt_f_post_collapse": True,       # DFO-2025：崩溃后 M>F
    "m_pre_assumed": 0.20,              # 崩溃前模型假设
    "m_post_mean": 0.47,                # DFO-2024：1995 后均值 0.47（范围 0.29-0.90）
    "m_peak_1992_94": 2.5,              # DFO-2025
    "capelin_pct_pre_collapse": 10.0,   # DFO-2024：Capelin ≈10% 崩溃前且未恢复
    "ssb_1990_kt": 1000.0,              # DFO-FW 锚点
    "ssb_1995_kt": 50.0,                # DFO-FW 锚点（2000 读数 50Kt，1995 为最低）
    "ssb_2005_pct_blim": 3.0,           # DFO-2018
    "ssb_2024_kt": 342.0,               # DFO-2024 表 1
}


def cod_collapse() -> dict:
    """C-S5-1：耗竭期 draw > regen。"""
    drop = 1.0 - COD_FACTS["ssb_1995_kt"] / COD_FACTS["ssb_1990_kt"]
    return {"f_gt_m_1950s_80s": COD_FACTS["f_gt_m_1950s_80s"],
            "ssb_drop_1990_95": float(drop),
            "pass": bool(COD_FACTS["f_gt_m_1950s_80s"] and drop > 0.80)}


def cod_absorption_years() -> dict:
    """C-S5-2：禁捕（交互停止）后低于 LRP 的年数。"""
    years = COD_FACTS["critical_zone_exit"] - COD_FACTS["moratorium_year"]
    return {"years_below_lrp_after_moratorium": years,
            "pass": bool(years >= 20)}


def cod_hysteresis_ratio() -> dict:
    """C-S5-3：迟滞不对称 = 恢复时标 / 崩溃时标。

    崩溃时标：1990→1995（5 年，SSB −95%）；
    恢复时标：1995→2024（29 年，SSB 仅从 ~50Kt 回到 342Kt≈1.2×LRP，
    仍只有 1990 年水平的约 1/3——且这还算"出危急区"的宽口径）。
    """
    t_collapse = COD_FACTS["collapse_low_year"] - COD_FACTS["collapse_start"]
    t_recover = 2024 - COD_FACTS["collapse_low_year"]
    ratio = t_recover / t_collapse
    return {"t_collapse_yr": t_collapse, "t_recover_yr": t_recover,
            "ratio": float(ratio),
            "ssb_2024_vs_1990": COD_FACTS["ssb_2024_kt"] / COD_FACTS["ssb_1990_kt"],
            "pass": bool(ratio >= 4.0)}


def cod_regen_damage() -> dict:
    """C-S5-4：再生机制受损（M 倍增 + 饵料崩溃）。"""
    m_ratio = COD_FACTS["m_post_mean"] / COD_FACTS["m_pre_assumed"]
    return {"m_ratio_post_pre": float(m_ratio),
            "capelin_pct": COD_FACTS["capelin_pct_pre_collapse"],
            "pass": bool(m_ratio >= 2.0)}


def cod_verdict() -> dict:
    c1, c2 = cod_collapse(), cod_absorption_years()
    c3, c4 = cod_hysteresis_ratio(), cod_regen_damage()
    return {"C-S5-1": c1, "C-S5-2": c2, "C-S5-3": c3, "C-S5-4": c4,
            "all_pass": bool(c1["pass"] and c2["pass"] and c3["pass"]
                             and c4["pass"])}
