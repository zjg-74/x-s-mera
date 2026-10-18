"""
R3-v2 案例定量升级（S3 流感 / S5 鳕鱼）的单元测试
====================================================

锚点：公开来源数据（CDC-VE / WHO-COMP / WHO-QA / DFO SAR 等）+
预注册判据 CASE_CRITERIA。数据完整性测试保证表中的数就是来源里的数。
"""

from __future__ import annotations

import pytest

from src.r3_case_upgrade import (
    FLU_VE, FLU_H3N2, DRIFT_SEASON, FLU_2014_TIMELINE,
    COD_SSB_ANCHORS, COD_SSB_RECENT, COD_FACTS, CASE_CRITERIA,
    flu_strain_update_rate, flu_ve_stats, flu_2014_window, flu_verdict,
    cod_collapse, cod_absorption_years, cod_hysteresis_ratio,
    cod_regen_damage, cod_verdict,
)


def test_criteria_preregistered():
    assert set(CASE_CRITERIA) >= {"C-S3-1", "C-S3-2", "C-S3-3",
                                  "C-S5-1", "C-S5-2", "C-S5-3", "C-S5-4"}
    assert "预注册" in CASE_CRITERIA["note"]


# ---------- 数据完整性（数=来源里的数；防回写漂移） ----------

def test_flu_ve_data_integrity():
    assert len(FLU_VE) == 22                      # 2004-05 … 2025-26
    assert FLU_VE["2020-21"] is None              # CDC 缺测季
    assert FLU_VE["2014-15"] == 19.0              # CDC 最终值
    assert FLU_VE["2010-11"] == 60.0
    assert FLU_VE["2024-25"] == 56.0


def test_flu_strain_data_integrity():
    assert FLU_H3N2["2014-15"] == "A/Texas/50/2012"
    assert FLU_H3N2["2015-16"] == "A/Switzerland/9715293/2013"
    assert FLU_H3N2["2006-07"] == FLU_H3N2["2007-08"] == "A/Wisconsin/67/2005"
    assert FLU_H3N2["2022-23"] == FLU_H3N2["2023-24"] == "A/Darwin/9/2021"


def test_flu_timeline_order():
    """2014-15 窗口：决策（2 月）→ 漂移检出（3 月）→ 季中失/改荐（9 月）→ 当季 VE。"""
    labels = [t[0] for t in FLU_2014_TIMELINE]
    assert labels == ["2014-02", "2014-03", "2014-09", "2014-09", "2014-15 季"]


def test_cod_facts_integrity():
    assert COD_FACTS["moratorium_year"] == 1992
    assert COD_FACTS["m_pre_assumed"] == 0.20
    assert COD_FACTS["m_post_mean"] == 0.47
    assert COD_SSB_RECENT[2024] == 342.0          # SAR 2024/049 表 1
    assert COD_SSB_ANCHORS[1990] == 1000.0


# ---------- 判据 ----------

def test_s3_c1_mismatch_penalty():
    r = flu_ve_stats()
    assert r["drift_ve"] < 0.5 * r["others_mean"]
    assert r["is_series_min_tight"] is True
    assert r["z_score"] < -2.0
    assert r["pass"] is True


def test_s3_c2_instantaneous_reading_is_operative():
    r = flu_strain_update_rate()
    assert r["n_transitions"] == 21               # 2004-05 → 2025-26
    assert r["update_rate"] >= 0.50
    assert r["pass"] is True


def test_s3_c3_window_breach():
    r = flu_2014_window()
    assert r["incremental_compliant"] is True     # 决策时点 >80% 匹配
    assert r["absolute_breach"] is True           # 季中 ≤50% 抑制
    assert r["correction_next_meeting"] is True   # 瞬时读法纠正
    assert r["pass"] is True


def test_s3_all_pass():
    assert flu_verdict()["all_pass"] is True


def test_s5_c1_draw_gt_regen():
    r = cod_collapse()
    assert r["f_gt_m_1950s_80s"] is True
    assert r["ssb_drop_1990_95"] > 0.80
    assert r["pass"] is True


def test_s5_c2_absorption():
    r = cod_absorption_years()
    assert r["years_below_lrp_after_moratorium"] == 24
    assert r["pass"] is True


def test_s5_c3_hysteresis():
    r = cod_hysteresis_ratio()
    assert r["ratio"] == pytest.approx(5.8)
    assert r["ssb_2024_vs_1990"] < 0.5            # 恢复 29 年仍不足崩溃前一半
    assert r["pass"] is True


def test_s5_c4_regen_damage():
    r = cod_regen_damage()
    assert r["m_ratio_post_pre"] >= 2.0
    assert r["capelin_pct"] <= 15.0
    assert r["pass"] is True


def test_s5_all_pass():
    assert cod_verdict()["all_pass"] is True
