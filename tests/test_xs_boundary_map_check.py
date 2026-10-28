# -*- coding: utf-8 -*-
"""xs_boundary_map_check 的结构性测试（仅标准库）。

锁：边界地图机检（A10 残留落地）的结构预言——必填字段、id 唯一、解决指向存在。
不锁具体状态词（枚举为软 warn），不锁文件内容数字。调底层函数、不写文件。
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import xs_boundary_map_check as bmc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DOCS = os.path.join(REPO, "docs", "xs-framework")


def _load():
    with open(os.path.join(DOCS, "boundary_map.json"), encoding="utf-8") as f:
        return json.load(f)


# ---- 顶层结构 ----
def test_top_level_required_fields():
    bm = _load()
    for f in bmc.TOP_LEVEL_REQUIRED:
        assert f in bm, f"缺顶层字段 {f}"


def test_version_nonempty():
    assert isinstance(_load()["version"], str) and _load()["version"].strip()


# ---- dimensions ----
def test_dimension_required_fields():
    for d in _load()["dimensions"]:
        for f in bmc.DIM_FIELDS:
            assert f in d, f"dimension 缺 {f}"


def test_dimension_id_unique():
    ids = [d["id"] for d in _load()["dimensions"]]
    assert len(ids) == len(set(ids))


# ---- r041_table ----
def test_r041_required_fields():
    for r in _load()["r041_table"]:
        for f in bmc.R041_FIELDS:
            assert f in r, f"r041 行缺 {f}"


def test_r041_level_unique():
    lvs = [r["level"] for r in _load()["r041_table"]]
    assert len(lvs) == len(set(lvs))


# ---- open_cells ----
def test_open_cell_required_fields():
    for c in _load()["open_cells"]:
        for f in bmc.OPEN_CELL_FIELDS:
            assert f in c, f"open_cell 缺 {f}"


def test_open_cell_id_unique():
    ids = [c["id"] for c in _load()["open_cells"]]
    assert len(ids) == len(set(ids))


# ---- gap_index ----
def test_gap_required_fields():
    for g in _load()["gap_index"]:
        for f in bmc.GAP_FIELDS:
            assert f in g, f"gap 缺 {f}"


def test_gap_resolved_in_points_to_existing_file():
    """已解决/已闭合格子的 resolved_in 必须指向仓库内真实文件。"""
    for g in _load()["gap_index"]:
        rp = g.get("resolved_in")
        if rp:
            assert os.path.exists(os.path.join(REPO, rp)), f"resolved_in 不存在: {rp}"


# ---- 校验器本身：干净登记件应零 error ----
def test_checker_zero_errors_on_current_map():
    report = bmc.run_all()
    assert report["error_count"] == 0, f"边界地图机检应零错误: {report['errors']}"
