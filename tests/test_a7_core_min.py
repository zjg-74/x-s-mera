#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A7 —— 核心最小化的结构性测试（锁结构，不锁数字）

对应：`experiments/a7_core_min/`（a7_core_min.py + core_table.json）、
18-v3.3-formal-framework-charter.md 待办 A7。

锁的结构预言：
  * 传递闭包自反、单调、传递（图论层面，与记号表内容无关）
  * 极小根子集 = 可达集 ∩ 根集（依赖结构上可计算）
  * 每个非根记号（约定输入除外）都能沿 derives_from 到达某个根 ⟹ 核心+规则⟹全部记号
  * 孤儿记号（非根、无依赖）必须被判为"不可达"（检查器有牙齿，不是一律通过）
  * η 被登记为约定输入（非根、无 derives_from）：与 B-B44 一致
  * 导出表覆盖全部非根记号

不锁：核心记号的具体条数、压缩率数值、根记号的具体名单。
（名单会随记号表修订而变，锁了就是锁数字。）

诚实边界：A7 是**依赖结构层面**的最小化，不是概念独立性主张；本测试只锁前者。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "experiments", "a7_core_min"))

import a7_core_min as R  # noqa: E402

# 假图（不依赖真实记号表，用于测图算法本身）
_FAKE = {
    "by_name": {"A": {}, "B": {}, "C": {}, "ORPHAN": {}},
    "edges": {"A": set(), "B": {"A"}, "C": {"B"}, "ORPHAN": set()},
}


def test_module_imports():
    for fn in ("load_registry", "build_graph", "transitive_deps",
               "minimal_roots_for", "main"):
        assert hasattr(R, fn), fn


def test_transitive_deps_is_reflexive():
    by_name, edges = _FAKE["by_name"], _FAKE["edges"]
    for n in by_name:
        assert n in R.transitive_deps(n, edges, by_name)


def test_transitive_deps_is_transitive():
    by_name, edges = _FAKE["by_name"], _FAKE["edges"]
    reach = R.transitive_deps("C", edges, by_name)
    assert {"C", "B", "A"} <= reach
    assert R.transitive_deps("A", edges, by_name) == {"A"}


def test_transitive_deps_tolerates_dangling_reference():
    """derives_from 指向不存在的记号时，不应崩溃也不应被计入可达集。"""
    by_name = {"X": {}}
    edges = {"X": {"GHOST"}}
    assert R.transitive_deps("X", edges, by_name) == {"X"}


def test_minimal_roots_equals_reachable_roots():
    by_name, edges = _FAKE["by_name"], _FAKE["edges"]
    roots = {"A"}
    assert R.minimal_roots_for("C", edges, by_name, roots) == ["A"]
    assert R.minimal_roots_for("B", edges, by_name, roots) == ["A"]
    assert R.minimal_roots_for("A", edges, by_name, roots) == ["A"]
    assert R.minimal_roots_for("ORPHAN", edges, by_name, roots) == []


def test_orphan_is_flagged_unreachable_has_teeth():
    """孤儿记号必须被判为不可达 —— 否则展开校验就是一律通过的空转。"""
    by_name, edges = _FAKE["by_name"], _FAKE["edges"]
    roots = {"A"}
    gaps = {s for s in by_name
            if s not in roots and not (R.transitive_deps(s, edges, by_name) & roots)}
    assert "ORPHAN" in gaps
    assert "B" not in gaps and "C" not in gaps


def test_registry_loads_with_required_fields():
    reg = R.load_registry()
    assert reg["symbols"]
    for s in reg["symbols"]:
        assert "symbol" in s and "kind" in s


def test_build_graph_covers_all_symbols():
    reg = R.load_registry()
    by_name, edges = R.build_graph(reg)
    assert len(by_name) == len(reg["symbols"])
    assert set(edges) == set(by_name)


def test_every_non_root_symbol_reaches_a_root():
    """核心 + 规则 ⟹ 全部记号（约定输入除外，它本就不该有依赖）。"""
    reg = R.load_registry()
    by_name, edges = R.build_graph(reg)
    roots = {s["symbol"] for s in reg["symbols"] if s.get("kind") in R.ROOT_KINDS}
    conventions = {s["symbol"] for s in reg["symbols"]
                   if s.get("kind") not in R.ROOT_KINDS
                   and not (s.get("derives_from") or [])}
    assert roots
    for sym in by_name:
        if sym in roots or sym in conventions:
            continue
        assert R.transitive_deps(sym, edges, by_name) & roots, sym


def test_eta_is_registered_as_convention_input():
    """η 是约定性输入（B-B44）：非根、无 derives_from，有意为之而非遗漏。"""
    reg = R.load_registry()
    by_name, edges = R.build_graph(reg)
    assert "eta" in by_name
    eta = by_name["eta"]
    assert eta.get("kind") not in R.ROOT_KINDS
    assert not (eta.get("derives_from") or [])
    assert not edges["eta"]


def test_derived_table_covers_all_non_root_symbols(tmp_path, monkeypatch):
    """导出表行数 = 非根记号数，且每个非约定输入都有非空最小根集。"""
    monkeypatch.setattr(R, "HERE", str(tmp_path))
    assert R.main() == 0
    import json
    with open(tmp_path / "core_table.json", "r", encoding="utf-8") as fh:
        table = json.load(fh)

    reg = R.load_registry()
    non_root = [s["symbol"] for s in reg["symbols"] if s.get("kind") not in R.ROOT_KINDS]
    assert len(table["derived_table"]) == len(non_root)
    assert table["gaps"] == []
    for row in table["derived_table"]:
        assert row["type"] != "不可达", row
        if row["type"] == "导出":
            assert row["min_roots"], row
    # 核心 = 根 ∪ 约定输入，是全部记号的子集
    core = set(table["core"]["roots"]) | set(table["core"]["conventions"])
    assert core <= {s["symbol"] for s in reg["symbols"]}
    assert table["core"]["core_count"] < table["core"]["total_symbols"]
