"""XS 推导链溯源与冻结窗闸门（C1.7 / C1.8）的测试（仅标准库）。

锁定的是**结构预言**，不锁具体数字：
- 真实登记表当前必须无 error（已登记的无根断言只允许 warn）；
- 闸门必须能抓住它该抓的每类违规（用合成坏样本验证，不改真实文件）。
"""

import copy
import json
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from src import xs_formal_checks as xfc  # noqa: E402
from src import xs_trace_gate as xtg  # noqa: E402


@pytest.fixture(scope="module")
def symbol_reg():
    return xfc.load_json(xfc.SYMBOL_REGISTRY)


@pytest.fixture(scope="module")
def trace_gaps():
    return xfc.load_json(xtg.TRACE_GAPS).get("gaps", {})


@pytest.fixture(scope="module")
def frozen_manifest():
    return xfc.load_json(xtg.FROZEN_MANIFEST)


def codes(issues):
    return [i["code"] for i in issues]


# --------------------------------------------------------------------------
# 真实登记表
# --------------------------------------------------------------------------

def test_real_registry_has_no_trace_errors(symbol_reg, trace_gaps):
    issues = xtg.check_derivation_trace(symbol_reg, trace_gaps)
    errs = [i for i in issues if i["level"] == "error"]
    assert not errs, "真实登记表不应有推导链 error：%s" % errs


def test_real_frozen_window_intact(frozen_manifest):
    issues = xtg.check_frozen_window(frozen_manifest)
    errs = [i for i in issues if i["level"] == "error"]
    assert not errs, "冻结窗被破坏：%s" % errs


def test_frozen_manifest_nonempty(frozen_manifest):
    assert frozen_manifest.get("frozen"), "冻结清单不得为空"


def test_every_gap_entry_has_reason_and_todo(trace_gaps):
    for name, entry in trace_gaps.items():
        assert str(entry.get("reason", "")).strip(), "无根断言 %s 缺 reason" % name
        assert str(entry.get("todo", "")).strip(), "无根断言 %s 缺 todo" % name


def test_whitelist_is_not_stale(symbol_reg, trace_gaps):
    """已登记为无根的记号，必须真的不可达根（否则 G004 要求移除）。"""
    issues = xtg.check_derivation_trace(symbol_reg, trace_gaps)
    stale = [i for i in issues if i["code"] == "G004"]
    assert not stale, "白名单过期（已登记但其实可达根）：%s" % stale


def test_main_checker_includes_trace_gate():
    report = xfc.run_all()
    assert report["trace_gate"]["enabled"], "主校验器必须包含推导链闸门"
    assert not report["errors"], "合并后不应有 error：%s" % report["errors"]


# --------------------------------------------------------------------------
# 合成坏样本：推导链
# --------------------------------------------------------------------------

def _mini_reg(symbols):
    return {"version": "test", "symbols": symbols}


def _sym(name, kind, derives):
    return {"symbol": name, "name_zh": name, "layer": "L2", "kind": kind,
            "unit": "bit", "definition": "t", "derives_from": derives,
            "defined_in": ["test"], "identifiability": None,
            "allowed_readings": ["a"], "forbidden_readings": ["b"]}


def test_orphan_symbol_is_caught():
    reg = _mini_reg([_sym("Ax", "primitive", []),
                     _sym("orphan", "observable", [])])
    issues = xtg.check_derivation_trace(reg, {})
    assert "G001" in codes(issues)


def test_symbol_derived_from_assumption_has_root():
    """根包括「已声明假设」——采自 C1.7 原文「回溯到公理与已声明假设」。"""
    reg = _mini_reg([_sym("A_D", "assumption", []),
                     _sym("a", "parameter", ["A_D"])])
    issues = xtg.check_derivation_trace(reg, {})
    assert not [i for i in issues if i["level"] == "error"]


def test_chain_to_primitive_has_root():
    reg = _mini_reg([_sym("phi", "primitive", []),
                     _sym("P", "observable", ["phi"]),
                     _sym("U", "observable", ["P"])])
    issues = xtg.check_derivation_trace(reg, {})
    assert not [i for i in issues if i["level"] == "error"]


def test_declared_gap_becomes_warning_not_error():
    reg = _mini_reg([_sym("Ax", "primitive", []),
                     _sym("eta", "observable", [])])
    gaps = {"eta": {"reason": "约定性输入", "todo": "B-B44 解决后移除"}}
    issues = xtg.check_derivation_trace(reg, gaps)
    errs = [i for i in issues if i["level"] == "error"]
    warns = [i for i in issues if i["level"] == "warn" and i["code"] == "G001"]
    assert not errs
    assert warns, "已登记的无根断言应报 warn 而非 error"


def test_stale_whitelist_entry_is_caught():
    reg = _mini_reg([_sym("Ax", "primitive", []),
                     _sym("good", "observable", ["Ax"])])
    gaps = {"good": {"reason": "x", "todo": "y"}}
    issues = xtg.check_derivation_trace(reg, gaps)
    assert "G004" in codes(issues)


def test_gap_entry_without_reason_is_error():
    reg = _mini_reg([_sym("Ax", "primitive", []),
                     _sym("bad", "observable", [])])
    gaps = {"bad": {"reason": "  ", "todo": "y"}}
    issues = xtg.check_derivation_trace(reg, gaps)
    assert "G005" in codes(issues)


def test_no_roots_at_all_is_error():
    reg = _mini_reg([_sym("x", "observable", [])])
    issues = xtg.check_derivation_trace(reg, {})
    assert "G000" in codes(issues)


def test_reachable_roots_includes_self_and_transitive():
    reg = _mini_reg([_sym("phi", "primitive", []),
                     _sym("P", "observable", ["phi"]),
                     _sym("U", "observable", ["P"])])
    assert xtg.reachable_roots(reg, "U") == {"U", "P", "phi"}


# --------------------------------------------------------------------------
# 合成坏样本：冻结窗
# --------------------------------------------------------------------------

def _write(tmp_path, rel, text):
    path = tmp_path / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_frozen_window_detects_modification(tmp_path):
    _write(tmp_path, "frozen.md", "original text")
    manifest = {"frozen": [{"path": "frozen.md",
                            "sha256": xtg.sha256_of(str(tmp_path / "frozen.md")),
                            "reason": "test"}]}
    assert not [i for i in xtg.check_frozen_window(manifest, str(tmp_path))
                if i["level"] == "error"]
    _write(tmp_path, "frozen.md", "tampered text")
    issues = xtg.check_frozen_window(manifest, str(tmp_path))
    assert "G002" in codes(issues)


def test_frozen_window_detects_missing_file(tmp_path):
    manifest = {"frozen": [{"path": "gone.md", "sha256": "0" * 64,
                            "reason": "test"}]}
    issues = xtg.check_frozen_window(manifest, str(tmp_path))
    assert "G003" in codes(issues)


def test_missing_manifest_is_error(tmp_path):
    """冻结清单不存在 ⟹ C1.8 只是口头声明，必须是 error。"""
    issues = xtg.check_frozen_window(None, str(tmp_path),
                                     manifest_path=str(tmp_path / "nope.json"))
    assert "G006" in codes(issues)


def test_empty_manifest_is_error(tmp_path):
    issues = xtg.check_frozen_window({"frozen": []}, str(tmp_path))
    assert "G006" in codes(issues)


def test_hash_ignores_crlf_vs_lf(tmp_path):
    """换行符差异不得被判为"冻结窗被破坏"——否则跨平台克隆会全库假阳性。"""
    _write(tmp_path, "lf.md", "line1\nline2\n")
    (tmp_path / "crlf.md").write_bytes(b"line1\r\nline2\r\n")
    assert xtg.sha256_of(str(tmp_path / "lf.md")) == \
        xtg.sha256_of(str(tmp_path / "crlf.md"))


def test_hash_still_catches_content_change(tmp_path):
    """归一化只忽略换行，文字改动必须仍被抓到。"""
    _write(tmp_path, "a.md", "line1\nline2\n")
    _write(tmp_path, "b.md", "line1\nline2 changed\n")
    assert xtg.sha256_of(str(tmp_path / "a.md")) != \
        xtg.sha256_of(str(tmp_path / "b.md"))


def test_frozen_entry_without_sha_is_error(tmp_path):
    _write(tmp_path, "f.md", "x")
    manifest = {"frozen": [{"path": "f.md"}]}
    issues = xtg.check_frozen_window(manifest, str(tmp_path))
    assert "G007" in codes(issues)


def test_init_frozen_manifest_refreshes_hashes(tmp_path):
    _write(tmp_path, "f.md", "content")
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps({"frozen": [{"path": "f.md", "sha256": "",
                                            "reason": "r"}]}),
                    encoding="utf-8")
    manifest = xtg.init_frozen_manifest(str(path), str(tmp_path))
    assert manifest["frozen"][0]["sha256"] == xtg.sha256_of(str(tmp_path / "f.md"))
