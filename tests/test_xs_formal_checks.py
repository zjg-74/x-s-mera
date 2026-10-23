"""XS 框架形式校验器的测试（仅标准库）。

锁定的是**结构预言**，不锁具体数字：
- 登记表本身必须自洽（无错误）；
- 校验器必须能抓住它该抓的每类违规（用合成坏样本验证，不改真实登记表）。
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


@pytest.fixture(scope="module")
def symbol_reg():
    return xfc.load_json(xfc.SYMBOL_REGISTRY)


@pytest.fixture(scope="module")
def distinction_reg():
    return xfc.load_json(xfc.DISTINCTION_REGISTRY)


@pytest.fixture(scope="module")
def report():
    return xfc.run_all()


# ---------------------------------------------------------------- 真实登记表

def test_真实登记表无错误(report):
    assert report["errors"] == [], report["errors"]


def test_真实登记表有实质内容(report):
    assert report["symbols"] >= 30
    assert report["distinctions"] >= 5
    assert report["legacy_conflicts"] >= 2


def test_三层各有记号(symbol_reg):
    layers = {}
    for sym in symbol_reg["symbols"]:
        layers.setdefault(sym["layer"], []).append(sym["symbol"])
    assert set(layers) == {"L0", "L1", "L2"}
    for layer, names in layers.items():
        assert len(names) >= 5, (layer, names)


def test_记号唯一(symbol_reg):
    names = [s["symbol"] for s in symbol_reg["symbols"]]
    assert len(names) == len(set(names))


def test_承重假设只有一个且被标注(symbol_reg):
    assumptions = [s for s in symbol_reg["symbols"] if s["kind"] == "assumption"]
    assert len(assumptions) == 1
    assert "未获现实检验" in assumptions[0]["definition"]
    assert "不能读作「已被证明」" in assumptions[0]["forbidden_readings"]


def test_两个历史冲突都被登记(symbol_reg):
    titles = [c["title"] for c in symbol_reg["legacy_conflicts"]]
    assert any("Ψ" in t for t in titles)
    assert any("Q" in t for t in titles)
    for con in symbol_reg["legacy_conflicts"]:
        assert con["status"] in ("resolved", "annotated", "open")
        assert con["resolution"].strip()


def test_区分登记含用户点名的两条区分(distinction_reg):
    names = " ".join(d["name"] for d in distinction_reg["distinctions"])
    assert "读得出来" in names
    assert "同一尺度下的耦合" in names


def test_每条区分都有非平凡实例与外部必要性(distinction_reg):
    for dist in distinction_reg["distinctions"]:
        inst = dist["nontrivial_instance"]["instances"]
        assert len(inst) >= 2
        assert len(set(inst)) == len(inst)
        assert dist["nontrivial_instance"]["why_differ"].strip()
        assert dist["external_necessity"].strip()
        assert dist["decision_procedure"]["executable"] is True
        assert dist["decision_procedure"]["steps"]
        assert dist["decision_procedure"]["procedure_ref"].strip()


def test_可辨识级别覆盖四个等级(symbol_reg):
    used = set(s["identifiability"] for s in symbol_reg["symbols"]
               if s["kind"] in ("observable", "parameter"))
    assert used == {"A", "B", "C", "D"}


def test_D级记号被标为伪目标(symbol_reg):
    for sym in symbol_reg["symbols"]:
        if sym["identifiability"] == "D":
            assert any("伪目标" in r for r in sym["forbidden_readings"])


# ---------------------------------------------------------------- 校验器行为

def _codes(issues):
    return set(i["code"] for i in issues)


def test_抓重复记号(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    bad["symbols"].append(copy.deepcopy(bad["symbols"][0]))
    assert "S012" in _codes(xfc.check_symbols(bad))


def test_抓非法层(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    bad["symbols"][0]["layer"] = "L9"
    assert "S002" in _codes(xfc.check_symbols(bad))


def test_抓缺量纲(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    for sym in bad["symbols"]:
        if sym["kind"] == "observable":
            sym["unit"] = None
            break
    assert "S004" in _codes(xfc.check_symbols(bad))


def test_抓未声明可辨识级别(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    for sym in bad["symbols"]:
        if sym["kind"] == "observable":
            sym["identifiability"] = "Z"
            break
    assert "S008" in _codes(xfc.check_symbols(bad))


def test_抓允许与禁止读法相交(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    bad["symbols"][0]["forbidden_readings"].append(
        bad["symbols"][0]["allowed_readings"][0])
    assert "S011" in _codes(xfc.check_symbols(bad))


def test_抓未登记记号引用(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    bad["symbols"][0]["derives_from"] = ["NOT_REGISTERED"]
    assert "S006" in _codes(xfc.check_symbols(bad))


def test_抓依赖环(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    by_name = {s["symbol"]: s for s in bad["symbols"]}
    a, b = by_name["E"], by_name["U"]
    a["derives_from"] = ["U"]
    b["derives_from"] = ["E"]
    assert "S013" in _codes(xfc.check_symbols(bad))


def test_抓自指(symbol_reg):
    bad = copy.deepcopy(symbol_reg)
    by_name = {s["symbol"]: s for s in bad["symbols"]}
    by_name["E"]["derives_from"] = ["E"]
    assert "S007" in _codes(xfc.check_symbols(bad))


def test_抓区分缺少外部必要性(distinction_reg, symbol_reg):
    bad = copy.deepcopy(distinction_reg)
    bad["distinctions"][0]["external_necessity"] = "  "
    names = set(s["symbol"] for s in symbol_reg["symbols"])
    assert "D005" in _codes(xfc.check_distinctions(bad, names))


def test_抓区分只给一个实例(distinction_reg, symbol_reg):
    bad = copy.deepcopy(distinction_reg)
    bad["distinctions"][0]["nontrivial_instance"]["instances"] = ["只有一个"]
    names = set(s["symbol"] for s in symbol_reg["symbols"])
    assert "D003" in _codes(xfc.check_distinctions(bad, names))


def test_抓不可执行却标为established(distinction_reg, symbol_reg):
    bad = copy.deepcopy(distinction_reg)
    bad["distinctions"][0]["decision_procedure"]["executable"] = False
    names = set(s["symbol"] for s in symbol_reg["symbols"])
    assert "D004" in _codes(xfc.check_distinctions(bad, names))


def test_抓区分侧使用未登记记号(distinction_reg, symbol_reg):
    bad = copy.deepcopy(distinction_reg)
    bad["distinctions"][0]["discriminates"][0]["left"] = "NOT_A_SYMBOL"
    names = set(s["symbol"] for s in symbol_reg["symbols"])
    assert "D006" in _codes(xfc.check_distinctions(bad, names))


def test_报告格式可读(report):
    text = xfc.format_report(report)
    assert "XS 框架形式校验" in text
    assert "记号登记" in text


def test_登记表是合法json():
    for path in (xfc.SYMBOL_REGISTRY, xfc.DISTINCTION_REGISTRY):
        with open(path, "r", encoding="utf-8") as fh:
            json.load(fh)
