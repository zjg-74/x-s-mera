"""XS 框架的推导链溯源与冻结窗闸门（仅标准库）。

对应合法性判据一的两条此前**未机检**的项（见 18-v3.3-formal-framework-charter.md §3.1）：

  - **C1.7 推导链可追**：每条结论能回溯到公理与已声明假设。
    机检形式：每个「非根记号」必须沿 `derives_from` 可达某个根记号。
    **根 = kind ∈ {primitive, assumption}**（公理级记号与显式声明的承重假设）。
    —— 采自 C1.7 原文「回溯到公理**与已声明假设**」，不是只有公理。

  - **C1.8 冻结窗未被破坏**：被冻结文本只带标注、不被改写。
    机检形式：`frozen_manifest.json` 登记冻结文件的 sha256，逐条比对。

白名单机制与其风险（必须写明）：
    `trace_gaps.json` 登记**已知的无根断言**。它会让校验通过，因此是一个滑坡点：
    往里加一条就能让一个洞消失。防线有三条：
      ① 每条必须带 `reason` 与 `todo`；
      ② 校验器对「已登记但其实可达」的条目报 G004（白名单过期，必须清理）；
      ③ 白名单条目数只减不增——新增需显式写入本文件并说明。

用法：
    python src/xs_trace_gate.py              # 打印报告，error > 0 时退出码 1
    python src/xs_trace_gate.py --json       # 输出 JSON
    python src/xs_trace_gate.py --init-frozen   # 生成/刷新冻结清单哈希
"""

from __future__ import annotations

import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DOCS = os.path.join(REPO, "docs", "xs-framework")

SYMBOL_REGISTRY = os.path.join(DOCS, "symbol_registry.json")
TRACE_GAPS = os.path.join(DOCS, "trace_gaps.json")
FROZEN_MANIFEST = os.path.join(DOCS, "frozen_manifest.json")

ROOT_KINDS = ("primitive", "assumption")


def _issue(level, code, where, message):
    return {"level": level, "code": code, "where": where, "message": message}


def load_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def sha256_of(path, normalize_newlines=True):
    """文件哈希。**默认先把 CRLF 归一化为 LF 再算**。

    为什么要归一化：git 在 Windows 上按 `core.autocrlf` 检出时会把 LF 写成 CRLF，
    同一份文本在不同平台的字节流不同 ⟹ 直接哈希会产生**假阳性**的"冻结窗被破坏"。
    （v3.4 首次生成清单后，在干净克隆里实测报出 9 条 G002，全部是换行差异。）
    归一化后，**文字内容的任何改动仍然会被抓到**，只有换行符差异被忽略——
    这是正确的取舍：我们要防的是**改写**，不是**搬运**。
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            if normalize_newlines:
                chunk = chunk.replace(b"\r\n", b"\n")
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# C1.7 推导链溯源
# --------------------------------------------------------------------------

def reachable_roots(reg, symbol):
    """从 symbol 出发沿 derives_from 可达的全部记号（含自身）。"""
    by_name = {s.get("symbol"): s for s in reg.get("symbols", [])}
    seen = set()
    stack = [symbol]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        for ref in by_name.get(node, {}).get("derives_from", []) or []:
            if ref in by_name:
                stack.append(ref)
    return seen


def check_derivation_trace(reg, gaps=None):
    """每个非根记号必须可达某个根记号（primitive 或 assumption）。"""
    issues = []
    if gaps is None:
        gaps = {}
    by_name = {s.get("symbol"): s for s in reg.get("symbols", [])}
    roots = {s.get("symbol") for s in reg.get("symbols", [])
             if s.get("kind") in ROOT_KINDS}

    if not roots:
        issues.append(_issue("error", "G000", "symbol_registry",
                             "没有任何根记号（primitive/assumption）⟹ 推导链无从追溯"))
        return issues

    for sym in reg.get("symbols", []):
        name = sym.get("symbol")
        where = "symbol(%s)" % name
        if sym.get("kind") in ROOT_KINDS:
            continue
        reach = reachable_roots(reg, name)
        if reach & roots:
            if name in gaps:
                issues.append(_issue("warn", "G004", where,
                                     "已登记为无根断言，但实际可达根 ⟹ 白名单过期，"
                                     "请从 trace_gaps.json 移除"))
            continue
        # 不可达
        if name in gaps:
            entry = gaps[name] or {}
            if not str(entry.get("reason", "")).strip():
                issues.append(_issue("error", "G005", where,
                                     "无根断言白名单条目必须给出 reason"))
            if not str(entry.get("todo", "")).strip():
                issues.append(_issue("error", "G005", where,
                                     "无根断言白名单条目必须给出 todo"))
            issues.append(_issue("warn", "G001", where,
                                 "无根断言（已登记）：%s"
                                 % str(entry.get("reason", "")).strip()))
        else:
            issues.append(_issue("error", "G001", where,
                                 "推导链断裂：无法回溯到公理或已声明假设，"
                                 "且未在 trace_gaps.json 登记"))
        _ = by_name  # noqa: F841  (保持可读性)
    return issues


# --------------------------------------------------------------------------
# C1.8 冻结窗
# --------------------------------------------------------------------------

def check_frozen_window(manifest=None, base=REPO, manifest_path=FROZEN_MANIFEST):
    issues = []
    if manifest is None:
        if not os.path.exists(manifest_path):
            issues.append(_issue("error", "G006", "frozen_manifest.json",
                                 "冻结清单不存在 ⟹ C1.8（冻结窗未被破坏）只是口头声明，"
                                 "不是可核验的框架性质"))
            return issues
        manifest = load_json(manifest_path)

    entries = manifest.get("frozen", [])
    if not entries:
        issues.append(_issue("error", "G006", "frozen_manifest.json",
                             "冻结清单为空 ⟹ 没有任何文本被冻结"))
        return issues

    for entry in entries:
        rel = entry.get("path", "")
        where = "frozen(%s)" % rel
        if not rel:
            issues.append(_issue("error", "G007", where, "冻结条目缺少 path"))
            continue
        path = os.path.join(base, rel)
        if not os.path.exists(path):
            issues.append(_issue("error", "G003", where,
                                 "被冻结文件不存在（被删除或移动）"))
            continue
        actual = sha256_of(path)
        expected = entry.get("sha256", "")
        if not expected:
            issues.append(_issue("error", "G007", where, "冻结条目缺少 sha256"))
            continue
        if actual != expected:
            issues.append(_issue("error", "G002", where,
                                 "冻结窗被破坏：内容哈希与冻结清单不一致 "
                                 "（期望 %s…，实为 %s…）"
                                 % (expected[:12], actual[:12])))
        if not str(entry.get("reason", "")).strip():
            issues.append(_issue("warn", "G008", where,
                                 "冻结条目未说明冻结理由"))
    return issues


def init_frozen_manifest(path=FROZEN_MANIFEST, base=REPO):
    """生成或刷新冻结清单。只刷新哈希，不改 path/reason 列表。"""
    if os.path.exists(path):
        manifest = load_json(path)
    else:
        manifest = {
            "version": "1.0",
            "framework": "XS v3.4",
            "purpose": "C1.8 冻结窗：被冻结文本只带标注、不被改写（U4）",
            "note": "新增冻结条目须说明 reason；刷新哈希用 --init-frozen",
            "frozen": [],
        }
    for entry in manifest.get("frozen", []):
        full = os.path.join(base, entry.get("path", ""))
        if os.path.exists(full):
            entry["sha256"] = sha256_of(full)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    return manifest


# --------------------------------------------------------------------------
# 总入口
# --------------------------------------------------------------------------

def run_all(symbol_path=SYMBOL_REGISTRY, gaps_path=TRACE_GAPS,
            manifest_path=FROZEN_MANIFEST):
    reg = load_json(symbol_path)
    gaps = {}
    if os.path.exists(gaps_path):
        gaps = load_json(gaps_path).get("gaps", {})
    issues = []
    issues += check_derivation_trace(reg, gaps)
    issues += check_frozen_window(
        load_json(manifest_path) if os.path.exists(manifest_path) else None)
    errors = [i for i in issues if i["level"] == "error"]
    warns = [i for i in issues if i["level"] == "warn"]
    return {
        "roots": sum(1 for s in reg.get("symbols", [])
                     if s.get("kind") in ROOT_KINDS),
        "trace_gaps_declared": len(gaps),
        "errors": errors,
        "warnings": warns,
    }


def format_report(report):
    lines = []
    lines.append("XS 推导链溯源与冻结窗闸门（C1.7 / C1.8）")
    lines.append("  根记号（公理+已声明假设）：%d 条" % report["roots"])
    lines.append("  已登记无根断言：%d 条" % report["trace_gaps_declared"])
    lines.append("  错误：%d    警告：%d"
                 % (len(report["errors"]), len(report["warnings"])))
    for issue in report["errors"] + report["warnings"]:
        lines.append("  [%s] %s %s — %s"
                     % (issue["level"], issue["code"], issue["where"],
                        issue["message"]))
    return "\n".join(lines)


def main(argv):
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8")
            except (ValueError, OSError):
                pass
    if "--init-frozen" in argv:
        manifest = init_frozen_manifest()
        print("冻结清单已刷新：%d 条" % len(manifest.get("frozen", [])))
        return 0
    report = run_all()
    if "--json" in argv:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        print(format_report(report))
    return 1 if report["errors"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
