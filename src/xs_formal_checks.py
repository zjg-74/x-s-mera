"""XS 框架的形式校验器（仅标准库）。

用途：把「内部一致性」与「术语明确性」从靠人读变成可机检。

对应合法性判据（见 docs/xs-framework/18-v3.3-formal-framework-charter.md）：
  - 判据一「内部一致性」→ 依赖无环、无自指、层归属唯一；
    以及 v3.4 新增的 C1.7（推导链可追）、C1.8（冻结窗未被破坏）→ 由 xs_trace_gate 提供；
  - 判据二「术语明确性」→ 记号唯一、量纲与层固定、允许读法与禁止读法不相交、
    遗留冲突必须被登记且给出处置；
  - 判据三「跨实例可复用性」→ 区分登记里的非平凡实例必须互不相同。

用法：
    python src/xs_formal_checks.py            # 打印报告，error > 0 时退出码 1
    python src/xs_formal_checks.py --json     # 输出 JSON
"""

from __future__ import annotations

import json
import os
import sys

try:  # 作为包被导入时（tests: from src import xs_formal_checks）
    from src import xs_trace_gate  # type: ignore
except ImportError:  # 作为脚本直接运行时（sys.path[0] 为 src/）
    import xs_trace_gate  # type: ignore

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DOCS = os.path.join(REPO, "docs", "xs-framework")
SYMBOL_REGISTRY = os.path.join(DOCS, "symbol_registry.json")
DISTINCTION_REGISTRY = os.path.join(DOCS, "distinction_registry.json")

LAYERS = ("L0", "L1", "L2")
KINDS = ("primitive", "derived", "observable", "parameter", "assumption",
         "procedure", "classification")
LEVELS = ("A", "B", "C", "D")
STATUSES = ("established", "provisional", "conjectured")
QUANTITY_KINDS = ("observable", "parameter")

REQUIRED_SYMBOL_FIELDS = (
    "symbol", "name_zh", "layer", "kind", "unit", "definition",
    "derives_from", "defined_in", "identifiability",
    "allowed_readings", "forbidden_readings",
)
REQUIRED_DISTINCTION_FIELDS = (
    "id", "name", "statement", "discriminates", "nontrivial_instance",
    "decision_procedure", "external_necessity", "status",
)


def _issue(level, code, where, message):
    return {"level": level, "code": code, "where": where, "message": message}


def load_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# 符号登记
# --------------------------------------------------------------------------

def check_symbols(reg):
    issues = []
    symbols = reg.get("symbols")
    if not isinstance(symbols, list) or not symbols:
        return [_issue("error", "S000", "symbol_registry", "symbols 必须是非空列表")]

    seen = {}
    for idx, sym in enumerate(symbols):
        where = "symbol[%d]" % idx
        name = sym.get("symbol", where)
        where = "symbol(%s)" % name

        for field in REQUIRED_SYMBOL_FIELDS:
            if field not in sym:
                issues.append(_issue("error", "S001", where,
                                     "缺少必填字段 %s" % field))

        if sym.get("symbol") in seen:
            issues.append(_issue("error", "S012", where,
                                 "记号重复登记（记号唯一性判据）"))
        seen[sym.get("symbol")] = idx

        if sym.get("layer") not in LAYERS:
            issues.append(_issue("error", "S002", where,
                                 "layer 必须是 %s 之一，实为 %r"
                                 % (list(LAYERS), sym.get("layer"))))
        if sym.get("kind") not in KINDS:
            issues.append(_issue("error", "S003", where,
                                 "kind 必须是 %s 之一，实为 %r"
                                 % (list(KINDS), sym.get("kind"))))

        kind = sym.get("kind")
        unit = sym.get("unit", None)
        if kind in QUANTITY_KINDS:
            if not isinstance(unit, str) or not unit.strip():
                issues.append(_issue("error", "S004", where,
                                     "观测量/参数必须声明量纲"))
        elif unit is not None and not isinstance(unit, str):
            issues.append(_issue("error", "S004", where,
                                 "unit 必须是字符串或 null"))

        derives = sym.get("derives_from", [])
        if not isinstance(derives, list):
            issues.append(_issue("error", "S005", where,
                                 "derives_from 必须是列表"))
            derives = []
        elif kind == "derived" and not derives:
            issues.append(_issue("error", "S005", where,
                                 "derived 记号必须声明 derives_from"))
        if sym.get("symbol") in derives:
            issues.append(_issue("error", "S007", where,
                                 "derives_from 含自指"))

        ident = sym.get("identifiability", None)
        if kind in QUANTITY_KINDS:
            if ident not in LEVELS:
                issues.append(_issue("error", "S008", where,
                                     "观测量/参数必须声明可辨识级别 %s，实为 %r"
                                     % (list(LEVELS), ident)))
        elif ident is not None:
            issues.append(_issue("warn", "S009", where,
                                 "非观测量/参数不应带可辨识级别"))

        allowed = sym.get("allowed_readings", [])
        forbidden = sym.get("forbidden_readings", [])
        if not isinstance(allowed, list) or not allowed:
            issues.append(_issue("error", "S010", where,
                                 "allowed_readings 必须是非空列表"))
            allowed = []
        if not isinstance(forbidden, list) or not forbidden:
            issues.append(_issue("error", "S010", where,
                                 "forbidden_readings 必须是非空列表"))
            forbidden = []
        overlap = sorted(set(allowed) & set(forbidden))
        if overlap:
            issues.append(_issue("error", "S011", where,
                                 "同一读法既被允许又被禁止：%s" % overlap))

    known = set(seen)
    for idx, sym in enumerate(symbols):
        where = "symbol(%s)" % sym.get("symbol", idx)
        for ref in sym.get("derives_from", []) or []:
            if ref not in known:
                issues.append(_issue("error", "S006", where,
                                     "derives_from 引用了未登记的记号 %r" % ref))

    for cycle in _find_cycles(symbols):
        issues.append(_issue("error", "S013", "symbol_registry",
                             "derives_from 存在环：%s" % " -> ".join(cycle)))

    for idx, con in enumerate(reg.get("legacy_conflicts", []) or []):
        where = "legacy_conflict(%s)" % con.get("id", idx)
        for field in ("id", "title", "symbols", "identified_in", "status",
                      "resolution"):
            if field not in con:
                issues.append(_issue("error", "S015", where,
                                     "遗留冲突缺少必填字段 %s" % field))
        if con.get("status") not in ("resolved", "annotated", "open"):
            issues.append(_issue("error", "S015", where,
                                 "遗留冲突 status 必须是 resolved/annotated/open"))
        for ref in con.get("symbols", []) or []:
            if ref not in known:
                issues.append(_issue("error", "S014", where,
                                     "遗留冲突引用了未登记的记号 %r" % ref))
    return issues


def _find_cycles(symbols):
    graph = {}
    for sym in symbols:
        graph[sym.get("symbol")] = list(sym.get("derives_from", []) or [])
    cycles = []
    state = {}
    stack = []

    def visit(node):
        if state.get(node) == "done":
            return
        if state.get(node) == "open":
            idx = stack.index(node)
            cycles.append(stack[idx:] + [node])
            return
        state[node] = "open"
        stack.append(node)
        for nxt in graph.get(node, []):
            if nxt in graph:
                visit(nxt)
        stack.pop()
        state[node] = "done"

    for node in graph:
        visit(node)
    return cycles


# --------------------------------------------------------------------------
# 区分登记
# --------------------------------------------------------------------------

def check_distinctions(reg, symbol_names):
    issues = []
    dists = reg.get("distinctions")
    if not isinstance(dists, list) or not dists:
        return [_issue("error", "D000", "distinction_registry",
                       "distinctions 必须是非空列表")]

    seen = set()
    for idx, dist in enumerate(dists):
        where = "distinction[%d]" % idx
        name = dist.get("id", where)
        where = "distinction(%s)" % name

        for field in REQUIRED_DISTINCTION_FIELDS:
            if field not in dist:
                issues.append(_issue("error", "D001", where,
                                     "缺少必填字段 %s" % field))
        if dist.get("id") in seen:
            issues.append(_issue("error", "D007", where, "区分 id 重复"))
        seen.add(dist.get("id"))

        if dist.get("status") not in STATUSES:
            issues.append(_issue("error", "D002", where,
                                 "status 必须是 %s 之一，实为 %r"
                                 % (list(STATUSES), dist.get("status"))))

        inst = dist.get("nontrivial_instance", {}) or {}
        cases = inst.get("instances", []) or []
        if len(cases) < 2:
            issues.append(_issue("error", "D003", where,
                                 "非平凡实例至少需要 2 个（否则该区分恒真或恒假）"))
        elif len(set(cases)) != len(cases):
            issues.append(_issue("error", "D003", where,
                                 "非平凡实例必须互不相同"))
        if not inst.get("why_differ", "").strip():
            issues.append(_issue("error", "D003", where,
                                 "非平凡实例必须说明「为什么它们给出不同答案」"))

        proc = dist.get("decision_procedure", {}) or {}
        executable = proc.get("executable")
        if not isinstance(executable, bool):
            issues.append(_issue("error", "D004", where,
                                 "decision_procedure.executable 必须是布尔值"))
        elif executable:
            steps = proc.get("steps", []) or []
            if not steps:
                issues.append(_issue("error", "D004", where,
                                     "可执行的判定程序必须给出步骤"))
            if not proc.get("procedure_ref", "").strip():
                issues.append(_issue("error", "D004", where,
                                     "可执行的判定程序必须给出 procedure_ref"))
        elif dist.get("status") == "established":
            issues.append(_issue("error", "D004", where,
                                 "不可执行的判定程序不能标为 established"))

        if not str(dist.get("external_necessity", "")).strip():
            issues.append(_issue("error", "D005", where,
                                 "外部必要性不得为空——这是框架意义的核心判据"))

        # v3.4：边界依赖与参数敏感度（礁石 R-041 / R-042）。
        # 这两项目前全部「待测绘/未测」，因此这里报的是 warn 而非 error：
        # 它是 A11/A12 待办的**显式回声**，让"还没画边界地图"这件事每次校验都可见。
        missing_meta = []
        bd = dist.get("boundary_dependence")
        if not isinstance(bd, dict) or not str(bd.get("status", "")).strip():
            missing_meta.append("boundary_dependence（R-041 边界地图）")
        ps = dist.get("parameter_sensitivity")
        if not isinstance(ps, dict) or not str(ps.get("status", "")).strip():
            missing_meta.append("parameter_sensitivity（R-042 参数敏感性）")
        if missing_meta:
            issues.append(_issue("warn", "D010", where,
                                 "未登记：%s" % "、".join(missing_meta)))

        pairs = dist.get("discriminates", []) or []
        if not pairs:
            issues.append(_issue("error", "D006", where,
                                 "discriminates 不得为空"))
        for pair in pairs:
            left, right = pair.get("left"), pair.get("right")
            if left == right:
                issues.append(_issue("error", "D006", where,
                                     "区分的两侧不得相同：%r" % left))
            for side in (left, right):
                if side not in symbol_names and side not in LEVELS:
                    issues.append(_issue("error", "D006", where,
                                         "区分侧 %r 既不是已登记记号也不是可辨识级别"
                                         % side))
    return issues


# --------------------------------------------------------------------------
# 总入口
# --------------------------------------------------------------------------

def run_all(symbol_path=SYMBOL_REGISTRY, distinction_path=DISTINCTION_REGISTRY,
            with_trace_gate=True):
    reg = load_json(symbol_path)
    dreg = load_json(distinction_path)
    names = set(s.get("symbol") for s in reg.get("symbols", []))
    issues = check_symbols(reg) + check_distinctions(dreg, names)

    trace = {"errors": [], "warnings": [], "roots": 0,
             "trace_gaps_declared": 0, "enabled": False}
    if with_trace_gate:
        try:
            trace = xs_trace_gate.run_all()
            trace["enabled"] = True
        except FileNotFoundError as exc:
            issues.append(_issue("warn", "G009", "xs_trace_gate",
                                 "推导链/冻结窗闸门未运行（缺文件：%s）" % exc))
        issues += trace.get("errors", []) + trace.get("warnings", [])

    errors = [i for i in issues if i["level"] == "error"]
    warns = [i for i in issues if i["level"] == "warn"]
    return {
        "symbols": len(reg.get("symbols", [])),
        "distinctions": len(dreg.get("distinctions", [])),
        "legacy_conflicts": len(reg.get("legacy_conflicts", [])),
        "trace_gate": trace,
        "errors": errors,
        "warnings": warns,
    }


def format_report(report):
    lines = []
    lines.append("XS 框架形式校验")
    lines.append("  记号登记：%d 条" % report["symbols"])
    lines.append("  区分登记：%d 条" % report["distinctions"])
    lines.append("  遗留冲突：%d 条" % report["legacy_conflicts"])
    tg = report.get("trace_gate") or {}
    if tg.get("enabled"):
        lines.append("  推导链闸门 C1.7/C1.8：根记号 %d 条，已登记无根断言 %d 条"
                     % (tg.get("roots", 0), tg.get("trace_gaps_declared", 0)))
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
    report = run_all()
    if "--json" in argv:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        print(format_report(report))
    return 1 if report["errors"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
