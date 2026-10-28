"""XS 边界地图机检校验器（仅标准库）。

用途：把边界地图 `docs/xs-framework/boundary_map.json` 从"仅登记"推进到"被机器看住"。
对应 A10 残留（27-v3.4-boundary-map.md §6.5：机检未接；逻辑桥"让地图被机器看住"）：
逐版本检查测绘状态的结构完整性、id 唯一性、gap_index 的解决指向文件存在性、
开放格状态非空等。error > 0 时退出码 1，接入 CI 门。

用法：
    python src/xs_boundary_map_check.py          # 报告，error>0 → 退出码 1
    python src/xs_boundary_map_check.py --json   # 输出 JSON
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DOCS = os.path.join(REPO, "docs", "xs-framework")
BOUNDARY_MAP = os.path.join(DOCS, "boundary_map.json")

TOP_LEVEL_REQUIRED = (
    "version", "framework", "name", "dimensions",
    "r041_table", "open_cells", "gap_index",
)
DIM_FIELDS = ("id", "name", "judge", "refs", "status")
R041_FIELDS = ("level", "quantity", "closed_instance",
               "open_instance", "instance_dependent", "evidence")
OPEN_CELL_FIELDS = ("id", "dim", "issue", "status", "next")
GAP_FIELDS = ("cell", "ref", "status")


def _issue(level, code, where, message):
    return {"level": level, "code": code, "where": where, "message": message}


def load_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def check_boundary_map(path=BOUNDARY_MAP):
    issues = []
    try:
        bm = load_json(path)
    except FileNotFoundError:
        return [_issue("error", "M001", "boundary_map",
                       "缺文件：%s（边界地图须登记）" % path)]
    except ValueError as exc:
        return [_issue("error", "M001", "boundary_map", "JSON 解析失败：%s" % exc)]

    # 顶层必填
    for field in TOP_LEVEL_REQUIRED:
        if field not in bm:
            issues.append(_issue("error", "M002", "boundary_map",
                                 "缺少必填顶层字段 %s" % field))
        elif not isinstance(bm[field], (list, dict, str, int)):
            issues.append(_issue("error", "M002", "boundary_map",
                                 "顶层字段 %s 类型非法" % field))

    # 版本
    ver = bm.get("version")
    if not isinstance(ver, str) or not ver.strip():
        issues.append(_issue("error", "M003", "boundary_map",
                             "version 必须是非空字符串"))

    # dimensions
    dims = bm.get("dimensions") or []
    if not isinstance(dims, list):
        issues.append(_issue("error", "M004", "boundary_map",
                             "dimensions 必须是列表"))
        dims = []
    seen_dim = set()
    for idx, dim in enumerate(dims):
        where = "dimensions[%d]" % idx
        for f in DIM_FIELDS:
            if f not in dim or (isinstance(dim.get(f), str) and not dim[f].strip()):
                issues.append(_issue("error", "M005", where, "缺少必填字段 %s" % f))
        did = dim.get("id")
        if did in seen_dim:
            issues.append(_issue("error", "M006", where, "dimension id 重复：%r" % did))
        seen_dim.add(did)

    # r041_table
    r041 = bm.get("r041_table") or []
    if not isinstance(r041, list):
        issues.append(_issue("error", "M007", "boundary_map",
                             "r041_table 必须是列表"))
        r041 = []
    seen_lv = set()
    for idx, row in enumerate(r041):
        where = "r041_table[%d]" % idx
        for f in R041_FIELDS:
            if f not in row:
                issues.append(_issue("error", "M008", where, "缺少必填字段 %s" % f))
        lv = row.get("level")
        if lv in seen_lv:
            issues.append(_issue("error", "M009", where, "level 重复：%r" % lv))
        seen_lv.add(lv)
        if lv not in ("A", "B", "C", "D"):
            issues.append(_issue("warn", "M010", where,
                                 "level 非 A/B/C/D：%r" % lv))
        for key in ("closed_instance", "open_instance"):
            cell = row.get(key) or {}
            if not isinstance(cell, dict):
                issues.append(_issue("error", "M011", where,
                                     "%s 必须是对象" % key))

    # open_cells
    oc = bm.get("open_cells") or []
    if not isinstance(oc, list):
        issues.append(_issue("error", "M012", "boundary_map",
                             "open_cells 必须是列表"))
        oc = []
    seen_oc = set()
    for idx, c in enumerate(oc):
        where = "open_cells[%d]" % idx
        for f in OPEN_CELL_FIELDS:
            if f not in c or (isinstance(c.get(f), str) and not c[f].strip()):
                issues.append(_issue("error", "M013", where, "缺少必填字段 %s" % f))
        cid = c.get("id")
        if cid in seen_oc:
            issues.append(_issue("error", "M014", where, "open_cell id 重复：%r" % cid))
        seen_oc.add(cid)
        st = str(c.get("status", "")).strip()
        if st and not any(seg in st for seg in
                          ("已", "待", "约定", "开放", "定稿", "解", "闭合", "验证")):
            issues.append(_issue("warn", "M015", where,
                                 "status 无法归入已知状态词：%r" % st))

    # gap_index：ref 必填；已解决/已做/已闭合 的应有 resolved_in 且文件存在
    gi = bm.get("gap_index") or []
    if not isinstance(gi, list):
        issues.append(_issue("error", "M016", "boundary_map",
                             "gap_index 必须是列表"))
        gi = []
    for idx, g in enumerate(gi):
        where = "gap_index[%d]" % idx
        for f in GAP_FIELDS:
            if f not in g or (isinstance(g.get(f), str) and not g[f].strip()):
                issues.append(_issue("error", "M017", where, "缺少必填字段 %s" % f))
        st = str(g.get("status", "")).strip()
        resolved = g.get("resolved_in")
        is_done = any(seg in st for seg in ("已", "解决", "闭合", "核验", "定稿"))
        if is_done and not resolved:
            issues.append(_issue("warn", "M018", where,
                                 "status 表明已解决/已闭合，但缺 resolved_in 指向"))
        if resolved:
            rp = os.path.join(REPO, str(resolved))
            if not os.path.exists(rp):
                issues.append(_issue("error", "M019", where,
                                     "resolved_in 指向不存在：%s" % resolved))

    return issues


def run_all(path=BOUNDARY_MAP):
    issues = check_boundary_map(path)
    errors = [i for i in issues if i["level"] == "error"]
    warns = [i for i in issues if i["level"] == "warn"]
    return {"path": path, "errors": errors, "warnings": warns,
            "error_count": len(errors), "warning_count": len(warns)}


def format_report(report):
    lines = []
    lines.append("XS 边界地图机检（A10 残留落地）")
    lines.append("  登记件：%s" % report["path"])
    lines.append("  错误：%d    警告：%d"
                 % (report["error_count"], report["warning_count"]))
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
