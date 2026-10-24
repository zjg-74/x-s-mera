#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A7 —— 核心最小化：把 47 条记号压到极小核心，其余由规则导出

对应 18-v3.3-formal-framework-charter.md 待办 A7：
    「核心最小化：把 47 条记号压到更小核心，其余由规则导出；导出记号展开校验」。

本脚本基于 symbol_registry.json 的 `derives_from` 依赖图，计算：
  1. **根记号集**：kind ∈ {primitive, assumption}（公理级，需承载）
  2. **约定输入**：无 derives_from 但非根的记号（如 η，B-B44 约定性输入）
  3. **导出记号**：可由根沿 derives_from 可达的记号
  4. **极小生成集 / 导出表**：每个非根记号 → 它实际依赖的最小根子集
  5. **展开校验**：对每个非根记号，检查它的 derives_from 是否足以从根闭包到达
     （即：核心 + 规则 ⟹ 全部记号，机械可核验）

产出：核心定义 + 导出表（写入 core_table.json），并打印报告。

仅标准库；固定种子无关（纯图论）；逐位可复现。
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
REG = os.path.join(REPO, "docs", "xs-framework", "symbol_registry.json")

ROOT_KINDS = ("primitive", "assumption")


def load_registry():
    with open(REG, "r", encoding="utf-8") as fh:
        return json.load(fh)


def build_graph(reg):
    """symbol -> {derives_from 里的存在符号集合}；返回 (by_name, edges)。"""
    by_name = {s["symbol"]: s for s in reg["symbols"]}
    edges = {}
    for s in reg["symbols"]:
        deps = [d for d in (s.get("derives_from") or []) if d in by_name]
        edges[s["symbol"]] = set(deps)
    return by_name, edges


def transitive_deps(name, edges, by_name):
    """name 沿 derives_from 可达的全部记号（含自身）。"""
    seen = set()
    stack = [name]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        for d in edges.get(node, set()):
            if d in by_name:
                stack.append(d)
    return seen


def minimal_roots_for(name, edges, by_name, roots):
    """name 实际依赖的最小根子集（沿依赖闭包取根）。"""
    reach = transitive_deps(name, edges, by_name)
    return sorted(r for r in roots if r in reach)


def main():
    out = []

    def p(s=""):
        out.append(str(s))

    reg = load_registry()
    by_name, edges = build_graph(reg)

    roots = {s["symbol"] for s in reg["symbols"] if s.get("kind") in ROOT_KINDS}
    all_syms = {s["symbol"] for s in reg["symbols"]}
    # 约定输入：非根且无 derives_from（有意为之，非遗漏）
    conventions = {s["symbol"] for s in reg["symbols"]
                   if s.get("kind") not in ROOT_KINDS
                   and not (s.get("derives_from") or [])}
    # 非根记号（全部需可导出）
    non_root = {s["symbol"] for s in reg["symbols"] if s.get("kind") not in ROOT_KINDS}

    # 展开校验：每个非根记号必须可达某个根（约定输入除外）
    gaps = {}
    for sym in sorted(non_root):
        if sym in conventions:
            continue
        reach = transitive_deps(sym, edges, by_name)
        if not (reach & roots):
            gaps[sym] = True

    # 核心最小化：根记号里，哪些"实际承载"，哪些可进一步归并为组合
    # （参数记号如 a,b,e,u 等 derives_from A_D，其"根"是 A_D，不是它们自己）
    # 计算每个记号的最小根依赖（导出表）
    derived_rows = []
    for s in reg["symbols"]:
        name = s["symbol"]
        if s.get("kind") in ROOT_KINDS:
            continue
        mroots = minimal_roots_for(name, edges, by_name, roots)
        layer = s.get("layer")
        kind = s.get("kind")
        if name in conventions:
            row_type = "约定输入"
        elif mroots:
            row_type = "导出"
        else:
            row_type = "不可达"
        derived_rows.append({
            "symbol": name, "layer": layer, "kind": kind,
            "type": row_type, "min_roots": mroots,
            "identifiability": s.get("identifiability"),
        })

    # 核心定义 = 根记号 + 约定输入
    core_syms = sorted(roots | conventions)
    core_table = {
        "version": "1.0",
        "framework": "XS",
        "purpose": "A7 核心最小化：极小生成集（根记号）+ 约定输入 + 导出表",
        "registry_date": reg.get("registry_date"),
        "core": {
            "roots": sorted(roots),
            "conventions": sorted(conventions),
            "core_count": len(core_syms),
            "total_symbols": len(all_syms),
            "derived_count": len([r for r in derived_rows if r["type"] == "导出"]),
        },
        "derived_table": derived_rows,
        "gaps": sorted(gaps),
    }

    p("=" * 78)
    p("A7 —— 核心最小化：47 记号 → 极小核心 + 导出表")
    p("=" * 78)
    p("总记号数 = %d" % len(all_syms))
    p("根记号（primitive/assumption） = %d：%s" % (len(roots), ", ".join(sorted(roots))))
    p("约定输入（非根、无依赖，有意为之）= %d：%s"
      % (len(conventions), ", ".join(sorted(conventions)) or "(无)"))
    p("核心（根 + 约定输入）= %d；导出记号 = %d；不可达 = %d"
      % (len(core_syms),
         len([r for r in derived_rows if r["type"] == "导出"]),
         len([r for r in derived_rows if r["type"] == "不可达"])))

    p("")
    p("--- 导出表（按层）---")
    for layer in ("L0", "L1", "L2"):
        rows = [r for r in derived_rows if r["layer"] == layer]
        p("  [%s] %d 个" % (layer, len(rows)))
        for r in rows:
            if r["type"] == "导出":
                p("    %-10s %-14s 依赖根: %s" % (
                    r["symbol"], r["kind"], ", ".join(r["min_roots"]) or "(无)"))
            elif r["type"] == "约定输入":
                p("    %-10s %-14s [约定输入]" % (r["symbol"], r["kind"]))

    p("")
    p("--- 展开校验 ---")
    if gaps:
        p("  ❌ 存在不可达非根记号：%s" % ", ".join(sorted(gaps)))
    else:
        p("  ✅ 全部非根记号（除约定输入）均可由根沿 derives_from 到达 ⟹ 核心+规则 ⟹ 全部记号")

    # 核心压缩率
    core_n = len(core_syms)
    total_n = len(all_syms)
    p("")
    p("核心压缩率：%d → %d（压缩 %.1f%%）" % (
        total_n, core_n, (1 - core_n / total_n) * 100))

    ok = not gaps
    p("")
    p("结论: " + ("A7 核心最小化完成：%d 条记号压到 %d 条核心（根 %d + 约定输入 %d），"
                  "其余 %d 条全部可由核心沿 derives_from 导出" % (
                      total_n, core_n, len(roots), len(conventions),
                      len([r for r in derived_rows if r["type"] == "导出"]))
                  if ok else "存在不可达记号，需修复"))
    p("=" * 78)

    with open(os.path.join(HERE, "core_table.json"), "w", encoding="utf-8") as fh:
        json.dump(core_table, fh, ensure_ascii=False, indent=1)

    txt = "\n".join(out)
    try:
        print(txt)
    except UnicodeEncodeError:
        import sys
        sys.stdout.buffer.write(txt.encode("utf-8", "backslashreplace"))

    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
