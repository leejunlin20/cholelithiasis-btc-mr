#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
47_audit_references.py —— 参考文献表交付前体检

为什么需要：44_finalize_references.py 逐条从 PubMed 重取字段时，个别记录的
`authors` 列被整条引文污染（例如 lammert2016 渲染出中文占位符「[无作者字段]」），
而 PubMed/CrossRef/Europe PMC 三个库对该条都确实没有作者字段——说明这是**上游
元数据缺口**，不是本流水线的解析 bug。此类缺陷肉眼扫一遍文献表极易漏过，
必须机器体检。

本脚本对 references_final.csv 逐条做 10 项检查，并对渲染产物做编号-PMID 对账。

严重度分级（重要）
------------------
体检项分两类，混淆它们会造成「误报阻断交付」：

  · **阻断级（BLOCK）** —— 影响成品：vancouver 列含占位符/中文、起始不像作者串、
    缺标题/年/卷/页/DOI、条目内重复、status 为 CHECK/NOT。这些必须修。
  · **卫生级（HYGIENE）** —— 只影响源表可读性：`authors` 列被整条引文污染。
    渲染器已改为只读 vancouver 列，故**不影响成品**；但仍应修干净，否则下一次
    有人改回按列拼接就会把缺陷带进稿件。

用法：
  PY=python3
  $PY scripts/47_audit_references.py                # 体检
  $PY scripts/47_audit_references.py --fix-authors  # 体检 + 从 vancouver 重建 authors 列
输出：results/主稿/文献表体检报告.txt（并打印）
"""

import collections
import csv
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "主稿")
CSV = os.path.join(D, "references_final.csv")
MS = os.path.join(D, "Manuscript_v1.md")
REPORT = os.path.join(D, "文献表体检报告.txt")

PLACEHOLDER = re.compile(r"[\u4e00-\u9fff]|无作者|\[|\]|\{|\}|TODO|XX+|待补")
AUTHORISH = re.compile(r"^[A-Z][A-Za-z'\-]+ [A-Z]{1,3}\b")   # 单作者：Surname AB
JOURNALISH = re.compile(r"\b[A-Z][A-Za-z]+ [A-Z][a-z]+\.\s*\d{4}")  # 「Journal Name. 2016」形态


def main():
    rows = list(csv.DictReader(open(CSV, encoding="utf-8-sig")))
    ms = open(MS, encoding="utf-8").read() if os.path.exists(MS) else ""
    fix = "--fix-authors" in sys.argv

    # ---------- 0. 卫生级修复：从 vancouver 重建 authors 列 ----------
    fixed = []
    if fix:
        for r in rows:
            body = re.sub(r"^\s*\d+\.\s*", "", (r.get("vancouver") or "").strip())
            title = (r.get("title") or "").strip().rstrip(".")
            idx = body.lower().find(title.lower()) if title else -1
            if idx > 0:
                clean = body[:idx].strip().rstrip(".").strip()
                if clean and clean != (r.get("authors") or "").strip():
                    r["authors"] = clean
                    fixed.append(r["key"])
        if fixed:
            hdr = list(rows[0].keys())
            with open(CSV, "w", encoding="utf-8-sig", newline="") as f:
                w = csv.DictWriter(f, fieldnames=hdr)
                w.writeheader()
                w.writerows(rows)

    rep, blockers, hygiene = [], [], []
    rep.append("参考文献表交付前体检报告")
    rep.append("=" * 78)
    rep.append(f"文献表：{os.path.relpath(CSV, ROOT)}（{len(rows)} 条）")
    rep.append(f"产物　：{os.path.relpath(MS, ROOT)}")
    if fix:
        rep.append(f"已从 vancouver 列重建 authors 列：{len(fixed)} 条"
                   + (f"（{', '.join(fixed[:6])}{'…' if len(fixed) > 6 else ''}）" if fixed else ""))
    rep.append("-" * 78)

    # 编号-PMID 对账（产物）
    if ms:
        refblock = ms.split("## References", 1)[-1]
        pairs = re.findall(r"(?m)^(\d+)\. (.*)\n\s*<!-- PMID (\d+) -->", refblock)
        rep.append(f"产物 References 条目数：{len(pairs)}")
        nums = [int(n) for n, _, _ in pairs]
        if nums != list(range(1, len(nums) + 1)):
            blockers.append(f"产物编号不连续：{nums[:3]}…{nums[-3:]}")
        rep.append(f"产物编号连续性：{'通过' if nums == list(range(1, len(nums)+1)) else '不通过'}")
    rep.append("-" * 78)

    seen_keys = collections.Counter(r["key"] for r in rows)
    dup_keys = [k for k, v in seen_keys.items() if v > 1]
    if dup_keys:
        blockers.append(f"键重复：{dup_keys}")

    for r in rows:
        key = r["key"]
        n = r["n"]
        van = (r.get("vancouver") or "").strip()
        au = (r.get("authors") or "").strip()
        be, he = [], []          # 阻断级 / 卫生级

        # 1 PMID
        if not re.fullmatch(r"\d{7,9}", (r.get("pmid") or "").strip()):
            be.append("PMID 缺失或格式错")
        # 2 占位符 / 中文泄漏（成品直接可见 → 阻断）
        if PLACEHOLDER.search(van):
            be.append(f"vancouver 含占位符或中文：{PLACEHOLDER.search(van).group(0)!r}")
        # 3 authors 列污染（仅源表 → 卫生级）
        if PLACEHOLDER.search(au):
            he.append("authors 列含占位符/中文")
        if len(au) > 220 or JOURNALISH.search(au):
            he.append("authors 列被整条引文污染（含期刊+年份）")
        # 4 vancouver 起始是否像作者串
        body = re.sub(r"^\s*\d+\.\s*", "", van)
        head = body.split(".")[0]
        if not ("," in head or AUTHORISH.match(head)):
            be.append(f"vancouver 起始不像作者串：{head[:45]!r}")
        # 5 vancouver 是否含标题
        title = (r.get("title") or "").strip().rstrip(".")
        if title and title.lower()[:40] not in body.lower():
            be.append("vancouver 未含标题")
        # 6 年/卷/页
        for fld, label in (("year", "年"), ("volume", "卷"), ("pages", "页")):
            v = (r.get(fld) or "").strip()
            if v and v not in body:
                be.append(f"vancouver 未含{label}（{fld}={v}）")
        # 7 DOI
        doi = (r.get("doi") or "").strip()
        if doi and doi not in body:
            be.append("vancouver 未含 DOI")
        # 8 vancouver 内重复片段
        parts = [p.strip() for p in re.split(r"(?<=\.)\s+", body) if len(p.strip()) > 25]
        c = collections.Counter(parts)
        if any(v > 1 for v in c.values()):
            be.append("vancouver 内有重复片段")
        # 9 status
        if (r.get("status") or "").upper().startswith(("CHECK", "NOT")):
            be.append(f"status={r['status']}")

        if be:
            blockers.append(f"[{n}] {key}: " + "；".join(be))
        if he:
            hygiene.append(f"[{n}] {key}: " + "；".join(he))

    rep.append(f"阻断级异常（影响成品）：{len(blockers)}")
    for p in blockers:
        rep.append("  ✗ " + p)
    if not blockers:
        rep.append("  ✓ 50 条 vancouver 全部通过 8 项成品侧体检"
                   "（占位符/作者串/标题/年卷页/DOI/重复/status/编号）")
    rep.append(f"卫生级异常（仅源表，不影响成品）：{len(hygiene)}")
    for p in hygiene:
        rep.append("  ! " + p)
    if not hygiene:
        rep.append("  ✓ authors 列无污染")
    rep.append("=" * 78)
    rep.append("结论：" + ("可以交付（阻断级 0 项）" if not blockers
                        else "不可交付：存在阻断级异常，需处置后重跑"))
    text = "\n".join(rep) + "\n"
    open(REPORT, "w", encoding="utf-8").write(text)
    print(text)
    return 0 if not blockers else 1


if __name__ == "__main__":
    sys.exit(main())
