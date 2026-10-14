#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
51_render_reference_list.py —— 从定稿表 references_final.csv 重新渲染人读清单 .txt

为什么需要单独的渲染脚本
------------------------
`44_finalize_references.py` 会**重新联网抓 PubMed**，而它的 FINAL 字典里若干
条目在定稿后被人工定点替换过（如 portincasa2006 替换掉无作者字段的 27121544）。
直接重跑 44 有把这些已核验的修正**改回去**的风险。
因此：**要重抓 PubMed 才跑 44；只是想让 .txt 与 .csv 同步，就跑本脚本。**
本脚本只读 csv、不联网，因此不可能改变任何文献内容。

用法
----
  PY=python3
  $PY scripts/51_render_reference_list.py
  $PY scripts/51_render_reference_list.py --check   # 只报差异，不写文件
"""

import csv
import datetime
import difflib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "主稿")
CSV_PATH = os.path.join(D, "references_final.csv")
TXT_PATH = os.path.join(D, "references_final.txt")

HEADER = "主稿参考文献（定稿，全部经 PubMed E-utilities 逐字段复核）"
SOURCE_LINE = ("复核时间：{ts} ｜ 内容来源：scripts/44_finalize_references.py（PubMed 逐字段核验）"
               " ｜ 渲染：scripts/51_render_reference_list.py（不联网，仅同步 .txt 与 .csv）")
RULE = "=" * 100


def render() -> str:
    with open(CSV_PATH, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    lines = [HEADER,
             SOURCE_LINE.format(ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
             RULE, ""]
    for r in rows:
        lines.append(r["vancouver"].strip())
        note = (r.get("note") or "").strip()
        note = note.split("；")[0][:60]          # 用途摘要，过长截断
        lines.append(f"   PMID {r['pmid']} | {note}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    check_only = "--check" in sys.argv
    new = render()
    old = open(TXT_PATH, encoding="utf-8").read() if os.path.exists(TXT_PATH) else ""

    if new == old:
        print("references_final.txt 已与 references_final.csv 一致，无需重写。")
        return 0

    diff = list(difflib.unified_diff(old.splitlines(), new.splitlines(),
                                     fromfile="old(.txt)", tofile="new(.txt)", lineterm="", n=1))
    print(f"检测到 .txt 与 .csv 不一致，差异 {len(diff)} 行：")
    for line in diff[:40]:
        print("  " + line)
    if len(diff) > 40:
        print(f"  …（共 {len(diff)} 行）")
    if check_only:
        print("\n--check 模式：未写文件。")
        return 1
    open(TXT_PATH, "w", encoding="utf-8").write(new)
    print(f"\n已重写：{TXT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
