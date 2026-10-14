#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
50_word_count.py —— 稿件字数口径的显式定义与统计（可复现）

为什么需要这个脚本
------------------
稿件的「字数」此前是用一次性命令算的，口径没有落盘，导致两件事：
① 同一个数字无法被第三方复算；② 换口径（例如某刊要求"含参考文献"）时
无法快速给出对照值。本脚本把口径写死，任何人跑一次都能得到同一个数。

口径定义（primary = whitespace 分词，与 Word 的计数方式最接近）
--------------------------------------------------------------
主文（Main text）
  = 从 "## 1. Introduction" 起，到 "## Acknowledgments" 前止
  排除：① 各级标题行；② 文内引用括号数字 (1,2)；③ Markdown 强调符号 * _ `
       与表格竖线 |；④ HTML 注释（<!-- ... -->）
  保留：章节体内的全部空白分词，含数字、缩写、统计学符号

摘要（Abstract）
  = "## Abstract" 与 "## 1. Introduction" 之间的正文
  排除：① "**Keywords**" 行；② "**Word count**" 行；③ 同上②③④
  用途：多数期刊的摘要上限（AJCN = 300 词）按此口径核算

次要口径（secondary = 只计含字母或数字的分词）
  排除纯标点分词（如孤立的破折号 "—"），给出一个更保守的下限值。

用法
----
  PY=python3
  $PY scripts/50_word_count.py                       # 默认统计 Manuscript_v1.md
  $PY scripts/50_word_count.py --file <path.md>
  $PY scripts/50_word_count.py --json
"""

import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT = os.path.join(ROOT, "results", "主稿", "Manuscript_v1.md")

CITE_RE = re.compile(r"\(\d{1,2}(?:,\d{1,2})*\)")
MD_RE = re.compile(r"[*_`|]")
HEAD_RE = re.compile(r"(?m)^#+\s*")
COMMENT_RE = re.compile(r"<!--.*?-->", re.S)
ALNUM_RE = re.compile(r"[A-Za-z0-9]")

# 章节锚点：主文 = Introduction .. Acknowledgments 之前
MAIN_START = "## 1. Introduction"
MAIN_END = "## Acknowledgments"
CHAPTERS = [
    ("Introduction", "## 1. Introduction", "## 2. Methods"),
    ("Methods", "## 2. Methods", "## 3. Results"),
    ("Results", "## 3. Results", "## 4. Discussion"),
    ("Discussion (incl. 4.7 Conclusions)", "## 4. Discussion", MAIN_END),
]


def clean(s: str) -> str:
    s = COMMENT_RE.sub(" ", s)
    s = CITE_RE.sub(" ", s)
    s = HEAD_RE.sub("", s)
    s = MD_RE.sub(" ", s)
    return s


def count(s: str) -> tuple:
    c = clean(s)
    toks = c.split()
    return len(toks), len([w for w in toks if ALNUM_RE.search(w)])


def slice_between(text: str, a: str, b: str) -> str:
    if a not in text or b not in text:
        return ""
    i = text.index(a)
    j = text.index(b, i + len(a))
    return text[i:j]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default=DEFAULT)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    text = open(args.file, encoding="utf-8").read()
    result = {"file": os.path.relpath(args.file, ROOT), "chapters": {}}

    for name, a, b in CHAPTERS:
        result["chapters"][name] = count(slice_between(text, a, b))

    main_txt = slice_between(text, MAIN_START, MAIN_END)
    if not main_txt:
        main_txt = text[text.index(MAIN_START):]
    result["main_text"] = count(main_txt)

    ab = slice_between(text, "## Abstract", "## 1. Introduction")
    ab = re.sub(r"(?s)\*\*Keywords\*\*.*", "", ab)          # 去关键词行
    ab = re.sub(r"(?s)\*\*Word count\*\*.*", "", ab)         # 去字数声明行
    result["abstract"] = count(ab)

    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0

    print("字数统计（口径见本脚本 docstring）")
    print("=" * 66)
    print(f"文件：{result['file']}")
    print("-" * 66)
    for name, (a, b) in result["chapters"].items():
        print(f"  {name:<36} {a:>6}  （保守口径 {b}）")
    a, b = result["main_text"]
    print("-" * 66)
    print(f"  {'主文 Main text（Introduction–Discussion）':<36} {a:>6}  （保守口径 {b}）")
    a, b = result["abstract"]
    print(f"  {'摘要 Abstract':<36} {a:>6}  （保守口径 {b}）")
    print("=" * 66)
    print("口径：空白分词；排除标题行、文内引用括号数字、Markdown 标记与注释。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
