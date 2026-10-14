#!/usr/bin/env python3
"""49_parse_author_guide.py —— 出版社作者指南 HTML/PDF 关键条款抽取器

用途（journal-due-diligence Step 3）：把抓下来的作者指南页面转成纯文本，
再按关键词定位「稿件类型 / 字数 / 摘要 / 参考文献 / 图表 / APC / 排除条款」相关段落，
供人工逐字核对。**只做定位与摘录，不做判断，不替用户下结论。**

支持 .html/.htm 与 .pdf（PDF 需 pypdf）。

用法：
  PY=python3   # 有 pypdf
  $PY scripts/49_parse_author_guide.py /tmp/eje.out --kw "word" --kw "APC"
  $PY scripts/49_parse_author_guide.py /tmp/ajcn.pdf --list-sections
"""
from __future__ import annotations

import argparse
import html
import re
import sys


def to_text(path: str) -> str:
    raw = open(path, "rb").read()
    if raw[:4] == b"%PDF":
        try:
            import pypdf  # noqa: PLC0415
        except ImportError:
            sys.exit("读 PDF 需要 pypdf：pip install pypdf")
        txt = "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(path).pages)
    else:
        s = raw.decode("utf-8", "ignore")
        s = re.sub(r"(?is)<(script|style|noscript)[^>]*>.*?</\1>", " ", s)
        s = re.sub(r"(?i)<br\s*/?>|</(p|div|li|tr|h[1-6])>", "\n", s)
        txt = html.unescape(re.sub(r"<[^>]+>", "\n", s))
    txt = re.sub(r"[ \t\u00a0]+", " ", txt)
    txt = re.sub(r"\n\s*\n+", "\n", txt)
    return "\n".join(x.strip() for x in txt.split("\n") if x.strip())


# 常用检索组：每组给若干正则，命中即摘录上下文
GROUPS = {
    "字数": r"(word limit|words? \(|not exceed|no more than|maximum of|limited to|fewer than)[^.\n]{0,140}\bword",
    "字数2": r"\b\d[\d,]{2,6}\s*word",
    "摘要": r"abstract[^.\n]{0,160}(word|structured|unstructured|sections?)",
    "参考文献": r"(references?|reference list)[^.\n]{0,130}(exceed|maximum|limit|up to|no more than|around)",
    "图表": r"(tables? (and|&|\+) figures?|figures? (and|&|\+) tables?|displays?|illustrations?)[^.\n]{0,130}(up to|maximum|no more than|limit)",
    "排除条款": r"(will not be considered|not considered for|not accepted|does not accept|are not considered|desk reject|rejected without review)",
    "描述性/方法学": r"(purely descriptive|descriptive research|methodolog\w+ papers?|preliminary|confirmatory|no novelty)",
    "生信/数据库": r"(bioinformatic\w*|in silico|public\w* (data|database)|re-?analysis of|computational)",
    "MR": r"Mendelian randomi\w*",
    "报告清单": r"(STROBE|RECORD|PRISMA|CONSORT|MOOSE|ARRIVE|reporting guideline|checklist extension)",
    "APC": r"(\bAPC\b|article processing charge|open access (charge|licence|license)|open choice)",
    "费用": r"(page charge|colour charge|color charge|publication charge|fee of|USD|GBP|EUR)\s?\$?\d",
    "稿酬/减免": r"(waiver|discount|Research4Life|free of charge|no charge)",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--kw", action="append", default=[], help="自定义正则（可多次）")
    ap.add_argument("--group", action="append", default=[], help="从内置组里选（可多次）")
    ap.add_argument("--list-groups", action="store_true")
    ap.add_argument("--list-sections", action="store_true", help="列出疑似章节标题行")
    ap.add_argument("--width", type=int, default=300)
    ap.add_argument("--max", type=int, default=12, help="每组最多打印多少条")
    a = ap.parse_args()

    if a.list_groups:
        for k, v in GROUPS.items():
            print(f"{k:<14} {v}")
        return

    txt = to_text(a.path)
    print(f"文件: {a.path}   纯文本长度: {len(txt):,} 字符\n")

    if a.list_sections:
        print("=== 疑似章节标题（短行且首字母大写） ===")
        for i, ln in enumerate(txt.split("\n")):
            if 3 <= len(ln) <= 70 and re.match(r"^[A-Z][A-Za-z0-9 ,/&()\-']+$", ln) and not ln.isupper() or (ln.isupper() and 3 <= len(ln) <= 60):
                print(f"  L{i:>5}  {ln}")
        print()

    chosen = {}
    if a.group:
        for g in a.group:
            if g not in GROUPS:
                sys.exit(f"未知组 {g!r}；用 --list-groups 查看")
            chosen[g] = GROUPS[g]
    if a.kw:
        for i, k in enumerate(a.kw, 1):
            chosen[f"自定义{i}: {k}"] = k
    if not chosen:
        chosen = GROUPS

    for label, pat in chosen.items():
        print("=" * 100)
        print(f"### {label}")
        print("=" * 100)
        hits = 0
        seen = set()
        for m in re.finditer(pat, txt, re.I):
            frag = re.sub(r"\s+", " ", txt[max(0, m.start() - 90): m.end() + 220]).strip()
            key = frag[:90]
            if key in seen:
                continue
            seen.add(key)
            hits += 1
            print(f"  [{hits}] …{frag[:a.width + 130]}…")
            if hits >= a.max:
                print(f"  （已达上限 {a.max}，如需更多请调 --max）")
                break
        if hits == 0:
            print("  （未命中）")
        print()


if __name__ == "__main__":
    main()
