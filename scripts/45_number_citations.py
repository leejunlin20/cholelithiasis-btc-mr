#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
45_number_citations.py —— 正文引用自动编号 + 参考文献渲染 + 三对齐机器校验

为什么这样做
------------
Vancouver 编号的定义就是「按正文首次出现顺序」。手工编号在反复修改中必然与
参考文献表脱节，所以编号必须由脚本生成，人手不介入。

上一版的三个缺陷（已修）
------------------------
① **写回源文件**：旧版把渲染结果写回 Manuscript_v1.md 这个「源」，导致多键组
   `[@a; @b]` 的一个正则失误就把整组引用键烧毁（只剩一个数字），且不可逆。
   → 现在严格区分：源（带 [@key]）与产物（带 [n]）是两个文件，本脚本只读源、
     只写产物。
② **条目重复**：旧版用 `authors` 列重新拼接，而该列已被污染（含整条引文），
   于是标题/期刊被打印两遍；且首条编号丢失。
   → 现在直接取 references_final.csv 的 `vancouver` 列（该列已由 PubMed
     esummary 重建、逐条核验），仅剥掉其自带的行首序号。
③ **补充项校验逻辑错误**：旧版把「补充表编号」与「主表编号」取交集判冲突，
   而两者是**独立命名空间**（Table 1–3 vs Table S1–S5），必然误报。
   → 现在按各自命名空间分别校验「从 1 起连续」，并额外检查不存在漏掉 S 前缀
     的越界编号。

用法
----
  PY=python3
  $PY scripts/45_number_citations.py            # 渲染 + 校验 + 写出产物
  $PY scripts/45_number_citations.py --check    # 只校验，不写产物

输入：results/主稿/Manuscript_source.md（含 [@key]）、results/主稿/references_final.csv
输出：results/主稿/Manuscript_v1.md（含 [n] 与 References 段）
      results/主稿/三对齐校验报告.txt
"""

import csv
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "主稿")
SRC = os.path.join(D, "Manuscript_source.md")
OUT = os.path.join(D, "Manuscript_v1.md")
REFS = os.path.join(D, "references_final.csv")
REPORT = os.path.join(D, "三对齐校验报告.txt")

REF_HEADING = "## References"
TABLE_ANCHOR = "## Tables"
CITE_RE = re.compile(r"\[@([^\]]+)\]")

# 文内引用渲染格式。AJCN / ASN 系 = 括号数字、逗号无空格、句读前，如 (1,2)。
# 换刊时只改这两个常量（以及 CITE_SEP），不要改上面的 CITE_RE。
CITE_OPEN = "("
CITE_CLOSE = ")"
CITE_SEP = ","


def load_refs():
    with open(REFS, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    return {r["key"]: r for r in rows}


def strip_leading_number(vancouver: str) -> str:
    return re.sub(r"^\s*\d+\.\s*", "", vancouver).strip()


def first_seen_order(text, pattern):
    """返回按首次出现顺序去重的捕获组列表。"""
    seen, order = set(), []
    for m in re.finditer(pattern, text):
        g = m.group(1)
        if g not in seen:
            seen.add(g)
            order.append(g)
    return order


def main():
    dry = "--check" in sys.argv
    refs = load_refs()
    src = open(SRC, encoding="utf-8").read()

    # ---------- 1. 引用键 → 编号（按首次出现顺序） ----------
    keys_in_order = []
    seen = set()
    for grp in CITE_RE.findall(src):
        for k in (x.strip().lstrip("@") for x in grp.split(";")):
            if k and k not in seen:
                seen.add(k)
                keys_in_order.append(k)
    number = {k: i for i, k in enumerate(keys_in_order, 1)}

    missing = [k for k in keys_in_order if k not in refs]
    uncited = sorted(k for k in refs if k not in number)

    # ---------- 2. 正文替换为 (n) / (n,m) ----------
    # 括号数字是 AJCN / ASN 系期刊的强制文内引用格式（官方 Guide for Authors:
    # "identify references by Arabic numerals in parentheses before the punctuation,
    #  not as superscripts"）。切换目标期刊时只需改 CITE_OPEN / CITE_CLOSE。
    def repl(m):
        ks = [x.strip().lstrip("@") for x in m.group(1).split(";")]
        nums = sorted(number[k] for k in ks if k in number)
        return CITE_OPEN + CITE_SEP.join(str(n) for n in nums) + CITE_CLOSE

    body = CITE_RE.sub(repl, src)
    # 源文件不应再含 [@ 或 ?? 残留
    leftover_at = body.count("[@")
    leftover_q = body.count("[@??") + len(re.findall(r"@\?\?", body))

    # ---------- 3. 渲染 References 段（取自 vancouver 列） ----------
    lines = [REF_HEADING, "",
             "*References are numbered in order of first citation in the text. "
             "Every entry was retrieved from and verified against PubMed "
             "(see Supplementary Methods).*", ""]
    ref_gaps = []
    for k in keys_in_order:
        if k not in refs:
            continue
        r = refs[k]
        van = strip_leading_number(r["vancouver"])
        lines.append(f"{number[k]}. {van}")
        lines.append(f"   <!-- PMID {r['pmid']} -->")
    # 编号连续性
    emitted = sorted(number[k] for k in keys_in_order if k in refs)
    if emitted and emitted != list(range(1, len(emitted) + 1)):
        ref_gaps.append(f"编号不连续：{emitted[:5]}…{emitted[-5:]}")

    # References 段插在正文（含 Acknowledgments / AI 声明）之后、Tables 之前——这是
    # Elsevier / ASN 系期刊的稿件文件排布惯例：正文 → 致谢 → 参考文献 → 表 → 图注。
    refs_block = "\n".join(lines)
    if "\n" + TABLE_ANCHOR in body:
        out = body.replace(TABLE_ANCHOR, refs_block + "\n\n---\n\n" + TABLE_ANCHOR, 1)
    else:
        out = body.rstrip() + "\n\n---\n\n" + refs_block + "\n"

    # ---------- 4. 三对齐校验 ----------
    head, _, tail = out.partition("## Tables")
    tab_sec, _, leg_sec = tail.partition("## Figure legends")

    fig_cites = first_seen_order(head, r"\*\*Figure (\d)\*\*")
    tab_cites = first_seen_order(head, r"\*\*Table (\d)\*\*")
    fig_seq = [int(x) for x in fig_cites]
    tab_seq = [int(x) for x in tab_cites]
    leg_def = [int(m.group(1)) for m in re.finditer(r"(?m)^\*\*Figure (\d)\.", leg_sec)]
    tab_def = [int(m.group(1)) for m in re.finditer(r"(?m)^### Table (\d)\.", tab_sec)]

    supp_t = sorted(set(int(m) for m in re.findall(r"Table S(\d)", out)))
    supp_f = sorted(set(int(m) for m in re.findall(r"Figure S(\d)", out)))

    checks = []

    def chk(name, ok, detail):
        checks.append((name, ok, detail))

    # —— 图 ——
    chk("①-1 Figure 正文首次引用顺序递增", fig_seq == sorted(fig_seq), f"顺序 {fig_seq}")
    chk("①-2 Figure 图注物理顺序 = 1..N", leg_def == list(range(1, len(leg_def) + 1)),
        f"物理顺序 {leg_def}")
    chk("①-3 正文引用图集 == 图注集（无多引/漏引）", set(fig_seq) == set(leg_def),
        f"正文 {sorted(set(fig_seq))} vs 图注 {sorted(set(leg_def))}")
    chk("①-4 Figure 编号从 1 起且连续", fig_seq == list(range(1, len(fig_seq) + 1)),
        f"{fig_seq}")

    # —— 表 ——
    chk("②-1 Table 正文首次引用顺序递增", tab_seq == sorted(tab_seq), f"顺序 {tab_seq}")
    chk("②-2 Table 表体物理顺序 = 1..N", tab_def == list(range(1, len(tab_def) + 1)),
        f"物理顺序 {tab_def}")
    chk("②-3 正文引用表集 == 表体集（无多引/漏引）", set(tab_seq) == set(tab_def),
        f"正文 {sorted(set(tab_seq))} vs 表体 {sorted(set(tab_def))}")
    chk("②-4 Table 编号从 1 起且连续", tab_seq == list(range(1, len(tab_seq) + 1)),
        f"{tab_seq}")

    # —— 补充项：独立命名空间 ——
    chk("③-1 补充表编号从 S1 起连续（独立命名空间）",
        supp_t == list(range(1, len(supp_t) + 1)) if supp_t else False,
        f"Table S{supp_t}")
    chk("③-2 补充图编号从 S1 起连续（独立命名空间）",
        supp_f == list(range(1, len(supp_f) + 1)) if supp_f else False,
        f"Figure S{supp_f}")
    chk("③-3 主文编号不带 S 前缀、补充编号必带 S 前缀",
        not re.search(r"Table S0|Figure S0", out) and not re.search(r"\*\*Table S\d+\*\*", head),
        "主文正文不出现裸 S 编号")

    # —— 引用完整性 ——
    chk("④-1 引用键全部解析到已核验文献", not missing, f"缺 {missing or '无'}")
    chk("④-2 已核验文献全部被引用（无僵尸文献）", not uncited, f"未被引用 {uncited or '无'}")
    chk("④-3 References 编号 1..N 连续无缺号", not ref_gaps, ref_gaps[0] if ref_gaps else "连续")
    chk("④-4 正文无 [@key] / ?? 残留", leftover_at == 0 and leftover_q == 0,
        f"[@ 残留 {leftover_at} 处；?? 残留 {leftover_q} 处")
    chk("④-5 正文引用键数 == 文献表条数",
        len(keys_in_order) == len(refs) and not missing and not uncited,
        f"引用 {len(keys_in_order)} 条 / 文献 {len(refs)} 条")

    # ---------- 5. 报告 ----------
    rep = ["主稿「三对齐」机器校验报告", "=" * 78,
           f"源文件：{os.path.relpath(SRC, ROOT)}",
           f"产物　：{os.path.relpath(OUT, ROOT)}",
           f"文献表：{os.path.relpath(REFS, ROOT)}（{len(refs)} 条，全部经 PubMed 逐字段核验）",
           f"正文引用：{len(keys_in_order)} 个键 → 编号 1..{len(number)}",
           "-" * 78]
    allok = True
    for name, ok, detail in checks:
        rep.append(f"[{'通过' if ok else '不通过'}] {name}")
        rep.append(f"          {detail}")
        allok &= ok
    rep.append("-" * 78)
    rep.append("三对齐明细")
    rep.append(f"  图：正文首次引用顺序 {fig_seq}")
    rep.append(f"      图注物理顺序     {leg_def}")
    rep.append(f"  表：正文首次引用顺序 {tab_seq}")
    rep.append(f"      表体物理顺序     {tab_def}")
    rep.append(f"  补充：Table S{supp_t} ｜ Figure S{supp_f}")
    rep.append("=" * 78)
    rep.append("结论：" + ("全部通过（三对齐达成，可交付）" if allok
                        else "存在未通过项，需修正后重跑"))
    text = "\n".join(rep) + "\n"

    if not dry:
        open(OUT, "w", encoding="utf-8").write(out)
        open(REPORT, "w", encoding="utf-8").write(text)

    print(text)
    if not dry:
        print(f"已写出：{OUT}")
        print(f"已写出：{REPORT}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
