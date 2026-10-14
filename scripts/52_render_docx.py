#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
52_render_docx.py —— 把主稿 Markdown 渲染为投稿用 Word（AJCN / Elsevier 口径）

用法
----
    PY=python3
    $PY scripts/52_render_docx.py            # 全部产出
    $PY scripts/52_render_docx.py --check    # 只做产出后自检

产出（写入 results/主稿/，并同步拷贝到 ~/Desktop/投稿/AJCN/）
--------------------------------------------------------------
    AJCN_Manuscript.docx                投稿主件：
                                        标题页 → 摘要 → 正文 → 致谢 → 生成式 AI 声明
                                        → 参考文献 → 表 1–3（横向节）→ 图注 1–5
                                        → 补充材料索引
    AJCN_Manuscript_with_figures.docx   同上 + 文末「现有图件附录」（内部存档，非投稿件）
    AJCN_STROBE-MR_Checklist.docx       AJCN 强制随稿上传的报告清单

排版约定（依 ASN 官方 Guide for Authors 与 Elsevier 稿件惯例）
--------------------------------------------------------------
  * 正文 Times New Roman 12 pt、双倍行距、左对齐、连续行号
  * 页脚居中页码；Letter 纸、四边 1 英寸
  * 投稿稿件惯例：各级标题不放大字号，仅加粗
  * 表格区切换为横向（landscape）独立节，表格字号 8.5 pt、表头跨页重复
  * 文内引用已是 AJCN 强制的括号数字（由 45_number_citations.py 生成）

行内标记解析
------------
  **粗体** / *斜体* / `代码`（去反引号，不作等宽，避免稿件里出现等宽字体）
  [[PLACEHOLDER]] → 黄底加粗（提醒待填）
  Unicode 上下标（10⁻⁵、*p*₁₂）→ Word 真上下标
  HTML 注释（<!-- PMID ... -->）不进正文

⚠️ 本脚本只做「Markdown → Word」的忠实搬运，不改写任何文字。
   要改文字请改 results/主稿/Manuscript_source.md，然后重跑 45 号脚本再跑本脚本。
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import zipfile
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_COLOR_INDEX
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "results" / "主稿"
MD_MANUSCRIPT = MAIN / "Manuscript_v1.md"
MD_CHECKLIST = MAIN / "STROBE-MR_checklist.md"
DESKTOP_OUT = Path.home() / "Desktop" / "投稿" / "AJCN"

FONT = "Times New Roman"
BODY_PT = 12
TABLE_PT = 8.5

# ---------------------------------------------------------------- 行内标记

SUP_CHARS = set("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺ⁿ")
SUB_CHARS = set("₀₁₂₃₄₅₆₇₈₉₊₋")
SUP_MAP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺ⁿ", "0123456789-+n")
SUB_MAP = str.maketrans("₀₁₂₃₄₅₆₇₈₉₊₋", "0123456789+-")

PLACEHOLDER = re.compile(r"\[\[[A-Z_]+\]\]")


def tokenize(s: str, bold: bool = False, italic: bool = False):
    """把一段行内 Markdown 拆成 (text, bold, italic, highlight) 序列。"""
    out, buf, i = [], [], 0

    def flush():
        if buf:
            out.append(("".join(buf), bold, italic, False))
            buf.clear()

    while i < len(s):
        if s.startswith("**", i):
            j = s.find("**", i + 2)
            if j != -1:
                flush()
                out.extend(tokenize(s[i + 2:j], True, italic))
                i = j + 2
                continue
        if s.startswith("[[", i):
            j = s.find("]]", i + 2)
            if j != -1:
                flush()
                out.append((s[i:j + 2], bold, italic, True))
                i = j + 2
                continue
        if s[i] == "*":
            j = s.find("*", i + 1)
            if j != -1 and j > i + 1:
                flush()
                out.extend(tokenize(s[i + 1:j], bold, True))
                i = j + 1
                continue
        if s[i] == "`":
            j = s.find("`", i + 1)
            if j != -1 and j > i + 1:
                flush()
                out.append((s[i + 1:j], bold, italic, False))
                i = j + 1
                continue
        buf.append(s[i])
        i += 1
    flush()
    return out


def split_vertical(text: str):
    """把一段文字按「上标 / 正常 / 下标」切段。"""
    segs, cur, mode = [], [], None
    for ch in text:
        m = "sup" if ch in SUP_CHARS else ("sub" if ch in SUB_CHARS else None)
        if m != mode:
            if cur:
                segs.append(("".join(cur), mode))
            cur, mode = [], m
        cur.append(ch)
    if cur:
        segs.append(("".join(cur), mode))
    return segs


def style_run(run, bold=False, italic=False, highlight=False, vert=None, size=None):
    run.font.name = FONT
    run._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    run.bold = bold
    run.italic = italic
    if size:
        run.font.size = Pt(size)
    if highlight:
        run.font.highlight_color = WD_COLOR_INDEX.YELLOW
    if vert == "sup":
        run.font.superscript = True
    elif vert == "sub":
        run.font.subscript = True
    return run


def add_runs(par, text: str, size=None, base_bold=False, colour=None):
    for seg_text, b, it, hl in tokenize(text, base_bold):
        for chunk, mode in split_vertical(seg_text):
            if not chunk:
                continue
            if mode == "sup":
                chunk = chunk.translate(SUP_MAP)
            elif mode == "sub":
                chunk = chunk.translate(SUB_MAP)
            r = par.add_run(chunk)
            style_run(r, bold=b or hl, italic=it, highlight=hl, vert=mode, size=size)
            if colour is not None:
                r.font.color.rgb = colour
    return par


# ---------------------------------------------------------------- 段落 / 表格

def new_par(doc, style=None, page_break=False, space_before=0, space_after=0,
            line_spacing=2.0, align=WD_ALIGN_PARAGRAPH.LEFT, indent_left=None,
            first_line=None):
    p = doc.add_paragraph(style=style)
    pf = p.paragraph_format
    pf.page_break_before = page_break
    pf.space_before = Pt(space_before)
    pf.space_after = Pt(space_after)
    pf.line_spacing = line_spacing
    pf.alignment = align
    if indent_left is not None:
        pf.left_indent = Inches(indent_left)
    if first_line is not None:
        pf.first_line_indent = Inches(first_line)
    return p


def add_table(doc, rows, size=TABLE_PT):
    if not rows:
        return None
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    t.autofit = True
    for ri, row in enumerate(rows):
        for ci, cell in enumerate(row):
            c = t.cell(ri, ci)
            c.text = ""
            p = c.paragraphs[0]
            pf = p.paragraph_format
            pf.space_before = Pt(1)
            pf.space_after = Pt(1)
            pf.line_spacing = 1.0
            add_runs(p, cell, size=size, base_bold=(ri == 0))
    # 表头跨页重复
    trPr = t.rows[0]._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader")
    th.set(qn("w:val"), "true")
    trPr.append(th)
    return t


# ---------------------------------------------------------------- 节 / 页码 / 行号

def configure_section(sec, landscape=False):
    sec.orientation = WD_ORIENT.LANDSCAPE if landscape else WD_ORIENT.PORTRAIT
    w, h = sec.page_width, sec.page_height
    long_side, short_side = max(w, h), min(w, h)
    if landscape:
        sec.page_width, sec.page_height = long_side, short_side
    else:
        sec.page_width, sec.page_height = short_side, long_side
    for attr in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, attr, Inches(1))


def add_line_numbers(sec):
    sectPr = sec._sectPr
    drop_line_numbers(sec)
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")
    sectPr.insert_element_before(
        ln, "w:pgNumType", "w:cols", "w:formProt", "w:vAlign", "w:noEndnote",
        "w:titlePg", "w:textDirection", "w:bidi", "w:rtlGutter", "w:docGrid",
        "w:printerSettings", "w:sectPrChange",
    )


def drop_line_numbers(sec):
    """关掉本节行号。ASN 官方：行号「止于参考文献之前」。

    行号是节级属性：sectPr 中不存在 w:lnNumType 即为关闭，无需写 countBy=0。
    """
    sectPr = sec._sectPr
    for el in sectPr.findall(qn("w:lnNumType")):
        sectPr.remove(el)


def add_footer_page_numbers(doc):
    for sec in doc.sections:
        try:
            sec.footer.is_linked_to_previous = False
        except Exception:
            pass
        p = sec.footer.paragraphs[0]
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.line_spacing = 1.0
        run = p.add_run()
        style_run(run, size=10)
        f1 = OxmlElement("w:fldChar")
        f1.set(qn("w:fldCharType"), "begin")
        it = OxmlElement("w:instrText")
        it.set(qn("xml:space"), "preserve")
        it.text = "PAGE"
        f2 = OxmlElement("w:fldChar")
        f2.set(qn("w:fldCharType"), "end")
        for el in (f1, it, f2):
            run._r.append(el)


# ---------------------------------------------------------------- Markdown 解析

def parse_blocks(md: str):
    """返回 [(kind, value)]，kind ∈ {h, p, table, ref, rule}。"""
    lines = md.split("\n")
    blocks, i = [], 0
    while i < len(lines):
        s = lines[i].strip()
        if not s or s.startswith("<!--"):
            i += 1
            continue
        if set(s) <= set("-—=*_ ") and len(s) >= 3:
            blocks.append(("rule", None))
            i += 1
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            blocks.append(("h", (len(m.group(1)), m.group(2).strip())))
            i += 1
            continue
        if s.startswith("|"):
            raw = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                raw.append(lines[i].strip())
                i += 1
            raw = [r for r in raw if not set(r) <= set("|-: ")]
            rows = [[c.strip() for c in r.strip("|").split("|")] for r in raw]
            blocks.append(("table", rows))
            continue
        if re.match(r"^\d+\.\s", s):
            blocks.append(("ref", s))
            i += 1
            continue
        buf = [s]
        i += 1
        while i < len(lines):
            s2 = lines[i].strip()
            if not s2 or s2.startswith(("#", "|", "<!--")):
                break
            if set(s2) <= set("-—=*_ ") and len(s2) >= 3:
                break
            if re.match(r"^\d+\.\s", s2):
                break
            buf.append(s2)
            i += 1
        blocks.append(("p", " ".join(buf)))
    return blocks


# 「本节之前分页」的 h2 标题
BREAK_BEFORE_H2 = {"Abstract", "Supplementary material"}
# 需要另起一页的 h3（不含首个表格）
BREAK_BEFORE_H3 = {"Table 2.", "Table 3."}


def render_markdown(md: str, title_for_docprops: str) -> Document:
    doc = Document()
    doc.core_properties.title = title_for_docprops
    doc.core_properties.author = "Jun-Lin Li"

    normal = doc.styles["Normal"]
    normal.font.name = FONT
    normal.font.size = Pt(BODY_PT)
    normal.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    normal.paragraph_format.line_spacing = 2.0
    normal.paragraph_format.space_after = Pt(0)
    normal.paragraph_format.space_before = Pt(0)

    configure_section(doc.sections[0], landscape=False)

    blocks = parse_blocks(md)
    first_table_seen = False
    # ASN 官方：行号须「自摘要起，止于参考文献之前」。标题页独立成节、不加行号。
    ln_stop = None
    ln_start = None

    for kind, val in blocks:
        if kind == "rule":
            continue

        if kind == "h":
            level, text = val
            if level == 2 and text.startswith("Abstract"):
                configure_section(doc.add_section(WD_SECTION.NEW_PAGE),
                                  landscape=False)
                ln_start = len(doc.sections) - 1
                p = new_par(doc, line_spacing=1.5, space_before=0)
                add_runs(p, text, base_bold=True)
                continue
            if level == 2 and text.startswith("References"):
                configure_section(doc.add_section(WD_SECTION.NEW_PAGE),
                                  landscape=False)
                ln_stop = len(doc.sections) - 1
                p = new_par(doc, line_spacing=1.5, space_before=0)
                add_runs(p, text, base_bold=True)
                continue
            if level == 2 and text.startswith("Tables"):
                configure_section(doc.add_section(WD_SECTION.NEW_PAGE), landscape=True)
                p = new_par(doc, line_spacing=1.5)
                add_runs(p, text, base_bold=True)
                continue
            if level == 2 and text.startswith("Figure legends"):
                configure_section(doc.add_section(WD_SECTION.NEW_PAGE), landscape=False)
                p = new_par(doc, line_spacing=1.5)
                add_runs(p, text, base_bold=True)
                continue

            brk = False
            if level == 2 and text in BREAK_BEFORE_H2:
                brk = True
            if level == 3:
                head = text.split(" ")[0] + " " if text else ""
                key = (text[:8])
                if any(text.startswith(k) for k in BREAK_BEFORE_H3) and first_table_seen:
                    brk = True
                if "Table" in text[:8]:
                    first_table_seen = True
            p = new_par(doc, line_spacing=1.5, page_break=brk,
                        space_before=10 if level == 2 else 8)
            add_runs(p, text, base_bold=True)
            continue

        if kind == "table":
            add_table(doc, val)
            p = new_par(doc, line_spacing=1.0, space_after=6)
            continue

        if kind == "ref":
            p = new_par(doc, line_spacing=2.0, indent_left=0.3, first_line=-0.3)
            add_runs(p, val)
            continue

        # 普通段落
        s = val
        if PLACEHOLDER.search(s):
            p = new_par(doc, line_spacing=2.0, space_after=4)
            add_runs(p, s)
            continue
        st = s.strip()
        if st.startswith("**") and st.endswith("**") and st.count("**") == 2:
            # 纯粗体行（面板标题 / 图注标题）
            p = new_par(doc, line_spacing=1.5, space_before=6, space_after=2)
            add_runs(p, st)
            continue
        p = new_par(doc, line_spacing=2.0)
        add_runs(p, s)

    add_footer_page_numbers(doc)
    for i, sec in enumerate(doc.sections):
        configure_section(sec, landscape=(sec.page_width > sec.page_height))
        in_window = ((ln_start is None or i >= ln_start)
                     and (ln_stop is None or i < ln_stop))
        if in_window:
            add_line_numbers(sec)
        else:
            drop_line_numbers(sec)   # 标题页、参考文献及其后：无行号
    return doc


# ---------------------------------------------------------------- 图件附录

FIGURES_SUBMISSION = [
    ("Figure1_study_design", "Figure 1. Study design and instrument selection.",
     "A analytical framework; B screening funnel (379 screened \u2192 153 with a genome-wide "
     "significant variant \u2192 86 analysable in the full chain \u2192 38 with \u2265 5 instruments); "
     "C distribution of the number of instruments per dietary exposure."),
    ("Figure2_diet_cholelithiasis", "Figure 2. Causal effects of dietary intake on cholelithiasis.",
     "A forest plot of the three exposures reaching FDR < 0.05 under MR-PRESSO; "
     "B plain IVW versus MR-PRESSO outlier-corrected estimates; "
     "C concordance of IVW, MR-PRESSO and weighted median across all exposures with \u2265 5 instruments."),
    ("Figure3_cholelithiasis_btc_chain", "Figure 3. Cholelithiasis, biliary tract cancer and the dietary chain.",
     "A variant-level Wald ratios with the IVW fit and 95% band; B IVW, MR-Egger, weighted median "
     "and weighted mode; C funnel plot with the MR-Egger intercept; "
     "D two-step mediation (paths a, b, c and the indirect effect a\u00d7b)."),
    ("Figure4_mvmr_and_failure_modes", "Figure 4. Multivariable MR and the two failure modes of standard analyses.",
     "A multivariable MR estimates by adjustment set; B conditional F-statistics with the "
     "identification threshold; C effect of region definition on colocalisation, with SuSiE "
     "independent-signal counts and conditional residual |z|."),
    ("Figure5_cross_ancestry_replication", "Figure 5. Cross-ancestry replication in BioBank Japan.",
     "A variant-level transferability of the instrument (slope 0.915, z = 11.17); "
     "B aggregate direction tests; C outcome panel in BioBank Japan; D power curves."),
]


def append_figure_appendix(doc: Document, root: Path):
    """把投稿图 1\u20135（results/submission/）按稿内顺序嵌进同一个 Word 文件。

    这取代了 2026-09-27 之前的旧版附录：旧版嵌入的是 5 张过程稿图，并在正文里声明
    「Figure 1 在本工作区内没有任何图片文件」「这些图不是投稿件」，且其中
    results/Figure_MR_diet_cholelithiasis_BTC.png 的图内脚注仍是旧口径
    （81 个 IV / 416 个联合 IV / 条件 F 22.2, 44.5），与现稿的 84 / 401 / 23.75, 45.11 冲突。
    现四张新制图 + 新制 Figure 5 已全部就绪，故改为嵌入真图。
    """
    configure_section(doc.add_section(WD_SECTION.NEW_PAGE), landscape=False)
    p = new_par(doc, line_spacing=1.5)
    add_runs(p, "Appendix. Figures 1\u20135 \u2014 embedded for internal review",
             base_bold=True)
    p = new_par(doc, line_spacing=1.5, space_after=6)
    add_runs(
        p,
        "*The American Journal of Clinical Nutrition requires figure files to be uploaded "
        "separately; these pages are included only so that the text, the tables and the figures "
        "of this submission can be reviewed in a single file. Delete this appendix before "
        "submission.*",
    )

    sub = root / "results" / "submission"
    for stem, legend, panels in FIGURES_SUBMISSION:
        path = sub / (stem + ".png")
        p = new_par(doc, line_spacing=1.3, page_break=True, space_after=4)
        add_runs(p, f"**{legend}**", base_bold=True)
        if path.exists():
            doc.add_picture(str(path), width=Inches(6.4))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        else:
            p = new_par(doc, line_spacing=1.5)
            add_runs(p, f"\u26a0 File missing: results/submission/{stem}.png",
                     base_bold=True, colour=RGBColor(0xC0, 0x00, 0x00))
        p = new_par(doc, line_spacing=1.3, space_before=4)
        add_runs(p, f"Panels: {panels}", size=10)
        p = new_par(doc, line_spacing=1.3)
        add_runs(p, f"File: results/submission/{stem}.png (400 dpi) and {stem}.pdf (vector)",
                 size=10)


# ---------------------------------------------------------------- 自检

def docx_text(path: Path) -> str:
    r"""把 docx 抽成纯文本（读 document.xml 的 w:t），用于结构自检。

    必须**逐段落**抽取：直接全篇 findall(w:t) 会丢掉段落边界，
    导致「^\d+\.」这类按行匹配的引用计数恒为 0。
    """
    import html as _html
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "ignore")
    out = []
    for para in re.findall(r"<w:p(?:\s[^>]*)?>.*?</w:p>", xml, flags=re.S):
        parts = re.findall(r"<w:t(?:\s[^>]*)?>(.*?)</w:t>", para, flags=re.S)
        out.append(_html.unescape("".join(parts)))
    return "\n".join(out)


def self_check(path: Path, expect_refs: int | None, expect_tables: int | None,
               kind="manuscript"):
    problems = []
    txt = docx_text(path)
    if kind == "manuscript":
        for key in ["Title Page", "Abstract", "1. Introduction", "2. Methods",
                    "3. Results", "4. Discussion", "Acknowledgments",
                    "Declaration of Generative AI", "References", "Tables",
                    "Figure legends", "Supplementary material"]:
            if key not in txt:
                problems.append(f"缺章节：{key}")
        if expect_refs:
            got = len(re.findall(r"(?m)^\s*\d+\.\s\S", txt))
            if got < expect_refs:
                problems.append(f"参考文献条目只找到 {got} 条，应为 {expect_refs}")
        for bad in ["[@", "**", "WORKING NOTE", "<!-- "]:
            if bad in txt:
                problems.append(f"正文残留标记：{bad}")
        # 占位符核对改为「动态」：先看源文件里实际还剩哪些 token，再确认 docx 里逐个保真。
        # 硬编码清单会在用户填完某项后误报「占位符丢失」（2026-09-20 实际踩到）。
        expected_ph = sorted(set(PLACEHOLDER.findall(MD_MANUSCRIPT.read_text(encoding="utf-8"))))
        for ph in expected_ph:
            if ph not in txt:
                problems.append(f"占位符丢失（应保留待填）：{ph}")
    if expect_tables:
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8", "ignore")
        got = xml.count("<w:tbl>")
        if got != expect_tables:
            problems.append(f"表格数 {got}，应为 {expect_tables}")
    return problems


# ---------------------------------------------------------------- main

def build_all(root: Path = ROOT):
    md = MD_MANUSCRIPT.read_text(encoding="utf-8")

    doc = render_markdown(md, "Dietary intake, cholelithiasis and biliary tract cancer")
    out_main = MAIN / "AJCN_Manuscript.docx"
    doc.save(out_main)
    print(f"[写出] {out_main.relative_to(root)}")

    doc2 = render_markdown(md, "Dietary intake, cholelithiasis and biliary tract cancer")
    append_figure_appendix(doc2, root)
    add_footer_page_numbers(doc2)
    for sec in doc2.sections:
        configure_section(sec, landscape=(sec.page_width > sec.page_height))
        add_line_numbers(sec)
    out_fig = MAIN / "AJCN_Manuscript_with_figures.docx"
    doc2.save(out_fig)
    print(f"[写出] {out_fig.relative_to(root)}")

    if MD_CHECKLIST.exists():
        ck = MD_CHECKLIST.read_text(encoding="utf-8")
        doc3 = render_markdown(ck, "STROBE-MR Checklist")
        # 清单表格宽，整篇转横向
        for sec in doc3.sections:
            configure_section(sec, landscape=True)
            add_line_numbers(sec)
        out_ck = MAIN / "AJCN_STROBE-MR_Checklist.docx"
        doc3.save(out_ck)
        print(f"[写出] {out_ck.relative_to(root)}")

    DESKTOP_OUT.mkdir(parents=True, exist_ok=True)
    for f in ["AJCN_Manuscript.docx", "AJCN_Manuscript_with_figures.docx",
              "AJCN_STROBE-MR_Checklist.docx"]:
        src = MAIN / f
        if src.exists():
            shutil.copy2(src, DESKTOP_OUT / f)
            print(f"[拷贝] ~/Desktop/投稿/AJCN/{f}")


def run_checks(root: Path = ROOT):
    print("\n=== 产出后自检 ===")
    ok = True
    for f, refs, tabs in [("AJCN_Manuscript.docx", 50, 9),
                          ("AJCN_Manuscript_with_figures.docx", 50, 9)]:
        p = MAIN / f
        if not p.exists():
            print(f"  ✗ 缺失 {f}")
            ok = False
            continue
        probs = self_check(p, refs, tabs, "manuscript")
        size = p.stat().st_size / 1024
        if probs:
            ok = False
            print(f"  ✗ {f}（{size:.0f} KB）")
            for x in probs:
                print(f"      - {x}")
        else:
            nph = len(set(PLACEHOLDER.findall(MD_MANUSCRIPT.read_text(encoding="utf-8"))))
            phmsg = f"占位符 {nph} 个保留" if nph else "无占位符（可投稿状态）"
            print(f"  ✓ {f}（{size:.0f} KB）章节齐、文献 50 条、表格 9 个、{phmsg}")
    ck = MAIN / "AJCN_STROBE-MR_Checklist.docx"
    if ck.exists():
        probs = self_check(ck, None, 1, "checklist")
        print(f"  {'✓' if not probs else '✗'} AJCN_STROBE-MR_Checklist.docx"
              f"（{ck.stat().st_size/1024:.0f} KB）" + ("" if not probs else f" {probs}"))
    else:
        print("  ✗ 缺失 AJCN_STROBE-MR_Checklist.docx")
        ok = False
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="只自检，不重建")
    a = ap.parse_args()
    if not a.check:
        build_all()
    sys.exit(0 if run_checks() else 1)
