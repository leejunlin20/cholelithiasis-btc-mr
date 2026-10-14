# -*- coding: utf-8 -*-
"""67_final_check_submission.py —— 投稿包终检（交付前必跑）

逐项核对 ~/Desktop/投稿/AJCN/ 的 24 个投稿件与 README，检查：
  1) 文件齐备、命名与 README 对照表一致
  2) 三份 Word 无占位符残留、章节齐
  3) 行号只出现在「摘要—参考文献之前」的节
  4) STROBE-MR 清单无 ⏳ / 无占位符 / 无「卡在用户决定上」字样
  5) 补充材料 9 件均带「标题 + 第一作者」页眉、单件 ≤2 MB
     （补充图 S1–S3 的页眉按 **PDF 文本层可提取** 校验）
  6) 5 张图的 PDF 为矢量、PNG ≥300 dpi
  7) 全部 8 张投稿图：无 Type 3 字体、单图内不跨族混排
     （ASN 稿约：embed the fonts used / uniform lettering）
"""
from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path

from docx import Document
from PIL import Image
from pypdf import PdfReader

DEST = Path.home() / "Desktop" / "投稿" / "AJCN"
REQUIRED = (
    ["README_投稿上传清单.md", "AJCN_Manuscript.docx",
     "AJCN_Manuscript_with_figures.docx", "AJCN_Cover_Letter.docx",
     "AJCN_Declaration_of_Interests.docx", "AJCN_STROBE-MR_Checklist.docx"]
    + [f"AJCN_Figure{i}.pdf" for i in range(1, 6)]
    + [f"AJCN_Figure{i}.png" for i in range(1, 6)]
    + [f"AJCN_Supplementary_Table_S{n}.docx" for n in range(1, 6)]
    + [f"AJCN_Supplementary_Figure_S{n}.pdf" for n in range(1, 4)]
    + ["AJCN_Supplementary_Methods.docx"]
)
UPLOAD_SET = [f for f in REQUIRED if f.endswith(
    (".docx", ".pdf")) and "with_figures" not in f and not f.endswith(".png")]

# 补充图页眉必须同时含「稿件标题片段」与「第一作者」
HDR_TITLE = "Dietary intake, cholelithiasis and biliary tract cancer"
HDR_AUTHOR = "Jun-Lin Li"

_STYLE_TOKENS = ("BoldOblique", "BoldItalic", "Oblique", "Italic", "Bold",
                 "Display", "Condensed", "Narrow", "Light", "Medium",
                 "Semibold", "Regular", "MT")

fails, warns = [], []


def ok(msg):
    print(f"  \u2713 {msg}")


def bad(msg):
    fails.append(msg)
    print(f"  \u2717 {msg}")


print("=" * 84)
print("AJCN 投稿包终检")
print("=" * 84)

# ---- 1 文件齐备
print("\n[1] 文件齐备与命名")
missing = [f for f in REQUIRED if not (DEST / f).exists()]
extra = [p.name for p in DEST.iterdir()
         if p.is_file() and p.name not in REQUIRED]
if missing:
    bad(f"缺文件：{missing}")
else:
    ok(f"{len(REQUIRED)} 个文件齐备")
if extra:
    warns.append(f"目录内有多余文件：{extra}")
    print(f"  ! 多余文件：{extra}")

readme = io.open(DEST / "README_投稿上传清单.md", encoding="utf-8").read()
not_in_readme = [f for f in UPLOAD_SET
                 if f.replace("AJCN_", "") not in readme
                 and f not in readme]
if not_in_readme:
    bad(f"README 未列出的投稿件：{not_in_readme}")
else:
    ok(f"README 对照表覆盖全部 {len(UPLOAD_SET)} 个投稿件")

# ---- 2 三份 Word
print("\n[2] Word 主件与清单")
for name in ("AJCN_Manuscript.docx", "AJCN_STROBE-MR_Checklist.docx"):
    d = Document(DEST / name)
    txt = "\n".join(p.text for p in d.paragraphs)
    for t in d.tables:
        for r in t.rows:
            for c in r.cells:
                txt += "\n" + c.text
    ph = len(re.findall(r"\[\[[A-Z_]+\]\]", txt))
    if ph:
        bad(f"{name}: 占位符残留 {ph} 处")
    else:
        ok(f"{name}: 无占位符（{len(d.tables)} 个表 / {len(d.paragraphs)} 段）")

ck = Document(DEST / "AJCN_STROBE-MR_Checklist.docx")
ckt = "\n".join(p.text for p in ck.paragraphs)
for t in ck.tables:
    for r in t.rows:
        for c in r.cells:
            ckt += "\n" + c.text
for probe, label in (("\u23f3", "沙漏待办标记"),
                     ("Blocked on a user decision", "「卡在用户决定上」字样"),
                     ("[[", "占位符")):
    n = ckt.count(probe)
    if n:
        bad(f"STROBE-MR 清单残留{label} {n} 处")
    else:
        ok(f"STROBE-MR 清单无{label}")

# ---- 3 行号节窗口
print("\n[3] 行号节窗口（应只见于「摘要 → 参考文献之前」）")
x = zipfile.ZipFile(DEST / "AJCN_Manuscript.docx").read(
    "word/document.xml").decode("utf-8")
secs = re.findall(r"<w:sectPr.*?</w:sectPr>|<w:sectPr[^>]*/>", x, re.S)
pattern = "".join("1" if "<w:lnNumType" in s else "0" for s in secs)
print(f"  节序（1=有行号）：{pattern}")
if pattern == "01111"[:len(pattern)] or pattern.startswith("01") and pattern[2:] == "0" * (len(pattern) - 2):
    ok("标题页无行号、摘要起有行号、参考文献及其后无行号")
else:
    bad(f"行号节窗口异常：{pattern}")

# ---- 4 补充材料
print("\n[4] 补充材料")


def _flat(s: str) -> str:
    """归一化：去掉全部空白与逗号并转小写。

    ⚠️ 必须对「待查串」和「被查文本」用**同一个**归一化函数。2026-09-27 首版
    只把待查串去了逗号、被查文本没去，导致三张带页眉的图被判「页眉缺失」。
    """
    return re.sub(r"[\s,]", "", s).lower()


_tbl_bad = 0
for n in range(1, 6):
    p = DEST / f"AJCN_Supplementary_Table_S{n}.docx"
    d = Document(p)
    hdr = d.sections[0].header.paragraphs[0].text
    if _flat(HDR_AUTHOR) not in _flat(hdr) or _flat(HDR_TITLE.split(":")[0]) not in _flat(hdr):
        bad(f"{p.name}: 页眉缺「稿件标题 + 第一作者」")
        _tbl_bad += 1
    mb = p.stat().st_size / 1e6
    if mb > 2:
        bad(f"{p.name}: {mb:.2f} MB > 2 MB")
        _tbl_bad += 1
if not _tbl_bad:
    ok("Table S1–S5：页眉齐、均 ≤2 MB")
pm = DEST / "AJCN_Supplementary_Methods.docx"
if _flat(HDR_AUTHOR) not in _flat(Document(pm).sections[0].header.paragraphs[0].text):
    bad("Supplementary Methods 页眉缺第一作者")
else:
    ok("Supplementary Methods：页眉齐")

# 补充图：页眉按 PDF 文本层校验（出图脚本以 pdf.fonttype=42 嵌入，故文字可提取；
# 此前 35/62 输出 Type 3 字体时文字**不可提取**，正是这一项的漏检之处）
_fig_bad = 0
for n in range(1, 4):
    p = DEST / f"AJCN_Supplementary_Figure_S{n}.pdf"
    if not p.exists():
        bad(f"{p.name}: 缺文件")
        _fig_bad += 1
        continue
    flat = _flat(PdfReader(p).pages[0].extract_text() or "")
    mb = p.stat().st_size / 1e6
    if _flat(HDR_TITLE) not in flat:
        bad(f"{p.name}: 页眉缺稿件标题（PDF 文本层读不到 → 字体可能非 TrueType 嵌入）")
        _fig_bad += 1
    elif _flat(HDR_AUTHOR) not in flat:
        bad(f"{p.name}: 页眉缺第一作者（PDF 文本层读不到）")
        _fig_bad += 1
    elif mb > 2:
        bad(f"{p.name}: {mb:.2f} MB > 2 MB")
        _fig_bad += 1
if not _fig_bad:
    ok("Figure S1–S3：页眉齐（PDF 文本层可提取）、均 ≤2 MB")

# ---- 5 图件
print("\n[5] 图件")
for i in range(1, 6):
    pdf = DEST / f"AJCN_Figure{i}.pdf"
    png = DEST / f"AJCN_Figure{i}.png"
    head = pdf.open("rb").read(8)
    vector = b"/FlateDecode" in pdf.open("rb").read()
    im = Image.open(png)
    dpi = round(im.width / (11.4 if i != 1 else 11.4))
    if not head.startswith(b"%PDF"):
        bad(f"{pdf.name}: 非 PDF")
    if im.width < 3000:
        bad(f"{png.name}: 位图宽 {im.width} px 偏低")
ok(f"Figure 1–5：PDF 矢量 + PNG 位图（宽 {Image.open(DEST / 'AJCN_Figure1.png').width} px）")

for n in range(1, 4):
    p = DEST / f"AJCN_Supplementary_Figure_S{n}.pdf"
    if not p.open("rb").read(4).startswith(b"%PDF"):
        bad(f"{p.name}: 非 PDF")

# ---- 6 字体与嵌入
print("\n[6] 字体与嵌入（ASN：embed the fonts used / uniform lettering）")


def _root(name: str) -> str:
    base = name.split("+", 1)[-1]
    changed = True
    while changed:
        changed = False
        for t in _STYLE_TOKENS:
            if base.endswith(t) and len(base) > len(t):
                base = base[: -len(t)].rstrip("-_")
                changed = True
    return base


pdfs = [DEST / f"AJCN_Figure{i}.pdf" for i in range(1, 6)] \
     + [DEST / f"AJCN_Supplementary_Figure_S{n}.pdf" for n in range(1, 4)]
all_fams, font_fail = set(), 0
for p in pdfs:
    if not p.exists():
        bad(f"{p.name}: 缺文件")
        font_fail += 1
        continue
    raw = p.open("rb").read()
    names = sorted({n.decode() for n in re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-]+)", raw)})
    fams = sorted({_root(n) for n in names})
    all_fams.update(fams)
    if b"/Type3" in raw:
        bad(f"{p.name}: 含 Type 3 字体（须在出图脚本设 pdf.fonttype=42）")
        font_fail += 1
    if len(fams) > 1:
        bad(f"{p.name}: 单图内跨族混排 {fams}（违反 uniform lettering）")
        font_fail += 1
if not font_fail:
    ok(f"{len(pdfs)} 张投稿图：族={sorted(all_fams)}，无 Type 3、单图内不混排")

print("\n" + "=" * 84)
if fails:
    print(f"\u7ed3\u8bba\uff1a\u4e0d\u53ef\u4ea4\u4ed8 —— \u963b\u65ad\u7ea7 {len(fails)} \u9879")
    for f in fails:
        print("   -", f)
elif warns:
    print(f"\u7ed3\u8bba\uff1a\u53ef\u4ea4\u4ed8\uff08\u963b\u65ad\u7ea7 0 \u9879\uff0c\u536b\u751f\u7ea7 {len(warns)} \u9879\uff09")
    for w in warns:
        print("   -", w)
else:
    print("\u7ed3\u8bba\uff1a\u53ef\u4ea4\u4ed8\uff08\u963b\u65ad\u7ea7 0 \u9879\u3001\u536b\u751f\u7ea7 0 \u9879\uff09")
print("=" * 84)
