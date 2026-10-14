#!/usr/bin/env python3
"""57_audit_figures.py —— 投稿图交付前机器体检

对 results/submission/ 下的每张图 + 对应出图脚本做四项检查，任何一项 FAIL 都不应交付。

  1) 分辨率：PNG 像素宽 / figsize 宽 = 有效 dpi，须 ≥ --min-dpi（默认 300）。
     踩坑史：某脚本多写了一次 savefig(out.replace(".pdf",".pdf"))，这是**空操作**，
     会把 400 dpi 的 PNG 用默认 100 dpi 原地覆盖一遍，且不报任何错。
  2) 内容边界：非白像素的外接框到画布四边的距离，须 ≥ --edge-margin 像素。
     贴边即说明文字/图形被画布切掉（典型：数值标签右对齐后超出画布）。
  3) 最小字号：解析脚本里所有 fontsize=，取最小值，换算到目标印幅后的等效 pt，
     须 ≥ --min-pt（默认 6.0）。equiv = min_fontsize × print_width / canvas_width。
     **注意**：画布与字号同比例缩放不改变 equiv（比值不变），所以发现 equiv 偏小时
     只能减文字/减内容或加大画布占比，不能靠"改 figsize"绕过。
  4) 画布一致性：同一批图的 figsize 应一致，不一致通常意味着有脚本被改坏。

用法：
  PY=python3
  $PY scripts/57_audit_figures.py
  $PY scripts/57_audit_figures.py --print-width 7.5 --min-pt 6.0 --min-dpi 300
"""
import argparse
import os
import re
import sys

from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUB = os.path.join(ROOT, "results", "submission")
SCRIPTS = os.path.join(ROOT, "scripts")

# 图 → 出图脚本
PAIRS = [
    ("Figure1_study_design", "53_make_figure1.py"),
    ("Figure2_diet_cholelithiasis", "54_make_figure2.py"),
    ("Figure3_cholelithiasis_btc_chain", "56_make_figure3.py"),
    ("Figure4_mvmr_and_failure_modes", "55_make_figure4.py"),
    ("Figure5_cross_ancestry_replication", "61_make_figure5.py"),
]


def _consts(s):
    """收集脚本里的数值常量，支持 `A = 1.0` 与 `A, B = 1.0, 2.0` 两种写法。"""
    c = {m.group(1): float(m.group(2)) for m in
         re.finditer(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([\d.]+)\s*$", s, re.M)}
    for m in re.finditer(r"^\s*((?:[A-Za-z_][A-Za-z0-9_]*\s*,\s*)*[A-Za-z_][A-Za-z0-9_]*)"
                         r"\s*=\s*((?:[\d.]+\s*,\s*)*[\d.]+)\s*$", s, re.M):
        names = [x.strip() for x in m.group(1).split(",")]
        vals = [x.strip() for x in m.group(2).split(",")]
        if len(names) == len(vals):
            for k, v in zip(names, vals):
                c.setdefault(k, float(v))
    return c


def figsize_of(py_path):
    """支持 figsize=(11.4, 7.2) 与 figsize=(FIG_W, FIG_H) 两种写法。"""
    s = open(py_path).read()
    m = re.search(r"figsize\s*=\s*\(\s*([^,()]+?)\s*,\s*([^,()]+?)\s*\)", s)
    if not m:
        return None
    c = _consts(s)
    out = []
    for g in (m.group(1).strip(), m.group(2).strip()):
        try:
            out.append(float(g))
        except ValueError:
            if g not in c:
                return None
            out.append(c[g])
    return (out[0], out[1])


def min_fontsize(py_path):
    """解析所有 fontsize=<num> 与 rcParams 的 font.size，取最小值。"""
    s = open(py_path).read()
    vals = [float(v) for v in re.findall(r"fontsize\s*=\s*([\d.]+)", s)]
    for v in re.findall(r'["\']font\.size["\']\s*:\s*([\d.]+)', s):
        vals.append(float(v))
    return min(vals) if vals else None


def content_bbox(png, white_thresh=250):
    """返回非白内容的外接框 (l, t, r, b)，按「到四边的距离」使用。"""
    im = Image.open(png).convert("L")
    W, H = im.size
    # 缩小以加速判定（阈值不变，边界精度足够）
    small = im.resize((W // 4 or 1, H // 4 or 1))
    px = small.load()
    w, h = small.size
    l, t, r, b = w, h, -1, -1
    for y in range(h):
        for x in range(w):
            if px[x, y] < white_thresh:
                if x < l:
                    l = x
                if x > r:
                    r = x
                if y < t:
                    t = y
                if y > b:
                    b = y
    if r < 0:
        return None
    k = 4
    return (l * k, t * k, W - 1 - r * k, H - 1 - b * k)


def leaked_paths(py_path):
    """粗查：图内文字是否混入了内部文件路径（投稿图不能出现 'results/…csv'）。"""
    lines = open(py_path).read().splitlines()
    hits = []
    for i, ln in enumerate(lines):
        if re.search(r"\.text\s*\(", ln):
            # 只把紧跟其后的、以引号开头的续行算作同一段显示文字；
            # 同段内的 cf["xxx.csv"] 之类字典键不算（曾误报 Figure 4）。
            chunk = [ln]
            j = i + 1
            while j < len(lines) and lines[j].lstrip()[:1] in ('"', "'"):
                chunk.append(lines[j])
                j += 1
            text = "\n".join(chunk)
            for m in re.findall(r"[\"'][^\"']*(?:results/|\.csv|\.tsv)[^\"']*[\"']", text):
                hits.append((i + 1, m.strip()))
    return hits


def figure_inventory(manuscript, sub):
    """解析稿件中引用的 Figure 编号，核对 results/submission/ 是否真有对应图件文件。

    ★ 这是三对齐校验的固有盲区：45_number_citations.py 只比对
      「正文引用编号 == 图注编号 == 图注物理顺序」，**不校验图件文件是否存在**。
      2026-09-27 实测：稿件正文与图注都引用 Figure 5（BBJ 跨族裔复制），
      三对齐判 5/5 通过，但 submission/ 只有 Figure 1–4 —— 整整一张图是空的。
    """
    txt = open(manuscript, encoding="utf-8").read()
    cited = sorted({int(n) for n in re.findall(r"\*\*Figure\s+(\d+)", txt)})
    present = set()
    for f in os.listdir(sub):
        m = re.match(r"Figure(\d+)_", f)
        if m:
            present.add(int(m.group(1)))
    return cited, sorted(set(cited) - present), sorted(present)


_STYLE_TOKENS = ("BoldOblique", "BoldItalic", "Oblique", "Italic", "Bold",
                 "Display", "Condensed", "Narrow", "Light", "Medium",
                 "Semibold", "Regular", "MT")


def font_root(name):
    """把 BaseFont 名归到「字体族根名」：去掉子集前缀与样式/变体后缀。

    DejaVuSansDisplay → DejaVuSans（同一 typeface 的显示变体，matplotlib 字形
    回退产物，视觉一致，不算跨族）；Arial-BoldMT → Arial。
    """
    base = name.split("+", 1)[-1]
    changed = True
    while changed:
        changed = False
        for t in _STYLE_TOKENS:
            if base.endswith(t) and len(base) > len(t):
                base = base[: -len(t)].rstrip("-_")
                changed = True
    return base


def font_audit(pdf_path):
    """返回 (是否含 Type 3 字体, 字体族根名集合, 全部 BaseFont 名)。

    ★ 2026-09-27 新增。ASN/Elsevier 稿约原文两条：
        • "Embed the fonts used if the application provides that option."
        • "Aim to use the following fonts in your illustrations: Arial, Courier,
           Times New Roman, Symbol, or use fonts that look similar."
        • "Make sure you use uniform lettering and sizing of your original artwork."
      稿约**未**明文禁止 Type 3（实测全文 "Type 3" 命中 0），但 matplotlib 既已
      提供 `pdf.fonttype=42` 这个"option"，就应当用；Type 3 是程序化字体，
      生产系统常拒收、且文字不可检索。故本体检把 Type 3 与「单图内跨族混排」
      一律视为阻断级。

      实测教训（2026-09-27）：
        • 缺 `pdf.fonttype` 的脚本（35、62）输出 Type 3 字体；
        • 把 font.family 写成 sans-serif 列表 ["Arial","Helvetica","DejaVu Sans"]
          时，上标负号等 Arial 缺失字形会静默回退 DejaVu → 同一张图里
          Arial 与 DejaVu **两族**并存（真实混排，须修）；
        • 而 DejaVuSans + DejaVuSansDisplay 属同一族的不同变体，**不是**混排。
    """
    raw = open(pdf_path, "rb").read()
    names = sorted({n.decode() for n in re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-]+)", raw)})
    fams = {font_root(n) for n in names}
    return (b"/Type3" in raw), fams, names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--print-width", type=float, default=6.9,
                    help="目标印幅（英寸）。ASN/AJCN 全文宽约 17.5 cm ≈ 6.9 in；"
                         "单栏 < 8.5 cm ≈ 3.35 in")
    ap.add_argument("--min-pt", type=float, default=8.0,
                    help="印幅缩小后的最小字号（pt）。ASN Guide for Authors 原文："
                         "reduced text size should be approximately 8 points")
    ap.add_argument("--min-dpi", type=float, default=300.0)
    ap.add_argument("--manuscript", default=os.path.join(
        ROOT, "results", "主稿", "Manuscript_v1.md"),
        help="用于核对「稿件引用的图 == submission/ 实有图件」的稿件文件")
    ap.add_argument("--edge-margin", type=int, default=6,
                    help="非白内容到画布边缘的最小像素距离")
    a = ap.parse_args()

    print("=" * 100)
    print(f"投稿图体检  目标印幅 {a.print_width} in | 字号下限 {a.min_pt} pt | "
          f"dpi 下限 {a.min_dpi} | 贴边容差 {a.edge_margin} px")
    print("=" * 100)
    fails, warns = [], []
    sizes = {}

    for stem, script in PAIRS:
        png = os.path.join(SUB, stem + ".png")
        py = os.path.join(SCRIPTS, script)
        if not os.path.exists(png):
            print(f"✗ {stem:<40} 缺文件")
            fails.append(f"{stem}: 缺 PNG")
            continue
        fs = figsize_of(py)
        mf = min_fontsize(py)
        im = Image.open(png)
        W, H = im.size
        sizes[stem] = (W, H)
        dpi = W / fs[0] if fs else float("nan")
        equiv = mf * a.print_width / fs[0] if (fs and mf) else float("nan")
        bb = content_bbox(png)
        margins = (min(bb),) if bb else (0,)

        ok_dpi = dpi >= a.min_dpi
        ok_pt = equiv >= a.min_pt
        # 同时报单栏场景（官方：most figures can be composed in a single column < 8.5 cm）
        equiv_narrow = mf * 3.35 / fs[0] if (fs and mf) else float("nan")
        ok_edge = bool(bb) and all(m >= a.edge_margin for m in bb)
        mark = lambda f: "✓" if f else "✗"
        print(f"{mark(ok_dpi and ok_pt and ok_edge)} {stem:<40} "
              f"{W}×{H}  {dpi:.0f}dpi (需≥{a.min_dpi:.0f}) | "
              f"最小字号 {mf:.1f}pt → 全宽({a.print_width}in)等效 {equiv:.2f}pt / "
              f"单栏(3.35in)等效 {equiv_narrow:.2f}pt (需≥{a.min_pt}) | "
              f"留边 左{bb[0]} 上{bb[1]} 右{bb[2]} 下{bb[3]}")
        if not ok_dpi:
            fails.append(f"{stem}: 有效 dpi {dpi:.0f} < {a.min_dpi:.0f}（多半是二次 savefig 覆盖）")
        if not ok_edge:
            fails.append(f"{stem}: 内容贴边 {bb} —— 文字或图形被画布切掉")
        if not ok_pt:
            warns.append(f"{stem}: 印幅等效最小字号 {equiv:.2f}pt < {a.min_pt}pt"
                         f"（画布宽 {fs[0] if fs else float('nan')} in）"
                         "→ 需减少文字量或加大面板占比")
        ls = leaked_paths(py)
        for ln, txt in ls:
            fails.append(f"{stem}: L{ln} 图内疑似印出内部路径 {txt}")
        if ls:
            print(f"    ↳ 图内疑似内部路径：{ls}")

    # ---- 字体与嵌入（ASN：embed the fonts used / uniform lettering）----
    SUPP = os.path.join(SUB, "supplementary")
    pdf_list = [(f"Figure{i}", os.path.join(SUB, s + ".pdf"))
                for i, (s, _) in enumerate(PAIRS, 1)]
    pdf_list += [(f"Supplementary Figure S{n}",
                  os.path.join(SUPP, f"Supplementary_Figure_S{n}.pdf"))
                 for n in (1, 2, 3)]
    print("-" * 100)
    print("字体与嵌入（ASN：embed the fonts used / uniform lettering）")
    for label, p in pdf_list:
        if not os.path.exists(p):
            fails.append(f"{label}: 缺 PDF {p}")
            print(f"  ✗ {label:<26} 缺 PDF")
            continue
        t3, fams, _names = font_audit(p)
        isbad = t3 or len(fams) > 1
        extra = f"  原始: {_names}" if (t3 or len(fams) > 1) else ""
        print(f"  {'✗' if isbad else '✓'} {label:<26} 族={sorted(fams)}  Type3={t3}{extra}")
        if t3:
            fails.append(f"{label}: 含 Type 3 字体（应在出图脚本设 pdf.fonttype=42）")
        if len(fams) > 1:
            fails.append(f"{label}: 单图内多字体族混排 {sorted(fams)}（违反 uniform lettering）")

    # ---- 图件清单核对：稿件引用的每张图都必须有实体文件 ----
    ms = a.manuscript
    if os.path.exists(ms):
        cited, missing, present = figure_inventory(ms, SUB)
        print(f"图件清单：稿件引用 Figure {cited} | submission/ 实有 Figure {present}")
        for k in missing:
            print(f"✗ Figure {k}: 稿件引用了，但 submission/ 下无 Figure{k}_* 文件")
            fails.append(f"Figure {k}: 稿件引用但无图件文件（三对齐校验查不出此项）")
        extra = [k for k in present if k not in cited]
        if extra:
            print(f"    ↳ 提示：submission/ 有未被稿件引用的图 Figure {extra}")
    else:
        print(f"! 未找到稿件 {ms}，跳过图件清单核对")
        warns.append(f"未找到稿件 {ms}，跳过图件清单核对")

    # 只比宽度：同一批图的图幅宽度应一致（高度按内容长短可以不同）
    widths = {k: v[0] for k, v in sizes.items()}
    if len(set(widths.values())) > 1:
        warns.append(f"图幅宽度不一致：{widths}")

    print("-" * 100)
    if fails:
        print(f"阻断级 {len(fails)} 项：")
        for f in fails:
            print("  ✗", f)
    else:
        print("✓ 无阻断级问题")
    if warns:
        print(f"卫生级 {len(warns)} 项：")
        for w in warns:
            print("  !", w)
    else:
        print("✓ 无卫生级问题")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
