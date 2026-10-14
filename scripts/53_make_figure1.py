#!/usr/bin/env python3
"""53_make_figure1.py —— Figure 1 研究设计与工具变量筛选

面板：
  (A) 分析框架：饮食暴露 → 胆石症 → 胆道癌 + 五类分析
  (B) 筛选漏斗：379 → 153 → 86 → 38
  (C) 每个膳食暴露的工具变量数分布（86 个可分析暴露）

同时承担 STROBE-MR 第 10a 项（筛选流程）的流程图职责。

数据来源（脚本内硬断言，任一数字与稿件不符即失败）：
  results/full/exposure_log.csv           379 行，status ∈ {ok, no_iv}
  results/medi2/mediation_paths_all.csv   86 个可分析暴露的 nsnp
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FULL = os.path.join(ROOT, "results", "full")
MEDI2 = os.path.join(ROOT, "results", "medi2")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

RED = "#C0392B"
GREEN = "#1E8449"
GREY = "#566573"
BLUE = "#2471A3"
PURP = "#7D3C98"
ORNG = "#D68910"
SAND, MIST, LEAF = "#FDF2E3", "#EAF2F8", "#EAF7EF"
INK, LINE = "#212F3D", "#95A5A6"
GRID = "#BDC3C7"


def load(p, **kw):
    with open(p) as f:
        return list(csv.DictReader(f, **kw))


# ---------------------------------------------------------------- 数据 + 断言
log = load(os.path.join(FULL, "exposure_log.csv"))
medi = load(os.path.join(MEDI2, "mediation_paths_all.csv"))

n_screened = len(log)
n_iv = sum(1 for r in log if r["status"] == "ok")
n_analysed = len(medi)
nsnp = np.array([int(r["nsnp"]) for r in medi])
n_ge5 = int((nsnp >= 5).sum())
median_nsnp = float(np.median(nsnp))

assert n_screened == 379, f"筛选起始数应为 379，实为 {n_screened}"
assert n_iv == 153, f"有工具变量的暴露应为 153，实为 {n_iv}"
assert n_analysed == 86, f"可分析暴露应为 86，实为 {n_analysed}"
assert n_ge5 == 38, f">=5 工具的暴露应为 38，实为 {n_ge5}"
assert median_nsnp == 3.0, f"工具数中位数应为 3，实为 {median_nsnp}"

# ---------------------------------------------------------------- 画布
FIG_W, FIG_H = 11.4, 7.05
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = fig.add_gridspec(2, 2, height_ratios=[1.16, 1.0],
                      hspace=0.40, wspace=0.24,
                      left=0.052, right=0.975, top=0.945, bottom=0.072)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])
for ax in (axA, axB):
    ax.set_axis_off()

# ================================================================ Panel A
axA.set_xlim(0, 100)
axA.set_ylim(0, 50)

axA.text(0, 50, "A", fontsize=11.5, fontweight="bold", color=INK, va="top",
         ha="left")
axA.text(3.4, 50, "Analytical framework", fontsize=9.2, fontweight="bold",
         color=INK, va="top", ha="left")

BOX_Y, BOX_H, BOX_W = 30.5, 16.0, 27.2
X1, X2, X3 = 1.0, 36.4, 71.8


def stage_box(ax, x, y, w, h, fc, ec, title, lines, tsize=8.0, lsize=6.7):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.0,rounding_size=1.0",
                                linewidth=1.0, edgecolor=ec, facecolor=fc, zorder=2))
    ax.text(x + w / 2, y + h - 1.35, title, ha="center", va="top",
            fontsize=tsize, fontweight="bold", color=INK, zorder=3)
    ax.text(x + w / 2, y + h - 4.35, "\n".join(lines), ha="center", va="top",
            fontsize=lsize, color=GREY, linespacing=1.55, zorder=3)


stage_box(axA, X1, BOX_Y, BOX_W, BOX_H, SAND, ORNG,
          "UK Biobank — dietary exposures",
          ["379 phenotypes screened", "n \u2248 450,000 · European ancestry",
           "median 3 instruments per exposure"])
stage_box(axA, X2, BOX_Y, BOX_W, BOX_H, MIST, BLUE,
          "FinnGen R12 — cholelithiasis",
          ["49,834 cases / 437,418 controls", "ICD-10 K80",
           "84 instruments · mean F 122.6, min 31.2"])
stage_box(axA, X3, BOX_Y, BOX_W, BOX_H, LEAF, GREEN,
          "FinnGen R12 — biliary tract cancer",
          ["2,298 cases / 378,749 controls", "ICD-10 C23\u2013C24",
           "primary outcome"])


def harrow(ax, x1, x2, y, label, color=RED):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                 mutation_scale=11, linewidth=1.15, color=color,
                                 zorder=2, shrinkA=0, shrinkB=0))
    ax.text((x1 + x2) / 2, y + 0.9, label, ha="center", va="bottom",
            fontsize=7.6, fontweight="bold", color=color, zorder=3)


harrow(axA, X1 + BOX_W + 0.3, X2 - 0.3, BOX_Y + BOX_H / 2, "MR  a")
harrow(axA, X2 + BOX_W + 0.3, X3 - 0.3, BOX_Y + BOX_H / 2, "MR  b")

# 直接效应 c′（虚线，绕行下方）
DY = 24.6
axA.add_patch(FancyArrowPatch((X1 + BOX_W / 2, BOX_Y - 0.3), (X1 + BOX_W / 2, DY),
                              arrowstyle="-", linewidth=1.0, ls=(0, (3.5, 2.2)),
                              color=GREY, zorder=1))
axA.add_patch(FancyArrowPatch((X1 + BOX_W / 2, DY), (X3 + BOX_W / 2, DY),
                              arrowstyle="-", linewidth=1.0, ls=(0, (3.5, 2.2)),
                              color=GREY, zorder=1))
axA.add_patch(FancyArrowPatch((X3 + BOX_W / 2, DY), (X3 + BOX_W / 2, BOX_Y - 0.3),
                              arrowstyle="-|>", mutation_scale=9, linewidth=1.0,
                              ls=(0, (3.5, 2.2)), color=GREY, zorder=1))
axA.text((X1 + X3 + BOX_W) / 2, DY + 0.75,
         "direct effect  c\u2032  \u2014  estimable only in the multivariate model",
         ha="center", va="bottom", fontsize=6.7, color=GREY, style="italic")

# 方法条
methods = [
    ("Univariable MR", "IVW · MR-Egger\nweighted median / mode", "#F8F9F9", GRID),
    ("Multivariable MR", "adjustment for\nBMI / cholecystitis", "#F8F9F9", GRID),
    ("Colocalisation", "coloc.abf\n+ coloc.susie", "#F8F9F9", GRID),
    ("Two-step mediation", "diet \u2192 cholelithiasis\n\u2192 BTC", "#F8F9F9", GRID),
    ("Replication", "BioBank Japan\n+ FinnGen panel", "#F5EEF8", PURP),
]
M_W, M_GAP, M_Y, M_H = 18.0, 1.75, 2.0, 12.2
for i, (t, s, fc, ec) in enumerate(methods):
    x = 1.0 + i * (M_W + M_GAP)
    axA.add_patch(FancyBboxPatch((x, M_Y), M_W, M_H,
                                 boxstyle="round,pad=0.0,rounding_size=0.8",
                                 linewidth=0.75, edgecolor=ec, facecolor=fc,
                                 zorder=2))
    axA.text(x + M_W / 2, M_Y + M_H - 1.8, t, ha="center", va="top", fontsize=6.9,
             fontweight="bold", color=PURP if ec == PURP else INK, zorder=3)
    axA.text(x + M_W / 2, M_Y + M_H - 4.9, s, ha="center", va="top", fontsize=6.0,
             color=GREY, linespacing=1.5, zorder=3)

# ================================================================ Panel B
axB.set_xlim(0, 100)
axB.set_ylim(0, 100)

axB.text(0, 100, "B", fontsize=11.5, fontweight="bold", color=INK, va="top",
         ha="left")
axB.text(6.6, 100, "Screening funnel for dietary exposures", fontsize=9.2,
         fontweight="bold", color=INK, va="top", ha="left")

stages = [
    (379, "diet-related UK Biobank phenotypes screened", None, SAND, ORNG, 1.00),
    (153, "with \u2265 1 genome-wide significant variant",
     "\u2212226  no instrument", MIST, BLUE, 0.63),
    (86, "estimable on both cholelithiasis and BTC",
     "\u221267  not estimable on both", LEAF, GREEN, 0.41),
    (38, "with \u2265 5 instruments \u2014 outlier testing",
     "\u221248  fewer than 5", "#F5EEF8", PURP, 0.25),
]
BAR_H, VGAP, BAR_X, BAR_MAXW = 15.5, 8.0, 5.0, 30.0
TOP = 95.0
TXT_X = 38.5
for i, (n, lab, note, fc, ec, wf) in enumerate(stages):
    y = TOP - i * (BAR_H + VGAP) - BAR_H
    w = BAR_MAXW * wf
    axB.add_patch(FancyBboxPatch((BAR_X, y), w, BAR_H,
                                 boxstyle="round,pad=0.0,rounding_size=0.9",
                                 linewidth=1.0, edgecolor=ec, facecolor=fc,
                                 zorder=2))
    axB.text(BAR_X + w / 2, y + BAR_H / 2, str(n), ha="center", va="center",
             fontsize=13.0, fontweight="bold", color=ec, zorder=3)
    ytxt = y + BAR_H / 2 + (1.9 if note else 0.0)
    axB.text(TXT_X, ytxt, lab, ha="left", va="center", fontsize=7.3, color=INK,
             zorder=3)
    if note:
        axB.text(TXT_X, y + BAR_H / 2 - 2.9, note, ha="left", va="center",
                 fontsize=6.4, color=RED, zorder=3)
    if i < len(stages) - 1:
        axB.add_patch(FancyArrowPatch((BAR_X + 3.6, y - 0.9),
                                      (BAR_X + 3.6, y - VGAP + 0.9),
                                      arrowstyle="-|>", mutation_scale=9,
                                      linewidth=1.0, color=GREY, zorder=2))

# ================================================================ Panel C
axC.set_axis_on()
bins = [(2, 4), (5, 9), (10, 19), (20, 49), (50, 99)]
labels = ["2\u20134", "5\u20139", "10\u201319", "20\u201349", "50\u201399"]
counts = [int(((nsnp >= a) & (nsnp <= b)).sum()) for a, b in bins]
colors = [GREY, BLUE, "#5DADE2", GREEN, PURP]
assert sum(counts) == 86, f"分组计数合计应为 86，实为 {sum(counts)}"

xs = np.arange(len(bins))
axC.bar(xs, counts, width=0.66, color=colors, edgecolor="white", linewidth=0.8,
        zorder=3)
for x, c in zip(xs, counts):
    axC.text(x, c + 1.3, str(c), ha="center", va="bottom", fontsize=8.2,
             fontweight="bold", color=INK, zorder=4)
axC.set_xticks(xs)
axC.set_xticklabels(labels, fontsize=7.6)
axC.set_xlabel("Number of instruments per exposure", fontsize=7.8)
axC.set_ylabel("Dietary exposures (n)", fontsize=7.8)
axC.set_ylim(0, max(counts) * 1.26)
axC.tick_params(labelsize=7.4, length=2.5, color=LINE)
axC.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axC.spines[s].set_color(LINE)
axC.text(0.0, 1.115, "C", transform=axC.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axC.text(0.0, 1.045, "Instruments per dietary exposure", transform=axC.transAxes,
         fontsize=9.2, fontweight="bold", color=INK, va="top", ha="left")
axC.text(0.985, 0.955,
         "median 3  ·  IQR 2\u201318  ·  range 2\u201399\n"
         "38 of 86 eligible for outlier testing",
         transform=axC.transAxes, ha="right", va="top", fontsize=6.6,
         color=GREY, linespacing=1.5,
         bbox=dict(boxstyle="round,pad=0.40", facecolor="#F8F9F9",
                   edgecolor=GRID, linewidth=0.6))

# ================================================================ 总标题
# ⚠ 官方口径（ASN Guide for Authors 原文）："Compile figure titles and legends on one
# or more pages in the manuscript's text file rather than on the figure itself."
# 故此处**不再在图内印主标题与副标题**；图题与图注已写入稿件文本（Figure 1–5 图注）。
# 面板字母 A/B/C/D 与图例键保留在图内（官方允许 "figure keys on figures"）。

out = os.path.join(ROOT, "results", "submission", "Figure1_study_design.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=400, facecolor="white")
fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
print("Figure 1 已写出:", out)
print(f"  漏斗: {n_screened} → {n_iv} → {n_analysed} → {n_ge5}")
print(f"  分布: {dict(zip(labels, counts))}  median={median_nsnp}")
