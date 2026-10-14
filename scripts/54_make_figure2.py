#!/usr/bin/env python3
"""54_make_figure2.py —— Figure 2 膳食摄入对胆石症的因果效应

(A) 三个 FDR<0.05 暴露的森林图（MR-PRESSO 离群校正，主估计量）
(B) 同一三暴露：普通 IVW vs MR-PRESSO，展示剔除离群点的影响
(C) 全部 ≥5 工具暴露（38 个）的三估计量一致性散点

数据来源与自校验断言：
  results/full/mr_presso_results.csv          38 行 = ≥5 工具暴露全集
  results/full/mr_results_Cholelithiasis_v2.csv  IVW / 加权中位数（κ 修正版）
  results/medi2/mediation_paths_all.csv       86 暴露的 nsnp（定 ≥5 工具者）
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FULL = os.path.join(ROOT, "results", "full")
MEDI2 = os.path.join(ROOT, "results", "medi2")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
RED, GREEN, GREY, BLUE = "#C0392B", "#1E8449", "#566573", "#2471A3"
PURP, INK, LINE, GRID = "#7D3C98", "#212F3D", "#95A5A6", "#BDC3C7"


def L(p, **k):
    with open(p) as f:
        return list(csv.DictReader(f, **k))


def f(v, d=None):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return round(x, d) if d is not None else x


# ------------------------------------------------------------------ 数据
mp = {r["id"]: r for r in L(os.path.join(FULL, "mr_presso_results.csv"))}
v2 = L(os.path.join(FULL, "mr_results_Cholelithiasis_v2.csv"))
ivw = {r["id.exposure"]: r for r in v2 if r["method"] == "Inverse variance weighted"}
wmed = {r["id.exposure"]: r for r in v2 if r["method"] == "Weighted median"}
medi = L(os.path.join(MEDI2, "mediation_paths_all.csv"))
ge5 = {r["id"] for r in medi if int(r["nsnp"]) >= 5}

assert len(mp) == 38, f"MR-PRESSO 结果应为 38 个暴露，实为 {len(mp)}"
assert set(mp) == ge5, "MR-PRESSO 暴露集与 ≥5 工具暴露集不一致"

KEY = [("ukb-b-5779", "Alcohol intake frequency"),
       ("ukb-b-1489", "Cheese intake"),
       ("ukb-b-16576", "Dried fruit intake")]
shown = {
    "ukb-b-5779": ("Alcohol intake\nfrequency", 1.308, 1.190, 1.437, "3.1 × 10⁻⁷"),
    "ukb-b-1489": ("Cheese intake", 0.687, 0.576, 0.819, "9.8 × 10⁻⁵"),
    "ukb-b-16576": ("Dried fruit intake", 0.539, 0.407, 0.713, "1.3 × 10⁻⁴"),
}
# 与稿件 Table 2A 逐格核对
for eid, (nm, o, lo, hi, p) in shown.items():
    r = mp[eid]
    assert abs(f(r["or_corr"]) - o) < 0.001, f"{eid} MR-PRESSO OR 不符"
    assert abs(f(r["or_lci_corr"]) - lo) < 0.001, f"{eid} 下界不符"
    assert abs(f(r["or_uci_corr"]) - hi) < 0.001, f"{eid} 上界不符"


# MR-PRESSO 全局检验不显著时不出校正列（离群点数 0），此时校正估计 == 普通 IVW。
# 用 b_raw 回填，并断言回填值与 v2 的普通 IVW 完全一致——这同时是一次内部一致性检验。
def presso_est(eid):
    r = mp[eid]
    b, se = f(r["b_corr"]), f(r["se_corr"])
    if b is None:
        b, se = f(r["b_raw"]), f(r["se_raw"])
        assert int(r["nsnp_outlier"]) == 0, f"{eid} 有离群点却缺校正估计"
        assert abs(math.exp(b) - f(ivw[eid]["or"])) < 1e-6, \
            f"{eid} 无离群点时校正值应等于普通 IVW"
        return b, se, False
    return b, se, True


est = {eid: presso_est(eid) for eid in mp}
n_corrected = sum(1 for v in est.values() if v[2])
n_uncorrected = len(est) - n_corrected
assert n_corrected == 30 and n_uncorrected == 8, \
    f"离群校正分布异常: {n_corrected} 校正 / {n_uncorrected} 未校正"

# ------------------------------------------------------------------ 画布
fig = plt.figure(figsize=(11.4, 6.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.06], hspace=0.42, wspace=0.24,
                      left=0.052, right=0.975, top=0.945, bottom=0.085)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, :])

for ax in (axA, axB):
    ax.set_axis_off()
axA.set_xlim(0, 100)
axA.set_ylim(0, 100)
axB.set_xlim(0, 100)
axB.set_ylim(0, 100)

# ================================================================ Panel A
axA.text(0, 100, "A", fontsize=11.5, fontweight="bold", color=INK, va="top")
axA.text(6.8, 100, "Primary estimates (MR-PRESSO outlier-corrected)",
         fontsize=8.6, fontweight="bold", color=INK, va="top")

rows = list(shown.items())
y0, dy = 66.0, 22.0
xL, xR = 30.0, 84.0
lo_x, hi_x = 0.35, 2.60
def tx(v):
    return xL + (math.log(v) - math.log(lo_x)) / (math.log(hi_x) - math.log(lo_x)) * (xR - xL)

# 网格
for v in (0.5, 0.75, 1.0, 1.5, 2.0, 2.5):
    x = tx(v)
    axA.plot([x, x], [17, 84], color="#E5E8E8", lw=0.7, zorder=1)
    axA.text(x, 13.5, f"{v:g}", ha="center", va="top", fontsize=6.8, color=GREY)
axA.plot([tx(1.0), tx(1.0)], [17, 84], color=GREY, lw=0.9, ls="--", zorder=2)

for i, (eid, (nm, o, lo, hi, p)) in enumerate(rows):
    y = y0 - i * dy
    c = RED if o > 1 else BLUE
    axA.plot([tx(lo), tx(hi)], [y, y], color=c, lw=1.6, solid_capstyle="round", zorder=3)
    for xv in (lo, hi):
        axA.plot([tx(xv), tx(xv)], [y - 1.6, y + 1.6], color=c, lw=1.2, zorder=3)
    axA.plot([tx(o)], [y], "o", color=c, ms=6.2, zorder=4,
             markeredgecolor="white", markeredgewidth=0.9)
    axA.text(xR + 2.0, y + 1.8, f"{o:.3f}  ({lo:.3f}–{hi:.3f})", ha="left",
             va="center", fontsize=7.4, color=INK, zorder=4)
    axA.text(xR + 2.0, y - 2.6, f"P = {p}", ha="left", va="center",
             fontsize=6.6, color=GREY, zorder=4)
    axA.text(xL - 3.0, y, nm, ha="right", va="center", fontsize=7.6, color=INK,
             linespacing=1.4, zorder=4)
    axA.text(xL - 3.0, y - 7.4,
             f"{mp[eid]['nsnp']} IV · {mp[eid]['nsnp_outlier']} out",
             ha="right", va="center", fontsize=6.2, color=GREY, zorder=4)

axA.text(xL, 6.0, "Odds ratio per category increase (log scale)", ha="left",
         va="top", fontsize=6.9, color=GREY, style="italic")
axA.text(xR + 2.0, 88.5, "OR (95% CI)", ha="left", va="center", fontsize=6.9,
         color=GREY, style="italic")

# ================================================================ Panel B
axB.text(0, 100, "B", fontsize=11.5, fontweight="bold", color=INK, va="top")
axB.text(6.8, 100, "Effect of removing outlying variants", fontsize=8.6,
         fontweight="bold", color=INK, va="top")

bxL, bxW = 34.0, 44.0
def bx(v):
    return bxL + (math.log(v) - math.log(lo_x)) / (math.log(hi_x) - math.log(lo_x)) * bxW

for v in (0.5, 0.75, 1.0, 1.5, 2.0, 2.5):
    x = bx(v)
    axB.plot([x, x], [14, 80], color="#E5E8E8", lw=0.7, zorder=1)
    axB.text(x, 10.5, f"{v:g}", ha="center", va="top", fontsize=6.8, color=GREY)
axB.plot([bx(1.0), bx(1.0)], [14, 80], color=GREY, lw=0.9, ls="--", zorder=2)

for i, (eid, (nm, o, lo, hi, p)) in enumerate(rows):
    y = y0 - i * dy
    o_ivw = f(ivw[eid]["or"])
    lo_i, hi_i = f(ivw[eid]["or_lci95"]), f(ivw[eid]["or_uci95"])
    axB.plot([bx(lo_i), bx(hi_i)], [y, y], color=GREY, lw=1.5,
             solid_capstyle="round", zorder=3)
    axB.plot([bx(o_ivw)], [y], "s", color=GREY, ms=5.4, zorder=4,
             markeredgecolor="white", markeredgewidth=0.9)
    c = RED if o > 1 else BLUE
    axB.plot([bx(lo), bx(hi)], [y - 6.4, y - 6.4], color=c, lw=1.5,
             solid_capstyle="round", zorder=3)
    axB.plot([bx(o)], [y - 6.4], "o", color=c, ms=5.4, zorder=4,
             markeredgecolor="white", markeredgewidth=0.9)
    axB.add_patch(plt.matplotlib.patches.FancyArrowPatch(
        (bx(o_ivw), y - 1.0), (bx(o), y - 5.4), arrowstyle="-|>",
        mutation_scale=8, lw=0.9, color="#A6ACAF", zorder=2))
    axB.text(bxL - 3.0, y - 3.2, nm.replace("\n", " "), ha="right", va="center",
             fontsize=7.4, color=INK, wrap=True)
    axB.text(bxL + bxW + 3.0, y - 3.2,
             f"−{mp[eid]['nsnp_outlier']} out of {mp[eid]['nsnp']}", ha="left",
             va="center", fontsize=6.4, color=GREY)

axB.legend(handles=[
    Line2D([], [], marker="s", ls="-", color=GREY, ms=5.4,
           markerfacecolor=GREY, label="Plain IVW"),
    Line2D([], [], marker="o", ls="-", color=RED, ms=5.4,
           markerfacecolor=RED, label="MR-PRESSO outlier-corrected"),
], loc="upper right", bbox_to_anchor=(1.005, 0.985), ncol=2, fontsize=6.8,
    frameon=True, framealpha=0.95, edgecolor=GRID, borderpad=0.55,
    handlelength=1.5)
axB.text(bxL, 5.0, "Odds ratio (log scale)", ha="left", va="top",
         fontsize=6.9, color=GREY, style="italic")

# ================================================================ Panel C
axC.set_axis_on()
axC.text(-0.052, 1.135, "C", transform=axC.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axC.text(-0.052, 1.055, "Concordance of estimators across all exposures with ≥ 5 instruments",
         transform=axC.transAxes, fontsize=8.6, fontweight="bold", color=INK,
         va="top", ha="left")

xr, yr, xw, yw = [], [], [], []
for eid in mp:
    b, se, _ = est[eid]
    m = math.exp(b)
    i_ = f(ivw[eid]["or"]) if eid in ivw else None
    w_ = f(wmed[eid]["or"]) if eid in wmed else None
    if m <= 0:
        continue
    if i_ and i_ > 0:
        xr.append(math.log(m)); yr.append(math.log(i_))
    if w_ and w_ > 0:
        xw.append(math.log(m)); yw.append(math.log(w_))

lim = (-1.85, 1.85)
axC.plot(lim, lim, color=GREY, lw=0.9, ls="--", zorder=2)
axC.axhline(0, color="#E5E8E8", lw=0.8, zorder=1)
axC.axvline(0, color="#E5E8E8", lw=0.8, zorder=1)
axC.scatter(xr, yr, s=26, marker="o", facecolor="none", edgecolor=BLUE,
            linewidth=0.9, zorder=3, label="Plain IVW")
axC.scatter(xw, yw, s=30, marker="^", facecolor="none", edgecolor=GREEN,
            linewidth=0.9, zorder=3, label="Weighted median")

for eid, (nm, o, lo, hi, p) in shown.items():
    b, _, _ = est[eid]
    m = b
    i_ = math.log(f(ivw[eid]["or"]))
    axC.plot([m], [i_], "o", color=RED, ms=6.0, zorder=5,
             markeredgecolor="white", markeredgewidth=0.9)
    axC.annotate(nm.replace("\n", " "), (m, i_), xytext=(6, 7),
                 textcoords="offset points", fontsize=6.6, color=RED,
                 fontweight="bold", zorder=6)

axC.set_xlim(*lim)
axC.set_ylim(*lim)
axC.set_xticks([math.log(v) for v in (0.2, 0.35, 0.6, 1.0, 1.7, 3.0, 5.0)])
axC.set_xticklabels(["0.2", "0.35", "0.6", "1.0", "1.7", "3.0", "5.0"], fontsize=7.2)
axC.set_yticks([math.log(v) for v in (0.2, 0.35, 0.6, 1.0, 1.7, 3.0, 5.0)])
axC.set_yticklabels(["0.2", "0.35", "0.6", "1.0", "1.7", "3.0", "5.0"], fontsize=7.2)
axC.set_xlabel("OR under MR-PRESSO outlier correction", fontsize=7.8)
axC.set_ylabel("OR under the comparator estimator", fontsize=7.8)
axC.tick_params(length=2.5, color=LINE)
axC.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axC.spines[s].set_color(LINE)
axC.legend(loc="lower right", fontsize=7.0, frameon=True, framealpha=0.95,
           edgecolor=GRID, borderpad=0.55, handletextpad=0.5)
axC.text(0.015, 0.955,
         f"38 exposures with ≥ 5 instruments\n"
         f"{n_corrected} with outlying variants removed · "
         f"{n_uncorrected} with none detected (corrected = plain IVW)\n"
         "dashed line = identity",
         transform=axC.transAxes, ha="left", va="top", fontsize=6.6,
         color=GREY, linespacing=1.6,
         bbox=dict(boxstyle="round,pad=0.42", facecolor="#F8F9F9",
                   edgecolor=GRID, linewidth=0.6))

# ------------------------------------------------------------------ 标题
# ⚠ 官方口径（ASN Guide for Authors 原文）："Compile figure titles and legends on one
# or more pages in the manuscript's text file rather than on the figure itself."
# 故此处**不再在图内印主标题与副标题**；图题与图注已写入稿件文本（Figure 1–5 图注）。
# 面板字母 A/B/C/D 与图例键保留在图内（官方允许 "figure keys on figures"）。

out = os.path.join(ROOT, "results", "submission", "Figure2_diet_cholelithiasis.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=400, facecolor="white")
fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
print("Figure 2 已写出:", out)
print(f"  38 暴露 · IVW 配对 {len(xr)} · 加权中位数配对 {len(xw)}")
print(f"  三估计量与稿件 Table 2A 断言全部通过")
