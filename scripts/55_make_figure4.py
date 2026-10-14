#!/usr/bin/env python3
"""55_make_figure4.py —— Figure 4 多变量 MR 与标准分析的两类失效模式

(A) 四个校正集下胆石症→胆道癌的 MVMR 估计，按模型是否可识别着色
(B) 各校正集的条件 F 统计量，标示识别阈值 F = 10
(C) 区域定义对共定位的影响：69 个区域的 PP.H4 + R07 合并/拆分对照 + SuSiE 诊断

数据来源与自校验断言：
  results/mvmr/mvmr_ivw_chole_bmi.csv            主模型（401 IV）
  results/mvmr/mvmr_ivw_chole_cholecyst.csv      胆囊炎模型（86 IV）
  results/mvmr/mvmr_robust_chole_cholecyst.csv   同模型 median / Lasso
  results/mvmr/mvmr_ivw_with_bmi.csv             四暴露模型（401 IV）
  results/mvmr/mvmr_condF_*.csv                  条件 F
  results/coloc/coloc_main_p12_1e-5.csv          69 区域 PP.H4
  results/coloc/coloc_R07_split.csv              R07 合并 vs 拆分
  results/coloc_susie/independent_signals.tsv    SuSiE 独立信号数与条件残差 |z|
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
RED, GREEN, GREY, BLUE = "#C0392B", "#1E8449", "#566573", "#2471A3"
PURP, ORNG, INK, LINE, GRID = "#7D3C98", "#D68910", "#212F3D", "#95A5A6", "#BDC3C7"


def L(p, **k):
    with open(p) as f:
        return list(csv.DictReader(f, **k))


def orci(b, se):
    return math.exp(b), math.exp(b - 1.959964 * se), math.exp(b + 1.959964 * se)


# ------------------------------------------------------------------ 数据
def one(fn, expo):
    for r in L(os.path.join(R, "mvmr", fn)):
        if r["exposure"] == expo:
            return float(r["b"]), float(r["se"])
    raise KeyError(fn)


sets = []
b, se = one("mvmr_ivw_chole_bmi.csv", "cholelith")
sets.append(("+ BMI\n(primary model)", "401", (b, se), True))
b, se = one("mvmr_ivw_chole_cholecyst.csv", "cholelith")
sets.append(("+ cholecystitis\n(IVW)", "86", (b, se), False))
# 每个方法有两行（cholelith / cholecyst），必须按暴露筛选，否则会取到胆囊炎那一行
rb = {r["method"].split(" (")[0]: r
      for r in L(os.path.join(R, "mvmr", "mvmr_robust_chole_cholecyst.csv"))
      if r["exposure"] == "cholelith"}
assert set(rb) >= {"MVMR-median", "MVMR-Lasso"}, f"稳健估计方法缺失: {sorted(rb)}"
for meth, name in (("MVMR-median", "+ cholecystitis\n(weighted median)"),
                   ("MVMR-Lasso", "+ cholecystitis\n(Lasso)")):
    r = rb[meth]
    sets.append((name, "86", (float(r["b"]), float(r["se"])), False))
b, se = one("mvmr_ivw_with_bmi.csv", "cholelith")
sets.append(("+ cholecystitis\n+ cholecystectomy\n+ BMI", "401", (b, se), False))

TABLE3A = [1.244, 1.13, 1.24, 1.35, 1.79]
for (nm, n_iv, (bb, ss), ok), want in zip(sets, TABLE3A):
    o, lo, hi = orci(bb, ss)
    assert abs(o - want) < 0.012, f"{nm}: OR {o:.3f} 与 Table 3A 的 {want} 不符"

condF = {
    "mvmr_condF_chole_bmi.csv": ["cholelithiasis", "BMI"],
    "mvmr_condF_chole_cholecyst.csv": ["cholelithiasis", "cholecystitis"],
    "mvmr_condF_with_bmi.csv": ["cholelithiasis", "cholecystitis", "cholecystectomy", "BMI"],
}
cf = {}
for fn, names in condF.items():
    r = L(os.path.join(R, "mvmr", fn))[0]
    cf[fn] = (names, [float(r[f"exposure{i+1}"]) for i in range(len(names))
                      if r.get(f"exposure{i+1}")])
assert abs(cf["mvmr_condF_chole_bmi.csv"][1][0] - 23.7529) < 1e-3

coloc = [r for r in L(os.path.join(R, "coloc", "coloc_main_p12_1e-5.csv"))
         if r["pair"] == "CHOLE_vs_BTC"]
coloc.sort(key=lambda r: -float(r["PP.H4"]))
assert len(coloc) == 69, f"共定位区域应为 69（R01–R69），实为 {len(coloc)}"

split = {(r["region"], r["p12"]): r
         for r in L(os.path.join(R, "coloc", "coloc_R07_split.csv"))}
R07 = {"merged": float(split[("R07", "1e-05")]["PP.H4"]),
       "R07a": float(split[("R07a", "1e-05")]["PP.H4"]),
       "R07b": float(split[("R07b", "1e-05")]["PP.H4"])}
assert abs(R07["R07a"] - 0.1877) < 1e-3 and abs(R07["R07b"] - 0.1158) < 1e-3

sig = {}
for r in L(os.path.join(R, "coloc_susie", "independent_signals.tsv"), delimiter="\t"):
    if r["r2_threshold"] == "0.1" and r["trait"] == "Cholelithiasis":
        sig[r["region"]] = (int(r["n_independent"]),
                            float(r["max_abs_z_cond"]) if r["max_abs_z_cond"] else None)
assert sig["R07"][0] == 19, f"R07 独立信号数应为 19，实为 {sig['R07'][0]}"

# ------------------------------------------------------------------ 画布
fig = plt.figure(figsize=(11.4, 7.35))
gs = fig.add_gridspec(2, 2, height_ratios=[0.94, 1.0], hspace=0.50, wspace=0.20,
                      left=0.105, right=0.975, top=0.938, bottom=0.072)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
sub = gs[1, 1].subgridspec(2, 1, height_ratios=[1.30, 1.0], hspace=0.95)
axC = fig.add_subplot(sub[0])
axC2 = fig.add_subplot(sub[1])
axA.set_axis_off()

# ================================================================ Panel A
axA.set_xlim(0, 100)
axA.set_ylim(0, 116)
axA.text(0, 116, "A", fontsize=11.5, fontweight="bold", color=INK, va="top")
axA.text(6.0, 116, "Multivariable MR of cholelithiasis on biliary tract cancer",
         fontsize=8.6, fontweight="bold", color=INK, va="top")
axA.scatter([70.0], [112.0], s=38, marker="o", facecolor=GREEN, edgecolor=GREEN, zorder=3)
axA.text(73.0, 112.0, "identified", ha="left", va="center", fontsize=6.4, color=GREY)
axA.scatter([86.0], [112.0], s=38, marker="o", facecolor="white", edgecolor=RED,
            linewidth=1.1, zorder=3)
axA.text(89.0, 112.0, "not identified", ha="left", va="center", fontsize=6.4, color=GREY)

xL, xR = 26.0, 76.0
lo_x, hi_x = 0.72, 3.5
def tx(v):
    return xL + (math.log(v) - math.log(lo_x)) / (math.log(hi_x) - math.log(lo_x)) * (xR - xL)

for v in (0.8, 1.0, 1.5, 2.0, 3.0):
    x = tx(v)
    axA.plot([x, x], [12, 100], color="#E5E8E8", lw=0.7, zorder=1)
    axA.text(x, 8.5, f"{v:g}", ha="center", va="top", fontsize=6.8, color=GREY)
axA.plot([tx(1.0), tx(1.0)], [12, 100], color=GREY, lw=0.9, ls="--", zorder=2)

y0, dy = 90.0, 19.0
for i, (nm, n_iv, (bb, ss), ok) in enumerate(sets):
    y = y0 - i * dy
    o, lo, hi = orci(bb, ss)
    c = GREEN if ok else RED
    axA.plot([tx(lo), tx(hi)], [y, y], color=c, lw=1.6, solid_capstyle="round", zorder=3)
    axA.plot([tx(o)], [y], "o", color=c if ok else "white", ms=6.2, zorder=4,
             markeredgecolor=c, markeredgewidth=1.2)
    axA.text(xR + 2.4, y, f"{o:.2f} ({lo:.2f}–{hi:.2f})", ha="left", va="center",
             fontsize=7.2, color=INK, zorder=4)
    nm1 = nm.replace("\n", "  ·  ") if "\n" in nm else nm
    axA.text(xL - 2.4, y + 2.9, nm1.replace("  ·  ", "  "), ha="right", va="center",
             fontsize=7.1, color=INK, zorder=4)
    axA.text(xL - 2.4, y - 3.4, f"{n_iv} instruments", ha="right", va="center",
             fontsize=6.2, color=GREY, zorder=4)
    if not ok:
        axA.text(tx(o), y - 3.4, "not identified", ha="center", va="center",
                 fontsize=6.0, color=RED, style="italic", zorder=4)

axA.text(xL, 3.0, "Odds ratio per category increase (log scale)", ha="left",
         va="top", fontsize=6.8, color=GREY, style="italic")

# ================================================================ Panel B
axB.set_axis_on()
axB.text(-0.155, 1.125, "B", transform=axB.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axB.text(-0.055, 1.125, "Conditional $F$-statistics", transform=axB.transAxes,
         fontsize=8.6, fontweight="bold", color=INK, va="top", ha="left")

groups = [
    ("+ BMI\n401 IV", cf["mvmr_condF_chole_bmi.csv"]),
    ("+ cholecystitis\n86 IV", cf["mvmr_condF_chole_cholecyst.csv"]),
    ("+ cholecystitis + cholecystectomy + BMI\n401 IV",
     cf["mvmr_condF_with_bmi.csv"]),
]
palette = {"cholelithiasis": BLUE, "BMI": ORNG, "cholecystitis": RED,
           "cholecystectomy": PURP}
xs, vals, cols = [], [], []
centres, gname_list = [], []
cur = 0.0
for gname, (names, v) in groups:
    start = cur
    for nm, val in zip(names, v):
        xs.append(cur); vals.append(val); cols.append(palette[nm])
        cur += 1.0
    centres.append((start + cur - 1.0) / 2)
    gname_list.append(gname)
    cur += 0.9

axB.bar(xs, vals, width=0.72, color=cols, edgecolor="white", linewidth=0.8, zorder=3)
axB.axhline(10, color=RED, lw=1.0, ls="--", zorder=4)
axB.text(-0.6, 11.5, "identification threshold  $F$ = 10", ha="left",
         va="bottom", fontsize=6.4, color=RED, zorder=5)
for x, v in zip(xs, vals):
    axB.text(x, v * 1.13, f"{v:.2f}", ha="center", va="bottom", fontsize=6.2,
             color=INK, zorder=5)
axB.set_yscale("log")
axB.set_ylim(0.18, 900)
axB.set_xlim(-0.8, max(xs) + 0.8)
axB.set_xticks(centres)
axB.set_xticklabels(gname_list, fontsize=6.2, linespacing=1.45)
axB.set_yticks([0.2, 1, 10, 100])
axB.set_yticklabels(["0.2", "1", "10", "100"])
axB.set_ylabel("Conditional $F$-statistic (log scale)", fontsize=7.4)
axB.tick_params(axis="y", labelsize=6.8, length=2.5, color=LINE)
axB.tick_params(axis="x", length=0)
axB.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axB.spines[s].set_color(LINE)
handles = [plt.Rectangle((0, 0), 1, 1, facecolor=palette[k], edgecolor="white")
           for k in ("cholelithiasis", "BMI", "cholecystitis", "cholecystectomy")]
axB.legend(handles, ["cholelithiasis", "BMI", "cholecystitis", "cholecystectomy"],
           loc="upper center", fontsize=6.1, frameon=True, framealpha=0.95,
           edgecolor=GRID, ncol=2, borderpad=0.45, handlelength=1.0,
           columnspacing=0.9, bbox_to_anchor=(0.62, 1.02))

# ================================================================ Panel C
axC.text(-0.10, 1.30, "C", transform=axC.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axC.text(-0.005, 1.30, "Effect of region definition on colocalisation",
         transform=axC.transAxes, fontsize=8.6, fontweight="bold", color=INK,
         va="top", ha="left")

HILITE = {"R69": ("PNPLA3", GREEN), "R51": ("SERPINA1", GREEN),
          "R67": ("HNF4A", GREEN), "R07": ("ABCG5/ABCG8–CAMKMT (merged)", ORNG)}
xs = np.arange(1, len(coloc) + 1)
ys = [float(r["PP.H4"]) for r in coloc]
for x, y, r in zip(xs, ys, coloc):
    rid = r["region_id"]
    if rid in HILITE:
        _, c = HILITE[rid]
        axC.vlines(x, 0, y, color=c, lw=1.2, zorder=3)
        axC.plot([x], [y], "o", color=c, ms=4.8, zorder=4,
                 markeredgecolor="white", markeredgewidth=0.7)
    else:
        axC.vlines(x, 0, y, color="#D5D8DC", lw=0.8, zorder=2)
        axC.plot([x], [y], "o", color="#AEB6BF", ms=2.2, zorder=3)

axC.axhline(0.8, color=GREY, lw=0.9, ls="--", zorder=2)
axC.text(len(coloc) + 0.4, 0.815, "PP.H4 = 0.8", ha="right", va="bottom",
         fontsize=6.3, color=GREY, zorder=5)
# 三个共定位位点紧挨在最左，用错开的引出线避免标签互相覆盖
for k, rid in enumerate(("R69", "R51", "R67")):
    gene, c = HILITE[rid]
    idx = [i for i, r in enumerate(coloc) if r["region_id"] == rid][0]
    axC.annotate(gene, (xs[idx], ys[idx]), xytext=(16 + k * 34, 46 + k * 20),
                 textcoords="offset points", fontsize=6.5, color=c,
                 fontweight="bold", ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=c, lw=0.6,
                                 connectionstyle="angle,angleA=0,angleB=90,rad=2"))
idx07 = [i for i, r in enumerate(coloc) if r["region_id"] == "R07"][0]
axC.annotate(HILITE["R07"][0], (xs[idx07], ys[idx07]), xytext=(14, 34),
             textcoords="offset points", fontsize=6.5, color=ORNG,
             fontweight="bold", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=ORNG, lw=0.6))
axC.set_xlim(0, len(coloc) + 1)
axC.set_ylim(0, 1.55)
axC.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
axC.set_ylabel("PP.H4", fontsize=7.0)
axC.tick_params(labelsize=6.8, length=2.2, color=LINE)
axC.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axC.spines[s].set_color(LINE)

# R07 合并 vs 拆分
labs = ["R07\nmerged\n2.39 Mb", "R07a\nsplit 1", "R07b\nsplit 2"]
v07 = [R07["merged"], R07["R07a"], R07["R07b"]]
xx = np.arange(3)
axC2.bar(xx, v07, width=0.58, color=[ORNG, BLUE, BLUE], edgecolor="white",
         linewidth=0.8, zorder=3)
for x, v in zip(xx, v07):
    axC2.text(x, v + 0.03, f"{v:.3f}", ha="center", va="bottom", fontsize=6.2,
              color=INK, zorder=4)
axC2.axhline(0.8, color=GREY, lw=0.8, ls="--", zorder=2)
axC2.set_xticks(xx)
axC2.set_xticklabels(labs, fontsize=5.9, linespacing=1.35)
axC2.set_ylim(0, 1.02)
axC2.set_yticks([0, 0.4, 0.8])
axC2.set_ylabel("PP.H4", fontsize=7.0)
axC2.tick_params(labelsize=6.6, length=2.2, color=LINE)
axC2.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axC2.spines[s].set_color(LINE)
axC2.set_title(f"merged region carries {sig['R07'][0]} independent SuSiE signals "
               f"(max conditional |z| = {sig['R07'][1]:.1f})",
               fontsize=6.0, color=ORNG, pad=4.0)

# ------------------------------------------------------------------ 标题
# ⚠ 官方口径（ASN Guide for Authors 原文）："Compile figure titles and legends on one
# or more pages in the manuscript's text file rather than on the figure itself."
# 故此处**不再在图内印主标题与副标题**；图题与图注已写入稿件文本（Figure 1–5 图注）。
# 面板字母 A/B/C/D 与图例键保留在图内（官方允许 "figure keys on figures"）。

out = os.path.join(R, "submission", "Figure4_mvmr_and_failure_modes.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=400, facecolor="white")
fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
print("Figure 4 已写出:", out)
print("  Table 3A 断言通过：OR", [f"{orci(*s[2])[0]:.2f}" for s in sets])
print("  R07 合并/拆分 PP.H4:", {k: round(v, 4) for k, v in R07.items()})
