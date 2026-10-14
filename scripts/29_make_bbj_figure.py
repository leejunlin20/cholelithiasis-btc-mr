#!/usr/bin/env python3
"""
29_make_bbj_figure.py —— BBJ 跨队列复制汇总图
Panel A: 6 个结局的 IVW OR 森林图（阳性/阴性/主结局分色）
Panel B: FinnGen 胆石症 β vs BBJ 胆道癌 β（61 位点，加权回归）
Panel C: FinnGen 胆道癌 β vs BBJ 胆道癌 β（61 位点，加权回归）

配色遵循中文习惯：风险升高（OR>1 / β 同向为正）用红，降低用绿。
"""
import csv, math, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.environ.get("MR_BASE", ".")
RES = os.path.join(ROOT, "results", "bbj_replication")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})
RED, GREEN, GREY, BLUE, DARK = "#C0392B", "#1E8449", "#7F8C8D", "#2471A3", "#212F3D"


def load(p, delim=None):
    if delim is None:
        delim = "\t" if p.endswith(".tsv") else ","
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


# ---------------- 读数据 ----------------
mr = load(os.path.join(RES, "bbj_mr_results.csv"))
ivw = {}
for r in mr:
    if r["method"] == "Inverse variance weighted":
        ivw[r["pheno"]] = r

ORDER = [
    ("BtC",           "Biliary tract cancer (main)", 418,    DARK),
    ("Cholelithiasis", "Cholelithiasis (pos. ctrl)", 9305,   RED),
    ("Cholecystitis",  "Cholecystitis (pos. ctrl)",  1860,   RED),
    ("GB_polyp",       "Gallbladder polyp",          843,    RED),
    ("Urolithiasis",   "Urolithiasis (neg. ctrl)",   11699,  GREEN),
    ("PaC",            "Pancreatic cancer (neg. ctrl)", 499, GREEN),
]

loci = load(os.path.join(RES, "locus_replication_v2.tsv"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def pairs(xk, yk, yse):
    out = []
    for r in loci:
        x, y, s = num(r[xk]), num(r[yk]), num(r[yse])
        if x is None or y is None or s is None or s <= 0:
            continue
        out.append((x, y, s))
    return out


def wls(x, y, w):
    S = sum(w)
    mx = sum(wi * xi for wi, xi in zip(w, x)) / S
    my = sum(wi * yi for wi, yi in zip(w, y)) / S
    Sxx = sum(wi * (xi - mx) ** 2 for wi, xi in zip(w, x))
    Sxy = sum(wi * (xi - mx) * (yi - my) for wi, xi, yi in zip(w, x, y))
    b = Sxy / Sxx
    resid = [yi - my - b * (xi - mx) for xi, yi in zip(x, y)]
    s2 = sum(wi * ri ** 2 for wi, ri in zip(w, resid)) / (len(x) - 2)
    return b, math.sqrt(s2 / Sxx), mx, my


fig = plt.figure(figsize=(7.2, 2.55))
gs = fig.add_gridspec(1, 3, width_ratios=[1.30, 1.0, 1.0], wspace=0.42)

# ================= Panel A: forest =================
axA = fig.add_subplot(gs[0, 0])
ys = np.arange(len(ORDER))[::-1]
for y, (ph, lab, n, col) in zip(ys, ORDER):
    r = ivw[ph]
    orv, lo, hi = math.exp(float(r["b"])), math.exp(float(r["lo_ci"])), math.exp(float(r["hi_ci"]))
    p = float(r["pval"])
    axA.plot([lo, hi], [y, y], color=col, lw=1.5, solid_capstyle="butt", zorder=2)
    axA.plot([orv], [y], "o", ms=5.0, color=col,
             markerfacecolor=col if p < 0.05 else "white",
             markeredgecolor=col, markeredgewidth=1.1, zorder=3)
    # P 值放在绘图区最右端（右对齐），与置信区间不重叠
    axA.text(0.995, y, "P = %.2g" % p, transform=axA.get_yaxis_transform(),
             va="center", ha="right", fontsize=6.3, color=DARK if col == DARK else GREY)
axA.axvline(1.0, color=GREY, lw=0.8, ls=(0, (4, 2)), zorder=1)
axA.set_yticks(ys)
axA.set_yticklabels(["%s\n(%s cases)" % (lab, format(n, ",")) for _, lab, n, _ in ORDER],
                    fontsize=6.6)
axA.set_xscale("log")
axA.set_xlim(0.55, 4.2)
axA.set_xticks([0.6, 1, 1.5, 2.0, 2.5])
axA.set_xticklabels(["0.6", "1.0", "1.5", "2.0", "2.5"])
axA.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axA.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(
    ["0.6", "1.0", "1.5", "2.0", "2.5"]))
axA.set_ylim(-0.7, len(ORDER) - 0.25)
axA.set_xlabel("IVW odds ratio per log-odds cholelithiasis\n(FinnGen R12 $\\rightarrow$ BioBank Japan)",
               fontsize=6.8)
axA.tick_params(axis="y", length=0)
axA.tick_params(axis="x", labelsize=6.4)
axA.text(-0.02, 1.16, "A", transform=axA.transAxes, fontsize=10,
         fontweight="bold", va="top", ha="left")
for sp in ("right", "top"):
    axA.spines[sp].set_visible(False)

# ================= Panel B & C: scatter =================
for k, (ax_idx, xk, xlab, tag) in enumerate([
        (1, "fg_chole_beta", "FinnGen cholelithiasis $\\beta$", "B"),
        (2, "fg_btc_beta",   "FinnGen biliary tract cancer $\\beta$", "C")]):
    ax = fig.add_subplot(gs[0, ax_idx])
    d = pairs(xk, "bbj_btc_beta", "bbj_btc_se")
    x = [a for a, _, _ in d]; y = [b for _, b, _ in d]
    w = [1.0 / s ** 2 for _, _, s in d]
    b, se, mx, my = wls(x, y, w)
    z = b / se
    xs = np.array([min(x), max(x)])
    ax.plot(xs, my + b * (xs - mx), color=BLUE, lw=1.1, zorder=3)
    ax.scatter(x, y, s=8, color=GREY, alpha=0.75, linewidths=0, zorder=2)
    ax.axhline(0, color=GREY, lw=0.6, ls=(0, (3, 2)))
    ax.axvline(0, color=GREY, lw=0.6, ls=(0, (3, 2)))
    ax.set_xlabel(xlab, fontsize=6.8)
    ax.set_ylabel("BBJ biliary tract cancer $\\beta$", fontsize=6.8)
    ax.set_xlim(min(x) - 0.03, max(x) + 0.03)
    ax.set_ylim(min(y) - 0.35, max(y) + 0.35)
    ax.text(0.05, 0.95, "slope $=%.3f$  ($z=%.2f$)\n$n=%d$ lead SNPs" % (b, z, len(x)),
            transform=ax.transAxes, fontsize=6.4, va="top", ha="left",
            bbox=dict(fc="white", ec="#D5D8DC", lw=0.5, pad=2.0, alpha=0.95))
    ax.text(-0.16, 1.15, tag, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")
    ax.tick_params(labelsize=6.4)

OUT = os.path.join(RES, "fig_bbj_cross_cohort")
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf")
print("已输出:", OUT + ".png / .pdf")

# ---------------- 交付前自检：像素边界是否被裁 ----------------
from PIL import Image
im = Image.open(OUT + ".png").convert("L")
a = np.array(im)
h, w = a.shape
for name, strip in [("顶部", a[:12, :]), ("底部", a[-12:, :]),
                    ("左侧", a[:, :12]), ("右侧", a[:, -12:])]:
    nz = int((strip < 250).sum())
    print("  %s 12px 内非白像素 = %d %s" % (name, nz, "OK" if nz <= 50 else "⚠️ 可能被裁"))
print("  画布 %dx%d px" % (w, h))
