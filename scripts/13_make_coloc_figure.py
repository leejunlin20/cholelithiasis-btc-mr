#!/usr/bin/env python3
"""
13_make_coloc_figure.py  (v2)
-----------------------------
COLOC 结果图（Nature 风格，上下双面板 + 独立纵轴）。

A  胆石病 69 个位点中「胆石 ↔ 胆道癌」的 PP.H4 排序（H4>0.8 高亮）
B  阳性对照：BMI ↔ 胆石病的 PP.H4 排序
C–F  区域图：上=胆石病，下=胆道癌，各自独立纵轴；虚线=共定位首选 SNP
     C PNPLA3 / D SERPINA1 / E HNF4A（三个共定位位点）  F ABCG5-ABCG8（最强胆石位点，未共定位）

输出: results/coloc/fig_coloc.png / .pdf / .svg
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

BASE = os.environ.get("MR_BASE", ".")
RES = os.path.join(BASE, "results/coloc")
DAT = os.path.join(BASE, "data/coloc")

plt.rcParams.update({
    # ASN 稿约：图内字体「aim to use Arial, Courier, Times New Roman, Symbol,
    # or fonts that look similar」＋「use uniform lettering」。本稿 8 张投稿图
    # 统一为单一无衬线族 DejaVu Sans（含全部数学/希腊字符，不会跨族回退），
    # 且一律以 TrueType（pdf.fonttype=42）嵌入，不用 Type 3。
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "savefig.bbox": "tight", "savefig.dpi": 600,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

C_HIT, C_MID, C_LOW = "#B2182B", "#4D4D4D", "#C9C9C9"
C_CTRL = "#1B7837"
C_CHOLE, C_BTC = "#2166AC", "#B2182B"

summ = pd.read_csv(os.path.join(RES, "coloc_all_summary.csv"))
main = summ[summ.p12 == 1e-5].copy()
REG = pd.read_csv(os.path.join(RES, "regions.tsv"), sep="\t")

CHOLE = pd.read_csv(os.path.join(DAT, "CHOLE_regions.tsv"), sep="\t")
BTC = pd.read_csv(os.path.join(DAT, "BTC_regions.tsv"), sep="\t")

PANELS = [("R69", "PNPLA3", "chr22"), ("R51", "SERPINA1", "chr14"),
          ("R67", "HNF4A", "chr20"), ("R07", "ABCG5 / ABCG8", "chr2")]


def stat(pair, rid, col):
    r = summ[(summ.pair == pair) & (summ.region_id == rid) & (summ.p12 == 1e-5)]
    return float(r[col].iloc[0]) if len(r) else np.nan


def top_shared(rid):
    r = summ[(summ.pair == "CHOLE_vs_BTC") & (summ.region_id == rid) & (summ.p12 == 1e-5)]
    return str(r.top_shared_snp.iloc[0]) if len(r) else ""


# ================= 画布 =================
fig = plt.figure(figsize=(7.2, 7.4))
gs = GridSpec(2, 4, figure=fig, height_ratios=[0.85, 1.35],
              hspace=0.42, wspace=0.50,
              left=0.085, right=0.985, top=0.945, bottom=0.075)

# ---------- Panel A ----------
axA = fig.add_subplot(gs[0, 0:2])
d = main[main.pair == "CHOLE_vs_BTC"].sort_values("PP.H4", ascending=False).reset_index(drop=True)
x = np.arange(len(d))
v = d["PP.H4"].values
axA.bar(x, v, width=0.85,
        color=[C_HIT if a > 0.8 else (C_MID if a > 0.5 else C_LOW) for a in v],
        edgecolor="none")
axA.axhline(0.8, ls="--", lw=0.7, color="#666666", zorder=3)
axA.set_xlim(-1.5, len(d) + 1)
axA.set_ylim(0, 1.22)
axA.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
axA.set_ylabel("PP.H4 (shared causal variant)")
axA.set_xlabel("Cholelithiasis loci (n = %d), ranked by PP.H4" % len(d))
axA.set_title("A   Cholelithiasis vs biliary tract cancer", loc="left", fontweight="bold")
axA.spines[["top", "right"]].set_visible(False)

annot = ("PNPLA3 (chr22)    PP.H4 = 1.00\n"
         "SERPINA1 (chr14)  PP.H4 = 0.99\n"
         "HNF4A (chr20)     PP.H4 = 0.94")
axA.annotate(annot, xy=(9, 0.16), xytext=(28, 0.60),
             fontsize=6.3, color=C_HIT, ha="left", va="center",
             arrowprops=dict(arrowstyle="-|>", color=C_HIT, lw=0.7,
                             connectionstyle="arc3,rad=-0.25"))
abcg = int(d.index[d.region_id == "R07"][0])
axA.annotate("ABCG5/ABCG8 (chr2)\nstrongest gallstone locus,\nbut PP.H4 = 0.31",
             xy=(abcg, v[abcg] + 0.02), xytext=(17, 0.36),
             fontsize=6.3, color="#333333", ha="left", va="center",
             arrowprops=dict(arrowstyle="-|>", color="#333333", lw=0.7,
                             connectionstyle="arc3,rad=0.25"))
axA.legend(handles=[Patch(facecolor=C_HIT, label="PP.H4 > 0.8"),
                    Patch(facecolor=C_MID, label="0.5 – 0.8"),
                    Patch(facecolor=C_LOW, label="< 0.5")],
           frameon=False, loc="upper center", bbox_to_anchor=(0.60, 1.03),
           ncol=3, handlelength=1.1, columnspacing=1.1)

# ---------- Panel B ----------
axB = fig.add_subplot(gs[0, 2:4])
d2 = main[main.pair == "CHOLE_vs_BMI"].sort_values("PP.H4", ascending=False).reset_index(drop=True)
v2 = d2["PP.H4"].values
axB.bar(np.arange(len(d2)), v2, width=0.85,
        color=[C_CTRL if a > 0.8 else (C_MID if a > 0.5 else C_LOW) for a in v2],
        edgecolor="none")
axB.axhline(0.8, ls="--", lw=0.7, color="#666666", zorder=3)
axB.set_xlim(-1.5, len(d2) + 1)
axB.set_ylim(0, 1.22)
axB.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
axB.set_ylabel("PP.H4 (shared causal variant)")
axB.set_xlabel("Cholelithiasis loci (n = %d), ranked by PP.H4" % len(d2))
axB.set_title("B   Positive control: BMI vs cholelithiasis", loc="left", fontweight="bold")
axB.spines[["top", "right"]].set_visible(False)
axB.annotate("FTO (0.99)\nARNTL (0.98)\nSKIDA1 (0.98)", xy=(8, 0.12), xytext=(24, 0.60),
             fontsize=6.3, color=C_CTRL, ha="left", va="center",
             arrowprops=dict(arrowstyle="-|>", color=C_CTRL, lw=0.7,
                             connectionstyle="arc3,rad=-0.25"))

# ---------- Panels C–F ----------
def nice_max(m):
    for step in (5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 800, 1200):
        if m <= step:
            return step
    return m


for k, (rid, gene, chrom_lab) in enumerate(PANELS):
    sub = gs[1, k].subgridspec(2, 1, hspace=0.10, height_ratios=[1, 1])
    ax1 = fig.add_subplot(sub[0])
    ax2 = fig.add_subplot(sub[1], sharex=ax1)

    c = CHOLE[CHOLE.region_id == rid]
    b = BTC[BTC.region_id == rid]
    lo = REG.loc[REG.region_id == rid, "start"].iloc[0] / 1e6
    hi = REG.loc[REG.region_id == rid, "end"].iloc[0] / 1e6
    sp = int(top_shared(rid).split(":")[1])

    ax1.scatter(c.pos / 1e6, c.mlogp, s=1.5, c=C_CHOLE, alpha=0.8, linewidths=0, rasterized=True)
    ax2.scatter(b.pos / 1e6, b.mlogp, s=1.5, c=C_BTC, alpha=0.8, linewidths=0, rasterized=True)

    m1, m2 = nice_max(c.mlogp.max() * 1.06), nice_max(b.mlogp.max() * 1.30)
    ax1.set_ylim(0, m1)
    ax2.set_ylim(m2, 0)
    ax1.set_xlim(lo, hi)

    for ax in (ax1, ax2):
        ax.axvline(sp / 1e6, ls="--", lw=0.7, color="#333333", zorder=4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", labelbottom=False)
    ax2.tick_params(axis="x", labelbottom=True)

    ax1.set_yticks(np.linspace(0, m1, 3)[1:].round(0))
    ax2.set_yticks(np.linspace(0, m2, 3)[1:].round(0))
    if k == 0:
        ax1.set_ylabel(r"$-\log_{10}P$")
        ax2.set_ylabel(r"$-\log_{10}P$")
    ax2.set_xlabel("Position (Mb)")

    h3, h4 = stat("CHOLE_vs_BTC", rid, "PP.H3"), stat("CHOLE_vs_BTC", rid, "PP.H4")
    col = C_HIT if h4 > 0.8 else "#444444"
    rs = c.loc[(c.pos - sp).abs().idxmin(), "rsids"]
    rs_txt = rs if (isinstance(rs, str) and 0 < len(rs) < 20 and rs != ".") else ("chr%d:%d" % (
        int(top_shared(rid).split(":")[0]), sp))
    ax1.set_title("%s   %s (%s)\n$PP.H3$ = %.2f   $PP.H4$ = %.2f\n%s"
                  % ("CDEF"[k], gene, chrom_lab, h3, h4, rs_txt),
                  loc="left", fontweight="bold", color=col, fontsize=7.2, pad=4)

    if k == 0:
        ax1.text(0.98, 0.93, "Cholelithiasis", transform=ax1.transAxes,
                 ha="right", va="top", fontsize=6.4, color=C_CHOLE, fontweight="bold")
        ax2.text(0.98, 0.07, "Biliary tract cancer", transform=ax2.transAxes,
                 ha="right", va="bottom", fontsize=6.4, color=C_BTC, fontweight="bold")

fig.text(0.085, 0.012,
         "Regions = lead SNPs of cholelithiasis ± 500 kb (merged)  ·  FinnGen R12  ·  "
         "coloc.abf, p12 = 1×10$^{-5}$",
         fontsize=6.2, color="#555555")

# ASN Guide for Authors：「Add the manuscript title and first author's name as a header
# to all supplementary data files.」本图作为 Supplementary Figure S2 投稿，故加页眉。
fig.text(0.085, 1.000,
         "Dietary intake, cholelithiasis and biliary tract cancer: a two-sample "
         "Mendelian randomisation study with cross-ancestry replication and "
         "region-definition-aware colocalisation",
         fontsize=6.2, color="#555555", ha="left", va="top")
fig.text(0.085, 0.982,
         "Jun-Lin Li et al. \u2014 Supplementary material, Supplementary Figure S2",
         fontsize=6.2, color="#555555", ha="left", va="top")

fig.savefig(os.path.join(RES, "fig_coloc.png"))
fig.savefig(os.path.join(RES, "fig_coloc.pdf"))
fig.savefig(os.path.join(RES, "fig_coloc.svg"))
SUPP = os.path.join(BASE, "results", "submission", "supplementary")
os.makedirs(SUPP, exist_ok=True)
fig.savefig(os.path.join(SUPP, "Supplementary_Figure_S2.png"), dpi=400)
fig.savefig(os.path.join(SUPP, "Supplementary_Figure_S2.pdf"))
print("已输出: fig_coloc.png / .pdf / .svg  ｜ 补充图 S2（带页眉）")
