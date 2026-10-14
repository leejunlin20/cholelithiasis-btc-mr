#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
60_make_figureS_qtl.py —— 候选补充图：四共定位区的 QTL 层证据

面板
  A  eQTL × GWAS 共定位 PP.H4（肝/血 × 胆石症，三档先验）
  B  PP.H3 vs PP.H4（主口径 p12=1e-5）—— 区分"可判阴性"与"不可判定"
  C  cis-MR 工具强度 vs 效应量（标出不可解释者）
  D  R07 先验敏感性：现稿区域级 vs 本轮 eQTL 层

口径：图内不印主标题（ASN：图题图注写在稿件文本里）；不出现仓库路径；
      字体 ≈8 pt；400 dpi PNG + 矢量 PDF；每面板内置数字闭合断言。
"""
import os, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJ = os.environ.get("MR_BASE", ".")
QT = os.path.join(PROJ, "results", "qtl")
OUT = os.path.join(QT, "FigureS_qtl_layer_coloc")

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.facecolor": "white", "figure.facecolor": "white",
    "savefig.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

C_H3 = "#2166ac"     # 可判阴性（H3）
C_H2 = "#b2182b"     # 不可判定（H2）
C_GREY = "#666666"
C_LIVER = "#1b7837"
C_BLOOD = "#762a83"
P12C = {1e-6: "#a6cee3", 1e-5: "#1f78b4", 1e-4: "#08306b"}

# ============================ 数据 ============================
co = pd.read_csv(os.path.join(QT, "coloc_eqtl_gwas.csv")).rename(
    columns={"PP.H0": "PPH0", "PP.H1": "PPH1", "PP.H2": "PPH2",
             "PP.H3": "PPH3", "PP.H4": "PPH4"})
mr = pd.read_csv(os.path.join(QT, "qtl_mr_by_gene.csv"))
lead = pd.read_csv(os.path.join(QT, "lead_snp_eqtl_check.csv"))

chole = co[co.outcome == "CHOLE"].copy()
main = chole[chole.p12 == 1e-5].copy()
main["label"] = main.apply(
    lambda r: f"{r['region']}  {r['gene']}  ({'liver' if r['tissue']=='Liver' else 'blood'})", axis=1)

# ---- 闭合断言 ----
assert abs(co["PPH4"].max() - 0.435519192719236) < 1e-9, co["PPH4"].max()
assert (co["PPH4"] > 0.5).sum() == 0
assert len(co) == 48, len(co)
_m = mr[(mr.region == "R51") & (mr.tissue == "Whole_Blood")]
assert abs(float(_m.OR_CHOLE.iloc[0]) - 1.051104) < 1e-4
_m2 = mr[(mr.region == "R07") & (mr.gene == "ABCG5") & (mr.tissue == "Liver")]
assert abs(float(_m2.OR_CHOLE.iloc[0]) - 3.975143) < 1e-4
_leads = lead[(lead.region == "R69") & (lead.lead_snp == "rs738408") & (lead.tissue == "Liver")]
assert float(_leads.eqtl_p.iloc[0]) == 0.546615, float(_leads.eqtl_p.iloc[0])

fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.6))
fig.subplots_adjust(hspace=0.45, wspace=0.32)

# ==================== A: PP.H4 三档先验 ====================
ax = axes[0, 0]
order = main.sort_values("PPH4").reset_index(drop=True)
ypos = np.arange(len(order))
for p12 in (1e-6, 1e-5, 1e-4):
    v = [float(chole[(chole.region == r.region) & (chole.gene == r.gene) &
                     (chole.tissue == r.tissue) & (chole.p12 == p12)]["PPH4"].iloc[0])
         for r in order.itertuples()]
    ax.barh(ypos + {1e-6: -0.26, 1e-5: 0.0, 1e-4: 0.26}[p12], v, height=0.25,
            color=P12C[p12], label=f"p12 = {p12:.0e}", edgecolor="none")
ax.axvline(0.5, color=C_GREY, lw=0.8, ls="--")
ax.axvline(0.8, color=C_GREY, lw=0.8, ls=":")
ax.text(0.5, 1.0, "0.5", transform=ax.get_xaxis_transform(), fontsize=6,
        color=C_GREY, ha="center", va="bottom")
ax.text(0.8, 1.0, "0.8", transform=ax.get_xaxis_transform(), fontsize=6,
        color=C_GREY, ha="center", va="bottom")
ax.set_yticks(ypos); ax.set_yticklabels(order["label"], fontsize=7)
ax.set_xlim(0, 1.0); ax.set_xlabel("PP.H4 (eQTL shared with cholelithiasis)")
ax.legend(frameon=False, loc="lower right", handlelength=1.1, borderpad=0.2)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_title("A   Expression-level colocalisation (liver / blood)", loc="left", fontsize=8)

# ==================== B: H3 vs H4 ====================
ax = axes[0, 1]
sub = main.copy()
col = [C_LIVER if t == "Liver" else C_BLOOD for t in sub.tissue]
ax.axvspan(0.5, 1.02, color=C_H3, alpha=0.06, lw=0)
ax.scatter(sub["PPH3"], sub["PPH4"], s=np.sqrt(sub["n_snp"]) * 1.4,
           c=col, alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
ax.axhline(0.5, color=C_GREY, lw=0.8, ls="--")
ax.axvline(0.5, color=C_GREY, lw=0.8, ls="--")
_lab = [r for r in sub.itertuples() if r.PPH3 > 0.5 or r.PPH4 > 0.09]
_hi = [r for r in _lab if r.PPH3 > 0.9]
for r in _lab:
    if r.PPH3 <= 0.9:
        ax.annotate(f"{r.region}/{r.gene}", (r.PPH3, r.PPH4), fontsize=6,
                    xytext=(-2, 8), textcoords="offset points", ha="right")
if _hi:
    ax.annotate("R07/LRPPRC ·\nR51/SERPINA1\n(whole blood)",
                (0.985, 0.02), fontsize=6, ha="right", va="bottom")

ax.set_xlim(0, 1.05); ax.set_ylim(-0.02, 0.56)
ax.set_xlabel("PP.H3 (two distinct signals)")
ax.set_ylabel("PP.H4 (shared causal variant)")
ax.text(0.98, 0.30, "H3 > 0.5\n\"not mediated\"\n(decidable)",
        fontsize=6, color=C_H3, ha="right", va="top")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C_LIVER, label="liver"),
                   Line2D([], [], marker="o", ls="", color=C_BLOOD, label="whole blood")],
          frameon=False, loc="upper left", handletextpad=0.2, borderpad=0.2)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_title("B   Distinguishing 'no mediation' from 'not estimable'", loc="left", fontsize=8)

# ==================== C: cis-MR ====================
ax = axes[1, 0]
mrp = mr[mr.OR_CHOLE.notna() & (mr.n_shared > 0)].copy()
bad = mrp.status.str.startswith("no_cis") | (mrp.top_dist_kb > 100)
colr = np.where(bad, "#d6604d", "#4393c3")
ax.errorbar(mrp.F, mrp.OR_CHOLE,
            yerr=[mrp.OR_CHOLE - mrp.ORl_CHOLE, mrp.ORu_CHOLE - mrp.OR_CHOLE],
            fmt="none", ecolor="#bbbbbb", elinewidth=0.8, capsize=1.5, zorder=2)
ax.scatter(mrp.F, mrp.OR_CHOLE, s=26, c=colr, edgecolor="white", linewidth=0.5, zorder=3)
ax.axhline(1.0, color=C_GREY, lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([10, 20, 40, 80]); ax.set_xticklabels(["10", "20", "40", "80"])
ax.set_yticks([0.7, 1.0, 1.5, 2.0, 3.0, 4.0])
ax.set_yticklabels(["0.7", "1.0", "1.5", "2.0", "3.0", "4.0"])
ax.minorticks_off()
ax.set_xlabel("Instrument strength  F")
ax.set_ylabel("Wald-ratio OR for cholelithiasis")
for r in mrp.itertuples():
    if r.OR_CHOLE > 2:
        ax.annotate(f"{r.gene} (liver)", (r.F, r.OR_CHOLE), fontsize=6,
                    xytext=(-5, 6), textcoords="offset points", ha="right")
    elif r.OR_CHOLE < 0.92:
        ax.annotate(f"{r.gene} (blood)", (r.F, r.OR_CHOLE), fontsize=6,
                    xytext=(-5, -13), textcoords="offset points", ha="right")
ax.text(0.98, 0.97, "OR = 3.98 implausible\n(region harbours 5\nindependent signals)",
        transform=ax.transAxes, fontsize=6, color="#d6604d", va="top", ha="right")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#4393c3",
                          label="proximal & cis-significant"),
                   Line2D([], [], marker="o", ls="", color="#d6604d",
                          label="distal or weak instrument")],
          frameon=False, loc="lower right", handletextpad=0.2, borderpad=0.2)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_title("C   Single-instrument cis-MR and its failure modes", loc="left", fontsize=8)

# ==================== D: 先验敏感性对比 ====================
ax = axes[1, 1]
xs = np.array([1e-6, 1e-5, 1e-4])
legacy = np.array([0.0436930162883554, 0.313607890768151, 0.705633976561156])
abcg5 = np.array([float(chole[(chole.region == "R07") & (chole.gene == "ABCG5") &
                              (chole.tissue == "Liver") & (chole.p12 == p)]["PPH4"].iloc[0])
                  for p in xs])
abcg8 = np.array([float(chole[(chole.region == "R07") & (chole.gene == "ABCG8") &
                              (chole.tissue == "Liver") & (chole.p12 == p)]["PPH4"].iloc[0])
                  for p in xs])
ax.plot(xs, legacy, "-o", color="#b2182b", ms=4, lw=1.2,
        label="manuscript: regional GWAS–GWAS (R07)")
ax.plot(xs, abcg8, "-s", color="#1b7837", ms=4, lw=1.2,
        label="this work: ABCG8 eQTL (liver)")
ax.plot(xs, abcg5, "-^", color="#2166ac", ms=4, lw=1.2,
        label="this work: ABCG5 eQTL (liver)")
ax.axhline(0.5, color=C_GREY, lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_ylim(-0.02, 1.02)
ax.set_xlabel("coloc prior  p12  (log scale)")
ax.set_ylabel("PP.H4")
ax.text(1.05e-6, 0.62, "crosses 0.5\n(prior-dependent)", fontsize=6, color="#b2182b")
ax.text(1.05e-6, 0.40, "never crosses 0.5\n(prior-robust)", fontsize=6,
        color="#2166ac", va="center")
ax.legend(frameon=False, loc="upper left", handlelength=1.4, borderpad=0.2,
          labelspacing=0.25)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_title("D   Prior sensitivity: regional vs expression level", loc="left", fontsize=8)

for ext, kw in ((".png", {"dpi": 400}), (".pdf", {})):
    fig.savefig(OUT + ext, **kw)

# 只允许写一次 PNG：落盘后断言像素宽度
import struct
_b = open(OUT + ".png", "rb").read(24)
w = struct.unpack(">I", _b[16:20])[0]
_h = struct.unpack(">I", _b[20:24])[0]
assert w >= 2500, f"PNG 宽度异常 {w}（疑似被低 dpi 覆盖）"
print(f"OK  {OUT}.png  {w}x{_h}")
print(f"OK  {OUT}.pdf")
print(f"闭合断言：max PP.H4 = {co['PPH4'].max():.6f} | H4>0.5 = {(co['PPH4']>0.5).sum()} | 行数 = {len(co)}")
