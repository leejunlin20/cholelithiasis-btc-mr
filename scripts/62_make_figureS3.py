#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""62_make_figureS3.py —— Supplementary Figure S3：敏感性分析三件套

面板
  A  Steiger 方向性检验（含患病率敏感性）：胆石症 → 胆道癌
  B  MR-PRESSO 多效性/离群诊断：38 个可分析饮食暴露 → 胆石症
  C  留一法（leave-one-out）：胆石症 → 胆道癌，84 个工具变量

口径（与 53–56 / 61 一致）：
  - 图内不印主标题（ASN：图题图注写在稿件文本里，不写在图上）
  - 图内不出现任何仓库路径
  - 单栏/全宽皆可，字体基线 8 pt
  - 每个面板内置数字闭合断言，不闭合即拒绝出图

留一法口径：与主分析同一套随机效应约定（单 κ，下限钳 1）——
  w_k = β_exp,k² / se_out,k²；r_k = β_out,k / β_exp,k
  b = Σw r / Σw；se_FE = 1/√Σw；Q = Σw(r−b)²；κ = max(1, √(Q/df))；se_RE = se_FE × κ
  **绝不用双 κ**（TwoSampleMR::mr_ivw() 默认已返回随机效应）。
"""
import os
import math
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

PROJ = os.environ.get("MR_BASE", ".")
SUB = os.path.join(PROJ, "results", "submission")
OUT = os.path.join(SUB, "FigureS3_sensitivity")

plt.rcParams.update({
    # 与全稿其它投稿图统一为单一无衬线族 + TrueType 嵌入（见 13_make_coloc_figure.py 注）
    "font.family": "DejaVu Sans",
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.facecolor": "white", "figure.facecolor": "white",
    "savefig.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

C_OK = "#1b7837"      # 通过 / 正确方向
C_BAD = "#b2182b"     # 失败 / 反向
C_GREY = "#666666"
C_BLUE = "#2166ac"
C_ORANGE = "#d95f02"

# ==================================================================
# 1. 留一法：胆石症 → 胆道癌（84 IV），自造、单 κ 随机效应
# ==================================================================
IV_TSV = os.path.join(PROJ, "results", "mvmr", "iv_cholelith_clumped.tsv")
OUT_TSV = "/tmp/outcome_btc_84iv.tsv"
if not os.path.exists(OUT_TSV):  # 沙箱重启后 /tmp 会被清空 → 从仓库内的权威副本取
    OUT_TSV = os.path.join(PROJ, "results", "mvmr", "outcome_btc_84iv.tsv")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _read(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


ivs = _read(IV_TSV)
out = _read(OUT_TSV)
idx = {}
for r in out:
    for k in ("rsids", "rsid"):
        if r.get(k):
            idx[r[k]] = r
            break

recs = []
for v in ivs:
    o = idx.get(v["rsid"])
    if o is None:
        continue
    be, ae = v["ref"], v["alt"]
    bo, ao = o["ref"], o["alt"]
    if ae == ao and be == bo:
        sgn = 1.0
    elif ae == bo and be == ao:
        sgn = -1.0
    elif COMP.get(ae) == ao and COMP.get(be) == bo:
        sgn = 1.0
    elif COMP.get(ae) == bo and COMP.get(be) == ao:
        sgn = -1.0
    else:
        continue
    bx, sx = float(v["beta"]), float(v["se"])
    by, sy = sgn * float(o["beta"]), float(o["sebeta"])
    if sx <= 0 or sy <= 0:
        continue
    recs.append({"rsid": v["rsid"], "bx": bx, "sx": sx, "by": by, "sy": sy})


def ivw_re(sub):
    """单 κ 随机效应 IVW；返回 (b, se, kappa, Q, df)。"""
    w = [r["bx"] ** 2 / r["sy"] ** 2 for r in sub]
    ratio = [r["by"] / r["bx"] for r in sub]
    b = sum(wi * ri for wi, ri in zip(w, ratio)) / sum(w)
    se_fe = 1.0 / math.sqrt(sum(w))
    Q = sum(wi * (ri - b) ** 2 for wi, ri in zip(w, ratio))
    df = len(sub) - 1
    kappa = max(1.0, math.sqrt(Q / df)) if df > 0 else 1.0
    return b, se_fe * kappa, kappa, Q, df


b_full, se_full, kap_full, Q_full, df_full = ivw_re(recs)

# ---- 闭合断言：必须复现项目对数锚点 ----
assert len(recs) == 84, f"配对 IV 数 {len(recs)} ≠ 84"
assert abs(b_full - 0.2023964295) < 1e-9, b_full
assert abs(se_full - 0.0732635699) < 1e-9, se_full
assert abs(Q_full - 232.394325) < 1e-3, Q_full

loo = []
for k in range(len(recs)):
    sub = recs[:k] + recs[k + 1:]
    bk, sek, kappak, _, _ = ivw_re(sub)
    loo.append({"rsid": recs[k]["rsid"], "b": bk, "se": sek,
                "or": math.exp(bk), "lci": math.exp(bk - 1.959964 * sek),
                "uci": math.exp(bk + 1.959964 * sek)})
loo = pd.DataFrame(loo).sort_values("or").reset_index(drop=True)

OR_FULL = math.exp(b_full)
LCI_FULL = math.exp(b_full - 1.959964 * se_full)
UCI_FULL = math.exp(b_full + 1.959964 * se_full)
assert abs(OR_FULL - 1.224333) < 1e-6, OR_FULL
assert loo["lci"].min() > 1.0 or loo["or"].min() > 1.0, "留一法出现方向翻转，须核查"

loo.to_csv(os.path.join(SUB, "FigureS3_leave_one_out.csv"), index=False)

# ==================================================================
# 2. Steiger 与 MR-PRESSO 数据
# ==================================================================
st = pd.read_csv(os.path.join(PROJ, "results", "steiger", "steiger_lor_results.csv"))
ch = st[st.pair == "Cholelithiasis → BTC"]
ch_corr = ch[ch.method.str.contains("lor")]
ch_pn = ch[ch.method.str.contains("pn")]
# 闭合断言：修正法比值区间与反向的原方法
assert abs(ch_corr.ratio.min() - 2.203802) < 1e-5, ch_corr.ratio.min()
assert abs(ch_corr.ratio.max() - 2.261307) < 1e-5, ch_corr.ratio.max()
assert bool(ch_corr.correct_direction.all())
assert len(ch_pn) == 1 and not bool(ch_pn.correct_direction.iloc[0])
assert abs(float(ch_pn.ratio.iloc[0]) - 0.946285) < 1e-5

pr = pd.read_csv(os.path.join(PROJ, "results", "full", "mr_presso_results.csv"))
FDR_HITS = ["Alcohol intake frequency.", "Cheese intake", "Dried fruit intake"]
assert set(FDR_HITS).issubset(set(pr.exposure)), "三个 FDR<0.05 暴露未全部出现"
assert (pr[pr.exposure.isin(FDR_HITS)].fdr_corr < 0.05).all()

# ==================================================================
# 3. 出图
# ==================================================================
fig = plt.figure(figsize=(11.4, 7.2))
gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.25],
                       hspace=0.30, wspace=0.24,
                       left=0.06, right=0.985, top=0.930, bottom=0.065)

# ---------------- A: Steiger 患病率敏感性 ----------------
axA = fig.add_subplot(gs[0, 0])
prev_out_grid = [0.0005, 0.001, 0.002, 0.005]
xpos = np.arange(len(prev_out_grid))
lo, hi, mid = [], [], []
for po in prev_out_grid:
    s = ch_corr[np.isclose(ch_corr.prev_out, po)]
    assert len(s) == 4, f"prev_out={po} 应有 4 档暴露患病率，实得 {len(s)}"
    lo.append(s.ratio.min()); hi.append(s.ratio.max()); mid.append(s.ratio.mean())
axA.vlines(xpos, lo, hi, color=C_BLUE, lw=2.4, alpha=0.85)
axA.plot(xpos, mid, "o", color=C_BLUE, ms=5, zorder=3,
         label="Variance explained ratio,\ncorrected $r^2$ ($L_oR$)")
axA.axhline(1.0, color=C_GREY, lw=0.9, ls="--")
axA.text(len(xpos) - 0.45, 1.03, "ratio = 1", fontsize=6.5, color=C_GREY, ha="right")
axA.plot([len(xpos) - 0.5], [float(ch_pn.ratio.iloc[0])], "s", color=C_BAD, ms=5.5,
         zorder=3, label="Uncorrected $r^2$ ($P_N$)\n- direction reversed")
axA.axhline(float(ch_pn.ratio.iloc[0]), color=C_BAD, lw=0.8, ls=":")
axA.set_xticks(xpos)
axA.set_xticklabels([f"{p*100:g}%" for p in prev_out_grid])
axA.set_xlabel("Assumed outcome (BTC) prevalence")
axA.set_ylabel("$r^2$ exposure / $r^2$ outcome")
axA.set_ylim(0.75, 2.45)
axA.set_xlim(-0.6, len(xpos) - 0.15)
axA.legend(loc="upper right", frameon=False, handletextpad=0.5, borderaxespad=0.1)
axA.set_title("Assumptions yield ratio > 1 at every prevalence",
              fontsize=7.5, loc="left", pad=4)
for sp in ("top", "right"):
    axA.spines[sp].set_visible(False)

# ---------------- B: MR-PRESSO 离群诊断 ----------------
axB = fig.add_subplot(gs[0, 1])
pr2 = pr.copy()
pr2["is_hit"] = pr2.exposure.isin(FDR_HITS)
axB.scatter(pr2.loc[~pr2.is_hit, "nsnp_outlier"],
            pr2.loc[~pr2.is_hit, "global_p"].clip(lower=1e-4),
            s=16, facecolor="white", edgecolor=C_GREY, lw=0.7, label="Other exposures")
axB.scatter(pr2.loc[pr2.is_hit, "nsnp_outlier"],
            pr2.loc[pr2.is_hit, "global_p"].clip(lower=1e-4),
            s=26, facecolor=C_ORANGE, edgecolor="none",
            label="FDR < 0.05 after outlier correction")
axB.axhline(0.05, color=C_GREY, lw=0.9, ls="--")
axB.set_yscale("log")
axB.set_xlabel("Outlying variants detected (n)")
axB.set_ylabel("MR-PRESSO global test $P$")
axB.set_ylim(8e-4, 2)
axB.set_xlim(-1.2, pr2.nsnp_outlier.max() + 1.2)
axB.text(pr2.nsnp_outlier.max(), 0.058, "$P$ = 0.05", fontsize=6.5,
         color=C_GREY, ha="right")
# 图例置于右上（面板 B 右上为空区）；原先置左下时，DejaVu Sans 比 Arial 宽，
# 第二行文字会压到 x≈4 处的数据点
axB.legend(loc="upper right", frameon=False, handletextpad=0.5, borderaxespad=0.1)
axB.set_title(f"Outlier detection, {len(pr2)} dietary exposures",
              fontsize=7.5, loc="left", pad=4)
for sp in ("top", "right"):
    axB.spines[sp].set_visible(False)

# ---------------- C: 留一法 ----------------
axC = fig.add_subplot(gs[1, :])
n = len(loo)
y = np.arange(n)
axC.vlines(y, loo.lci, loo.uci, color=C_GREY, lw=0.6, alpha=0.9)
axC.plot(y, loo["or"], "o", color=C_BLUE, ms=2.6, zorder=3,
         label=f"Leave-one-out IVW ({n} variants)")
axC.axhline(OR_FULL, color=C_BAD, lw=1.1, ls="--", zorder=2,
            label=f"All variants: OR {OR_FULL:.3f} ({LCI_FULL:.3f}-{UCI_FULL:.3f})")
axC.axhline(1.0, color="#999999", lw=0.8, ls=":", zorder=1)
axC.set_xlim(-1, n)
axC.set_ylim(min(1.0, loo.lci.min() - 0.03), loo.uci.max() + 0.03)
axC.set_xlabel(f"Variant omitted (one of 84 per refit, ordered by leave-one-out "
               f"estimate; lowest {loo.rsid.iloc[0]} to highest {loo.rsid.iloc[-1]})")
axC.set_ylabel("OR for biliary tract cancer\nper unit log-odds of cholelithiasis")
axC.legend(loc="upper left", frameon=False, handletextpad=0.5, borderaxespad=0.2, ncol=2)
for sp in ("top", "right"):
    axC.spines[sp].set_visible(False)
axC.text(0.995, 0.06,
         f"Direction preserved in {n}/{n} refits; "
         f"estimates span OR {loo["or"].min():.3f}-{loo["or"].max():.3f}",
         transform=axC.transAxes, ha="right", va="bottom", fontsize=6.8, color=C_GREY)

for lab, ax in (("A", axA), ("B", axB), ("C", axC)):
    ax.text(-0.085 if lab != "C" else -0.045, 1.045, lab,
            transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")

# ASN Guide for Authors：「Add the manuscript title and first author's name as a header
# to all supplementary data files.」本图作为 Supplementary Figure S3 投稿，故加页眉。
fig.text(0.06, 1.000,
         "Dietary intake, cholelithiasis and biliary tract cancer: a two-sample "
         "Mendelian randomisation study with cross-ancestry replication and "
         "region-definition-aware colocalisation",
         fontsize=7.6, color="#555555", ha="left", va="top")
fig.text(0.06, 0.973,
         "Jun-Lin Li et al. \u2014 Supplementary material, Supplementary Figure S3",
         fontsize=7.6, color="#555555", ha="left", va="top")

fig.savefig(OUT + ".png", dpi=400)
fig.savefig(OUT + ".pdf")
SUPP = os.path.join(SUB, "supplementary")
os.makedirs(SUPP, exist_ok=True)
fig.savefig(os.path.join(SUPP, "Supplementary_Figure_S3.png"), dpi=400)
fig.savefig(os.path.join(SUPP, "Supplementary_Figure_S3.pdf"))
plt.close(fig)

print(f"配对 IV: {len(recs)}  →  b = {b_full:.10f}  se = {se_full:.10f}  "
      f"kappa = {kap_full:.4f}  Q = {Q_full:.3f} on {df_full} df")
print(f"全样本 OR {OR_FULL:.4f} ({LCI_FULL:.4f}-{UCI_FULL:.4f})")
print(f"留一法: OR 区间 {loo["or"].min():.4f}–{loo["or"].max():.4f}；"
      f"全部 CI 下界 > 1: {bool((loo.lci > 1).all())}")
print(f"Steiger: 修正法比值 {ch_corr.ratio.min():.3f}–{ch_corr.ratio.max():.3f}（方向全部正确）；"
      f"原方法 {float(ch_pn.ratio.iloc[0]):.3f}（方向错误）")
print(f"MR-PRESSO: {len(pr)} 个暴露，离群数 {pr.nsnp_outlier.min()}–{pr.nsnp_outlier.max()}")
print(f"输出: {OUT}.png / .pdf")
