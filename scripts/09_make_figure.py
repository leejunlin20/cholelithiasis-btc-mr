#!/usr/bin/env python3
"""09_make_figure.py —— 【已作废，勿用于投稿】

⚠️ 作废日期 2026-09-27。本脚本产出的 results/Figure_MR_diet_cholelithiasis_BTC.{png,pdf}
   图内脚注仍是旧口径（81 个 IV / 416 个联合 IV / 条件 F 22.2, 44.5），
   而现稿权威值为 84 / 401 / 23.75, 45.11 —— 即该图与现稿表格自相矛盾。
   旧文件已重命名为 results/_DEPRECATED_Figure_MR_diet_cholelithiasis_BTC.{png,pdf}。

✅ 投稿图件的唯一来源是 results/submission/，由下列脚本生成：
     53_make_figure1.py    Figure1_study_design
     54_make_figure2.py    Figure2_diet_cholelithiasis
     56_make_figure3.py    Figure3_cholelithiasis_btc_chain
     55_make_figure4.py    Figure4_mvmr_and_failure_modes
     61_make_figure5.py    Figure5_cross_ancestry_replication
   交付前体检：57_audit_figures.py（含图件存在性 + 印幅字号 + dpi + 贴边 + 内部路径检查）

如需重新启用本脚本，必须先按现稿权威值更新第 100–140 行的硬编码数字，并重跑 57 号体检。
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FULL = os.path.join(ROOT, "results", "full")
MEDI = os.path.join(ROOT, "results", "medi")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

RED, GREEN, GREY, BLUE = "#C0392B", "#1E8449", "#566573", "#2471A3"


def load(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def fnum(v, d=None):
    try:
        x = float(v)
        return round(x, d) if d is not None else x
    except (TypeError, ValueError):
        return None


def forest(ax, labels, ors, los, his, ps, color, title, xlim, xticks):
    y = np.arange(len(labels))[::-1]
    for i, (o, lo, hi, p) in enumerate(zip(ors, los, his, ps)):
        c = color(o, p)
        ax.plot([lo, hi], [y[i], y[i]], color=c, lw=1.2, solid_capstyle="round", zorder=2)
        ax.plot([o], [y[i]], "o", color=c, ms=5, zorder=3)
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        suffix = f"  {star}" if star else ""
        ax.text(xlim[1] * 1.05, y[i], f"{o:.2f} ({lo:.2f}–{hi:.2f}){suffix}", va="center",
                ha="left", fontsize=7, color="#212F3D")
    ax.axvline(1.0, color="#7F8C8D", lw=0.8, ls="--", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(t) for t in xticks], fontsize=7)
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_title(title, fontsize=9, fontweight="bold", pad=8, loc="left")
    ax.set_xlabel("Odds ratio (95% CI), per 1-SD / per category increase", fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#95A5A6")
    ax.tick_params(length=2.5, color="#95A5A6")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.6),
                               gridspec_kw={"width_ratios": [1.45, 1]})

# ---------- Panel A: 饮食 → 胆囊结石病（MR-PRESSO 校正）----------
mp = {r["id"]: r for r in load(os.path.join(FULL, "mr_presso_results.csv"))}
allr = load(os.path.join(FULL, "mr_results_Cholelithiasis.csv"))
show = [
    ("ukb-b-5779", "Alcohol intake frequency", True),
    ("ukb-b-1489", "Cheese intake", True),
    ("ukb-b-16576", "Dried fruit intake", True),
    ("ukb-b-5237", "Coffee intake", False),
    ("ukb-b-6066", "Tea intake", False),
    ("ukb-b-3881", "Fresh fruit intake", False),
    ("ukb-b-8089", "Cooked vegetable intake", False),
    ("ukb-b-2209", "Oily fish intake", False),
    ("ukb-b-6324", "Processed meat intake", False),
]
labels, ors, los, his, ps, sigs = [], [], [], [], [], []
for eid, name, _ in show:
    r = mp.get(eid)
    if r and fnum(r.get("or_corr")):
        o, lo, hi, p = fnum(r["or_corr"]), fnum(r["or_lci_corr"]), fnum(r["or_uci_corr"]), fnum(r["p_corr"])
    else:
        iv = [x for x in allr if x["id.exposure"] == eid and x["method"] == "Inverse variance weighted"]
        if not iv:
            continue
        o, lo, hi, p = fnum(iv[0]["or"]), fnum(iv[0]["or_lci95"]), fnum(iv[0]["or_uci95"]), fnum(iv[0]["pval"])
    labels.append(name)
    ors.append(o); los.append(lo); his.append(hi); ps.append(p)
    sigs.append(fnum(mp.get(eid, {}).get("fdr_corr")) is not None and fnum(mp[eid]["fdr_corr"]) < 0.05)

forest(ax1, labels, ors, los, his, ps,
       lambda o, p: RED if o > 1 else BLUE,
       "A   Diet → cholelithiasis  (MR-PRESSO outlier-corrected)",
       (0.45, 2.4), [0.5, 0.75, 1.0, 1.5, 2.0])

# ---------- Panel B: 胆囊结石病 → 胆道癌 ----------
s2 = load(os.path.join(MEDI, "step2_chole_to_btc.csv"))
mv = load(os.path.join(MEDI, "mvmr_finngen_chole_bmi_to_btc.csv"))
NAME = {"Inverse variance weighted": "IVW (2-sample)",
        "MR Egger": "MR-Egger",
        "Weighted median": "Weighted median",
        "Weighted mode": "Weighted mode"}
lbl2, o2, l2, h2, p2 = [], [], [], [], []
for r in s2:
    if r["method"] in NAME:
        lbl2.append(NAME[r["method"]])
        o2.append(fnum(r["or"], 3)); l2.append(fnum(r["or_lci95"], 3)); h2.append(fnum(r["or_uci95"], 3))
        p2.append(fnum(r["pval"]))
for r in mv:
    if "Chole" in r["exposure"]:
        lbl2.append("MVMR (adj. BMI)")
        o2.append(fnum(r["or"], 3)); l2.append(fnum(r["or_lci"], 3)); h2.append(fnum(r["or_uci"], 3))
        p2.append(fnum(r["pval"]))

forest(ax2, lbl2, o2, l2, h2, p2,
       lambda o, p: RED,
       "B   Cholelithiasis → biliary tract cancer",
       (0.9, 2.2), [1.0, 1.25, 1.5, 2.0])

fig.text(0.008, 0.985,
         "Genetic evidence for a cholelithiasis → biliary tract cancer causal pathway",
         fontsize=10.5, fontweight="bold", va="top")
fig.text(0.008, 0.935,
         "Two-sample Mendelian randomization  |  379 UK Biobank dietary traits → FinnGen R12 cholelithiasis (49,834 cases) → biliary tract cancer (2,298 cases)  |  FDR < 0.05",
         fontsize=6.8, color="#566573", va="top")
fig.text(0.008, 0.012,
         "*P<0.05  **P<0.01  ***P<0.001.  Red = risk-increasing, blue = protective. Panel A: MR-PRESSO outlier-corrected IVW (or IVW if <5 SNPs). Panel B: 81 independent instruments; MVMR adjusted for BMI using 416 union instruments (conditional F = 22.2 / 44.5).",
         fontsize=6.2, color="#7F8C8D", va="bottom")

fig.subplots_adjust(left=0.128, right=0.868, top=0.845, bottom=0.155, wspace=1.05)
out = os.path.join(ROOT, "results", "Figure_MR_diet_cholelithiasis_BTC.png")
fig.savefig(out, dpi=320, facecolor="white")
fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
print("图已写出:", out)
