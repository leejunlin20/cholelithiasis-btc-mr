#!/usr/bin/env python3
"""
16c_make_susie_figure.py
------------------------
coloc.susie 补充图（Nature 风格，与 fig_coloc 同配色/字号）。

A  PNPLA3 区域图（唯一能跑通 coloc.susie 的区域）
   上=胆石病，下=胆道癌，各自独立纵轴；实心点=susie credible set 成员；
   虚线=两侧共享 lead SNP (rs738408)
B  PNPLA3 区域两侧的 PIP 曲线（细定位后验概率），阴影=credible set
C  5 个区域：coloc.abf 与 coloc.susie 的 PP.H4 对照
   灰色=abf，深红=susie；斜纹且标 n/a=胆道癌侧建不出 credible set，susie 不可用

输出: results/coloc_susie/fig_coloc_susie.png / .pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

BASE = os.environ.get("MR_BASE", ".")
SUS = os.path.join(BASE, "results/coloc_susie")
COL = os.path.join(BASE, "results/coloc")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
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

C_CHOLE, C_BTC = "#2166AC", "#B2182B"
C_ABF = "#BDBDBD"
REGIONS = [("R69", "PNPLA3"), ("R51", "SERPINA1"), ("R67", "HNF4A"),
           ("R07", "ABCG5/8"), ("R54", "FTO")]

pip = pd.read_csv(os.path.join(SUS, "R69_pip.csv"))
feas = pd.read_csv(os.path.join(SUS, "coloc_susie_feasibility.csv"))
abf = pd.read_csv(os.path.join(COL, "coloc_CHOLE_vs_BTC.csv"))
abf = abf[(abf.p12 == 1e-5) & (abf.region_id.isin([r for r, _ in REGIONS]))]
ABF = dict(zip(abf.region_id, abf["PP.H4"]))
abf3 = dict(zip(abf.region_id, abf["PP.H3"]))

sus_path = os.path.join(SUS, "coloc_susie_all.csv")
SUSR = {}
if os.path.exists(sus_path):
    _s = pd.read_csv(sus_path)
    for _, r in _s.iterrows():
        SUSR[r["region"]] = (r["PP.H3"], r["PP.H4"])

FEAS = {r["region"]: r for _, r in feas.iterrows()}

LEAD_POS = 43928850
lead_row = pip.loc[pip.pos == LEAD_POS].iloc[0]

# ================= 画布 =================
fig = plt.figure(figsize=(7.2, 3.75))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.42, 1.0],
              height_ratios=[1, 1], wspace=0.30, hspace=0.62)

# ---------------- A 区域图 ----------------
axA1 = fig.add_subplot(gs[0, 0])
axA2 = fig.add_subplot(gs[1, 0], sharex=axA1)

x = pip.pos.values / 1e6
for ax, pcol, ccol, lab, cscol in (
        (axA1, "p_1", C_CHOLE, "Cholelithiasis", "in_cs_1"),
        (axA2, "p_2", C_BTC, "Biliary tract cancer", "in_cs_2")):
    y = -np.log10(np.clip(pip[pcol].values, 1e-300, None))
    incs = pip[cscol].values.astype(bool)
    ax.vlines(x[~incs], 0, y[~incs], color="#D9D9D9", lw=0.35)
    ax.vlines(x[incs], 0, y[incs], color=ccol, lw=0.5)
    ax.scatter(x[incs], y[incs], s=2.6, color=ccol, zorder=3, linewidths=0)
    ax.axvline(LEAD_POS / 1e6, color="#333333", ls="--", lw=0.7, zorder=2)
    ymax = y.max()
    ax.set_ylim(0, ymax * 1.30)
    ax.set_ylabel(r"$-\log_{10}P$", fontsize=7)
    n_cs = int(FEAS["R69"][f"n_cs_{'chole' if pcol == 'p_1' else 'btc'}"])
    ax.text(0.012, 0.93, f"{lab}   ({n_cs} credible set{'s' if n_cs != 1 else ''})",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.6,
            color=ccol, fontweight="bold")

axA1.set_title("A   PNPLA3 (chr22) \u2014 susie fine-mapping", loc="left", fontsize=8)
axA2.set_xlabel("Position on chr22 (Mb)", fontsize=7.5)
plt.setp(axA1.get_xticklabels(), visible=False)

# 共享 lead 标注（置于坐标轴内侧，避免溢出被裁）
axA1.annotate("rs738408 (shared lead)",
              xy=(LEAD_POS / 1e6, axA1.get_ylim()[1] * 0.97),
              xytext=(0.42, 0.72), textcoords="axes fraction",
              fontsize=6.2, color="#333333", ha="left", va="center",
              arrowprops=dict(arrowstyle="-", lw=0.55, color="#333333",
                              shrinkA=0, shrinkB=1))

n_cs1 = int(FEAS["R69"]["n_cs_chole"]); n_cs2 = int(FEAS["R69"]["n_cs_btc"])
p3, p4 = SUSR.get("R69", (np.nan, np.nan))
axA2.text(0.985, 0.10,
          f"PP.H4 = {p4:.3f}   PP.H3 = {p3:.3f}",
          transform=axA2.transAxes, ha="right", va="bottom", fontsize=6.9,
          fontweight="bold", bbox=dict(fc="white", ec="#999999", lw=0.5,
                                       boxstyle="round,pad=0.30"))

# ---------------- B PIP ----------------
axB = fig.add_subplot(gs[0, 1])
for pcol, ccol, lab in (("pip_1", C_CHOLE, "Cholelithiasis"),
                        ("pip_2", C_BTC, "Biliary tract cancer")):
    axB.plot(x, pip[pcol].values, color=ccol, lw=0.8, label=lab)
axB.axvline(LEAD_POS / 1e6, color="#333333", ls="--", lw=0.7)
cs1 = pip[pip.in_cs_1].pos.values / 1e6
cs2 = pip[pip.in_cs_2].pos.values / 1e6
axB.axvspan(cs1.min(), cs1.max(), color=C_CHOLE, alpha=0.13, lw=0)
axB.axvspan(cs2.min(), cs2.max(), color=C_BTC, alpha=0.16, lw=0)
axB.set_ylabel("PIP", fontsize=7)
axB.set_xlim(x.min(), x.max())
axB.set_ylim(0, max(pip.pip_1.max(), pip.pip_2.max()) * 1.42)
axB.set_title("B   Posterior inclusion probability", loc="left", fontsize=8)
axB.set_xlabel("Position on chr22 (Mb)", fontsize=7.5)
axB.legend(frameon=False, loc="upper right", handlelength=1.2,
           borderaxespad=0.2)

# ---------------- C abf vs susie ----------------
axC = fig.add_subplot(gs[1, 1])
labels = [g for _, g in REGIONS]
xs = np.arange(len(REGIONS))
w = 0.38
v_abf, v_sus, ok = [], [], []
for rid, _ in REGIONS:
    v_abf.append(ABF.get(rid, np.nan))
    if rid in SUSR:
        v_sus.append(SUSR[rid][1]); ok.append(True)
    else:
        v_sus.append(np.nan); ok.append(False)

axC.bar(xs - w / 2, v_abf, w, color=C_ABF, edgecolor="#7A7A7A",
        linewidth=0.5, label="coloc.abf")
sus_drawn = False
for i, (v, o) in enumerate(zip(v_sus, ok)):
    if o:
        axC.bar(i + w / 2, v, w, color=C_BTC, edgecolor="#7A2E2E",
                linewidth=0.5, label=None if sus_drawn else "coloc.susie")
        sus_drawn = True
    else:
        # 该区域胆道癌侧建不出 credible set ⇒ susie 不可用
        axC.bar(i + w / 2, 0.012, w, color="none", edgecolor="#7A7A7A",
                linewidth=0.5, hatch="////")
        axC.text(i + w / 2, 0.035, "n/a", ha="center", va="bottom",
                 fontsize=6.0, color="#555555")
axC.axhline(0.8, color="#333333", ls=":", lw=0.7)
axC.text(3.95, 0.815, "PP.H4 = 0.8", fontsize=5.9,
         color="#333333", ha="center", va="bottom")
axC.set_xticks(xs)
axC.set_xticklabels(labels, fontsize=5.9, rotation=30, ha="right")
axC.set_ylim(0, 1.42)
axC.set_ylabel("PP.H4", fontsize=7)
axC.set_title("C   coloc.abf vs coloc.susie (cholelithiasis vs BTC)",
              loc="left", fontsize=8)
hh = [Line2D([], [], color=C_ABF, lw=3),
      Line2D([], [], color=C_BTC, lw=3),
      Line2D([], [], color="none", marker="s", markersize=4.5,
             markerfacecolor="none", markeredgecolor="#7A7A7A",
             linestyle="none")]
ll = ["coloc.abf", "coloc.susie", "susie n/a (BTC: no CS)"]
axC.legend(hh, ll, frameon=False, loc="upper right", handlelength=1.1,
           borderaxespad=0.15, labelspacing=0.35)

out_png = os.path.join(SUS, "fig_coloc_susie.png")
out_pdf = os.path.join(SUS, "fig_coloc_susie.pdf")
fig.savefig(out_png)
fig.savefig(out_pdf, bbox_inches="tight")

# ---- 交付前裁切自检：任何文本/补丁不得越出画布 ----
fig.canvas.draw()
W, H = fig.canvas.get_width_height()
bad = []
for ax_i in fig.get_axes():
    for t in ax_i.texts + [ax_i.title, ax_i.xaxis.label, ax_i.yaxis.label]:
        if not t.get_text():
            continue
        bb = t.get_window_extent(fig.canvas.get_renderer())
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
            bad.append((ax_i.get_title()[:14], t.get_text()[:28],
                        round(bb.x0), round(bb.y0), round(bb.x1), round(bb.y1)))
if bad:
    print("\n[裁切自检] ✗ 以下文本越出画布:")
    for b in bad:
        print("   ", b)
else:
    print("\n[裁切自检] ✓ 全部文本落在画布内 (%dx%d px)" % (W, H))
print("已输出:", out_png)
print("已输出:", out_pdf)

# ---- 自检：报告数值，供人工核对 ----
print("\n[自检] R69 lead:", lead_row.snp,
      "| PIP1=%.4f  PIP2=%.4f" % (lead_row.pip_1, lead_row.pip_2))
print("[自检] CS1 =", pip.in_cs_1.sum(), "SNP | CS2 =", pip.in_cs_2.sum(), "SNP")
print("\n[自检] PP.H4 对照")
for (rid, g), a, o in zip(REGIONS, v_abf, ok):
    s = SUSR[rid][1] if o else None
    print(f"  {rid} {g:<10} abf={a:.3f}  susie={'%.3f' % s if s else 'n/a'}")
