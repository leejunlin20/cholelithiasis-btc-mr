#!/usr/bin/env python3
"""61_make_figure5.py —— Figure 5 跨族裔复制（BioBank Japan，日本 / GRCh37）

(A) 位点级一致性：FinnGen 胆石症 β vs BBJ 胆石症 β（60 位点），含截距回归线 + y=x 参考线
(B) 聚合方向检验：以 1/SE_BBJ² 为权的加权回归（两张小图）
    B1  FinnGen 胆石症 β → BBJ 胆道癌 β
    B2  FinnGen 胆道癌 β → BBJ 胆道癌 β
(C) BBJ 中主结局 / 阳性对照 / 阴性对照的 IVW 森林图
(D) 可检出最小 OR 随例数变化的效能曲线，标注 FinnGen R12 与 BBJ 两个队列的位置

★ 本脚本自带闭合校验（不闭合即拒绝出图）：
  1) 三条跨血统回归斜率与 SE / z 必须复现 `BBJ跨队列复制报告.md` §6.3 的印刷值
     （0.915 / z 11.17；0.229 / z 1.43；0.304 / z 1.69，EIV 校正 0.920 / 0.232 / 0.430）；
  2) 森林图 OR 必须等于 `bbj_mr_results.csv` 中 IVW 行的 exp(b)；
  3) 效能曲线的隐含 R²（= 1 / (SE² · n_eff)）必须在全表恒定，
     且 or_min80 = exp(2.8015852 × SE)。

★ 出图硬规矩（踩过坑，勿删）：
  1) 只允许一次写 PNG（`out.replace(".pdf",".pdf")` 是空操作，会把 400 dpi 覆盖成 100 dpi）；
  2) 图内不得出现仓库内部路径；
  3) 图内不印主标题与副标题（ASN 官方：图题图注写在稿件文本里，非图内）。
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
BBJ = os.path.join(R, "bbj_replication")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
RED, GREEN, GREY, BLUE = "#C0392B", "#1E8449", "#566573", "#2471A3"
PURP, ORNG, INK, LINE, GRID = "#7D3C98", "#D68910", "#212F3D", "#95A5A6", "#BDC3C7"
Z80 = 1.959963985 + 0.841621234      # 80% 功效对应的 z 之和 = 2.8015852


def L(p, **k):
    with open(p) as f:
        return list(csv.DictReader(f, **k))


def num(r, k):
    try:
        return float(r.get(k, ""))
    except (TypeError, ValueError):
        return None


def wls(x, y, sx, sy):
    """误差变量校正的加权最小二乘（照抄 scripts/28_sign_concordance.py）。
    sy 全为 0 时退化为等权 OLS。返回 (b, b_corr, se_b, n)。"""
    use_w = any(s > 0 for s in sy)
    w = [1.0 / (sy[i] ** 2) if (use_w and sy[i] > 0) else 1.0 for i in range(len(x))]
    S = sum(w)
    mx = sum(w[i] * x[i] for i in range(len(x))) / S
    my = sum(w[i] * y[i] for i in range(len(x))) / S
    Sxx = sum(w[i] * (x[i] - mx) ** 2 for i in range(len(x)))
    Sxy = sum(w[i] * (x[i] - mx) * (y[i] - my) for i in range(len(x)))
    b = Sxy / Sxx
    Sxx_corr = Sxx - sum(w[i] * sx[i] ** 2 for i in range(len(x)))
    b_corr = Sxy / Sxx_corr if Sxx_corr > 0 else float("nan")
    resid = [y[i] - my - b * (x[i] - mx) for i in range(len(x))]
    s2 = sum(w[i] * resid[i] ** 2 for i in range(len(x))) / (len(x) - 2)
    return b, b_corr, math.sqrt(s2 / Sxx), len(x), mx, my


# ======================================================= 数据
loci = L(os.path.join(BBJ, "locus_replication_v2.tsv"), delimiter="\t")
mr = L(os.path.join(BBJ, "bbj_mr_results.csv"))
pw = L(os.path.join(R, "power", "power_curve_by_cases.csv"))


def gather(xk, yk, ys):
    x, y, sx, sy = [], [], [], []
    for r in loci:
        a, b = num(r, xk), num(r, yk)
        sa = num(r, xk.replace("_beta", "_se"))
        sb = num(r, ys) if ys else 0.0
        if a is None or b is None or sa is None:
            continue
        x.append(a); y.append(b); sx.append(sa); sy.append(sb if sb else 0.0)
    return x, y, sx, sy


PAIRS = {
    "A": ("fg_chole_beta", "bbj_chole_beta", None),        # 胆石症 → 胆石症（等权）
    "B1": ("fg_chole_beta", "bbj_btc_beta", "bbj_btc_se"),  # 胆石症 → 胆道癌
    "B2": ("fg_btc_beta", "bbj_btc_beta", "bbj_btc_se"),    # 胆道癌 → 胆道癌
}
FIT = {}
for tag, (xk, yk, ys) in PAIRS.items():
    x, y, sx, sy = gather(xk, yk, ys)
    FIT[tag] = dict(zip(("x", "y", "wls"),
                        (x, y, wls(x, y, sx, sy))))

# ------------------------------------------------- 闭合校验 1：斜率 / SE / z
AUTH = {           # 报告 §6.3 印刷值：b, se, b_corr
    "A":  (0.915, 0.0819, 0.920),
    "B1": (0.229, 0.1604, 0.2321),
    "B2": (0.304, 0.1799, 0.4295),
}
for tag, (eb, ese, ebc) in AUTH.items():
    b, bc, se, n, mx, my = FIT[tag]["wls"]
    assert abs(b - eb) < 5e-4, f"[{tag}] 斜率不闭合：{b:.4f} vs {eb}"
    assert abs(se - ese) < 5e-4, f"[{tag}] SE 不闭合：{se:.4f} vs {ese}"
    assert abs(bc - ebc) < 5e-4, f"[{tag}] EIV 校正斜率不闭合：{bc:.4f} vs {ebc}"
assert FIT["A"]["wls"][3] == 60 and FIT["B1"]["wls"][3] == 61 and FIT["B2"]["wls"][3] == 61

# ------------------------------------------------- 闭合校验 2：森林图用的 OR
IVW = {}
for r in mr:
    if r["method"] == "Inverse variance weighted":
        IVW[r["pheno"]] = r
ORDER = [
    ("BtC",            "Biliary tract cancer",   418,   INK,  "main outcome"),
    ("Cholelithiasis", "Cholelithiasis",         9305,  RED,  "positive ctrl"),
    ("Cholecystitis",  "Cholecystitis",          1860,  RED,  "positive ctrl"),
    ("GB_polyp",       "Gallbladder polyp",      843,   ORNG, "bile-duct"),
    ("Urolithiasis",   "Urolithiasis",           11699, GREEN, "negative ctrl"),
    ("PaC",            "Pancreatic cancer",      499,   GREEN, "negative ctrl"),
]
AUTH_OR = {"BtC": 1.2478, "Cholelithiasis": 1.9481, "Cholecystitis": 1.3675,
           "GB_polyp": 1.4705, "Urolithiasis": 1.0028, "PaC": 0.8230}
for ph, _lab, ncase, _c, _ctrl in ORDER:
    r = IVW[ph]
    assert int(r["ncase"]) == ncase, f"{ph} 例数不符：{r['ncase']} vs {ncase}"
    orv = math.exp(float(r["b"]))
    assert abs(orv - AUTH_OR[ph]) < 5e-4, f"{ph} OR 不闭合：{orv:.4f} vs {AUTH_OR[ph]}"

# ------------------------------------------------- 闭合校验 3：效能曲线
NW = [int(float(p["ncase"])) for p in pw]
R2S = []
for p in pw:
    se, ne = float(p["se"]), float(p["n_eff"])
    R2S.append(1.0 / (se ** 2 * ne))
    assert abs(math.exp(Z80 * se) - float(p["or_min80"])) < 1e-4, \
        f"or_min80 不闭合 @ {p['ncase']} 例"
R2_CURVE = sum(R2S) / len(R2S)
assert max(R2S) - min(R2S) < 1e-4, f"效能曲线的隐含 R² 不恒定：{min(R2S)}–{max(R2S)}"
MDE_FG, MDE_BBJ = 1.22783706102307, 1.55063609792366
assert abs(math.exp(Z80 * 0.0732635699) - MDE_FG) < 1e-6
assert abs(math.exp(Z80 * 0.15657750821466) - MDE_BBJ) < 1e-6

print(f"闭合校验通过：A {FIT['A']['wls'][0]:.4f} (z {FIT['A']['wls'][0]/FIT['A']['wls'][2]:.2f}) | "
      f"B1 {FIT['B1']['wls'][0]:.4f} (z {FIT['B1']['wls'][0]/FIT['B1']['wls'][2]:.2f}) | "
      f"B2 {FIT['B2']['wls'][0]:.4f} (z {FIT['B2']['wls'][0]/FIT['B2']['wls'][2]:.2f}) | "
      f"效能曲线 R² = {R2_CURVE:.5f}")

# ======================================================= 画布
fig = plt.figure(figsize=(11.4, 7.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.06, 1.0],
                      hspace=0.40, wspace=0.30,
                      left=0.150, right=0.978, top=0.945, bottom=0.088)

# ======================================================== Panel A
axA = fig.add_subplot(gs[0, 0])
axA.text(-0.145, 1.09, "A", transform=axA.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axA.text(-0.055, 1.09, "Variant-level transferability of the instrument",
         transform=axA.transAxes, fontsize=8.6, fontweight="bold",
         color=INK, va="top", ha="left")

xa, ya = np.array(FIT["A"]["x"]), np.array(FIT["A"]["y"])
bA, bcA, seA, nA, mxA, myA = FIT["A"]["wls"]
lim = (min(xa.min(), ya.min()) - 0.05, max(xa.max(), ya.max()) + 0.05)
axA.plot(lim, lim, color=GRID, lw=1.0, ls="--", zorder=2,
         label="identity ($y=x$)")
xs = np.linspace(lim[0], lim[1], 50)
axA.plot(xs, myA + bA * (xs - mxA), color=BLUE, lw=1.5, zorder=3,
         label="weighted fit")
axA.scatter(xa, ya, s=13, facecolor="white", edgecolor=GREY, linewidth=0.6,
            zorder=4, label=None)
axA.axhline(0, color="#E5E8E8", lw=0.8, zorder=0)
axA.axvline(0, color="#E5E8E8", lw=0.8, zorder=0)
axA.set_xlim(*lim); axA.set_ylim(*lim)
axA.set_xlabel("Variant effect on cholelithiasis, FinnGen (log OR)", fontsize=7.4)
axA.set_ylabel("Variant effect on cholelithiasis,\nBioBank Japan (log OR)", fontsize=7.4)
axA.tick_params(labelsize=7.0, length=2.4, color=LINE)
axA.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axA.spines[s].set_color(LINE)
axA.legend(loc="upper left", fontsize=6.3, frameon=True, framealpha=0.95,
           edgecolor=GRID, borderpad=0.5, handlelength=1.5, handletextpad=0.5)
axA.text(0.975, 0.06,
         f"slope {bA:.3f} (SE {seA:.3f}), $z$ = {bA/seA:.2f}\n"
         f"error-in-variables slope {bcA:.3f};  {nA} lead SNPs\n"
         f"equal weights (no BBJ standard error for this trait)",
         transform=axA.transAxes, ha="right", va="bottom", fontsize=6.3,
         color=INK, linespacing=1.5,
         bbox=dict(boxstyle="round,pad=0.38", facecolor="#F4F6F7",
                   edgecolor=GRID, linewidth=0.6))

# ======================================================== Panel B（两小图）
sub = gs[0, 1].subgridspec(2, 1, hspace=0.80)
BLAB = [("B1", "FinnGen cholelithiasis $\\beta$", FIT["B1"]),
        ("B2", "FinnGen biliary tract cancer $\\beta$", FIT["B2"])]
for i, (tag, xlab, fit) in enumerate(BLAB):
    ax = fig.add_subplot(sub[i])
    if i == 0:
        ax.text(-0.235, 1.25, "B", transform=ax.transAxes, fontsize=11.5,
                fontweight="bold", color=INK, va="top", ha="left")
        ax.text(-0.055, 1.25, "Aggregate direction tests", transform=ax.transAxes,
                fontsize=8.6, fontweight="bold", color=INK, va="top", ha="left")
    x, y = np.array(fit["x"]), np.array(fit["y"])
    b, bc, se, n, mx, my = fit["wls"]
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, my + b * (xs - mx), color=BLUE, lw=1.3, zorder=3)
    ax.scatter(x, y, s=9, color=GREY, alpha=0.75, linewidths=0, zorder=2)
    ax.axhline(0, color="#E5E8E8", lw=0.8, zorder=0)
    ax.axvline(0, color="#E5E8E8", lw=0.8, zorder=0)
    ax.set_xlabel(xlab, fontsize=6.9)
    ax.set_ylabel("BioBank Japan\nBTC $\\beta$", fontsize=6.9, linespacing=1.3)
    ax.tick_params(labelsize=6.5, length=2.2, color=LINE)
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(LINE)
    ax.text(0.975, 0.94, f"slope {b:.3f} (SE {se:.3f}), $z$ = {b/se:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.3,
            color=INK,
            bbox=dict(boxstyle="round,pad=0.30", facecolor="#F4F6F7",
                      edgecolor=GRID, linewidth=0.6))

# ======================================================== Panel C
axC = fig.add_subplot(gs[1, 0])
axC.text(-0.145, 1.10, "C", transform=axC.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axC.text(-0.055, 1.10, "Outcome panel in BioBank Japan (IVW)",
         transform=axC.transAxes, fontsize=8.6, fontweight="bold",
         color=INK, va="top", ha="left")

ys = np.arange(len(ORDER))[::-1]
for y, (ph, lab, ncase, col, ctrl) in zip(ys, ORDER):
    r = IVW[ph]
    orv = math.exp(float(r["b"]))
    lo, hi = math.exp(float(r["lo_ci"])), math.exp(float(r["hi_ci"]))
    p = float(r["pval"])
    axC.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="round", zorder=3)
    axC.plot([orv], [y], "o", ms=5.6, color=col, markeredgecolor="white",
             markeredgewidth=1.0, zorder=4)
    axC.text(2.38, y, f"{orv:.3f} ({lo:.3f}–{hi:.3f})", ha="left", va="center",
             fontsize=6.9, color=INK, zorder=4)
    axC.text(4.20, y, "P = %.2g" % p if p >= 1e-3 else "P = %.1e" % p,
             ha="left", va="center", fontsize=6.2, color=GREY, zorder=4)
axC.axvline(1.0, color=GREY, lw=0.8, ls=(0, (4, 2)), zorder=2)
axC.set_xscale("log")
axC.set_xlim(0.60, 5.60)
axC.set_xticks([0.7, 1.0, 1.5, 2.0])
axC.set_xticklabels(["0.7", "1.0", "1.5", "2.0"])
axC.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axC.set_yticks(ys)
axC.set_yticklabels([f"{lab}\n{ncase:,} cases · {ctrl}"
                     for _p, lab, ncase, _c, ctrl in ORDER],
                    fontsize=6.0, linespacing=1.35)
axC.set_ylim(-0.75, len(ORDER) - 0.25)
axC.set_xlabel("IVW odds ratio per log-odds cholelithiasis\n"
               "(FinnGen R12 $\\rightarrow$ BioBank Japan)", fontsize=7.4)
axC.tick_params(axis="y", length=0)
axC.tick_params(axis="x", labelsize=6.6, length=2.4, color=LINE)
axC.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axC.spines[s].set_color(LINE)


# ======================================================== Panel D
axD = fig.add_subplot(gs[1, 1])
axD.text(-0.155, 1.10, "D", transform=axD.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axD.text(-0.062, 1.10, "Smallest detectable effect vs. sample size",
         transform=axD.transAxes, fontsize=8.6, fontweight="bold",
         color=INK, va="top", ha="left")

nc = np.array([float(p["ncase"]) for p in pw])
mde = np.array([float(p["or_min80"]) for p in pw])
grid = np.logspace(math.log10(350), math.log10(160000), 200)
se_of = lambda x: math.sqrt(1.0 / (R2_CURVE * (x * 378749.0 / (x + 378749.0))))
axD.plot(grid, [math.exp(Z80 * se_of(g)) for g in grid], color=BLUE, lw=1.6,
         zorder=3, label="80% power")
axD.axhline(1.22, color=RED, lw=1.0, ls=":", zorder=2)
axD.text(360, 1.250, "European estimate\n(OR 1.22)", fontsize=6.0,
         color=RED, va="bottom", ha="left", linespacing=1.5)
for ncase, col, lab, dy in ((418, ORNG, "BioBank Japan\n418 cases", 0.10),
                            (2298, INK, "FinnGen R12\n2,298 cases", 0.085)):
    axD.axvline(ncase, color=col, lw=0.9, ls="--", zorder=2)
    axD.scatter([ncase], [math.exp(Z80 * se_of(ncase))], s=26, color=col,
                zorder=5, edgecolor="white", linewidth=0.8)
axD.annotate("BioBank Japan\n418 cases  ·  MDE 1.55", xy=(418, MDE_BBJ),
             xytext=(500, 1.700), fontsize=6.3, color=ORNG, ha="left", va="center",
             linespacing=1.5,
             arrowprops=dict(arrowstyle="-", color=ORNG, lw=0.7,
                             shrinkA=1, shrinkB=3))
axD.annotate("FinnGen R12\n2,298 cases  ·  MDE 1.23", xy=(2298, MDE_FG),
             xytext=(3000, 1.440), fontsize=6.3, color=INK, ha="left", va="center",
             linespacing=1.5,
             arrowprops=dict(arrowstyle="-", color=INK, lw=0.7,
                             shrinkA=1, shrinkB=3))
axD.set_xscale("log"); axD.set_yscale("log")
axD.set_xticks([500, 1000, 2000, 5000, 10000, 50000])
axD.set_xticklabels(["500", "1,000", "2,000", "5,000", "10,000", "50,000"])
axD.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axD.set_yticks([1.05, 1.10, 1.20, 1.30, 1.50, 1.70])
axD.set_yticklabels(["1.05", "1.10", "1.20", "1.30", "1.50", "1.70"])
axD.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axD.set_xlim(350, 160000); axD.set_ylim(1.01, 1.82)
axD.set_xlabel("Number of cases", fontsize=7.4)
axD.set_ylabel("Smallest detectable odds ratio\n(at 80% power)", fontsize=7.4,
               linespacing=1.4)
axD.tick_params(labelsize=6.6, length=2.4, color=LINE)
axD.spines[["top", "right"]].set_visible(False)
for s in ("left", "right", "bottom"):
    pass
for s in ("left", "bottom"):
    axD.spines[s].set_color(LINE)
axD.text(0.015, 0.035,
         "80%% power; 84-variant cholelithiasis\n"
         "instrument set (liability-scale $R^2$ = %.3f)" % R2_CURVE,
         transform=axD.transAxes, ha="left", va="bottom", fontsize=6.0,
         color=GREY, linespacing=1.5)

# --------------------------------------------- 不印主标题（见文件头硬规矩 3）

# ======================================================= 落盘（只写一次 PNG！）
stem = os.path.join(R, "submission", "Figure5_cross_ancestry_replication")
os.makedirs(os.path.dirname(stem), exist_ok=True)
png, pdf = stem + ".png", stem + ".pdf"
fig.savefig(png, dpi=400, facecolor="white")
fig.savefig(pdf, facecolor="white")
px = Image.open(png).size
assert px[0] == int(round(fig.get_figwidth() * 400)), \
    f"PNG 分辨率不符：{px}，期望宽 {int(round(fig.get_figwidth()*400))}——" \
    "检查是否有第二次 savefig 把高分辨率版本覆盖了"
print("Figure 5 已写出:", png, f"({px[0]}×{px[1]} @400dpi)")
