#!/usr/bin/env python3
"""56_make_figure3.py —— Figure 3 胆石症、胆道癌与饮食链条

(A) 变异层面 Wald 比：β_out vs β_exp，附 IVW 拟合（过原点 + 95% 楔形带）与 MR-Egger 拟合线
(B) 四估计量森林图（IVW / MR-Egger / 加权中位数 / 加权众数）
(C) 漏斗图（Wald 比 vs 精度）+ IVW / Cochran's Q 标注
(D) 两步中介：a / b / c / c′ 与经胆石症的中介比例（紧凑表）

★ 本脚本自带闭合校验：从 SNP 级数据重算 IVW 与 MR-Egger，须与权威
  results/mvmr/mvmr_univariable.csv 逐位吻合，否则拒绝出图。

TwoSampleMR 的 mr_egger_regression 内部会先按 sign(b_exp) 翻转 b_out 并把 b_exp 取绝对值，
再以 1/se_out^2 为权重做 WLS；SE 除以 min(1, sigma)，P 用 t 分布（df = n-2）。
不照此实现会得到完全不同的斜率（实测 0.185 vs 0.361），故此处逐步照抄。

★ 出图硬规矩（踩过坑，勿删）：
  1) 四个估计量必须同用同一套工具变量（曾出现加权众数残留 81-IV 旧值 OR 1.226）；
  2) 只允许一次写 PNG（曾因 out.replace(".pdf",".pdf") 空操作把 400 dpi 覆盖成 100 dpi）；
  3) 图内不得出现内部文件路径。
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
import numpy as np
from scipy import stats

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


# =============================================================== 数据 + 校验
ivs = L(os.path.join(R, "mvmr", "iv_cholelith_clumped.tsv"), delimiter="\t")
outp = "/tmp/outcome_btc_84iv.tsv"
if not os.path.exists(outp) or os.path.getsize(outp) == 0:
    raise SystemExit(
        "缺少 84 工具变量的结局提取件。先运行：\n"
        "  python3 scripts/02_extract_outcome_snps.py "
        "data/gwas/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz "
        "/tmp/iv84_rsids.txt /tmp/outcome_btc_84iv.tsv")
ou = {r["rsids"]: r for r in L(outp, delimiter="\t")}

SNP = []
for v in ivs:
    o = ou.get(v["rsid"])
    if o is None:
        raise SystemExit(f"结局缺 {v['rsid']}")
    assert v["alt"] == o["alt"] and v["ref"] == o["ref"], f"{v['rsid']} 等位不一致"
    SNP.append((v["rsid"], float(v["beta"]), float(v["se"]),
                float(o["beta"]), float(o["sebeta"])))
SNP.sort(key=lambda x: x[0])
bx = np.array([s[1] for s in SNP]); sx = np.array([s[2] for s in SNP])
by = np.array([s[3] for s in SNP]); sy = np.array([s[4] for s in SNP])
n = len(bx)
assert n == 84, f"应为 84 个工具变量，实为 {n}"

# ---- IVW（自造单 κ，全项目唯一口径）----
w = bx ** 2 / sy ** 2
ratio = by / bx
b_ivw = float((w * ratio).sum() / w.sum())
se_fe = float(1.0 / math.sqrt(w.sum()))
Q = float((w * (ratio - b_ivw) ** 2).sum())
df = n - 1
kappa = max(1.0, math.sqrt(Q / df))
se_ivw = se_fe * kappa

# ---- MR-Egger（照抄 TwoSampleMR 内部算法）----
sgn = np.where(bx < 0, -1.0, 1.0)
by2, bx2 = by * sgn, np.abs(bx)
We = 1.0 / sy ** 2
X = np.column_stack([np.ones(n), bx2])
XtWX = X.T @ (We[:, None] * X)
coef = np.linalg.solve(XtWX, X.T @ (We * by2))
resid = by2 - X @ coef
sigma = math.sqrt(float((We * resid ** 2).sum() / (n - 2)))
se_raw = np.sqrt(np.diag(np.linalg.inv(XtWX)) * sigma ** 2)
sc = min(1.0, sigma)
b_eg, a_eg = float(coef[1]), float(coef[0])
se_eg = float(se_raw[1] / sc)
se_a = float(se_raw[0] / sc)
p_eg = float(2 * stats.t.sf(abs(b_eg / se_eg), n - 2))
p_a = float(2 * stats.t.sf(abs(a_eg / se_a), n - 2))
Q_eg = sigma ** 2 * (n - 2)

# ---- 闭环校验：IVW 与 MR-Egger 必须逐位吻合；加权中位数仅记录偏差 ----
auth = {r["method"]: r for r in L(os.path.join(R, "mvmr", "mvmr_univariable.csv"))
        if r["exposure"] == "cholelith"}
for meth, b, se in (("Inverse variance weighted", b_ivw, se_ivw),
                    ("MR Egger", b_eg, se_eg)):
    a = auth[meth]
    assert abs(b - float(a["b"])) < 5e-7, \
        f"{meth} 点估计不闭合：本算 {b:.8f} vs 权威 {float(a['b']):.8f}"
    assert abs(se - float(a["se"])) < 5e-7, \
        f"{meth} SE 不闭合：本算 {se:.8f} vs 权威 {float(a['se']):.8f}"
assert abs(a_eg - (-0.018304525)) < 5e-7, f"Egger 截距不闭合：{a_eg:.8f}"
assert abs(Q_eg - 224.6625) < 0.01, f"Egger Q 不闭合：{Q_eg:.4f}"

# ---- 自造式加权中位数（仅作记录，森林图取权威值）----
order = np.argsort(ratio)
cw = np.cumsum(w[order])
b_wmed = float(ratio[order[min(int(np.searchsorted(cw, w.sum() / 2)), n - 1)]])
wm_auth = float(auth["Weighted median"]["b"])
print(f"加权中位数：本算 {b_wmed:.6f} vs 权威 {wm_auth:.6f}"
      f"（差 {abs(b_wmed-wm_auth)/wm_auth*100:.2f}%，森林图取权威值）")
print(f"闭合校验通过：IVW OR {math.exp(b_ivw):.4f} | Egger OR {math.exp(b_eg):.4f} | "
      f"截距 {a_eg:+.5f} | Q {Q:.1f}/{df} | κ {kappa:.4f}")

# ---- 四估计量（森林图直接取权威值，避免混口径）----
FOREST = [
    ("IVW (random effects)", float(auth["Inverse variance weighted"]["or"]),
     float(auth["Inverse variance weighted"]["lci"]),
     float(auth["Inverse variance weighted"]["uci"]),
     float(auth["Inverse variance weighted"]["pval"]), BLUE),
    ("MR-Egger", float(auth["MR Egger"]["or"]), float(auth["MR Egger"]["lci"]),
     float(auth["MR Egger"]["uci"]), p_eg, RED),
    ("Weighted median", float(auth["Weighted median"]["or"]),
     float(auth["Weighted median"]["lci"]), float(auth["Weighted median"]["uci"]),
     float(auth["Weighted median"]["pval"]), GREEN),
    ("Weighted mode", float(auth["Weighted mode"]["or"]),
     float(auth["Weighted mode"]["lci"]), float(auth["Weighted mode"]["uci"]),
     float(auth["Weighted mode"]["pval"]), PURP),
]
# ★ 防回归：四个估计量必须同用 84 个工具变量。历史稿件的加权众数是 81-IV 遗留值
#   （OR 1.226, P 0.006），若此处不拦，图会再次与 Table 2A 自相矛盾。
for _m, _r in auth.items():
    assert int(_r["nsnp"]) == n, f"{_m} 用了 {_r['nsnp']} 个 IV，应为 {n}"
assert len(FOREST) == 4 and all(x[2] is not None and x[3] is not None for x in FOREST), \
    "森林图每个估计量都必须有 CI；出现 None 说明混入了未报 CI 的遗留值"

# ---- 中介（Table 2C）----
MED = [("Alcohol intake frequency", "ukb-b-5779"),
       ("Cheese intake", "ukb-b-1489"),
       ("Dried fruit intake", "ukb-b-16576")]
me = {r["id"]: r for r in L(os.path.join(R, "medi2", "mediation_estimators.csv"))
      if r["estimator"] == "IVW-RE"}
b_path = L(os.path.join(R, "medi2", "mediation_path_b.csv"))[0]
# 字段为 term, b, se_fe, se_re, kappa, Q, Q_df, Q_pval, nsnp, or, lci, uci, p
# 必须取 se_re（随机效应单 κ 口径），不是 se_fe
B_PATH = (float(b_path["b"]), float(b_path["se_re"]))
assert abs(B_PATH[0] - 0.2023964295) < 1e-8, f"b 路径点估计不闭合：{B_PATH[0]}"
assert abs(B_PATH[1] - 0.0732635699) < 1e-8, f"b 路径 SE 不闭合：{B_PATH[1]}"
assert abs(float(b_path["kappa"]) - kappa) < 1e-6, "b 路径 κ 与本脚本不一致"
A_AUTH = {"ukb-b-5779": 0.2275, "ukb-b-1489": -0.2833, "ukb-b-16576": -0.6403}
for _nm, _eid in MED:
    assert abs(float(me[_eid]["a"]) - A_AUTH[_eid]) < 1e-3, f"{_eid} 的 path a 不闭合"

# =============================================================== 画布
fig = plt.figure(figsize=(11.4, 7.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.02], hspace=0.46, wspace=0.26,
                      left=0.068, right=0.975, top=0.942, bottom=0.098)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
axB.set_axis_off()
axD.set_axis_off()


def panel_head(ax, letter, title, y=111.5):
    """四个面板的字母与标题必须同一基线，否则并排时会错位。"""
    ax.text(0, y, letter, fontsize=11.5, fontweight="bold", color=INK, va="top")
    ax.text(7.0, y, title, fontsize=8.6, fontweight="bold", color=INK, va="top")


# ================================================================ Panel A
axA.text(-0.14, 1.115, "A", transform=axA.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axA.text(-0.055, 1.115, "Variant-level Wald ratios", transform=axA.transAxes,
         fontsize=8.6, fontweight="bold", color=INK, va="top", ha="left")

xs = np.linspace(bx.min() * 1.05, bx.max() * 1.05, 200)
axA.fill_between(xs, (b_ivw - 1.96 * se_ivw) * xs, (b_ivw + 1.96 * se_ivw) * xs,
                 color=BLUE, alpha=0.13, lw=0, zorder=1,
                 label="IVW 95% CI")
axA.plot(xs, b_ivw * xs, color=BLUE, lw=1.4, zorder=3, label="IVW")
axA.plot(xs, a_eg + b_eg * xs, color=RED, lw=1.2, ls="--", zorder=3, label="MR-Egger")
axA.scatter(bx, by, s=13, facecolor="white", edgecolor=GREY, linewidth=0.6,
            zorder=4, label="Variant")
axA.axhline(0, color="#E5E8E8", lw=0.8, zorder=0)
axA.axvline(0, color="#E5E8E8", lw=0.8, zorder=0)
axA.set_xlabel("Variant effect on cholelithiasis (log OR)", fontsize=7.4)
axA.set_ylabel("Variant effect on BTC (log OR)", fontsize=7.4)
axA.tick_params(labelsize=7.0, length=2.4, color=LINE)
axA.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axA.spines[s].set_color(LINE)
axA.legend(loc="upper left", fontsize=6.3, frameon=True, framealpha=0.95,
           edgecolor=GRID, borderpad=0.5, handlelength=1.4, handletextpad=0.5)
axA.text(0.975, 0.06,
         f"MR-Egger intercept {a_eg:+.4f} (SE {se_a:.4f}), P = {p_a:.3f}",
         transform=axA.transAxes, ha="right", va="bottom", fontsize=6.3,
         color=RED,
         bbox=dict(boxstyle="round,pad=0.38", facecolor="#FDEDEC",
                   edgecolor=RED, linewidth=0.6))

# ================================================================ Panel B
axB.set_xlim(0, 100)
axB.set_ylim(0, 100)
panel_head(axB, "B", "Pleiotropy-robust estimators")
xL, xR = 40.0, 70.0          # 对数轴区段；右侧 78–99 留给数值标签，防止超出画布被裁
lo_x, hi_x = 0.92, 2.05


def tx(v):
    return xL + (math.log(v) - math.log(lo_x)) / (math.log(hi_x) - math.log(lo_x)) * (xR - xL)


for v in (1.0, 1.25, 1.5, 2.0):
    x = tx(v)
    axB.plot([x, x], [18, 88], color="#E5E8E8", lw=0.7, zorder=1)
    axB.text(x, 14.0, f"{v:g}", ha="center", va="top", fontsize=6.8, color=GREY)
axB.plot([tx(1.0), tx(1.0)], [18, 88], color=GREY, lw=0.9, ls="--", zorder=2)
y0, dy = 80.0, 19.0
for i, (nm, o, lo, hi, p, c) in enumerate(FOREST):
    y = y0 - i * dy
    axB.plot([tx(lo), tx(hi)], [y, y], color=c, lw=1.6, solid_capstyle="round",
             zorder=3)
    axB.plot([tx(o)], [y], "o", color=c, ms=6.0, zorder=4,
             markeredgecolor="white", markeredgewidth=1.0)
    axB.text(99.0, y, f"{o:.3f} ({lo:.3f}–{hi:.3f})", ha="right", va="center",
             fontsize=7.0, color=INK, zorder=4)
    axB.text(xL - 2.6, y + 2.7, nm, ha="right", va="center", fontsize=7.1,
             color=c if i == 0 else INK, fontweight="bold" if i == 0 else "normal",
             zorder=4)
    axB.text(xL - 2.6, y - 2.9, f"P = {p:.3g}" if p >= 1e-3 else f"P = {p:.1e}",
             ha="right", va="center", fontsize=6.3, color=GREY, zorder=4)
axB.text((xL + xR) / 2, 8.0, "Odds ratio for cholelithiasis on BTC (log scale)",
         ha="center", va="top", fontsize=6.8, color=GREY, style="italic")

# ================================================================ Panel C
axC.text(-0.14, 1.115, "C", transform=axC.transAxes, fontsize=11.5,
         fontweight="bold", color=INK, va="top", ha="left")
axC.text(-0.055, 1.115, "Funnel plot", transform=axC.transAxes, fontsize=8.6,
         fontweight="bold", color=INK, va="top", ha="left")

prec = np.abs(bx) / sy
axC.scatter(prec, ratio, s=12, facecolor="white", edgecolor=GREY, linewidth=0.6,
            zorder=3)
axC.axhline(b_ivw, color=BLUE, lw=1.3, zorder=4)
axC.axhline(b_ivw - 1.96 * se_ivw, color=BLUE, lw=0.8, ls=":", zorder=2)
axC.axhline(b_ivw + 1.96 * se_ivw, color=BLUE, lw=0.8, ls=":", zorder=2)
axC.axhline(0, color="#E5E8E8", lw=0.8, zorder=0)
axC.set_xlabel("Precision of the variant–exposure estimate  (|β$_{exp}$| / SE$_{out}$)",
               fontsize=7.4)
axC.set_ylabel("Wald ratio  (β$_{out}$ / β$_{exp}$)", fontsize=7.4)
axC.tick_params(labelsize=7.0, length=2.4, color=LINE)
axC.spines[["top", "right"]].set_visible(False)
for s in ("left", "bottom"):
    axC.spines[s].set_color(LINE)
axC.text(0.975, 0.94,
         f"IVW estimate {b_ivw:+.4f}  (95% CI {b_ivw-1.96*se_ivw:+.3f} to {b_ivw+1.96*se_ivw:+.3f})\n"
         f"Cochran's Q = {Q:.1f} on {df} df  (P = {stats.chi2.sf(Q, df):.1e})",
         transform=axC.transAxes, ha="right", va="top", fontsize=6.3,
         color=GREY, linespacing=1.6,
         bbox=dict(boxstyle="round,pad=0.38", facecolor="#F8F9F9",
                   edgecolor=GRID, linewidth=0.6))

# ================================================================ Panel D
# 紧凑表：行 = 暴露，列 = a / c / c′ / a×b / 中介比例。卡片式竖排会把文字挤破，
# 故改为表格（同样信息量下文字密度下降约 1.4 倍，为后续按印幅缩图留余量）。
axD.set_xlim(0, 100)
axD.set_ylim(0, 100)
panel_head(axD, "D", "Two-step mediation: diet → cholelithiasis → BTC")

axD.text(0.6, 94.0,
         f"path b (cholelithiasis → BTC) = {B_PATH[0]:+.4f} (SE {B_PATH[1]:.4f});  "
         f"{n} instruments",
         ha="left", va="center", fontsize=6.4, color=INK, zorder=3,
         bbox=dict(boxstyle="round,pad=0.32", facecolor="#F4F6F7",
                   edgecolor=GRID, linewidth=0.6))

COLS = [("path a", 35.0), ("path c", 48.0), ("path c′", 61.0),
        ("path a×b", 74.0), ("prop.", 90.0)]
x_name, x_vrule = 1.0, 29.5
y_top, y_bot = 66.0, 26.0
rows_y = [60.0, 47.0, 34.0]

for nm, _eid in MED:
    axD.add_patch(Rectangle((0.5, y_bot), 99.0, y_top - y_bot,
                            facecolor="#FBFCFC", edgecolor=GRID,
                            linewidth=0.7, zorder=1))
    break
axD.plot([0.5, 99.5], [y_top, y_top], color=GRID, lw=0.7, zorder=2)
axD.plot([0.5, 99.5], [y_bot, y_bot], color=GRID, lw=0.7, zorder=2)
axD.plot([x_vrule, x_vrule], [y_bot, y_top], color=GRID, lw=0.7, zorder=2)
for yy in ((rows_y[0] + rows_y[1]) / 2, (rows_y[1] + rows_y[2]) / 2):
    axD.plot([0.5, 99.5], [yy, yy], color=GRID, lw=0.5, zorder=2)

# 表头
axD.text(x_name, y_top + 8.0, "Exposure", ha="left", va="center", fontsize=6.4,
         fontweight="bold", color=INK, zorder=3)
for lab, xc in COLS:
    axD.text(xc, y_top + 8.0, lab, ha="center", va="center", fontsize=6.4,
             fontweight="bold", color=INK, zorder=3)

for (nm, eid), y in zip(MED, rows_y):
    r = me[eid]
    sig = float(r["p_indirect"]) < 0.05
    c = GREEN if sig else ORNG
    axD.text(x_name, y, nm, ha="left", va="center", fontsize=6.4,
             fontweight="bold", color=INK, zorder=3)
    cells = [
        (float(r["a"]), float(r["se_a"]), GREY),
        (float(r["c"]), float(r["se_c"]), GREY),
        (float(r["direct"]), float(r["se_direct"]), GREY),
        (float(r["indirect"]), float(r["se_indirect"]), c),
    ]
    for (val, se, cc), (_lab, xc) in zip(cells, COLS):
        if cc != GREY:
            axD.add_patch(Rectangle((xc - 6.5, y - 6.0), 13.0, 12.0,
                                    facecolor=cc, alpha=0.11, edgecolor=cc,
                                    linewidth=0.6, zorder=2))
        axD.text(xc, y + 1.6, f"{val:+.4f}", ha="center", va="center",
                 fontsize=6.8, color=INK if cc == GREY else cc,
                 fontweight="bold" if cc != GREY else "normal", zorder=3)
        axD.text(xc, y - 4.6, f"{se:.4f}", ha="center", va="center",
                 fontsize=5.6, color=GREY, zorder=3)
    axD.text(COLS[4][1], y + 1.6, f"{float(r['prop']):.3f}", ha="center",
             va="center", fontsize=6.8, color=INK, zorder=3)
    axD.text(COLS[4][1], y - 4.6, f"P {float(r['p_indirect']):.3f}", ha="center",
             va="center", fontsize=5.6, color=c if sig else GREY, zorder=3)

axD.text(0.6, 21.5,
         "a, diet → cholelithiasis;   c, diet → BTC (total);   "
         "c′, diet → BTC given cholelithiasis",
         ha="left", va="center", fontsize=5.8, color=GREY, zorder=3)
axD.text(0.6, 17.0,
         "a×b, indirect effect through cholelithiasis;   prop., a×b / c;   "
         "below each estimate: its SE (for prop., P)",
         ha="left", va="center", fontsize=5.8, color=GREY, zorder=3)
axD.text(0.6, 11.0,
         "Mediated proportions describe the pathway through cholelithiasis, "
         "not proof of mediation.",
         ha="left", va="center", fontsize=5.8, color=GREY, style="italic", zorder=3)
axD.text(0.6, 6.0,
         "Indirect effects use the corrected random-effects IVW convention.",
         ha="left", va="center", fontsize=5.8, color=GREY, style="italic", zorder=3)

# ------------------------------------------------------------------ 标题
# ⚠ 官方口径（ASN Guide for Authors 原文）："Compile figure titles and legends on one
# or more pages in the manuscript's text file rather than on the figure itself."
# 故此处**不再在图内印主标题与副标题**；图题与图注已写入稿件文本（Figure 1–5 图注）。
# 面板字母 A/B/C/D 与图例键保留在图内（官方允许 "figure keys on figures"）。

# =============================================================== 落盘（只写一次 PNG！）
stem = os.path.join(R, "submission", "Figure3_cholelithiasis_btc_chain")
os.makedirs(os.path.dirname(stem), exist_ok=True)
png, pdf = stem + ".png", stem + ".pdf"
fig.savefig(png, dpi=400, facecolor="white")
fig.savefig(pdf, facecolor="white")
px = Image.open(png).size
assert px[0] == int(round(fig.get_figwidth() * 400)), \
    f"PNG 分辨率不符：{px}，期望宽 {int(round(fig.get_figwidth()*400))}——" \
    "检查是否有第二次 savefig 把高分辨率版本覆盖了"
print("Figure 3 已写出:", png, f"({px[0]}×{px[1]} @400dpi)")
print(f"  IVW b={b_ivw:.8f} se={se_ivw:.8f} | Egger b={b_eg:.8f} se={se_eg:.8f} | 截距 {a_eg:+.8f}")
