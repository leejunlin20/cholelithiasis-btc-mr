#!/usr/bin/env python3
"""Figure_MVMR_identifiability — 交付级三面板图

Panel A: 86 个工具变量上 胆石症 β vs 胆囊炎 β 的共线性散点
         （证明两性状效应近乎成比例 → 中介/共线，非独立信号）
Panel B: 条件 F 统计量对照（同规模 401 IV 下 BMI vs 胆囊炎）
Panel C: 胆石症 → 胆道癌 的 OR 森林图（可识别 vs 不可识别模型 + Bai 2026 对照）

硬约束（figure-clipping-audit）：
  1) 出图后遍历 ax.patches / ax.texts，用 get_window_extent().transformed(ax.transData.inverted())
     断言全部落在 xlim/ylim 内，越界即中止；
  2) 用 PIL 检查四边各 10 行像素是否纯白，非白像素 > 50 即判为被裁。

输入：results/mvmr/identifiability_input86.tsv
输出：results/mvmr/Figure_MVMR_identifiability.{png,pdf}
"""
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# 与全稿其它投稿图统一：单一无衬线族 + TrueType 嵌入（此前为字体默认设置，
# 输出 Type 3 字体，且与主图字体族不一致 → 违反 ASN「embed the fonts used」
# 与「uniform lettering」两条要求）。字体族保持 matplotlib 默认 DejaVu Sans，
# 与 Figure 1–5 相同；仅显式声明并改为 TrueType 嵌入，不改变字号与版式。
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "mvmr"
SUPP = ROOT / "results" / "submission" / "supplementary"
INP = OUT / "identifiability_input86.tsv"

# ASN 补充材料页眉要素（官方要求「稿件标题 + 第一作者姓名」）
_SUPP_TITLE = ("Dietary intake, cholelithiasis and biliary tract cancer: a two-sample "
               "Mendelian randomisation study with cross-ancestry replication and "
               "region-definition-aware colocalisation")
_SUPP_STAMP = ("Jun-Lin Li et al. \u2014 Supplementary material, "
               "Supplementary Figure S1      ")
_SUPP_PROVENANCE = ("FinnGen R12, single harmonised freeze; outcome "
                    "C3_BILIARY_GALLBLADDER_EXALLC (2,298 cases); "
                    "505-SNP watch list, 0 missing, 0 allele mismatches")

# ── 配色（浅色主题，红=风险/正，蓝=对照） ──────────────────────────────
C_CHOLE = "#c0392b"   # 胆石症
C_CHOLEC = "#2c6fbb"  # 胆囊炎
C_BMI = "#1a9e6f"     # BMI（参照）
C_BAD = "#b0b0b0"     # 不可识别
C_TXT = "#1a1a1a"
C_GRID = "#dddddd"


def read_input():
    hdr = None
    rows = []
    with open(INP) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = c
                continue
            rows.append(dict(zip(hdr, c)))
    return rows


def ols_no_intercept(x, y, wx=None):
    if wx is None:
        wx = np.ones_like(x)
    sw = np.sum(wx * x * x)
    b = np.sum(wx * x * y) / sw
    resid = y - b * x
    r2 = 1.0 - np.sum(wx * resid ** 2) / np.sum(wx * (y - np.average(y, weights=wx)) ** 2)
    # 斜率 SE（加权最小二乘，过原点）
    m = len(x)
    s2 = np.sum(wx * resid ** 2) / (m - 1)
    se = np.sqrt(s2 / sw)
    return b, se, r2


def guard(ax, tol=1e-6):
    """断言所有 patch / text 落在数据范围内（越界即中止）。"""
    xl, yl = ax.get_xlim(), ax.get_ylim()
    bad = []
    for p in ax.patches:
        bb = p.get_window_extent().transformed(ax.transData.inverted())
        if (bb.x0 < xl[0] - tol or bb.x1 > xl[1] + tol
                or bb.y0 < yl[0] - tol or bb.y1 > yl[1] + tol):
            bad.append((type(p).__name__, bb))
    for t in ax.texts:
        bb = t.get_window_extent().transformed(ax.transData.inverted())
        if (bb.x0 < xl[0] - tol or bb.x1 > xl[1] + tol
                or bb.y0 < yl[0] - tol or bb.y1 > yl[1] + tol):
            bad.append((f"text:{t.get_text()[:20]}", bb))
    if bad:
        for n, bb in bad:
            print(f"  ✗ 越界 {n}: x[{bb.x0:.3f},{bb.x1:.3f}] y[{bb.y0:.3f},{bb.y1:.3f}]")
        raise SystemExit(f"裁切守卫中止：{len(bad)} 个元素越界（xlim={xl}, ylim={yl}）")
    print(f"  ✓ 裁切守卫通过（{len(ax.patches)} patches, {len(ax.texts)} texts）")


def check_white_margin(png, rows=11, thr=50):
    from PIL import Image
    im = Image.open(png).convert("L")
    a = np.asarray(im)
    h, w = a.shape
    edges = {
        "top": a[:rows, :], "bottom": a[h - rows:, :],
        "left": a[:, :rows], "right": a[:, w - rows:],
    }
    ok = True
    for name, blk in edges.items():
        n = int((blk < 250).sum())
        flag = "✓" if n <= thr else "✗ 被裁"
        if n > thr:
            ok = False
        print(f"  {flag} {name:6s} 非白像素 {n}")
    return ok


def main():
    rows = read_input()
    ch = np.array([float(r["ch_b"]) for r in rows])
    chs = np.array([float(r["ch_s"]) for r in rows])
    cc = np.array([float(r["cc_b"]) for r in rows])
    ccs = np.array([float(r["cc_s"]) for r in rows])
    btc = np.array([float(r["bt_b"]) for r in rows])
    btcs = np.array([float(r["bt_s"]) for r in rows])
    rsid = np.array([r["rsid"] for r in rows])
    n = len(ch)
    print(f"输入 {n} 个工具变量")

    # ── 共线性统计（方向与 identifiability_log.txt 一致：y = 胆石 β, x = 胆囊炎 β）──
    slope, se_slope, r2 = ols_no_intercept(cc, ch)          # ch = slope * cc
    r = float(np.corrcoef(ch, cc)[0, 1])
    wgt = 1.0 / np.sqrt(chs ** 2 + slope ** 2 * ccs ** 2)
    _, se_w, _ = ols_no_intercept(cc, ch, wx=wgt ** 2)
    r_w = float(np.average((ch - ch.mean()) * (cc - cc.mean()), weights=wgt ** 2) /
                np.sqrt(np.average((ch - ch.mean()) ** 2, weights=wgt ** 2) *
                        np.average((cc - cc.mean()) ** 2, weights=wgt ** 2)))
    delta = slope
    cond_res = (ch - delta * cc) / np.sqrt(chs ** 2 + delta ** 2 * ccs ** 2)
    n_gt3 = int((np.abs(cond_res) > 3).sum())
    max_z = float(np.abs(cond_res).max())
    print(f"  斜率 {slope:.4f} (SE {se_slope:.4f}), R2={r2:.4f}, r={r:.4f}, "
          f"|z|>3: {n_gt3}/{n} (max {max_z:.2f})")
    print(f"  加权 r = {r_w:.4f}")

    fig = plt.figure(figsize=(13.2, 5.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 0.95, 1.1], wspace=0.46,
                          left=0.062, right=0.985, top=0.800, bottom=0.120)

    # ══════════ Panel A ══════════
    axA = fig.add_subplot(gs[0, 0])
    p_cc = np.array([2 * (1 - 0.5 * (1 + math.erf(abs(b / s) / math.sqrt(2))))
                     for b, s in zip(cc, ccs)])
    sig = p_cc < 5e-8                 # 胆囊炎侧达全基因组显著（86 位点中仅 2 个）
    mid = (~sig) & (p_cc < 1e-3)
    axA.errorbar(cc[~sig], ch[~sig], xerr=ccs[~sig], yerr=chs[~sig],
                 fmt="o", ms=4.2, mfc="#9db8d6", mec="#5b7fa6", mew=0.6,
                 ecolor="#c7d4e2", elinewidth=0.8, capsize=0, alpha=0.9,
                 label=f"Cholecystitis P ≥ 5×10⁻⁸ (n={int((~sig).sum())})", zorder=3)
    axA.errorbar(cc[mid], ch[mid], xerr=ccs[mid], yerr=chs[mid],
                 fmt="^", ms=5.6, mfc="#f5c26b", mec="#a9761a", mew=0.7,
                 ecolor="#e8cfa0", elinewidth=0.9, capsize=0,
                 label=f"P < 1×10⁻³ (n={int(mid.sum())})", zorder=3)
    axA.errorbar(cc[sig], ch[sig], xerr=ccs[sig], yerr=chs[sig],
                 fmt="D", ms=7.2, mfc=C_CHOLEC, mec="#123a63", mew=0.9,
                 ecolor=C_CHOLEC, elinewidth=1.2, capsize=2, alpha=1.0,
                 label=f"P < 5×10⁻⁸ (n={int(sig.sum())})", zorder=4)
    xs = np.linspace(0, max(cc.max(), ch.max() / slope) * 1.02, 50)
    axA.plot(xs, slope * xs, color=C_CHOLE, lw=1.7, ls="-", zorder=2,
             label=f"Zero-intercept fit: slope {slope:.2f} (R²={r2:.2f})")
    axA.plot([0, xs[-1]], [0, xs[-1]], color="#888888", lw=1.1, ls="--", zorder=1,
             label="y = x (equal effects)")
    axA.axhline(0, color="#bbbbbb", lw=0.7); axA.axvline(0, color="#bbbbbb", lw=0.7)
    axA.set_xlabel("β for biliary tract cancer — cholecystitis (K81)", fontsize=8.6)
    axA.set_ylabel("β for biliary tract cancer — cholelithiasis (K80)", fontsize=8.6)
    axA.set_title("A  Near-proportional effects across 86 instruments\n"
                  f"r = {r:.2f};  |z| > 3 after conditioning: {n_gt3}/{n}",
                  fontsize=9.4, loc="left", color=C_TXT, pad=7)
    axA.text(0.97, 0.055, f"slope {slope:.2f} (SE {se_slope:.2f}),  R² = {r2:.2f}\n"
                          f"weighted r = {r_w:.2f}",
             transform=axA.transAxes, ha="right", va="bottom", fontsize=8.2, color="#444444",
             bbox=dict(fc="white", ec="#dddddd", lw=0.6, pad=2.6))
    axA.legend(fontsize=6.4, loc="upper left", framealpha=0.96, borderpad=0.42,
               handletextpad=0.5, labelspacing=0.30)
    axA.tick_params(labelsize=8)
    axA.set_ylim(-0.30, 0.66)
    axA.grid(alpha=0.28, lw=0.5, color=C_GRID)
    axA.set_axisbelow(True)

    # ══════════ Panel B ══════════
    axB = fig.add_subplot(gs[0, 1])
    groups = ["Cholelithiasis\n+ BMI\n(401 IVs)",
              "Cholelithiasis\n+ cholecystitis\n(401 IVs)",
              "Cholelithiasis\n+ cholecystitis\n(86 IVs)"]
    expo_names = ["Cholelithiasis", "Adjusting trait"]
    data = {
        groups[0]: (23.75295, 45.10888, C_BMI),
        groups[1]: (1.512824, 0.9773577, C_BAD),
        groups[2]: (0.960373, 0.863335, C_BAD),
    }
    x = np.arange(len(groups))
    bw = 0.34
    for i, (nm, ci) in enumerate(zip(expo_names, (C_CHOLE, None))):
        vals, cols = [], []
        for g in groups:
            v = data[g][i]
            vals.append(v)
            cols.append(ci if i == 0 else data[g][2])
        pos = x + (i - 0.5) * bw
        axB.bar(pos, vals, bw, color=cols, edgecolor="#555555", lw=0.6, zorder=3)
        for p, v in zip(pos, vals):
            axB.text(p, v + 1.6, f"{v:.2f}", ha="center", va="bottom", fontsize=7.4,
                     color=C_TXT, zorder=5)
    axB.axhline(10, color="#d35400", lw=1.3, ls="--", zorder=4)
    axB.text(0.40, 13.2, "F = 10  (identification threshold)", fontsize=7.0,
             color="#d35400", ha="left", va="bottom", zorder=5)
    axB.set_xticks(x); axB.set_xticklabels(groups, fontsize=7.6)
    axB.set_xlim(-0.58, 2.62)
    axB.set_ylabel("Conditional F-statistic", fontsize=8.6)
    axB.set_ylim(0, 52)
    axB.set_title("B  Conditional F on an identical IV set\n"
                  "diagnostic valid; model not identified",
                  fontsize=9.4, loc="left", color=C_TXT, pad=7)
    axB.legend(handles=[Patch(fc=C_CHOLE, ec="#555555", lw=0.6, label="Cholelithiasis"),
                        Patch(fc=C_BAD, ec="#555555", lw=0.6, label="Adjusting trait")],
               fontsize=7.2, loc="upper right", framealpha=0.94, borderpad=0.45)
    axB.tick_params(labelsize=8)
    axB.grid(alpha=0.28, lw=0.5, axis="y", color=C_GRID)
    axB.set_axisbelow(True)

    # ══════════ Panel C ══════════
    axC = fig.add_subplot(gs[0, 2])
    items = [
        ("BMI-adj. IVW (F 23.8/45.1)", 1.244, 1.131, 1.368, C_CHOLE),
        ("BMI-adj. Lasso", 1.153, 1.052, 1.264, C_CHOLE),
        ("Cholecyst.-adj. Lasso", 1.345, 1.089, 1.662, "#8e9aa6"),
        ("Cholecyst.-adj. median", 1.244, 0.960, 1.613, "#8e9aa6"),
        ("Cholecyst.-adj. IVW (F 1.0)", 1.133, 0.809, 1.586, "#8e9aa6"),
        ("Bai 2026 (published)", 1.45, 0.60, 3.52, "#7d3c98"),
    ]
    ys = np.arange(len(items))[::-1]
    for y, (lab, o, lo, hi, col) in zip(ys, items):
        axC.plot([lo, hi], [y, y], color=col, lw=2.1, solid_capstyle="round", zorder=3)
        axC.plot([o], [y], marker="s", ms=6.2, mfc=col, mec="white", mew=1.0, zorder=4)
        axC.text(hi * 1.055, y, f"{o:.2f} ({lo:.2f}–{hi:.2f})", fontsize=7.0,
                 va="center", ha="left", color=C_TXT, zorder=5)
    axC.axvline(1.0, color="#888888", lw=1.0, ls="--", zorder=1)
    axC.set_yticks(ys); axC.set_yticklabels([i[0] for i in items], fontsize=6.9)
    axC.set_xscale("log")
    axC.set_xlim(0.45, 12.0)
    axC.set_xticks([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0])
    axC.set_xticklabels(["0.5", "0.75", "1.0", "1.5", "2.0", "3.0", "5.0"], fontsize=8)
    axC.set_xlabel("OR for biliary tract cancer", fontsize=8.4)
    axC.set_ylim(-0.75, len(items) - 0.25)
    axC.set_title("C  Effect estimates and interval width\n"
                  "a 5.9-fold interval signals non-identification",
                  fontsize=9.4, loc="left", color=C_TXT, pad=7)
    axC.grid(alpha=0.28, lw=0.5, axis="x", color=C_GRID)
    axC.set_axisbelow(True)
    axC.text(0.985, 0.035, "red: identified (F > 10)   grey: not identified (F ≈ 1)",
             transform=axC.transAxes, ha="right", va="bottom", fontsize=6.5,
             style="italic", color="#666666")

    # ══════════ ASN 补充材料页眉 ══════════
    # ASN Guide for Authors：「Add the manuscript title and first author's name as a
    # header to all supplementary data files.」同时「Compile figure titles and legends
    # … in the manuscript's text file rather than on the figure itself」——
    # 因此原图内主标题（"Adjusting cholelithiasis → … the model is not identifiable"）
    # 已删除，此处只留「稿件标题 + 第一作者」页眉与数据来源行。
    fig.text(0.062, 0.958, _SUPP_TITLE, fontsize=7.2, color="#555555",
             ha="left", va="bottom")
    fig.text(0.062, 0.920, _SUPP_STAMP + _SUPP_PROVENANCE, fontsize=6.6,
             color="#555555", ha="left", va="bottom")

    for ax in (axA, axB, axC):
        for sp in ax.spines.values():
            sp.set_linewidth(0.7); sp.set_color("#999999")

    guard(axA); guard(axB); guard(axC)
    # ⚠️ 多行 ax.title 的 get_window_extent 在某些 matplotlib 版本返回零宽（实测 x1==x0），
    #    因此标题的越界**不能**靠上面的 guard 检出 —— 由下方「四边像素纯白」检查兜底（该方法为真值）。

    OUT.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    png = OUT / "Figure_MVMR_identifiability.png"
    fig.savefig(png, dpi=400, facecolor="white")
    fig.savefig(OUT / "Figure_MVMR_identifiability.pdf", facecolor="white")
    # ASN 投稿用补充图（带页眉）
    fig.savefig(SUPP / "Supplementary_Figure_S1.png", dpi=400, facecolor="white")
    fig.savefig(SUPP / "Supplementary_Figure_S1.pdf", facecolor="white")
    plt.close(fig)
    print(f"已输出 {png.name}")

    if not check_white_margin(png):
        raise SystemExit("边缘像素检查未通过 —— 图件可能被裁切")
    print("✓ 两项交付前检查全部通过")


if __name__ == "__main__":
    main()
