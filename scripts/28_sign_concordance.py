#!/usr/bin/env python3
"""
28_sign_concordance.py
----------------------
跨血统「效应方向一致性」的聚合检验。

逐个 SNP 数方向（51% vs 50%）没有信息量：FinnGen 胆道癌的逐 SNP 效应本身极小
（|β| 多为 0.01~0.09），BBJ 的 SE 却是 0.07~0.35，单个 SNP 的方向基本是噪声。
更有信息量的做法是：以 1/SE²(BBJ) 为权，把 BBJ 的 β 对 FinnGen 的 β 做加权回归，
检验斜率是否显著 > 0（即「FinnGen 中效应较大的等位基因，在 BBJ 中是否也较大」）。
按 (SE_BBJ² + SE_FinnGen²) 做误差变量(regression dilution)校正。
"""
import csv, math
import os

BASE = os.environ.get("MR_BASE", ".")
P = f"{BASE}/results/bbj_replication/locus_replication_v2.tsv"


def num(r, k):
    v = r.get(k, "")
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def wls(x, y, sx, sy):
    """误差变量校正的加权最小二乘：y = a + b*x；传入的 sy 全为 0 时退化为等权 OLS"""
    use_w = any(s > 0 for s in sy)
    w = [1.0 / (sy[i] ** 2) if (use_w and sy[i] > 0) else 1.0 for i in range(len(x))]
    S = sum(w)
    mx = sum(w[i] * x[i] for i in range(len(x))) / S
    my = sum(w[i] * y[i] for i in range(len(x))) / S
    Sxx = sum(w[i] * (x[i] - mx) ** 2 for i in range(len(x)))
    Sxy = sum(w[i] * (x[i] - mx) * (y[i] - my) for i in range(len(x)))
    b = Sxy / Sxx
    # 误差变量校正：减去 x 的测量误对方差
    Sxx_corr = Sxx - sum(w[i] * sx[i] ** 2 for i in range(len(x)))
    b_corr = Sxy / Sxx_corr if Sxx_corr > 0 else float("nan")
    # 斜率的标准误（近似）
    resid = [y[i] - my - b * (x[i] - mx) for i in range(len(x))]
    s2 = sum(w[i] * resid[i] ** 2 for i in range(len(x))) / (len(x) - 2)
    se_b = math.sqrt(s2 / Sxx)
    return b, b_corr, se_b, len(x)


rows = [r for r in csv.DictReader(open(P), delimiter="\t")]
pairs = {
    "FinnGen 胆石症 β → BBJ 胆道癌 β": ("fg_chole_beta", "fg_chole_se", "bbj_btc_beta", "bbj_btc_se"),
    "FinnGen 胆道癌 β → BBJ 胆道癌 β": ("fg_btc_beta", "fg_btc_se", "bbj_btc_beta", "bbj_btc_se"),
    "FinnGen 胆石症 β → BBJ 胆石症 β": ("fg_chole_beta", "fg_chole_se", "bbj_chole_beta", None),
}

print("=" * 78)
print("跨血统效应方向一致性 —— 误差变量校正的加权回归（权 = 1/SE_BBJ²）")
print("=" * 78)
for name, (xk, xs, yk, ys) in pairs.items():
    x, y, sx, sy = [], [], [], []
    for r in rows:
        a, b = num(r, xk), num(r, yk)
        sa, sb = num(r, xs), num(r, ys) if ys else 0.0
        if a is None or b is None or sa is None:
            continue
        sy.append(sb if sb else 0.0)
        sx.append(sa)
        x.append(a)
        y.append(b)
    if len(x) < 10:
        print(f"{name}: 可用 {len(x)} 个，样本太少")
        continue
    b, bc, se, n = wls(x, y, sx, sy)
    z = b / se
    print(f"\n{name}")
    print(f"  可用位点 n = {n}")
    print(f"  斜率（Naive）      b = {b:+.4f}  (SE {se:.4f},  z = {z:+.2f})")
    print(f"  斜率（误差变量校正）b = {bc:+.4f}")
    print(f"  ⇒ {'✅ 方向显著一致 (z>1.96)' if z > 1.96 else ('⚠️ 方向一致但未达显著 (0<z<1.96)' if z > 0 else '❌ 方向不一致 (z<0)')}")

print("\n" + "=" * 78)
print("说明：naive 斜率受回归稀释（regression dilution）向下偏，误差变量校正给出上界；")
print("      真实斜率介于两者之间。z 检验基于 naive 斜率（偏保守）。")
print("=" * 78)
