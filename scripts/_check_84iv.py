#!/usr/bin/env python3
"""自校验：从 84-IV 权威工具集 + BTC 结局 GWAS 独立重算 IVW/Egger/WMed，
与 results/mvmr/mvmr_univariable.csv 比对。通过后才用于出图。"""
import csv, math, os, statistics as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))

# --- 1. 84 IV 暴露侧 ---
ivs = load_tsv(os.path.join(ROOT, "results", "mvmr", "iv_cholelith_clumped.tsv"))
print(f"暴露侧 IV: {len(ivs)} 个")

# --- 2. 结局侧，按 rsid 建索引 ---
out = load_tsv("/tmp/outcome_btc_84iv.tsv")
print(f"结局侧记录: {len(out)} 个")
idx = {}
for r in out:
    for k in ("rsids", "rsid"):
        if k in r and r[k]:
            idx[r[k]] = r
            break

# --- 3. 等位基因对齐 ---
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
recs, n_missing, n_mismatch, n_palin = [], 0, 0, 0
for v in ivs:
    o = idx.get(v["rsid"])
    if o is None:
        n_missing += 1
        continue
    be, ae = v["ref"], v["alt"]
    bo, ao = o["ref"], o["alt"]
    if (ae == ao and be == bo):
        sgn = 1.0
    elif (ae == bo and be == ao):
        sgn = -1.0
    elif (COMP.get(ae) == ao and COMP.get(be) == bo):
        sgn = 1.0          # 链翻转但同向
        n_palin += 1
    elif (COMP.get(ae) == bo and COMP.get(be) == ao):
        sgn = -1.0
        n_palin += 1
    else:
        n_mismatch += 1
        continue
    bx, sx = float(v["beta"]), float(v["se"])
    by, sy = sgn * float(o["beta"]), float(o["sebeta"])
    if sx <= 0 or sy <= 0:
        continue
    recs.append({"rsid": v["rsid"], "bx": bx, "sx": sx, "by": by, "sy": sy,
                 "pexp": float(v["pval"]), "pout": float(o["pval"]),
                 "chrx": v["chr"], "posx": v["pos"]})

print(f"可配对: {len(recs)} | 缺 rsid: {n_missing} | 等位不匹配: {n_mismatch} | 回文: {n_palin}")

# --- 4. 各估计量 ---
# IVW 权重必须是 w_k = beta_exp^2 / se_out^2（Wald 比的一阶 delta）。
# 若误用 w = 1/se_out^2，点估计与 Q 都会系统性错（实测 Q 从 232 虚增到 6×10^4）。
w   = [r["bx"] ** 2 / r["sy"] ** 2 for r in recs]
ratio = [r["by"] / r["bx"] for r in recs]
b_ivw = sum(wi * ri for wi, ri in zip(w, ratio)) / sum(w)
se_ivw_fe = 1.0 / math.sqrt(sum(w))
Q = sum(wi * (ri - b_ivw) ** 2 for wi, ri in zip(w, ratio))
df = len(recs) - 1
kappa = max(1.0, math.sqrt(Q / df))
se_ivw = se_ivw_fe * kappa

# MR-Egger：把 by 直接回归到 bx（截距即多效性），权重是 1/se_out^2。
# 注意与 IVW 的权重不同——IVW 作用于 Wald 比，权重才是 beta_exp^2/se_out^2。
# 两者混用会让截距与斜率同时错（实测截距 +0.148 vs −0.018）。
we = [1.0 / r["sy"] ** 2 for r in recs]
sw = sum(we); swx = sum(wi * r["bx"] for wi, r in zip(we, recs))
swy = sum(wi * r["by"] for wi, r in zip(we, recs))
swxx = sum(wi * r["bx"] ** 2 for wi, r in zip(we, recs))
swxy = sum(wi * r["bx"] * r["by"] for wi, r in zip(we, recs))
den = sw * swxx - swx ** 2
b_eg = (sw * swxy - swx * swy) / den
b_int = (swy - b_eg * swx) / sw
resid = [r["by"] - (b_int + b_eg * r["bx"]) for r in recs]
s2 = sum(wi * resid[i] ** 2 for i, wi in enumerate(we)) / (len(recs) - 2)
se_eg = math.sqrt(sw * s2 / den)
se_int = math.sqrt(s2 * swxx / den)

# 加权中位数（标准误用 bootstrap 才准，此处只对点估计做核对）
order = sorted(range(len(recs)), key=lambda i: ratio[i])
cw, tot = 0.0, sum(w)
b_wmed = None
for i in order:
    cw += w[i]
    if cw >= tot / 2:
        b_wmed = ratio[i]
        break


def fmt(b, se=None):
    orv = math.exp(b)
    if se is None:
        return f"OR {orv:.4f}"
    return f"OR {orv:.4f} ({math.exp(b - 1.959964 * se):.4f}-{math.exp(b + 1.959964 * se):.4f})"


print("\n=== 本脚本重算（84 IV） vs 权威 mvmr_univariable.csv ===")
auth = {r["method"]: r for r in csv.DictReader(open(os.path.join(ROOT, "results", "mvmr", "mvmr_univariable.csv"))) if r["exposure"] == "cholelith"}
ok_all = True
for meth, b, se in (("Inverse variance weighted", b_ivw, se_ivw),
                    ("MR Egger", b_eg, se_eg),
                    ("Weighted median", b_wmed, None)):
    a = auth.get(meth)
    d_or = abs(math.exp(b) - float(a["or"]))
    d_se = abs(se - float(a["se"])) if se is not None else None
    flag = "✓" if (d_or < 2e-3 and (d_se is None or d_se < 2e-4)) else "✗"
    if flag == "✗":
        ok_all = False
    sese = f" | se 本算 {se:.6f} 权威 {float(a['se']):.6f}" if se is not None else ""
    print(f"  {flag} {meth:26s} 本算 {fmt(b, se):32s} 权威 OR {float(a['or']):.4f}{sese}")
print(f"  {'':2s} Egger intercept          {b_int:+.5f} (se {se_int:.4f})")
print(f"  {'':2s} Q (Cochran)              {Q:.1f} on {df} df  -> kappa {kappa:.4f}")
print(f"  {'':2s} se_FE {se_ivw_fe:.6f}  x  kappa {kappa:.4f}  =  se_RE {se_ivw:.6f}")

print("\n=== 与稿件 Table 2B 比对 ===")
print("  稿件: IVW 1.224 (1.061-1.413) P=0.0057 | Egger 1.435 P=0.0033 | WMed 1.242 P=0.0063")
print("         Q=232.4 on 83 df | Egger intercept -0.0183 P=0.097")
print(f"\n闭环结论: {'全部吻合' if ok_all else '存在不一致，需排查'}")
