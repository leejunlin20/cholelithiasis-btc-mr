#!/usr/bin/env python3
"""
27_direction_concordance_v2.py
------------------------------
重算「区域 top5 SNP 方向一致性」，区域坐标改用修正后的 regions_v2.tsv。

修正点：原 direction_concordance.tsv 的 R07 行取自 chr2:44.5-44.7Mb，
那实际是 **CAMKMT** 位点，却被贴上了 ABCG5/ABCG8 的标签。本脚本按
regions_v2.tsv 的 R07a（ABCG5/ABCG8）/ R07b（CAMKMT）分别重算。
"""
import os

BASE = os.environ.get("MR_BASE", ".")
DIR_DATA = f"{BASE}/data/coloc"
OUT = f"{BASE}/results/coloc_susie"

REG = {  # region_id -> (chrom, start, end, gene_label)
    "R69":  (22, None, None, "PNPLA3"),
    "R51":  (14, None, None, "SERPINA1"),
    "R67":  (20, None, None, "HNF4A"),
    "R07a": (2, 43306615, 44470983, "ABCG5/ABCG8"),
    "R07b": (2, 44470984, 45694235, "CAMKMT"),
    "R54":  (16, None, None, "FTO"),
}


def load_v2():
    spans = {}
    with open(f"{BASE}/results/coloc/regions_v2.tsv") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            spans[p[ix["region_id"]]] = (int(p[ix["chrom"]]), int(p[ix["start"]]),
                                         int(p[ix["end"]]))
    return spans


def read_trait(fn):
    rows = {}
    with open(f"{DIR_DATA}/{fn}") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            key = (int(p[ix["#chrom"]]), int(p[ix["pos"]]))
            rows[key] = dict(ref=p[ix["ref"]], alt=p[ix["alt"]],
                             beta=float(p[ix["beta"]]), p=float(p[ix["pval"]]),
                             af=float(p[ix["af_alt"]]),
                             snp=(p[ix["rsids"]] or "").split(",")[0])
    return rows


CHO = read_trait("CHOLE_regions.tsv")
BTC = read_trait("BTC_regions.tsv")
spans = load_v2()

header = ["region", "gene", "snp", "pos_hg38", "af_alt",
          "beta_chole", "p_chole", "beta_btc", "p_btc", "concordant"]
lines = ["\t".join(header)]
summary = []

for rid, gene in [(k, v[3]) for k, v in REG.items()]:
    c0, s0, e0 = spans[rid]
    cand = [k for k in CHO if k[0] == c0 and s0 <= k[1] <= e0 and k in BTC]
    cand.sort(key=lambda k: CHO[k]["p"])
    top = cand[:5]
    n_yes = 0
    for k in top:
        a, b = CHO[k], BTC[k]
        ok = "YES" if (a["beta"] > 0) == (b["beta"] > 0) else "NO_OPPOSITE"
        n_yes += ok == "YES"
        lines.append("\t".join(map(str, [rid, gene, a["snp"], k[1],
                                         "%.4f" % a["af"], "%.4f" % a["beta"],
                                         "%.3g" % a["p"], "%.4f" % b["beta"],
                                         "%.3g" % b["p"], ok])))
    summary.append((rid, gene, s0, e0, len(cand), n_yes, len(top)))

with open(f"{OUT}/direction_concordance_v2.tsv", "w") as o:
    o.write("\n".join(lines) + "\n")

print("区域 top5 SNP 方向一致性（regions_v2 坐标）")
print("%-6s %-12s %10s %10s %7s %9s" % ("region", "gene", "start", "end", "n_snp", "同向"))
for s in summary:
    print("%-6s %-12s %10d %10d %7d %5d/%d" % (s[0], s[1], s[2], s[3], s[4], s[5], s[6]))
print(f"\n输出: {OUT}/direction_concordance_v2.tsv")
