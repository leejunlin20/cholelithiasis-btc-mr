#!/usr/bin/env python3
"""
10_define_coloc_regions.py
--------------------------
以胆石病（K11_CHOLELITH）的独立先导 SNP 为锚，定义 COLOC 区域。

- 输入: results/medi/chole_iv_clumped.tsv  (clump r2<0.001, 10Mb)
- 窗口: 每个先导 SNP ±WINDOW bp（默认 500kb，即 1Mb 区域）
- 重叠区域合并（同染色体上区间相交即合并）
- 输出: results/coloc/regions.tsv  (region_id, chr, start, end, n_lead, lead_snps)

区域坐标统一 GRCh38（FinnGen R12 原生 build），1-based 闭区间。
"""
import os
import sys

BASE = os.environ.get("MR_BASE", ".")
IN_TSV = os.path.join(BASE, "results/medi/chole_iv_clumped.tsv")
OUT_TSV = os.path.join(BASE, "results/coloc/regions.tsv")

WINDOW = int(sys.argv[1]) if len(sys.argv) > 1 else 500_000


def load_leads(path):
    leads = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            leads.append(
                dict(
                    rsid=p[idx["rsid"]],
                    chrom=int(p[idx["chr"]]),
                    pos=int(p[idx["pos"]]),
                    pval=float(p[idx["pval"]]),
                    gene=p[idx.get("nearest_gene", idx["rsid"])] if "nearest_gene" in idx else "",
                )
            )
    return leads


def merge_intervals(items):
    """items: list of dict(chrom, start, end, lead...) -> merged regions"""
    out = []
    by_chr = {}
    for it in items:
        by_chr.setdefault(it["chrom"], []).append(it)
    for chrom in sorted(by_chr):
        arr = sorted(by_chr[chrom], key=lambda d: d["start"])
        cur = None
        for it in arr:
            if cur is None:
                cur = dict(chrom=chrom, start=it["start"], end=it["end"], leads=[it])
            elif it["start"] <= cur["end"]:
                cur["end"] = max(cur["end"], it["end"])
                cur["leads"].append(it)
            else:
                out.append(cur)
                cur = dict(chrom=chrom, start=it["start"], end=it["end"], leads=[it])
        if cur is not None:
            out.append(cur)
    out.sort(key=lambda d: (d["chrom"], d["start"]))
    return out


def main():
    leads = load_leads(IN_TSV)
    for l in leads:
        l["start"] = max(1, l["pos"] - WINDOW)
        l["end"] = l["pos"] + WINDOW

    regions = merge_intervals(leads)

    with open(OUT_TSV, "w") as out:
        out.write("region_id\tchrom\tstart\tend\tn_lead\tlead_snps\ttop_pval\tnearest_genes\n")
        for i, r in enumerate(regions, 1):
            rid = "R%02d" % i
            r["region_id"] = rid
            leads_sorted = sorted(r["leads"], key=lambda d: d["pval"])
            genes = ",".join(sorted({d["gene"] for d in r["leads"] if d["gene"] and d["gene"] != "NA"}))
            out.write(
                "\t".join(
                    [
                        rid,
                        str(r["chrom"]),
                        str(r["start"]),
                        str(r["end"]),
                        str(len(r["leads"])),
                        ",".join(d["rsid"] for d in leads_sorted),
                        "%.3g" % leads_sorted[0]["pval"],
                        genes if genes else "NA",
                    ]
                )
                + "\n"
            )

    total_bp = sum(r["end"] - r["start"] + 1 for r in regions)
    print("输入先导 SNP 数: %d" % len(leads))
    print("合并后区域数: %d  (窗口 ±%d bp)" % (len(regions), WINDOW))
    print("覆盖基因组总长: %.2f Mb" % (total_bp / 1e6))
    print("输出: %s" % OUT_TSV)


if __name__ == "__main__":
    main()
