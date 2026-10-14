#!/usr/bin/env python3
"""
24_split_regions.py
-------------------
修正 10_define_coloc_regions.py 的「过度合并」缺陷。

问题
----
10_*.py 用「先导 SNP ±500 kb，窗口相交即合并」定义区域。当两个物理上独立的
位点相距 600~900 kb 时，其 ±500 kb 窗口首尾相触，会被合并成一个跨越 2 Mb+
的「区域」。这类区域违反 COLOC / SuSiE 的「单一因果变异」前提，会导致：
  · SuSiE 报出大量伪「独立信号」（R07 报出 19 个）
  · 条件化后残余 z 值极大（R07 max|z|_cond = 20.2，其余区域仅 3.3~5.7）

修正规则
--------
合并后，若同一区域内相邻先导 SNP 的间距 > MAX_GAP（默认 500 kb），
则在该最大间隙处切开，迭代至无 >MAX_GAP 的间隙为止。

输出
----
results/coloc/regions_v2.tsv  （原 regions.tsv 保留不动，便于对照）
"""
import os

BASE = os.environ.get("MR_BASE", ".")
IN_TSV = os.path.join(BASE, "results/coloc/regions.tsv")
OUT_TSV = os.path.join(BASE, "results/coloc/regions_v2.tsv")

MAX_GAP = 500_000
WINDOW = 500_000


def load_leads(path):
    leads = []
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            leads.append(dict(rsid=p[ix["rsid"]], chrom=int(p[ix["chr"]]),
                              pos=int(p[ix["pos"]]), pval=float(p[ix["pval"]]),
                              gene=p[ix["nearest_gene"]]))
    return leads


def main():
    leads = load_leads(os.path.join(BASE, "results/medi/chole_iv_clumped.tsv"))
    by_rs = {l["rsid"]: l for l in leads}

    regions = []
    with open(IN_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            regions.append(dict(rid=p[ix["region_id"]], chrom=int(p[ix["chrom"]]),
                                start=int(p[ix["start"]]), end=int(p[ix["end"]]),
                                rsids=p[ix["lead_snps"]].split(",")))

    out_rows, split_log = [], []
    for r in regions:
        ls = sorted([by_rs[rs] for rs in r["rsids"] if rs in by_rs],
                    key=lambda d: d["pos"])
        if not ls:
            continue
        # 按 >MAX_GAP 的间隙切块
        chunks, cur = [], [ls[0]]
        for a, b in zip(ls, ls[1:]):
            if b["pos"] - a["pos"] > MAX_GAP:
                chunks.append(cur)
                cur = [b]
            else:
                cur.append(b)
        chunks.append(cur)

        # 切点处的边界：取两先导 SNP 的中点，避免相邻区域窗口重叠（同一 SNP 被计两次）
        bounds = []
        for k, ch in enumerate(chunks):
            start = max(1, ch[0]["pos"] - WINDOW)
            end = ch[-1]["pos"] + WINDOW
            if k > 0:
                gap_mid = (chunks[k - 1][-1]["pos"] + ch[0]["pos"]) // 2
                start = gap_mid + 1
            if k < len(chunks) - 1:
                gap_mid = (ch[-1]["pos"] + chunks[k + 1][0]["pos"]) // 2
                end = gap_mid
            bounds.append((start, end))

        for k, ch in enumerate(chunks):
            start, end = bounds[k]
            rid = r["rid"] if len(chunks) == 1 else "%s%c" % (r["rid"], ord("a") + k)
            genes = ",".join(sorted({d["gene"] for d in ch if d["gene"] and d["gene"] != "NA"}))
            top = min(ch, key=lambda d: d["pval"])
            out_rows.append(dict(rid=rid, chrom=r["chrom"], start=start, end=end,
                                 n_lead=len(ch),
                                 rsids=",".join(d["rsid"] for d in sorted(ch, key=lambda d: d["pval"])),
                                 top_p="%.3g" % top["pval"], genes=genes or "NA",
                                 span_mb=(end - start) / 1e6))
        if len(chunks) > 1:
            split_log.append((r["rid"], r["chrom"], r["start"], r["end"],
                              (r["end"] - r["start"]) / 1e6, " -> ".join(
                                  "%s%c[%.2f-%.2fMb]" % (r["rid"], ord("a") + i,
                                                         c[0]["pos"] / 1e6, c[-1]["pos"] / 1e6)
                                  for i, c in enumerate(chunks))))

    with open(OUT_TSV, "w") as out:
        out.write("region_id\tchrom\tstart\tend\tn_lead\tlead_snps\ttop_pval\tnearest_genes\tspan_mb\n")
        for r in out_rows:
            out.write("\t".join([r["rid"], str(r["chrom"]), str(r["start"]), str(r["end"]),
                                 str(r["n_lead"]), r["rsids"], r["top_p"], r["genes"],
                                 "%.2f" % r["span_mb"]]) + "\n")

    print("原区域数: %d  →  修正后区域数: %d" % (len(regions), len(out_rows)))
    print("被切分的区域: %d 个\n" % len(split_log))
    for s in split_log:
        print("  %s chr%d %d-%d (%.2fMb)  →  %s" % s)
    print("\n输出: %s" % OUT_TSV)
    print("修正后仍 >1.5Mb 的区域:")
    for r in out_rows:
        if r["span_mb"] > 1.5:
            print("  %s chr%d %d-%d (%.2fMb) %s" % (r["rid"], r["chrom"], r["start"],
                                                    r["end"], r["span_mb"], r["genes"]))


if __name__ == "__main__":
    main()
