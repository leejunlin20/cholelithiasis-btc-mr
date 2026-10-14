#!/usr/bin/env python3
"""
05_extract_iv_from_gwas.py —— 从 FinnGen 全量 GWAS .gz 中按 P 阈值抽取工具变量位点
（用于两步介导 MR 的中间表型：胆囊结石病 → 胆道癌）

为什么不用 R 读：文件解压后 ~2 GB，R 全读会爆内存。逐行流式过滤，只保留显著位点。

用法:
  python3 05_extract_iv_from_gwas.py <in.gz> <out.tsv> [p_threshold] [min_maf]
"""
import gzip
import io
import sys

COLS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp",
        "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]


def main():
    in_gz, out_tsv = sys.argv[1], sys.argv[2]
    p_thr = float(sys.argv[3]) if len(sys.argv) > 3 else 5e-8
    min_maf = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01

    n_total = n_hit = n_out = 0
    with gzip.open(in_gz, "rt") as fh, open(out_tsv, "w") as out:
        out.write("rsid\tchr\tpos\tref\talt\tbeta\tse\tpval\teaf\tnearest_gene\n")
        for line in fh:
            n_total += 1
            if n_total == 1 and line.startswith("#chrom"):
                continue
            parts = line.rstrip("\n").split("\t")
            # 连续变量端点（如 BMI_IRN）只有 11 列：无 af_alt_cases/af_alt_controls
            if len(parts) < 11:
                continue
            try:
                p = float(parts[6])
            except ValueError:
                continue
            if p >= p_thr:
                continue
            n_hit += 1
            try:
                af = float(parts[10])
                beta = float(parts[8])
                se = float(parts[9])
            except (ValueError, IndexError):
                continue
            maf = af if af <= 0.5 else 1 - af
            if maf < min_maf or se <= 0:
                continue
            # rsids 可能是逗号分隔的多个，取第一个
            rsid = parts[4].split(",")[0].strip()
            if not rsid.startswith("rs"):
                continue
            out.write("\t".join([rsid, parts[0], parts[1], parts[2], parts[3],
                                 f"{beta:.6g}", f"{se:.6g}", f"{p:.3g}", f"{af:.4g}",
                                 parts[5]]) + "\n")
            n_out += 1
            if n_out % 500 == 0:
                print(f"  已写出 {n_out} 个位点...", flush=True)

    print(f"总行数 {n_total}，P<{p_thr} 命中 {n_hit}，MAF>={min_maf} 后写出 {n_out}")
    print(f"输出: {out_tsv}")


if __name__ == "__main__":
    main()
