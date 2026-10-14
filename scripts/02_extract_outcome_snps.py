#!/usr/bin/env python3
"""
02_extract_outcome_snps.py —— 从 FinnGen 汇总统计（.gz）中流式抽取目标 rsID 的行。

FinnGen R12 列（tab 分隔）:
  #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls

用法:
  python3 02_extract_outcome_snps.py <finngen.gz> <rsids.txt> <out.tsv>

说明:
  - rsids 列可能是 "rs123,rs456" 形式，这里逐个匹配
  - 输出与原文件同格式（含表头），便于 R 端 fread 直接读
  - 先扫一遍 gzip 只保留命中的行，内存占用恒定
"""
import gzip
import sys


def load_rsids(path: str) -> set:
    out = set()
    with open(path) as f:
        for line in f:
            for tok in line.replace(",", " ").split():
                tok = tok.strip()
                if tok.startswith("rs"):
                    out.add(tok)
    return out


def main() -> int:
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    src, rsfile, out = sys.argv[1], sys.argv[2], sys.argv[3]
    targets = load_rsids(rsfile)
    print(f"目标 rsID: {len(targets)} 个")

    n_line = n_hit = 0
    header = None
    idx = {}

    with gzip.open(src, "rt") as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith("#chrom") or (header is None and "pval" in line):
                if header is None:
                    header = line.rstrip("\n")
                    cols = header.lstrip("#").split("\t")
                    idx = {c: i for i, c in enumerate(cols)}
                    fout.write(header + "\n")
                continue
            n_line += 1
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            rs_field = parts[idx.get("rsids", 4)]
            if not rs_field or rs_field == "NA":
                continue
            hit = False
            for r in rs_field.split(","):
                if r in targets:
                    hit = True
                    break
            if hit:
                n_hit += 1
                fout.write(line)
            if n_line % 5_000_000 == 0:
                print(f"  已扫描 {n_line/1e6:.0f}M 行，命中 {n_hit}")

    print(f"完成: 扫描 {n_line} 行，命中 {n_hit} 行 → {out}")
    if n_hit == 0:
        print("⚠️  0 命中：检查 rsID 与 FinnGen 版本是否一致，或该文件是否用 chr:pos:ref:alt 代替 rsid")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
