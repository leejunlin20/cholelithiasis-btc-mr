#!/usr/bin/env python3
"""
30_mvmr_extract.py —— MVMR 输入抽取：从 FinnGen 全量 GWAS 中一趟流式产出两类子集

为什么要一趟两出
--------------
MVMR 的识别前提是「所有暴露共享同一套工具变量」（取各暴露显著位点的并集）。
因此每份 GWAS 需要读两次逻辑：
  ① 抽自己的全基因组显著位点  → 构造各自的 IV 集，再 clump、取并集
  ② 在「并集观察名单」上取效应 → 补齐非自身 IV 的效应（缺失则 MVMR 无法识别）
两件事必须扫同一份大文件，合成一趟比读两遍省一半 I/O。

FinnGen R12 列（制表符分隔，首列名带 #）：
  #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt [af_alt_cases af_alt_controls]
注意：连续变量端点（如 BMI_IRN）无末两列，故长度判断不能用固定列数。

用法:
  python3 30_mvmr_extract.py iv     <in.gz> <out_iv.tsv>        [p_thr] [min_maf]
  python3 30_mvmr_extract.py watch  <in.gz> <watch.txt> <out_watch.tsv>
"""
from __future__ import annotations

import gzip
import sys

IV_HEADER = "rsid\tchr\tpos\tref\talt\tbeta\tse\tpval\teaf\tnearest_gene\n"
WATCH_HEADER = ("rsid\tchr\tpos\tref\talt\tbeta\tse\tpval\teaf\t"
                "af_cases\taf_controls\tnearest_gene\n")

NUM_KEYS = ("pos", "pval", "beta", "sebeta", "af_alt")


def parse(c: list[str], ix: dict):
    """从已切分的列里取数并做健壮性校验；不合法返回 None。

    注意 POS 在 FinnGen R12 里可能带小数（如 43928847.0），故先 float 再 int。
    """
    try:
        vals = {k: float(c[ix[k]]) for k in NUM_KEYS}
    except (ValueError, IndexError, KeyError):
        return None
    rsid = c[ix["rsids"]].split(",")[0].strip()
    if not rsid.startswith("rs"):
        return None
    return dict(rsid=rsid, chrom=c[ix["chrom"]], pos=str(int(vals["pos"])),
                ref=c[ix["ref"]], alt=c[ix["alt"]],
                beta=vals["beta"], se=vals["sebeta"], p=vals["pval"],
                af=vals["af_alt"], gene=c[ix["nearest_genes"]],
                afc=c[ix["af_alt_cases"]] if "af_alt_cases" in ix else "",
                afn=c[ix["af_alt_controls"]] if "af_alt_controls" in ix else "")


def header_index(fh) -> dict:
    cols = fh.readline().rstrip("\n").lstrip("#").split("\t")
    return {c: i for i, c in enumerate(cols)}


def mode_iv(in_gz: str, out_tsv: str, p_thr: str = "5e-8", min_maf: str = "0.01"):
    p_thr, min_maf = float(p_thr), float(min_maf)
    n_total = n_hit = n_out = 0
    with gzip.open(in_gz, "rt") as fh, open(out_tsv, "w") as out:
        ix = header_index(fh)
        out.write(IV_HEADER)
        for line in fh:
            n_total += 1
            cols = line.rstrip("\n").split("\t")
            try:
                p = float(cols[ix["pval"]])
            except (ValueError, IndexError):
                continue
            if p >= p_thr:
                continue
            n_hit += 1
            d = parse(cols, ix)
            if d is None:
                continue
            if min(d["af"], 1 - d["af"]) < min_maf or d["se"] <= 0:
                continue
            out.write("\t".join([d["rsid"], d["chrom"], d["pos"], d["ref"], d["alt"],
                                 f"{d['beta']:.6g}", f"{d['se']:.6g}", f"{d['p']:.3g}",
                                 f"{d['af']:.4g}", d["gene"]]) + "\n")
            n_out += 1
    print(f"[iv] {in_gz}")
    print(f"     总行 {n_total} | P<{p_thr} 命中 {n_hit} | MAF>={min_maf} 写出 {n_out}")
    print(f"     输出 {out_tsv}")


def mode_watch(in_gz: str, watch_f: str, out_tsv: str):
    watch = {w.strip() for w in open(watch_f) if w.strip()}
    print(f"[watch] {in_gz}  观察名单 {len(watch)} 个 rsID")
    n_total = n_hit = 0
    seen: set[str] = set()
    with gzip.open(in_gz, "rt") as fh, open(out_tsv, "w") as out:
        ix = header_index(fh)
        j = ix["rsids"]
        out.write(WATCH_HEADER)
        for line in fh:
            n_total += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= j:
                continue
            # 先用廉价的字符串比较筛掉绝大多数行，再做浮点解析
            rid = cols[j].split(",")[0].strip()
            if rid not in watch or rid in seen:
                continue
            d = parse(cols, ix)
            if d is None or d["rsid"] not in watch:
                continue
            seen.add(d["rsid"])
            n_hit += 1
            out.write("\t".join([d["rsid"], d["chrom"], d["pos"], d["ref"], d["alt"],
                                 f"{d['beta']:.6g}", f"{d['se']:.6g}", f"{d['p']:.3g}",
                                 f"{d['af']:.4g}", d["afc"], d["afn"],
                                 d["gene"]]) + "\n")
    miss = sorted(watch - seen)
    print(f"     总行 {n_total} | 名单命中 {n_hit} | 缺席 {len(miss)}")
    if miss:
        print(f"     缺席前 20: {', '.join(miss[:20])}")
    print(f"     输出 {out_tsv}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "iv":
        mode_iv(*sys.argv[2:])
    elif sys.argv[1] == "watch":
        mode_watch(*sys.argv[2:])
    else:
        sys.exit(__doc__)
