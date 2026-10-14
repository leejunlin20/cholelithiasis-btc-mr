#!/usr/bin/env python3
"""
从 BBJ.BtC 全量 summary stats 中抽取工具变量并构建 MR 结局表
=========================================================

为什么用文件而不是 API
----------------------
PheWeb 的 `/api/variant/{chrom}-{pos}-{ref}-{alt}` 需要猜 BBJ 侧的 ref/alt 顺序，
且响应体不含逐表型等位基因信息，实测在 ABCG8 上取到过 af≈1 的假记录。全量文件
里 `Allele1`=ref、`Allele2`=alt、`BETA` 以 Allele2 为参考，可直接按 `SNPID`
（含 rsID）精确匹配，无任何猜测。

BBJ 文件列
----------
v  CHR  POS  SNPID  Allele1  Allele2  AC_Allele2  AF_Allele2  imputationInfo  N
BETA  SE  Tstat  p.value  p.value.NA  Is.SPA.converge  varT  varTstar  AF.Cases  AF.Controls
（Allele1 = ref, Allele2 = alt, BETA 以 Allele2 计；N = 159,619）

输出
----
results/bbj_replication/bbj_btc_instruments.tsv
  rsid, chrom, pos_grch37, ref, alt, bbj_beta, bbj_se, bbj_p, bbj_af,
  effect_allele (对齐到 FinnGen 侧 alt), beta_aligned, se_aligned, sign, imputation_info
"""
import gzip
import io
import os
import sys
import zipfile

import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "results", "bbj_replication")
ZIP = os.path.join(BASE, "data/bbj/hum0197.v3.BBJ.BtC.v1.zip")
INNER = ("hum0197.v3.BBJ.BtC.v1/"
         "GWASsummary_BtC_Japanese_SakaueKanai2020.auto.txt.gz")
MIN_INFO = 0.3


def main():
    iv = pd.read_csv(os.path.join(BASE, "results/medi/chole_iv_clumped.tsv"), sep="\t")
    # FinnGen R12 为 GRCh38；BBJ 为 GRCh37，故按 rsID 匹配，不按位置
    eff = {r.rsid: dict(chr38=r.chr, pos38=r.pos, ref38=r.ref, alt38=r.alt,
                        beta_exp=r.beta, se_exp=r.se, p_exp=r.pval, eaf=r.eaf,
                        gene=r.nearest_gene)
           for r in iv.itertuples()}
    print(f"暴露工具变量 {len(eff)} 个（FinnGen R12, GRCh38）")

    keep = set(eff)
    rows = []
    n_line = 0
    with zipfile.ZipFile(ZIP) as z, z.open(INNER) as fh, \
            gzip.open(io.BufferedReader(fh), "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            n_line += 1
            c = line.rstrip("\n").split("\t")
            rs = c[idx["SNPID"]]
            if rs not in keep:
                continue
            a1, a2 = c[idx["Allele1"]], c[idx["Allele2"]]
            rows.append(dict(
                rsid=rs, chrom=c[idx["CHR"]], pos_grch37=int(c[idx["POS"]]),
                bbj_ref=a1, bbj_alt=a2,
                bbj_beta=float(c[idx["BETA"]]), bbj_se=float(c[idx["SE"]]),
                bbj_p=float(c[idx["p.value"]]), bbj_af=float(c[idx["AF_Allele2"]]),
                n=int(c[idx["N"]]), imputation_info=float(c[idx["imputationInfo"]]),
                chisq=float(c[idx["Tstat"]]) ** 2))
            if len(rows) == len(keep):
                break

    print(f"扫描 {n_line:,} 行，命中 {len(rows)} 个工具变量")
    d = pd.DataFrame(rows)
    d["chrom_grch38"] = d.rsid.map(lambda r: eff[r]["chr38"])
    d["pos_grch38"] = d.rsid.map(lambda r: eff[r]["pos38"])
    d["ref_grch38"] = d.rsid.map(lambda r: eff[r]["ref38"])
    d["effect_allele"] = d.rsid.map(lambda r: eff[r]["alt38"])
    d["beta_exposure"] = d.rsid.map(lambda r: eff[r]["beta_exp"])
    d["se_exposure"] = d.rsid.map(lambda r: eff[r]["se_exp"])
    d["p_exposure"] = d.rsid.map(lambda r: eff[r]["p_exp"])
    d["eaf_exposure"] = d.rsid.map(lambda r: eff[r]["eaf"])
    d["nearest_gene"] = d.rsid.map(lambda r: eff[r]["gene"])

    # 方向对齐：BBJ beta 以 bbj_alt 计；若 bbj_alt != 效应等位基因则反号
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    d["sign"] = [
        1 if (a == e or comp.get(a) == e) else -1
        for a, e in zip(d.bbj_alt, d.effect_allele)]
    d["beta_aligned"] = d.bbj_beta * d.sign
    d["allele_mismatch"] = [not (a == e or comp.get(a) == e)
                            for a, e in zip(d.bbj_alt, d.effect_allele)]

    print("\n== 等位基因对齐诊断 ==")
    print(d.groupby(["bbj_ref", "bbj_alt", "effect_allele", "sign"]).size()
          .rename("n").reset_index().to_string(index=False))
    bad = d[d.allele_mismatch]
    if len(bad):
        print(f"\n!! {len(bad)} 个位点效应等位基因与 BBJ 等位基因不一致（可能链向或位点错配）：")
        print(bad[["rsid", "bbj_ref", "bbj_alt", "ref_grch38", "effect_allele"]].to_string(index=False))
    lowq = d[d.imputation_info < MIN_INFO]
    print(f"\nimputationInfo < {MIN_INFO} 的位点: {len(lowq)} 个（建议剔除）")

    d.to_csv(os.path.join(OUT, "bbj_btc_instruments.tsv"), sep="\t", index=False)
    print(f"\n已写出 {os.path.join(OUT, 'bbj_btc_instruments.tsv')}  ({len(d)} 行)")
    print(f"可用于 MR 的位点（info>={MIN_INFO} 且等位基因匹配）: "
          f"{len(d[(d.imputation_info >= MIN_INFO) & (~d.allele_mismatch)])}")


if __name__ == "__main__":
    main()
