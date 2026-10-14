#!/usr/bin/env python3
"""
15_fetch_ld_regions.py —— 为 coloc.susie 抓取 FinnGen R12 官方同源 LD 矩阵。

对每个目标区域：
  1. 从已抽取的 CHOLE / BTC 区域 summary stats 中取 ±HALFWIN 窗口内的 SNP
  2. 与 LD 矩阵的 variant 命名 (chr{chrom}_{pos}_{ref}_{alt}) 对齐，取双方交集
  3. 用 14_tabix_fetch.RemoteTabix 并发抓取该区间的 LD 记录
  4. 只保留「两端都在目标 SNP 集内」的配对 → 落盘 long 格式

产出（data/ld/）：
  {rid}_snps.tsv.gz   目标 SNP 及两性状效应量（R 端 coloc.susie 输入）
  {rid}_ld_long.tsv.gz  snp1 / snp2 / r / r2（稀疏，r²>0.01 才有记录）
  ld_fetch_log.tsv     每区域抓取统计
"""
from __future__ import annotations

import gzip
import os
import sys
import time
import traceback
import pandas as pd
import concurrent.futures as cf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module

_m = import_module("14_tabix_fetch")
RemoteTabix = _m.RemoteTabix

BASE = os.environ.get("MR_BASE", ".")
LD_URL = "https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/finngen_r12_chr{chrom}_ld.tsv.gz"
OUT = os.path.join(BASE, "data/ld")

HALFWIN = 250_000          # ±250 kb，与 coloc.abf 窗宽敏感性分析一致
MAF_MIN = 0.01             # 常见变异阈值（LD 矩阵本身只含 r²>0.01 的配对）
THREADS = 48

# (region_id, chrom, center, 注释)
REGIONS = [
    ("R69", "22", 43928850, "PNPLA3"),
    ("R51", "14", 94378610, "SERPINA1"),
    ("R67", "20", 44413724, "HNF4A"),
    ("R07",  "2", 44500425, "ABCG5/ABCG8 (negative control)"),
    ("R54", "16", 53763996, "FTO (positive control)"),
]


def load_trait(name: str) -> pd.DataFrame:
    d = pd.read_csv(os.path.join(BASE, "data/coloc", f"{name}_regions.tsv"), sep="\t")
    d = d.rename(columns={"#chrom": "chrom"})
    return d


def vname(chrom, pos, ref, alt) -> str:
    return f"chr{int(chrom)}_{int(pos)}_{ref}_{alt}"


def build_region_targets(rid, chrom, ctr, chole, btc):
    """返回 (targets_df, vpos_map)。targets_df 含两性状效应量与 SNP 名。"""
    lo, hi = ctr - HALFWIN, ctr + HALFWIN
    c = chole[(chole.region_id == rid) & chole.pos.between(lo, hi)].copy()
    b = btc[(btc.region_id == rid) & btc.pos.between(lo, hi)].copy()

    c["snp"] = [vname(r.chrom, r.pos, r.ref, r.alt) for r in c.itertuples()]
    b["snp"] = [vname(r.chrom, r.pos, r.ref, r.alt) for r in b.itertuples()]
    # 位置键（备选匹配）
    c["pkey"] = c.chrom.astype(str) + "_" + c.pos.astype(str)
    b["pkey"] = b.chrom.astype(str) + "_" + b.pos.astype(str)

    cb = c.set_index("snp")
    bb = b.set_index("snp")
    common = cb.index.intersection(bb.index)
    exact_rate = len(common) / max(1, len(cb.index))

    if len(common) < 50:
        # fallback：按 染色体_位置 匹配（等位基因命名方向不一致时）
        cb2 = c.set_index("pkey")
        bb2 = b.set_index("pkey")
        pk = cb2.index.intersection(bb2.index)
        t = pd.DataFrame({
            "snp": cb2.loc[pk, "snp"], "chrom": cb2.loc[pk, "chrom"], "pos": cb2.loc[pk, "pos"],
            "ref": cb2.loc[pk, "ref"], "alt": cb2.loc[pk, "alt"],
            "beta_1": cb2.loc[pk, "beta"], "se_1": cb2.loc[pk, "sebeta"],
            "p_1": cb2.loc[pk, "pval"], "af": cb2.loc[pk, "af_alt"],
            "beta_2": bb2.loc[pk, "beta"], "se_2": bb2.loc[pk, "sebeta"],
            "p_2": bb2.loc[pk, "pval"],
        })
        match_mode = "pos"
    else:
        t = pd.DataFrame({
            "snp": common, "chrom": cb.loc[common, "chrom"].values, "pos": cb.loc[common, "pos"].values,
            "ref": cb.loc[common, "ref"].values, "alt": cb.loc[common, "alt"].values,
            "beta_1": cb.loc[common, "beta"].values, "se_1": cb.loc[common, "sebeta"].values,
            "p_1": cb.loc[common, "pval"].values, "af": cb.loc[common, "af_alt"].values,
            "beta_2": bb.loc[common, "beta"].values, "se_2": bb.loc[common, "sebeta"].values,
            "p_2": bb.loc[common, "pval"].values,
        })
        match_mode = "exact"

    t["maf"] = t.af.clip(0, 1).apply(lambda a: min(a, 1 - a))
    n_all = len(t)
    t = t[t.maf >= MAF_MIN].copy()
    t = t.sort_values("pos").reset_index(drop=True)
    t["idx"] = range(len(t))
    t.attrs["match_mode"] = match_mode
    t.attrs["exact_rate"] = exact_rate
    t.attrs["n_all"] = n_all
    return t


def fetch_one(rid, chrom, ctr, note):
    t_start = time.time()
    chole = load_trait("CHOLE")
    btc = load_trait("BTC")
    tg = build_region_targets(rid, chrom, ctr, chole, btc)
    if len(tg) < 20:
        return dict(region=rid, note=note, status="FAIL", detail=f"目标 SNP 仅 {len(tg)}")
    idx = {s: i for i, s in enumerate(tg.snp)}
    lo, hi = ctr - HALFWIN, ctr + HALFWIN

    url = LD_URL.format(chrom=chrom)
    tb = RemoteTabix(url, verbose=True)
    chunks = tb.chunks_for(chrom, lo, hi)
    nbytes = sum(b - a for a, b in chunks)

    pairs = 0
    nline = 0
    out_path = os.path.join(OUT, f"{rid}_ld_long.tsv.gz")
    with gzip.open(out_path, "wt") as fh:
        for line in tb.fetch(chrom, lo, hi, threads=THREADS):
            nline += 1
            # chr pos variant1 variant2 r r2
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            v1, v2, r, r2 = parts[2], parts[3], parts[4], parts[5]
            i1 = idx.get(v1)
            if i1 is None:
                continue
            i2 = idx.get(v2)
            if i2 is None:
                continue
            fh.write(f"{i1}\t{i2}\t{r}\t{r2}\n")
            pairs += 1

    tg.drop(columns=["idx"]).to_csv(
        os.path.join(OUT, f"{rid}_snps.tsv.gz"), sep="\t", index=False, compression="gzip")

    return dict(region=rid, note=note, status="OK", chrom=chrom, center=ctr,
                match_mode=tg.attrs["match_mode"],
                exact_rate=round(tg.attrs["exact_rate"], 3),
                n_snps_all=tg.attrs["n_all"], n_snps=len(tg),
                ld_bytes_mb=round(nbytes / 1e6, 1), ld_lines=nline, ld_pairs=pairs,
                secs=round(time.time() - t_start, 1))


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"=== 抓取 {len(REGIONS)} 个区域 LD（±{HALFWIN//1000}kb, MAF>={MAF_MIN}）===", flush=True)
    results = []
    with cf.ProcessPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(fetch_one, *r): r[0] for r in REGIONS}
        for f in cf.as_completed(futs):
            rid = futs[f]
            try:
                res = f.result()
            except Exception as e:
                traceback.print_exc()
                res = dict(region=rid, status="ERROR", detail=f"{type(e).__name__}: {e}")
            results.append(res)
            print(f"[{res['status']}] {res}", flush=True)

    df = pd.DataFrame([{k: v for k, v in r.items()} for r in results])
    df = df.sort_values("region")
    df.to_csv(os.path.join(OUT, "ld_fetch_log.tsv"), sep="\t", index=False)
    print("\n=== 汇总 ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
