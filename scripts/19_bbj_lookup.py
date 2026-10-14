#!/usr/bin/env python3
"""
BBJ（PheWeb.jp）东亚位点查证 —— 精确版
=====================================

要点（踩过的坑）
----------------
1. **PheWeb 表型码不是下载页的缩写**。真实 phenocode：`Cholelithiasis`、
   `BtC`、`Urolithiasis`、`GB_polyp`、`PaC`、`Cirrhosis`、`Cholecystitis`、
   `ALP/ALT/AST/GGT/TBil`、`BMI/HDLC/LDLC/TC/TG/HbA1c/Glucose/Ht/SBP/DBP/MAP`。
   用 `CL`/`Uro`/`GBP` 查会静默返回空集。
2. **坐标全为 GRCh37**，本项目 FinnGen R12 为 GRCh38 ⇒ 必须 liftover。
   用 Ensembl GRCh37 REST 批量映射（POST /variation/human，≤200 id/次）。
3. **等位基因顺序不能猜**。Ensembl 只给 `allele_string`；必须对全部排列试
   `chr-pos-ref-alt`，并用响应体自带的 **`rsids` 字段** 回验位点身份——否则
   会取到同位置的另一个罕见变异（实测 rs4299376 曾取到 af=0.9996 的假记录）。
4. 方向对齐：BBJ `beta` 以响应体的 `alt` 为参考；比对目标效应等位基因后定符号。
"""
import itertools
import json
import os
import sys
import time
import urllib.request

import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "results", "bbj_replication")
os.makedirs(OUT, exist_ok=True)
UA = "curl/8.4.0"

# BBJ/Sakaue-Kanai 2021 真实 phenocode
TARGET = {
    "Cholelithiasis": "胆石症(暴露)",
    "BtC": "肝内胆管癌(主结局)",
    "Urolithiasis": "尿路结石(阴性对照)",
    "PaC": "胰腺癌(阴性对照)",
    "GB_polyp": "胆囊息肉",
    "Cholecystitis": "胆囊炎",
    "Cirrhosis": "肝硬化",
    "ALP": "ALP", "ALT": "ALT", "AST": "AST", "GGT": "GGT", "TBil": "TBil",
    "BMI": "BMI", "HDLC": "HDL-C", "LDLC": "LDL-C", "TC": "TC", "TG": "TG",
    "HbA1c": "HbA1c", "Glucose": "血糖", "Ht": "红细胞压积",
    "SBP": "收缩压", "DBP": "舒张压", "MAP": "平均动脉压",
}
CIT = ("SakaueKanai2021", "Sato2023", "Kanai2021")


def _open(url, data=None, timeout=60, headers=None):
    h = {"User-Agent": UA}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, data=data, headers=h)
    return urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req, timeout=timeout)


def grch37_map(rsids):
    out = {}
    for i in range(0, len(rsids), 180):
        chunk = rsids[i:i + 180]
        raw = _open("https://grch37.rest.ensembl.org/variation/human",
                    data=json.dumps({"ids": chunk}).encode(),
                    headers={"Content-Type": "application/json",
                             "Accept": "application/json"}).read()
        for rs, v in json.loads(raw).items():
            if not isinstance(v, dict) or not v.get("mappings"):
                continue
            m = v["mappings"][0]
            alle = [a for a in m["allele_string"].split("/") if a not in ("-", "*")]
            out[rs] = dict(chrom=str(m["seq_region_name"]), pos=int(m["start"]), alleles=alle)
        time.sleep(0.4)
    return out


def fetch_variant(chrom, pos, ref, alt, tries=3):
    url = f"https://pheweb.jp/api/variant/{chrom}-{pos}-{ref}-{alt}"
    for k in range(tries):
        try:
            d = json.loads(_open(url, timeout=60).read())
            return d if isinstance(d, dict) and d.get("phenos") else None
        except Exception:                                        # noqa: BLE001
            time.sleep(1 + k)
    return None


def lookup(rsid, chrom, pos, alleles, effect_allele):
    """返回 (records, meta)；用响应体 rsids 回验位点身份。"""
    for ref, alt in itertools.permutations(alleles, 2):
        d = fetch_variant(chrom, pos, ref, alt)
        if not d:
            continue
        got_rs = d.get("rsids") or []
        if isinstance(got_rs, str):
            got_rs = [got_rs]
        if rsid not in got_rs:
            continue                                   # 同位置的其它变异，丢弃
        recs = {}
        bbj_alt = d["alt"]                    # ★ 等位基因在响应体顶层，不在 pheno 记录里
        for p in d["phenos"]:
            if p.get("citation") not in CIT:
                continue
            pc = p.get("phenocode")
            if pc in TARGET:
                sign = 1
                if effect_allele:
                    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
                    sign = 1 if (bbj_alt == effect_allele
                                 or comp.get(bbj_alt) == effect_allele) else -1
                recs[pc] = dict(ref=d["ref"], alt=bbj_alt, gene=d.get("nearest_genes"),
                                beta=float(p["beta"]), se=float(p["sebeta"]),
                                p=float(p["pval"]), af=float(p["af"]),
                                n_case=p.get("num_cases"), n_ctrl=p.get("num_controls"),
                                beta_aligned=float(p["beta"]) * sign, sign=sign,
                                citation=p.get("citation"))
        if recs:
            return recs, dict(ref=d["ref"], alt=d["alt"], rsids=got_rs)
    return {}, {}


def load_snps(mode):
    if mode == "leads":
        # 手工整理：rsID / hg38(chr,pos) / FinnGen 侧 ref/alt（= 效应等位基因）
        # 坐标与等位基因均取自 FinnGen R12 原始文件，非推测
        CUR = [
            ("R69 PNPLA3",   "rs738409",   "22", 43928847, "C", "G"),
            ("R69 PNPLA3",   "rs738408",   "22", 43928850, "C", "T"),
            ("R51 SERPINA1", "rs28929474", "14", 94378610, "C", "T"),
            ("R67 HNF4A",    "rs1800961",  "20", 44413724, "C", "T"),
            ("R07 ABCG8",    "rs4299376",  "2",  43845437, "G", "T"),
            ("R07 ABCG8",    "rs4245791",  "2",  43847292, "C", "T"),
            ("R07 ABCG8",    "rs11887534", "2",  43839108, "G", "C"),
            ("R54 FTO",      "rs7206790",  "16", 53763996, "C", "G"),
        ]
        return [dict(label=lab, snp=rs, ref=ref, alt=alt, chrom=c, pos=p,
                     beta1=None, p1=None, beta2=None, p2=None)
                for lab, rs, c, p, ref, alt in CUR]
    iv = pd.read_csv(os.path.join(BASE, "results/medi/chole_iv_clumped.tsv"), sep="\t")
    return [dict(label=r.nearest_gene or "", snp=r.rsid, ref=r.ref, alt=r.alt,
                 chrom=str(r.chr), pos=int(r.pos), beta1=r.beta, p1=r.pval, beta2=None,
                 p2=None) for r in iv.itertuples()]


def main(mode):
    snps = load_snps(mode)
    # 用 rsID 映射；leads 模式下 snp 字段是位置名，需要 dbSNP 逐步解析
    rsids = []
    for s in snps:
        if s["snp"].startswith("rs"):
            rsids.append(s["snp"])
    gmap = grch37_map(rsids) if rsids else {}
    rows = []
    for s in snps:
        rs = s["snp"] if s["snp"].startswith("rs") else None
        if rs is None:
            continue
        m = gmap.get(rs)
        if not m:
            print(f"  [无映射] {rs}"); continue
        recs, meta = lookup(rs, m["chrom"], m["pos"], m["alleles"], s["alt"])
        if not recs:
            print(f"  [无命中] {rs} hg19={m['chrom']}:{m['pos']} alleles={m['alleles']}")
            continue
        print(f"  {rs:<14} hg19={m['chrom']}:{m['pos']}  BBJ ref/alt={meta['ref']}/{meta['alt']}"
              f"  效应等位={s['alt']}  表型 {len(recs)} 个")
        for pc, r in recs.items():
            rows.append(dict(label=s["label"], rsid=rs, hg19=f"{m['chrom']}:{m['pos']}",
                             hg38=f"{s['chrom']}:{s['pos']}",
                             ref_grch38=s["ref"], alt_grch38=s["alt"],
                             bbj_ref=r["ref"], bbj_alt=r["alt"], sign=r["sign"],
                             pheno=pc, pheno_cn=TARGET[pc], citation=r["citation"],
                             gene=r["gene"], beta_bbj=r["beta"], beta_aligned=r["beta_aligned"],
                             se_bbj=r["se"], p_bbj=r["p"], af_bbj=r["af"],
                             n_case=r["n_case"], n_ctrl=r["n_ctrl"],
                             beta_finngen=s["beta1"], p_finngen=s["p1"],
                             beta_finngen_btc=s["beta2"], p_finngen_btc=s["p2"]))
        time.sleep(0.3)
    return rows


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "leads"
    rows = main(mode)
    if rows:
        path = os.path.join(OUT, f"bbj_lookup_{mode}.tsv")
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        print(f"\n已写出 {path}  ({len(rows)} 行)")
    else:
        print("无结果")
