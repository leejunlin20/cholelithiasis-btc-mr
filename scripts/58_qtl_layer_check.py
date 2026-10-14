#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
58_qtl_layer_check.py  v2 —— QTL 层因果检验（路 C / 执行顺序第 3 项）

目的
----
回答现稿四个共定位区在 **cis-eQTL 层** 的问题：
  R07 (ABCG5 / ABCG8 / LRPPRC)   R51 (SERPINA1)   R67 (HNF4A)   R69 (PNPLA3)

数据源
------
  · eQTL：eQTL Catalogue · GTEx v8 分组织 cis-eQTL（远程 tabix，按区域取数）
      https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/<Tissue>.tsv.gz
      列 variant r2 pvalue mt_object mt_id maf gene_id median_tpm beta se an ac
         chromosome position ref alt type rsid
  · 结局：本地 FinnGen R12 区域切片
      data/coloc/CHOLE_regions.tsv (K11_CHOLELITH) · BTC_regions.tsv (C3_BILIARY_…)

方法（v2 关键修正）
------------------
  v1 直接用 top eQTL SNP 会失配——GTEx 基于 WGS 含大量**罕见变异**，FinnGen 插补
  面板里没有（实测 R69 top SNP chr22:44372372_G_A 在 FinnGen 中不存在）。
  v2 因此：
    1. 基因 cis 窗口 = TSS ±1 Mb
    2. 在窗口内对 eQTL p 值做 BH-FDR，以 **q<0.05** 定义 cis-eQTL 显著集
    3. 与 FinnGen 按 **chrom:pos + 等位基因兼容** 取交集
    4. 工具 = 交集中 eQTL p 最小者（cis-MR 标准单工具 Wald ratio）
       若交集内存在 ≥2 个距离 >500 kb 的独立信号，另做 IVW 敏感性（标注未做 LD 校正）
    5. 等位基因以 eQTL alt 为准，必要时翻转 GWAS beta
    6. 工具强度 F = (b_exp/se_exp)^2；R2 = F/(F+n-2)
    7. 功效：现稿口径 SE ≈ 1/sqrt(R2·n_eff)，MDE = exp(2.8016·SE)
"""
import os, sys, math, json, time, urllib.request
import numpy as np
import pandas as pd
import pysam

PROJ = os.environ.get("MR_BASE", ".")
CACHE = "/tmp/GTEx_ge_cache"
EBI = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge"
OUT = os.path.join(PROJ, "results", "qtl")
os.makedirs(os.path.join(OUT, "raw"), exist_ok=True)

COLS = ["variant", "r2", "pvalue", "mt_object", "mt_id", "maf", "gene_id",
        "median_tpm", "beta", "se", "an", "ac", "chrom", "position",
        "ref", "alt", "type", "rsid"]

GENES = {
    "ABCG5":    dict(ensg="ENSG00000138075", chrom="2",  tss=43806213, end=43838865),
    "ABCG8":    dict(ensg="ENSG00000143921", chrom="2",  tss=43831942, end=43882988),
    "LRPPRC":   dict(ensg="ENSG00000138095", chrom="2",  tss=43886224, end=43996226),
    "SERPINA1": dict(ensg="ENSG00000197249", chrom="14", tss=94376745, end=94391024),
    "HNF4A":    dict(ensg="ENSG00000101076", chrom="20", tss=44355439, end=44432845),
    "PNPLA3":   dict(ensg="ENSG00000100344", chrom="22", tss=43923761, end=43964488),
}
SYM = {v["ensg"]: k for k, v in GENES.items()}

REGIONS = {
    "R07": dict(chrom="2",  start=43306636, end=45693031,
                genes=["ABCG5", "ABCG8", "LRPPRC"],
                gwas_leads=["rs56266464", "rs78451356", "rs11124977",
                            "rs72799962", "rs62134177"]),
    "R51": dict(chrom="14", start=93878610, end=94878610,
                genes=["SERPINA1"], gwas_leads=["rs28929474"]),
    "R67": dict(chrom="20", start=43913724, end=44913724,
                genes=["HNF4A"], gwas_leads=["rs1800961"]),
    "R69": dict(chrom="22", start=43428850, end=44428850,
                genes=["PNPLA3"], gwas_leads=["rs738408"]),
}
ORGANS = ["Liver", "Whole_Blood"]
N_EFF = {"CHOLE": 44737, "BTC": 2284}
Z80 = 2.8016
FDR_Q = 0.05


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def load_ensg_map(ensgs):
    """ENSG -> symbol，Ensembl REST 批量查询并本地缓存"""
    f = os.path.join(OUT, "ensg_symbol_map.json")
    try:
        m = json.load(open(f))
    except Exception:
        m = {}
    need = sorted({e for e in ensgs if e and e not in m})
    if need:
        for i in range(0, len(need), 300):
            chunk = need[i:i + 300]
            req = urllib.request.Request(
                "https://rest.ensembl.org/lookup/id?content-type=application/json",
                data=json.dumps({"ids": chunk}).encode(),
                headers={"Content-Type": "application/json"})
            try:
                d = json.load(urllib.request.urlopen(req, timeout=90))
                for k, v in d.items():
                    if isinstance(v, dict):
                        m[k] = v.get("display_name")
            except Exception as ex:
                log(f"  !! Ensembl lookup 失败: {ex}")
        try:
            json.dump(m, open(f, "w"))
        except Exception:
            pass
    return m


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    ranked = p[o] * n / (np.arange(1, n + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[o] = np.clip(ranked, 0, 1)
    return q


# ---------------------------------------------------------------- eQTL 取数
def fetch_eqtl(tissue, chrom, start, end, region, force=False):
    cache_f = os.path.join(OUT, "raw", f"{region}_{tissue}.tsv.gz")
    if os.path.exists(cache_f) and not force:
        df = pd.read_csv(cache_f, sep="\t", dtype={"chrom": str, "ref": str, "alt": str})
        log(f"  {tissue} chr{chrom}:{start}-{end} -> {len(df)} 行 (本地缓存)")
        return df
    url = f"{EBI}/{tissue}.tsv.gz"
    idx = f"{CACHE}/GTEx_{tissue}.tbi"
    t0 = time.time()
    tb = pysam.TabixFile(url, index=idx if os.path.exists(idx) else None)
    rows = [line.rstrip("\n").split("\t")
            for line in tb.fetch(str(chrom), int(start) - 1, int(end))]
    tb.close()
    log(f"  {tissue} chr{chrom}:{start}-{end} -> {len(rows)} 行 ({time.time()-t0:.1f}s)")
    if not rows:
        return pd.DataFrame(columns=COLS)
    df = pd.DataFrame(rows, columns=COLS)
    for c in ("pvalue", "maf", "beta", "se", "an", "ac", "position"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df.to_csv(cache_f, sep="\t", index=False, compression="gzip")
    return df


def load_gwas():
    out = {}
    for key, fn in [("CHOLE", "data/coloc/CHOLE_regions.tsv"),
                    ("BTC", "data/coloc/BTC_regions.tsv")]:
        df = pd.read_csv(os.path.join(PROJ, fn), sep="\t", dtype={"#chrom": str})
        df = df[df.region_id.isin(REGIONS.keys())].rename(columns={"#chrom": "chrom"})
        df["pos"] = pd.to_numeric(df["pos"])
        for c in ("beta", "sebeta", "pval"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        out[key] = df
    return out


def gwas_index(df):
    d = {}
    for r in df.itertuples(index=False):
        d[(str(r.chrom), int(r.pos))] = (r.ref, r.alt, float(r.beta),
                                         float(r.sebeta), float(r.pval), r.rsids)
    return d


def two_sided_p(z):
    return math.erfc(abs(z) / math.sqrt(2))


def wald(b_exp, se_exp, b_out, se_out):
    b = b_out / b_exp
    se = se_out / abs(b_exp)
    return b, se, two_sided_p(b / se)


# ---------------------------------------------------------------- 主流程
def main():
    global REGIONS, ORGANS
    force = "--force" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) > 0:
        REGIONS = {k: v for k, v in REGIONS.items() if k in set(args[0].split(","))}
    if len(args) > 1:
        ORGANS = args[1].split(",")
    log(f"regions={list(REGIONS)} organs={ORGANS} force={force}")

    gw = load_gwas()
    gi = {k: gwas_index(v) for k, v in gw.items()}
    log(f"GWAS 行数: " + ", ".join(f"{k}={len(v)}" for k, v in gw.items()))

    mr_rows, lead_rows, land_rows, ivw_rows = [], [], [], []
    qmap, pmap = {}, {}

    for rid, rg in REGIONS.items():
        chrom, start, end = rg["chrom"], rg["start"], rg["end"]
        log(f"== 区域 {rid} chr{chrom}:{start}-{end} ==")
        for tissue in ORGANS:
            qmap, pmap = {}, {}          # tissue-scoped
            e = fetch_eqtl(tissue, chrom, start, end, rid, force=force)
            if e.empty:
                log(f"  !! {rid}/{tissue} 无 eQTL 记录")
                continue

            # ---- 区域 cis 基因全景 ----
            ee = e.copy()
            ee["q"] = np.nan
            for gid, sub in ee.groupby("gene_id"):
                ee.loc[sub.index, "q"] = bh_fdr(sub.pvalue.values)
            sig_all = ee[ee.q < FDR_Q]
            if len(sig_all):
                g = (sig_all.groupby("gene_id")
                     .agg(n_cis_q05=("variant", "size"), min_p=("pvalue", "min"),
                          min_q=("q", "min")).reset_index())
                g["symbol"] = g.gene_id.map(SYM)
                g["region"], g["tissue"] = rid, tissue
                g["n_samples"] = e.an.max() / 2
                land_rows.extend(g.to_dict("records"))

            # ---- 目标基因 ----
            for gname in rg["genes"]:
                info = GENES[gname]
                ensg = info["ensg"]
                sub = e[e.gene_id == ensg].copy()
                sub = sub[(sub.position >= info["tss"] - 1_000_000) &
                          (sub.position <= info["end"] + 1_000_000)]
                if sub.empty:
                    log(f"  [{tissue}] {gname}: 窗口内无 eQTL 记录")
                    continue
                sub["pos_i"] = sub.position.astype(int)
                sub["q"] = bh_fdr(sub.pvalue.values)
                # 距基因体的距离（基因体内 = 0）
                sub["dist"] = np.where(
                    (sub.pos_i >= info["tss"]) & (sub.pos_i <= info["end"]), 0,
                    np.minimum(abs(sub.pos_i - info["tss"]),
                               abs(sub.pos_i - info["end"])))
                nsig_q05 = int((sub.q < FDR_Q).sum())
                nsig_5e8 = int((sub.pvalue < 5e-8).sum())
                nsig_near = int(((sub.q < FDR_Q) & (sub.dist <= 100_000)).sum())
                min_p = float(sub.pvalue.min())

                n = float(sub.an.max()) / 2.0

                # ---- 与 FinnGen 交集（等位基因兼容）----
                rows = []
                for r in sub.itertuples(index=False):
                    hitC = gi["CHOLE"].get((str(chrom), int(r.pos_i)))
                    if hitC is None:
                        continue
                    refC, altC = hitC[0], hitC[1]
                    if (refC, altC) == (r.ref, r.alt):
                        flip = False
                    elif (refC, altC) == (r.alt, r.ref):
                        flip = True
                    else:
                        continue
                    rows.append((r, flip))
                base = dict(region=rid, gene=gname, tissue=tissue, n_samples=n,
                            n_cis_q05=nsig_q05, n_cis_p5e8=nsig_5e8,
                            n_cis_near_q05=nsig_near, min_eqtl_p=min_p,
                            n_shared=len(rows))
                if not rows:
                    mr_rows.append(dict(**base,
                                        status="无可匹配工具（GTEx 与 FinnGen 变异集无交集）"))
                    log(f"  [{tissue}] {gname}: q<0.05={nsig_q05}, 无共享变异 -> 不可估计")
                    continue

                # 分层选择工具：优先近端显著 → 远端显著 → 最优可得
                def pick(pred):
                    c = [t for t in rows if pred(t[0])]
                    return min(c, key=lambda t: t[0].pvalue) if c else None

                chosen, tier = None, ""
                for label, pred in (
                        ("proximal_q05(<=100kb)", lambda r: r.q < FDR_Q and r.dist <= 100_000),
                        ("cis_q05", lambda r: r.q < FDR_Q),
                        ("nominal_only", lambda r: True)):
                    chosen = pick(pred)
                    if chosen:
                        tier = label
                        break
                best, flip = chosen
                F = (best.beta / best.se) ** 2
                R2 = F / (F + n - 2)

                rec = dict(**base,
                           tool_tier=tier,
                           top_variant=best.variant, top_rsid=best.rsid,
                           top_pos=int(best.pos_i), top_ref=best.ref, top_alt=best.alt,
                           top_dist_kb=float(best.dist) / 1000.0,
                           eqtl_beta=float(best.beta), eqtl_se=float(best.se),
                           eqtl_p=float(best.pvalue), eqtl_q=float(best.q),
                           eqtl_maf=float(best.maf),
                           eqtl_median_tpm=float(best.median_tpm),
                           F=float(F), R2=float(R2), flip=flip,
                           status=("ok" if (best.q < FDR_Q and best.dist <= 100_000)
                                   else ("cis_q05_distal（远距离工具，LD 混杂风险）"
                                         if best.q < FDR_Q
                                         else "no_cis_eqtl_q05（工具弱，数值不可解释）")))
                for okey in ("CHOLE", "BTC"):
                    hit = gi[okey].get((str(chrom), int(best.pos_i)))
                    if hit is None:
                        rec[f"OR_{okey}"] = np.nan
                        rec[f"p_{okey}"] = np.nan
                        continue
                    ref_o, alt_o, b_o, se_o, p_o, rs_o = hit
                    if flip:
                        b_o = -b_o
                    b_mr, se_mr, p_mr = wald(best.beta, best.se, b_o, se_o)
                    rec[f"OR_{okey}"] = math.exp(b_mr)
                    rec[f"ORl_{okey}"] = math.exp(b_mr - 1.96 * se_mr)
                    rec[f"ORu_{okey}"] = math.exp(b_mr + 1.96 * se_mr)
                    rec[f"p_{okey}"] = p_mr
                    rec[f"b_{okey}"], rec[f"se_{okey}"] = b_mr, se_mr
                    rec[f"mde_{okey}"] = math.exp(Z80 / math.sqrt(max(R2, 1e-12) * N_EFF[okey]))
                    rec[f"gwas_rsid_{okey}"] = rs_o
                    rec[f"p_gwas_{okey}"] = p_o
                    rec[f"gwas_b_{okey}"] = b_o
                mr_rows.append(rec)
                log(f"  [{tissue}] {gname}: tier={tier} tool={best.rsid} "
                    f"@{int(best.pos_i)} dist={best.dist/1000:.0f}kb "
                    f"p_eqtl={best.pvalue:.2e} q={best.q:.2e} F={F:.0f} R2={R2*100:.2f}% "
                    f"| OR_CHOLE={rec.get('OR_CHOLE', float('nan')):.3f} "
                    f"(P={rec.get('p_CHOLE', float('nan')):.3g}) "
                    f"| OR_BTC={rec.get('OR_BTC', float('nan')):.3f} "
                    f"(P={rec.get('p_BTC', float('nan')):.3g}) | shared={len(rows)}")

                # ---- 敏感性：距离 prune 后的 IVW（限 q<0.05）----
                sigrows = [t for t in rows if t[0].q < FDR_Q]
                pick2 = []
                for r, fl in sorted(sigrows, key=lambda t: t[0].pvalue):
                    if all(abs(int(r.pos_i) - int(p[0].pos_i)) > 500_000 for p in pick2):
                        pick2.append((r, fl))
                if len(pick2) >= 2:
                    for okey in ("CHOLE", "BTC"):
                        bx, by, sx, sy = [], [], [], []
                        for r, fl in pick2:
                            hit = gi[okey].get((str(chrom), int(r.pos_i)))
                            if hit is None:
                                continue
                            b_o, se_o = hit[2], hit[3]
                            if fl:
                                b_o = -b_o
                            bx.append(r.beta); sx.append(r.se)
                            by.append(b_o); sy.append(se_o)
                        if len(bx) >= 2:
                            bx, by = np.array(bx), np.array(by)
                            w = bx ** 2 / np.array(sy) ** 2
                            b_ivw = float((w * (by / bx)).sum() / w.sum())
                            se_fe = float(1 / math.sqrt(w.sum()))
                            Q = float((w * (by / bx - b_ivw) ** 2).sum())
                            kappa = max(1.0, math.sqrt(Q / max(len(bx) - 1, 1)))
                            se_ivw = se_fe * kappa
                            ivw_rows.append(dict(region=rid, gene=gname, tissue=tissue,
                                                 outcome=okey, n_snp=len(bx),
                                                 OR=math.exp(b_ivw),
                                                 ORl=math.exp(b_ivw - 1.96 * se_ivw),
                                                 ORu=math.exp(b_ivw + 1.96 * se_ivw),
                                                 P=two_sided_p(b_ivw / se_ivw),
                                                 Q=Q, kappa=kappa,
                                                 note="距离 prune >500kb；未做 LD 校正，仅供参考"))

                # 供 lead SNP 判定用的基因 cis q 映射
                qmap[ensg] = dict(zip(sub.variant, sub.q))
                pmap[ensg] = dict(zip(sub.variant, sub.pvalue))

            # ---- GWAS lead SNP 的 cis-eQTL 状态 ----
            tgt_ensg = [GENES[g]["ensg"] for g in rg["genes"]]
            for lead in rg["gwas_leads"]:
                ll = e[e.rsid == lead]
                if ll.empty:
                    lead_rows.append(dict(region=rid, tissue=tissue, lead_snp=lead,
                                          gene=None, ensg=None, target_gene=False,
                                          eqtl_beta=np.nan, eqtl_p=np.nan,
                                          q_within_cis=np.nan,
                                          is_sig_p5e8=False, is_cis_eqtl_q05=False,
                                          note="lead SNP 不在 eQTL 面板"))
                    continue
                recs = []
                for r in ll.itertuples(index=False):
                    gq = qmap.get(r.gene_id, {}).get(r.variant, np.nan)
                    gp = pmap.get(r.gene_id, {}).get(r.variant, np.nan)
                    recs.append(dict(
                        region=rid, tissue=tissue, lead_snp=lead,
                        gene=SYM.get(r.gene_id, r.gene_id), ensg=r.gene_id,
                        target_gene=r.gene_id in tgt_ensg,
                        eqtl_beta=float(r.beta),
                        eqtl_p=float(r.pvalue if (gp is None or
                                                  (isinstance(gp, float) and math.isnan(gp)))
                                     else gp),
                        q_within_cis=float(gq) if gq is not None else np.nan,
                        is_sig_p5e8=bool(r.pvalue < 5e-8),
                        is_cis_eqtl_q05=bool(gq is not None and not (isinstance(gq, float)
                                                                     and math.isnan(gq))
                                             and gq < FDR_Q),
                        note=""))
                keep = [x for x in recs if x["target_gene"]]
                if recs:
                    br = min(recs, key=lambda x: x["eqtl_p"])
                    if not any(x["gene"] == br["gene"] for x in keep):
                        keep.append(br)
                lead_rows.extend(keep or recs[:1])

    # ------------------------------------------------------------ ENSG→symbol 补齐
    ensgs = set()
    for rows in (mr_rows, lead_rows, land_rows):
        for r in rows:
            for k in ("best_ensg", "gene_id", "ensg"):
                if r.get(k):
                    ensgs.add(r[k])
    emap = load_ensg_map(ensgs)
    for r in land_rows:
        r["symbol"] = emap.get(r.get("gene_id"), r.get("gene_id"))

    # ------------------------------------------------------------ 落盘
    pd.DataFrame(mr_rows).to_csv(os.path.join(OUT, "qtl_mr_by_gene.csv"), index=False)
    pd.DataFrame(lead_rows).to_csv(os.path.join(OUT, "lead_snp_eqtl_check.csv"), index=False)
    if ivw_rows:
        pd.DataFrame(ivw_rows).to_csv(os.path.join(OUT, "qtl_ivw_sensitivity.csv"), index=False)
    if land_rows:
        dfa = pd.DataFrame(land_rows).sort_values(["region", "tissue", "min_p"])
        dfa.to_csv(os.path.join(OUT, "region_cis_gene_landscape.csv"), index=False)

    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 50)
    print("\n== 主结果（cis-eQTL → 疾病）==")
    d = pd.DataFrame(mr_rows)
    if not d.empty:
        cols = [c for c in ["region", "gene", "tissue", "n_cis_q05", "n_cis_near_q05",
                            "n_shared", "tool_tier", "top_rsid", "top_dist_kb",
                            "eqtl_p", "eqtl_q", "F", "R2",
                            "OR_CHOLE", "ORl_CHOLE", "ORu_CHOLE", "p_CHOLE",
                            "OR_BTC", "p_BTC", "mde_CHOLE", "mde_BTC", "status"]
                if c in d.columns]
        print(d[cols].to_string(index=False))
    print("\n== GWAS lead SNP 的 cis-eQTL 状态 ==")
    dl = pd.DataFrame(lead_rows)
    if not dl.empty:
        cols = [c for c in ["region", "tissue", "lead_snp", "gene", "target_gene",
                            "eqtl_beta", "eqtl_p", "q_within_cis", "is_sig_p5e8",
                            "is_cis_eqtl_q05", "note"] if c in dl.columns]
        print(dl[cols].to_string(index=False))
    if ivw_rows:
        print("\n== IVW 敏感性（距离 prune）==")
        print(pd.DataFrame(ivw_rows).to_string(index=False))
    log("done ->", OUT)


if __name__ == "__main__":
    main()
