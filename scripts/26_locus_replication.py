#!/usr/bin/env python3
"""
26_locus_replication.py
-----------------------
权威版「跨血统位点级复制表」。

三边数据全部来自官方原始文件，不做任何坐标/等位基因推测：

  ① FinnGen R12 胆石症      K11_CHOLELITH                      (hg38, 芬兰)
  ② FinnGen R12 胆道癌      C3_BILIARY_GALLBLADDER_EXALLC      (hg38, 芬兰)
  ③ BBJ 肝内胆管癌          BBJ.BtC                            (hg19, 日本)

匹配方式：三边一律按 rsID「token 级精确匹配」（不做子串匹配——子串匹配会让
rs738408 撞上 rs738409x 之类的近名位点，属已修复的真实缺陷）。

输出
----
results/bbj_replication/locus_replication_v2.tsv
results/bbj_replication/locus_replication_v2_key.tsv   （6 个关键区域，供正文用表）
"""
import gzip, io, os, zipfile

BASE = os.environ.get("MR_BASE", ".")
OUT = f"{BASE}/results/bbj_replication"
BBJ_ZIP = f"{BASE}/data/bbj/hum0197.v3.BBJ.BtC.v1.zip"
BBJ_INNER = "hum0197.v3.BBJ.BtC.v1/GWASsummary_BtC_Japanese_SakaueKanai2020.auto.txt.gz"
FG = f"{BASE}/data/gwas"

KEY_REGIONS = {"R69", "R51", "R67", "R07a", "R07b", "R54"}


def load_panel():
    panel = {}
    with open(f"{BASE}/results/coloc/regions_v2.tsv") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rid = p[ix["region_id"]]
            for rs in p[ix["lead_snps"]].split(","):
                panel[rs] = dict(region=rid, chrom=p[ix["chrom"]],
                                 start=p[ix["start"]], end=p[ix["end"]],
                                 genes=p[ix["nearest_genes"]])
    return panel


def scan_finngen(path, want):
    hit = {}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").lstrip("#").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            for w in set(c[ix["rsids"]].split(",")) & want:
                if w in hit:
                    continue
                hit[w] = dict(pos=int(c[ix["pos"]]), ref=c[ix["ref"]], alt=c[ix["alt"]],
                              beta=float(c[ix["beta"]]), se=float(c[ix["sebeta"]]),
                              p=float(c[ix["pval"]]), af=float(c[ix["af_alt"]]),
                              gene=c[ix["nearest_genes"]])
    return hit


def scan_bbj(want):
    hit = {}
    with zipfile.ZipFile(BBJ_ZIP) as z, z.open(BBJ_INNER) as fh, \
            gzip.open(io.BufferedReader(fh), "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            toks = set(c[ix["SNPID"]].replace(",", ";").split(";"))
            for w in toks & want:
                if w in hit:
                    continue
                hit[w] = dict(chrom=c[ix["CHR"]], pos=int(float(c[ix["POS"]])),
                              ref=c[ix["Allele1"]], alt=c[ix["Allele2"]],
                              af=float(c[ix["AF_Allele2"]]),
                              beta=float(c[ix["BETA"]]), se=float(c[ix["SE"]]),
                              p=float(c[ix["p.value"]]),
                              info=float(c[ix["imputationInfo"]]))
    return hit


def load_bbj_api(want):
    """BBJ 胆石症（PheWeb API，已 2 位小数舍入）——仅供方向核对"""
    out = {}
    p = f"{OUT}/bbj_lookup_instr.tsv"
    if not os.path.exists(p):
        return out
    with open(p) as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[ix["pheno"]] == "Cholelithiasis":
                out[c[ix["rsid"]]] = dict(
                    beta=float(c[ix["beta_aligned"]]), se=float(c[ix["se_bbj"]]),
                    p=float(c[ix["p_bbj"]]), af=float(c[ix["af_bbj"]]))
    return out


def main():
    panel = load_panel()
    want = set(panel)
    print(f"位点面板（regions_v2 先导 SNP）: {len(want)} 个，覆盖 {len(set(v['region'] for v in panel.values()))} 个区域")

    print("扫 FinnGen R12 胆石症 ...")
    chole = scan_finngen(f"{FG}/finngen_R12_K11_CHOLELITH.gz", want)
    print(f"  命中 {len(chole)}")
    print("扫 FinnGen R12 胆道癌 ...")
    btc = scan_finngen(f"{FG}/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz", want)
    print(f"  命中 {len(btc)}")
    print("扫 BBJ.BtC 全量（约 1.3 万行/秒，需 1~2 分钟）...")
    bbj = scan_bbj(want)
    print(f"  命中 {len(bbj)}")
    api = load_bbj_api(want)
    print(f"  BBJ 胆石症（API）命中 {len(api)}")

    hdr = ["region", "genes", "region_span", "rsid",
           "fg_chole_pos", "fg_chole_beta", "fg_chole_se", "fg_chole_p", "fg_chole_af",
           "fg_btc_pos", "fg_btc_beta", "fg_btc_se", "fg_btc_p", "fg_btc_af",
           "bbj_btc_pos", "bbj_btc_ref", "bbj_btc_alt", "bbj_btc_beta", "bbj_btc_se",
           "bbj_btc_p", "bbj_btc_af", "bbj_btc_info",
           "bbj_chole_beta", "bbj_chole_p", "bbj_chole_af",
           "dir_fgchole_vs_fgbtc", "dir_fgbtc_vs_bbjbtc", "n_source"]

    rows = []
    for rs in sorted(want, key=lambda r: (panel[r]["chrom"], int(panel[r]["start"]), r)):
        m = panel[rs]
        c, b, j, a = chole.get(rs), btc.get(rs), bbj.get(rs), api.get(rs)
        n_src = sum(x is not None for x in (c, b, j))

        # 方向一致性只在「同一效应等位基因」下比较：以 FinnGen 的 alt 为参照
        d1 = d2 = "NA"
        if c and b:
            d1 = "YES" if (c["beta"] > 0) == (b["beta"] > 0) else "NO_OPPOSITE"
        if b and j:
            # FinnGen alt 与 BBJ alt 是否同一等位基因
            same = (b["alt"] == j["alt"]) or (
                {"A": "T", "T": "A", "C": "G", "G": "C"}.get(b["alt"]) == j["alt"])
            if same and b["ref"] == j["ref"]:
                d2 = "YES" if (b["beta"] > 0) == (j["beta"] > 0) else "NO_OPPOSITE"
            else:
                d2 = "ALLELE_MISMATCH"

        rows.append(dict(
            region=m["region"], genes=m["genes"],
            region_span=f"{m['start']}-{m['end']}", rsid=rs,
            fg_chole_pos=c["pos"] if c else "", fg_chole_beta=c["beta"] if c else "",
            fg_chole_se=c["se"] if c else "", fg_chole_p=c["p"] if c else "",
            fg_chole_af=c["af"] if c else "",
            fg_btc_pos=b["pos"] if b else "", fg_btc_beta=b["beta"] if b else "",
            fg_btc_se=b["se"] if b else "", fg_btc_p=b["p"] if b else "",
            fg_btc_af=b["af"] if b else "",
            bbj_btc_pos=j["pos"] if j else "", bbj_btc_ref=j["ref"] if j else "",
            bbj_btc_alt=j["alt"] if j else "",
            bbj_btc_beta=j["beta"] if j else "", bbj_btc_se=j["se"] if j else "",
            bbj_btc_p=j["p"] if j else "", bbj_btc_af=j["af"] if j else "",
            bbj_btc_info=j["info"] if j else "",
            bbj_chole_beta=a["beta"] if a else "", bbj_chole_p=a["p"] if a else "",
            bbj_chole_af=a["af"] if a else "",
            dir_fgchole_vs_fgbtc=d1, dir_fgbtc_vs_bbjbtc=d2, n_source=n_src))

    def write(path, rs):
        with open(path, "w") as o:
            o.write("\t".join(hdr) + "\n")
            for r in rs:
                o.write("\t".join(str(r[k]) for k in hdr) + "\n")

    write(f"{OUT}/locus_replication_v2.tsv", rows)
    key = [r for r in rows if r["region"] in KEY_REGIONS]
    write(f"{OUT}/locus_replication_v2_key.tsv", key)

    print(f"\n输出: {OUT}/locus_replication_v2.tsv ({len(rows)} 行)")
    print(f"      {OUT}/locus_replication_v2_key.tsv ({len(key)} 行)")
    print("\n================ 关键区域先导位点跨血统对照 ================")
    print("%-6s %-9s %-30s %10s %10s %10s %6s" %
          ("region", "rsid", "genes", "β_chole", "β_BTC", "β_BBJ", "dir"))
    for r in key:
        print("%-6s %-9s %-30s %10s %10s %10s %6s" % (
            r["region"], r["rsid"], r["genes"][:30],
            ("%+.3f" % r["fg_chole_beta"]) if r["fg_chole_beta"] != "" else "-",
            ("%+.3f" % r["fg_btc_beta"]) if r["fg_btc_beta"] != "" else "-",
            ("%+.3f" % r["bbj_btc_beta"]) if r["bbj_btc_beta"] != "" else "-",
            r["dir_fgbtc_vs_bbjbtc"]))


if __name__ == "__main__":
    main()
