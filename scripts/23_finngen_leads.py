#!/usr/bin/env python3
"""确定 FinnGen R12 侧的先导位点坐标与效应（hg38），并用 Ensembl 交叉验证"""
import gzip, io, json
import os

BASE = os.environ.get("MR_BASE", ".")
FG = f"{BASE}/data/gwas"
FILES = {
    "chole": f"{FG}/finngen_R12_K11_CHOLELITH.gz",
    "btc":   f"{FG}/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz",
}
PANEL = {  # rsid -> (region, gene, hg37 from Ensembl, allele_maj)
    "rs11887534": ("R07", "ABCG8 D19H"),
    "rs4299376":  ("R07", "ABCG5/8"),
    "rs4245791":  ("R07", "ABCG5/8"),
    "rs738409":   ("R69", "PNPLA3"),
    "rs738408":   ("R69", "PNPLA3"),
    "rs1800961":  ("R67", "HNF4A"),
    "rs7206790":  ("R54", "FTO"),
    "rs28929474": ("R51", "SERPINA1 PiZ"),
}
want = set(PANEL)

print("== FinnGen R12 侧（按 rsids token 精确匹配）==")
out = {}
for tag, path in FILES.items():
    hit = {}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").lstrip("#").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            toks = set(c[ix["rsids"]].split(","))
            for w in toks & want:
                if w in hit:
                    continue
                hit[w] = dict(chrom=c[ix["chrom"]], pos=int(c[ix["pos"]]),
                              ref=c[ix["ref"]], alt=c[ix["alt"]],
                              pval=float(c[ix["pval"]]), beta=float(c[ix["beta"]]),
                              se=float(c[ix["sebeta"]]), af=float(c[ix["af_alt"]]),
                              gene=c[ix["nearest_genes"]])
    out[tag] = hit
    print(f"\n--- {tag} ({path.split('/')[-1]}) ---")
    for r in PANEL:
        if r in hit:
            d = hit[r]
            print("  %-14s %s:%d ref=%s alt=%s af_alt=%.4f beta=%+.4f se=%.4f p=%.3g  near=%s"
                  % (r, d["chrom"], d["pos"], d["ref"], d["alt"], d["af"],
                     d["beta"], d["se"], d["pval"], d["gene"]))
        else:
            print("  %-14s MISS" % r)

json.dump(out, open("/tmp/finngen_leads.json", "w"), indent=1)
print("\n已写 /tmp/finngen_leads.json")
