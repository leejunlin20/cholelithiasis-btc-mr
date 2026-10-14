#!/usr/bin/env python3
"""权威复核：Ensembl GRCh37 坐标 + BBJ 精确 rsID 匹配（token 级，非子串）"""
import gzip, io, json, sys, urllib.request, zipfile
import os

BASE = os.environ.get("MR_BASE", ".")
ZIP = f"{BASE}/data/bbj/hum0197.v3.BBJ.BtC.v1.zip"
INNER = "hum0197.v3.BBJ.BtC.v1/GWASsummary_BtC_Japanese_SakaueKanai2020.auto.txt.gz"

RSIDS = ["rs11887534", "rs4299376", "rs4245791", "rs738409", "rs738408",
         "rs28929474", "rs1800961", "rs7206790", "rs146855983", "rs11888048"]

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def ensembl_grch37(rsids):
    """POST to grch37.rest.ensembl.org → GRCh37 mapping per rsID"""
    out = {}
    url = "https://grch37.rest.ensembl.org/variation/human"
    for i in range(0, len(rsids), 180):
        chunk = rsids[i:i + 180]
        req = urllib.request.Request(
            url, data=json.dumps({"ids": chunk}).encode(),
            headers={"Content-Type": "application/json",
                     "Accept": "application/json"})
        d = json.load(opener.open(req, timeout=60))
        for rsid, info in d.items():
            maps = [m for m in (info.get("mappings") or [])
                    if m.get("assembly_name") == "GRCh37"
                    and m.get("seq_region_name") in [str(c) for c in range(1, 23)]
                    + ["X", "Y"]]
            # 优先选「蛋白编码基因内」的映射（排除 alt/patch scaffold）
            out[rsid] = [(str(m["seq_region_name"]), int(m["start"]),
                          m.get("allele_string"), m.get("assembly_name"))
                         for m in maps]
    return out


print("== ① Ensembl GRCh37 权威坐标 ==")
ens = ensembl_grch37(RSIDS)
for r in RSIDS:
    print(f"  {r:14s} {ens.get(r)}")

# ---- ② BBJ 精确匹配（token 级）----
print("\n== ② BBJ.BtC 精确 rsID 匹配（split 后逐 token 比对）==")
want = set(RSIDS)
found = {}
with zipfile.ZipFile(ZIP) as z, z.open(INNER) as fh, \
        gzip.open(io.BufferedReader(fh), "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(h)}
    for line in f:
        c = line.rstrip("\n").split("\t")
        toks = set(c[ix["SNPID"]].replace(",", ";").split(";"))
        for w in toks & want:
            found.setdefault(w, []).append(
                (c[ix["CHR"]], c[ix["POS"]], c[ix["Allele1"]], c[ix["Allele2"]],
                 c[ix["AF_Allele2"]], c[ix["BETA"]], c[ix["SE"]],
                 c[ix["p.value"]], c[ix["imputationInfo"]]))
            if len(found[w]) > 6:
                found[w] = found[w][:6]

for r in RSIDS:
    if r not in found:
        print(f"  {r:14s} MISS（该 rsID 不在 BBJ.BtC 面板/未通过质控）")
        continue
    for chr_, pos, a1, a2, af, beta, se, p, info in found[r]:
        print("  %-14s chr%s:%-10s ref=%-22s alt=%-22s AF=%-8s beta=%+9s p=%-10s info=%s"
              % (r, chr_, pos, a1, a2, af, beta, p, info))

json.dump(ens, open("/tmp/ens_grch37.json", "w"), ensure_ascii=False, indent=1)
print("\n已写 /tmp/ens_grch37.json")
