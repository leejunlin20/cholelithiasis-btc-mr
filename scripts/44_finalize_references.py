#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
44_finalize_references.py —— 参考文献定稿 + 逐条复核

在 42 号脚本（自动化检索）基础上，对需要人工判定的条目**显式指定 PMID**，
再逐条 fetch 回 PubMed 原始记录复核，输出最终可引用清单。
所有 PMID 均来自 PubMed 检索结果本身，未使用任何记忆值。

复核内容：标题关键词必须出现（防止命中同名近邻文献）；作者、期刊、年、卷、期、页逐字段回填。

输出：results/主稿/references_final.csv / .txt
"""

import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "主稿")
os.makedirs(OUT, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL, EMAIL = "wb_refcheck_final", "leejunlin@live.cn"

# key = (PMID, 标题必须包含的关键词片段列表, 用途)
FINAL = {
    # 流行病学与临床
    "globocan2022":     ("38572751", ["GLOBOCAN"], "胆道癌全球疾病负担"),
    "rizvi2013":        ("24140396", ["cholangiocarcinoma"], "胆管癌发病机制"),
    "banales2020":      ("32606456", ["Cholangiocarcinoma 2020"], "胆管癌机制与管理"),
    # ⚠️ 2026-09-17 定点替换（勿回改）：原 "lammert2016": ("27121544", ["Gallstones"] …)
    #    该 PMID 在 PubMed / CrossRef / Europe PMC 三库均无作者字段（出版社未存缴），
    #    成品侧会被判为「无作者」阻断项。已整条替换为 Portincasa 2006 (PMID 16844493)。
    #    若重跑本脚本把这一行改回去，会覆盖已核验的修正——务必保持现状。
    "portincasa2006":   ("16844493", ["Cholesterol gallstone disease"], "胆石症权威综述"),
    "stinton2012":      ("22570746", ["Epidemiology of gallbladder disease"], "胆石症与胆囊癌流行病学"),
    "hsing2007":        ("18000509", ["biliary tract cancer"], "胆石→胆道癌观察性证据（中国）"),
    "randi2006":        ("16397865", ["Gallbladder cancer worldwide"], "胆囊癌地理分布"),
    "nogueira2014":     ("24470530", ["digestive system cancers"], "胆石/胆囊切除与消化系肿瘤"),
    # 机制与基因
    "buch2007":         ("17632509", ["ABCG8"], "ABCG8 为胆石易感位点"),
    "joshi2016":        ("27094239", ["Susceptibility Loci for Gallstone"], "胆石症 GWAS meta"),
    "romeo2008":        ("18820647", ["PNPLA3"], "PNPLA3 148M 与肝病"),
    "silverman2009":    ("19553648", ["Alpha1-antitrypsin"], "SERPINA1 缺乏症（肝内炎症）"),
    "hayhurst2001":     ("11158324", ["Hepatocyte nuclear factor 4alpha"], "HNF4A 与肝脏脂质稳态"),
    "trauner2017":      ("27997980", ["bile acid"], "胆汁酸转运与信号（FXR 轴）"),
    # MR 方法学
    "davesmith2003":    ("12689998", ["Mendelian randomization"], "MR 原理"),
    "burgess2013":      ("24114802", ["summarized data"], "汇总数据 IVW"),
    "bowden2015":       ("26050253", ["Egger regression"], "MR-Egger"),
    "bowden2016":       ("27061298", ["weighted median"], "加权中位数"),
    "hartwig2017":      ("29040600", ["zero modal pleiotropy"], "加权众数"),
    "verbanck2018":     ("29686387", ["horizontal pleiotropy"], "MR-PRESSO"),
    "hemani2018":       ("29846171", ["MR-Base"], "TwoSampleMR 平台"),
    "strobemr2021":     ("34698778", ["STROBE-MR"], "STROBE-MR 报告规范"),
    "sanderson2019":    ("30535378", ["multivariable Mendelian randomization"], "MVMR 条件 F"),
    "burgess2015mvmr":  ("25632051", ["Multivariable Mendelian randomization"], "MVMR 原理"),
    "grant2021":        ("34342032", ["Pleiotropy robust methods"], "MVMR-median / Lasso"),
    "rees2017":         ("28960498", ["MR-Egger method for multivariable"], "MVMR-Egger"),
    "sanderson2021":    ("34338327", ["weak and pleiotropic instruments"], "条件 F 实现"),
    "giambartolomei2014": ("24830394", ["colocalisation"], "coloc.abf"),
    "wallace2021":      ("34587156", ["colocalisation"], "coloc.susie"),
    "zuber2022":        ("35452592", ["colocalization"], "MR 与共定位联合解读"),
    "burgess2016overlap": ("27625185", ["participant overlap"], "样本重叠偏倚"),
    "hemani2017steiger": ("29149188", ["imprecisely measured traits"], "Steiger 方向性检验"),
    "burgess2018binary": ("30039250", ["binary exposure"], "二分类暴露 r 换算"),
    "carter2021":       ("33961203", ["mediation analysis"], "MR 介导分析"),
    "burgess2015net":   ("25150977", ["Network Mendelian randomization"], "网络 MR"),
    "vanderweele2016":  ("26653405", ["Mediation Analysis"], "因果中介分析"),
    "burgess2011weak":  ("21414999", ["weak instruments"], "弱工具偏倚"),
    "lawlor2016":       ("27427429", ["Two-sample"], "两样本 MR 设计与局限"),
    "chang2015plink":   ("25722852", ["PLINK"], "LD clumping 工具"),
    "freedman2011":     ("21653922", ["dietary measurement error"], "膳食测量误差"),
    # 数据源
    "sudlow2015":       ("25826379", ["UK biobank"], "UK Biobank"),
    "bycroft2018":      ("30305743", ["UK Biobank resource"], "UK Biobank 基因型与表型"),
    "kurki2023":        ("36653562", ["FinnGen"], "FinnGen"),
    "sakaue2021":       ("34594039", ["cross-population atlas"], "BioBank Japan"),
    "1000g2015":        ("26432245", ["global reference for human genetic variation"], "1000 Genomes 参考面板"),
    # 胆石症危险因素 MR
    "chen2022":         ("34624136", ["modifiable risk factors of cholelithiasis"], "可干预危险因素与胆石症 MR"),
    # 竞品
    "han2026":          ("41627892", ["cholangiocarcinoma"], "竞品1 Hepatology 2026"),
    "bai2026":          ("40064441", ["cholelithiasis"], "竞品2 J Adv Res 2026"),
    "jiang2025":        ("40974978", ["cholelithiasis"], "竞品3 Transl Oncol 2025"),
    "an2026":           ("42350965", ["cholelithiasis"], "竞品4 BMC Gastroenterol 2026"),
}


def fetch(u):
    for i in range(4):
        try:
            r = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(r, timeout=45) as x:
                return json.loads(x.read().decode("utf-8", "replace"))
        except Exception:                                        # noqa: BLE001
            time.sleep(0.8 * (i + 1))
    raise RuntimeError("fetch failed")


def esummary(pmids):
    u = (f"{E}/esummary.fcgi?db=pubmed&retmode=json&tool={TOOL}"
         f"&email={urllib.parse.quote(EMAIL)}&id={','.join(pmids)}")
    r = fetch(u).get("result", {})
    return {k: v for k, v in r.items() if k.isdigit()}


def vancouver(rec, n):
    names = [a.get("name", "") for a in (rec.get("authors") or []) if a.get("name")]
    au = ", ".join(names[:6]) + (", et al" if len(names) > 6 else "") if names else "[无作者字段]"
    ti = (rec.get("title") or "").rstrip(".")
    so = rec.get("source") or ""
    yr = (rec.get("pubdate") or "").split(" ")[0]
    vol, iss, pg = rec.get("volume") or "", rec.get("issue") or "", rec.get("pages") or ""
    doi = next((a.get("value", "") for a in (rec.get("articleids") or []) if a.get("idtype") == "doi"), "")
    s = f"{n}. {au}. {ti}. {so}. {yr}"
    if vol:
        s += f";{vol}" + (f"({iss})" if iss else "")
    if pg:
        s += f":{pg}"
    s += "."
    if doi:
        s += f" doi:{doi}"
    return s


def main():
    keys = list(FINAL.keys())
    pmids = [FINAL[k][0] for k in keys]
    recs = {}
    for i in range(0, len(pmids), 20):
        recs.update(esummary(pmids[i:i + 20]))
        time.sleep(0.4)

    rows, bad = [], []
    for n, k in enumerate(keys, 1):
        pmid, kws, note = FINAL[k]
        rec = recs.get(pmid)
        if not rec:
            bad.append((k, pmid, "PubMed 无此记录"))
            print(f"[{n:02d}] {k:20s} ❌ PMID {pmid} 无记录")
            continue
        ti = rec.get("title") or ""
        miss = [w for w in kws if w.lower() not in ti.lower()]
        status = "OK" if not miss else "KEYWORD_MISS"
        if miss:
            bad.append((k, pmid, f"标题缺少关键词 {miss}"))
        row = dict(n=n, key=k, status=status, pmid=pmid, authors=vancouver(rec, n).split(". ", 1)[1],
                   title=ti, journal=rec.get("source") or "",
                   year=(rec.get("pubdate") or "").split(" ")[0], volume=rec.get("volume") or "",
                   issue=rec.get("issue") or "", pages=rec.get("pages") or "",
                   doi=next((a.get("value", "") for a in (rec.get("articleids") or []) if a.get("idtype") == "doi"), ""),
                   note=note, vancouver=vancouver(rec, n))
        rows.append(row)
        flag = "✅" if not miss else "⚠️"
        print(f"[{n:02d}] {k:20s} {flag} {pmid:>9s} {row['journal'][:26]:26s} {row['year']:5s} "
              f"{vol_pg(row):16s} {ti[:56]}")

    cols = ["n", "key", "status", "pmid", "authors", "title", "journal", "year", "volume",
            "issue", "pages", "doi", "note", "vancouver"]
    with open(os.path.join(OUT, "references_final.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(OUT, "references_final.txt"), "w", encoding="utf-8") as f:
        f.write("主稿参考文献（定稿，全部经 PubMed E-utilities 逐字段复核）\n")
        f.write(f"复核时间：{time.strftime('%Y-%m-%d %H:%M:%S')} ｜ 脚本：scripts/44_finalize_references.py\n")
        f.write("=" * 100 + "\n\n")
        for r in rows:
            f.write(r["vancouver"] + f"\n   PMID {r['pmid']} | {r['note']}\n\n")

    print(f"\n合计 {len(rows)} 条；其中需人工复核 {len(bad)} 条")
    for k, pid, why in bad:
        print(f"  ⚠️ {k} (PMID {pid}): {why}")
    print("已写出:", os.path.join(OUT, "references_final.csv"))
    return 0


def vol_pg(r):
    v = r.get("volume") or ""
    return f"{v}:{r.get('pages','')}" if v else ""


if __name__ == "__main__":
    sys.exit(main())
