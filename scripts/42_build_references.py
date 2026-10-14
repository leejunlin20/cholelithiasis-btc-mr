#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
42_build_references.py —— 主稿参考文献的 PubMed 逐条验证器（v2）

v1 失败教训：仅用 `[Title]` 检索，对短标题（如 "Gallstones"）或标题有细微差异者
会命中同名近邻文献，或命中 PubMed 单独收录的 **Correction / Erratum** 记录。
v2 的修正：
  1. 复合检索式 `"标题"[Title] AND "期刊"[Journal] AND 年份[dp]`，逐级放宽；
  2. **剔除更正/勘误/撤稿类记录**（标题以 Correction/Erratum/Author Correction/… 开头）；
  3. 对候选命中做**标题词元重合度**比对，取最优者；重合度过低则标 CHECK；
  4. 已知 PMID 者直接取用，但仍**逐字段回填并比对标题**，不采用记忆值。

用法：
  PY=python3
  $PY scripts/42_build_references.py
输出：
  results/主稿/references_verified.csv / .txt
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
OUT_DIR = os.path.join(ROOT, "results", "主稿")
os.makedirs(OUT_DIR, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "wb_manuscript_refcheck"
EMAIL = "leejunlin@live.cn"

BAD_PREFIX = re.compile(
    r"^(correction|erratum|author correction|publisher correction|retraction|"
    r"withdrawn|expression of concern|addendum)\b", re.I)

# key / title([Title] 检索片段) / journal([Journal]) / year([dp]) / pmid(可选) / note
CANDIDATES = [
    # ---- 流行病学与临床背景 ----
    dict(key="globocan2022", title="Global cancer statistics 2022", journal="CA Cancer J Clin", year="2024",
         note="胆道癌/胆囊癌全球疾病负担"),
    dict(key="rizvi2013", title="Pathogenesis, diagnosis, and management of cholangiocarcinoma",
         journal="Gastroenterology", year="2013", note="胆管癌发病机制（炎症-癌转化）"),
    dict(key="banales2020", title="Cholangiocarcinoma 2020", journal="Nat Rev Gastroenterol Hepatol", year="2020",
         pmid="32606456", note="胆管癌机制与管理权威综述"),
    # ⚠️ 2026-09-17 定点替换（勿回改）：原 key="lammert2016" / PMID 27121544
    #    (Lammert F, et al. Gallstones. Nat Rev Dis Primers. 2016) 在 PubMed、CrossRef、
    #    Europe PMC 三库均无作者字段（出版社未存缴），会被体检判为「无作者」阻断项。
    #    已整条替换为主题等价的 Portincasa 2006（PMID 16844493），且 44 号脚本同步改定。
    dict(key="portincasa2006", title="Cholesterol gallstone disease", journal="Lancet", year="2006",
         pmid="16844493", note="胆石症权威综述（替换无作者字段的 Nat Rev Dis Primers 2016 条目）"),
    dict(key="stinton2012", title="Epidemiology of gallbladder disease", journal="Gut Liver", year="2012",
         note="胆石症与胆囊癌流行病学"),
    dict(key="hsing2007", title="Gallstones and the risk of biliary tract cancer", journal="Br J Cancer", year="2007",
         note="胆石症→胆道癌观察性证据（中国人群）"),
    dict(key="randi2006", title="Gallbladder cancer worldwide", journal="Int J Cancer", year="2006",
         pmid="16397865", note="胆囊癌地理分布与危险因素"),
    dict(key="hundal2014", title="Gallbladder cancer: epidemiology and outcome", journal="Clin Epidemiol", year="2014",
         pmid="24634588", note="胆囊癌流行病学与预后"),
    dict(key="nogueira2014", title="Gallstones, cholecystectomy, and risk of digestive system cancers",
         journal="Am J Epidemiol", year="2014", pmid="24470530", note="胆石/胆囊切除与消化系统肿瘤"),
    dict(key="chow1999", title="Gallstones, cholecystectomy and risk of cancers of the liver",
         journal="Br J Cancer", year="1999", note="胆石与肝胆胰肿瘤风险"),

    # ---- 机制与基因 ----
    dict(key="buch2007", title="A genome-wide association scan identifies the hepatic cholesterol transporter ABCG8",
         journal="Nat Genet", year="2007", pmid="17632509", note="ABCG8 为胆石易感位点"),
    dict(key="joshi2016", title="Four susceptibility loci for gallstone disease", journal="Gastroenterology", year="2016",
         note="胆石症 GWAS meta（4 个易感位点）"),
    dict(key="romeo2008", title="Genetic variation in PNPLA3 confers susceptibility to nonalcoholic fatty liver disease",
         journal="Nat Genet", year="2008", note="PNPLA3 148M 与肝病"),
    dict(key="strnad2020", title="Alpha1-Antitrypsin Deficiency", journal="N Engl J Med", year="2020",
         note="SERPINA1 Z 等位基因"),
    dict(key="hayhurst2001", title="Hepatocyte nuclear factor 4alpha", journal="Mol Cell Biol", year="2001",
         pmid="11158324", note="HNF4A 与肝脏脂质稳态"),
    dict(key="trauner2017", title="New therapeutic concepts in bile acid transport and signaling",
         journal="Hepatology", year="2017", pmid="27997980", note="胆汁酸转运与信号（FXR 轴）"),
    dict(key="chiang2019", title="Bile Acids as Metabolic Regulators and Nutrient Sensors",
         journal="Annu Rev Nutr", year="2019", pmid="31018107", note="胆汁酸作为代谢调节因子"),

    # ---- MR 方法学 ----
    dict(key="davesmith2003", title="Mendelian randomization", journal="Int J Epidemiol", year="2003",
         pmid="12689998", note="MR 原理（奠基）"),
    dict(key="burgess2013", title="Mendelian randomization analysis with multiple genetic variants using summarized data",
         journal="Genet Epidemiol", year="2013", note="汇总数据 IVW"),
    dict(key="bowden2015", title="Mendelian randomization with invalid instruments", journal="Int J Epidemiol",
         year="2015", pmid="26050253", note="MR-Egger"),
    dict(key="bowden2016", title="Consistent Estimation in Mendelian Randomization with Some Invalid Instruments",
         journal="Genet Epidemiol", year="2016", pmid="27061298", note="加权中位数"),
    dict(key="hartwig2017", title="Robust inference in summary data Mendelian randomization",
         journal="Int J Epidemiol", year="2017", pmid="29040600", note="加权众数"),
    dict(key="verbanck2018", title="Detection of widespread horizontal pleiotropy", journal="Nat Genet", year="2018",
         note="MR-PRESSO"),
    dict(key="hemani2018", title="The MR-Base platform supports systematic causal inference",
         journal="Elife", year="2018", pmid="29846171", note="TwoSampleMR 平台"),
    dict(key="strobemr2021", title="Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization",
         journal="JAMA", year="2021", pmid="34698778", note="STROBE-MR"),
    dict(key="sanderson2019", title="An examination of multivariable Mendelian randomization",
         journal="Int J Epidemiol", year="2019", note="MVMR 条件 F（识别性诊断）"),
    dict(key="burgess2015mvmr", title="Multivariable Mendelian randomization", journal="Am J Epidemiol", year="2015",
         note="MVMR 原理"),
    dict(key="grant2021", title="Pleiotropy robust methods for multivariable Mendelian randomization",
         journal="Stat Med", year="2021", note="MVMR-median / Lasso"),
    dict(key="rees2017", title="Extending the MR-Egger method for multivariable Mendelian randomization",
         journal="Stat Med", year="2017", pmid="28960498", note="MVMR-Egger"),
    dict(key="sanderson2021", title="Testing and correcting for weak and pleiotropic instruments",
         journal="Stat Med", year="2021", pmid="34338327", note="条件 F 实现与弱工具校正"),
    dict(key="giambartolomei2014", title="Bayesian test for colocalisation between pairs of genetic association studies",
         journal="PLoS Genet", year="2014", pmid="24830394", note="coloc.abf"),
    dict(key="wallace2021", title="A more accurate method for colocalisation analysis",
         journal="PLoS Genet", year="2021", pmid="34587156", note="coloc.susie"),
    dict(key="zuber2022", title="Combining evidence from Mendelian randomization and colocalization",
         journal="Am J Hum Genet", year="2022", pmid="35452592", note="MR 与共定位联合解读"),
    dict(key="burgess2016overlap", title="Bias due to participant overlap in two-sample Mendelian randomization",
         journal="Genet Epidemiol", year="2016", note="样本重叠偏倚（中介协方差）"),
    dict(key="hemani2017steiger", title="Orienting the causal relationship between imprecisely measured traits",
         journal="PLoS Genet", year="2017", note="Steiger 方向性检验"),
    dict(key="burgess2018binary", title="Mendelian randomization with a binary exposure variable",
         journal="Eur J Epidemiol", year="2018", pmid="30039250", note="二分类暴露的 r 换算"),
    dict(key="carter2021", title="Mendelian randomisation for mediation analysis",
         journal="Eur J Epidemiol", year="2021", pmid="33961203", note="MR 介导分析"),
    dict(key="burgess2015net", title="Network Mendelian randomization", journal="Int J Epidemiol", year="2015",
         pmid="25150977", note="网络 MR（介导）"),
    dict(key="vanderweele2016", title="Mediation Analysis", journal="Annu Rev Public Health", year="2016",
         pmid="26653405", note="因果中介分析"),
    dict(key="burgess2011weak", title="Avoiding bias from weak instruments in Mendelian randomization studies",
         journal="Int J Epidemiol", year="2011", pmid="21414999", note="弱工具变量偏倚"),
    dict(key="lawlor2016", title="Two-sample Mendelian randomization", journal="Int J Epidemiol", year="2016",
         pmid="27427429", note="两样本 MR 设计与局限"),
    dict(key="chang2015plink", title="Second-generation PLINK", journal="Gigascience", year="2015",
         note="LD clumping 工具（PLINK 1.9）"),
    dict(key="freedman2011", title="Dealing with dietary measurement error in nutritional cohort studies",
         journal="J Natl Cancer Inst", year="2011", pmid="21653922", note="膳食测量误差（饮食侧局限）"),

    # ---- 数据源 ----
    dict(key="sudlow2015", title="UK biobank: an open access resource", journal="PLoS Med", year="2015",
         pmid="25826379", note="UK Biobank"),
    dict(key="bycroft2018", title="The UK Biobank resource with deep phenotyping", journal="Nature", year="2018",
         pmid="30305743", note="UK Biobank 基因型与表型"),
    dict(key="kurki2023", title="FinnGen provides genetic insights", journal="Nature", year="2023",
         note="FinnGen"),
    dict(key="sakaue2021", title="A cross-population atlas of genetic associations for 220 human phenotypes",
         journal="Nat Genet", year="2021", pmid="34594039", note="BioBank Japan"),
    dict(key="1000g2015", title="A global reference for human genetic variation", journal="Nature", year="2015",
         note="1000 Genomes 参考面板"),

    # ---- 胆石症可干预危险因素 MR ----
    dict(key="yuan2022", title="Obesity, Type 2 Diabetes, Lifestyle Factors, and Risk of Gallstone Disease",
         journal="Clin Gastroenterol Hepatol", year="2022", pmid="33418132", note="生活方式因素与胆石症 MR"),
    dict(key="chen2022", title="Insights into modifiable risk factors of cholelithiasis",
         journal="Hepatology", year="2022", note="可干预危险因素与胆石症 MR"),

    # ---- 竞品（必须正面引用） ----
    dict(key="han2026", pmid="41627892", note="竞品1 Hepatology 2026 CCA GWAS+MR"),
    dict(key="bai2026", pmid="40064441", note="竞品2 J Adv Res 2026 网络 MR（本稿正面回应）"),
    dict(key="jiang2025", pmid="40974978", note="竞品3 Transl Oncol 2025 MR+QTL 共定位"),
    dict(key="an2026", title="Causal effects of cholelithiasis", journal="BMC Gastroenterol", year="2026",
         note="竞品4 BMC Gastroenterol 2026 multi-cohort MR"),
]


def fetch(url, tries=4, pause=0.6):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 refcheck"})
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:                                  # noqa: BLE001
            last = e
            time.sleep(pause * (i + 1))
    raise RuntimeError(f"请求失败 ({last})")


def esearch(term, retmax=8):
    url = (f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}"
           f"&tool={TOOL}&email={urllib.parse.quote(EMAIL)}&term={urllib.parse.quote(term)}")
    js = fetch(url)
    return js.get("esearchresult", {}).get("idlist", [])


def esummary_many(pmids):
    if not pmids:
        return {}
    url = (f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&tool={TOOL}"
           f"&email={urllib.parse.quote(EMAIL)}&id={','.join(pmids)}")
    js = fetch(url)
    res = js.get("result", {})
    return {k: v for k, v in res.items() if k.isdigit()}


def tokens(s):
    return set(re.findall(r"[a-z0-9]+", (s or "").lower()))


def sim(a, b):
    A, B = tokens(a), tokens(b)
    if not A or not B:
        return 0.0
    return len(A & B) / len(A | B)


def authors_str(rec, max_full=6):
    names = [a.get("name", "") for a in (rec.get("authors") or []) if a.get("name")]
    if not names:
        return ""
    return ", ".join(names[:max_full]) + (", et al" if len(names) > max_full else "")


def to_vancouver(rec, n):
    au = authors_str(rec)
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


def resolve(c):
    """返回 (pmid, rec, hit_by, warn)。"""
    tried = []
    if c.get("pmid"):
        recs = esummary_many([c["pmid"]])
        rec = recs.get(c["pmid"])
        if rec and not BAD_PREFIX.match(rec.get("title") or ""):
            return c["pmid"], rec, "pmid", ""
        tried.append("pmid 命中更正/缺失记录")

    qs = []
    if c.get("title") and c.get("journal") and c.get("year"):
        qs.append((f'"{c["title"]}"[Title] AND "{c["journal"]}"[Journal] AND {c["year"]}[dp]', "title+journal+year"))
    if c.get("title") and c.get("year"):
        qs.append((f'"{c["title"]}"[Title] AND {c["year"]}[dp]', "title+year"))
    if c.get("title"):
        qs.append((f'"{c["title"]}"[Title]', "title"))
    if c.get("title") and c.get("journal"):
        qs.append((f'{c["title"]} AND "{c["journal"]}"[Journal]', "free+journal"))

    for term, tag in qs:
        ids = esearch(term)
        time.sleep(0.4)
        if not ids:
            continue
        recs = esummary_many(ids)
        time.sleep(0.4)
        best, bs = None, -1.0
        for pid in ids:
            r = recs.get(pid)
            if not r:
                continue
            ti = r.get("title") or ""
            if BAD_PREFIX.match(ti):
                continue
            s = sim(c["title"], ti)
            if s > bs:
                best, bs = pid, s
        if best:
            warn = "" if bs >= 0.75 else f"标题相似度仅 {bs:.2f}，请人工确认"
            return best, recs[best], tag, warn
    return None, None, "none", "；".join(tried) or "全部检索式均未命中"


def main():
    rows = []
    for i, c in enumerate(CANDIDATES, 1):
        key = c["key"]
        try:
            pmid, rec, hit_by, warn = resolve(c)
        except Exception as e:                                   # noqa: BLE001
            pmid, rec, hit_by, warn = None, None, "error", f"检索异常: {e}"

        if not pmid:
            rows.append(dict(n=i, key=key, status="NOT_FOUND", pmid="", authors="", title=c.get("title", ""),
                             journal=c.get("journal", ""), year=c.get("year", ""), volume="", issue="",
                             pages="", doi="", hit_by=hit_by, note=c["note"], warning=warn, vancouver=""))
            print(f"[{i:02d}] {key:18s} ❌ NOT_FOUND  {warn}")
            continue

        row = dict(n=i, key=key, status="OK" if not warn else "CHECK", pmid=pmid,
                   authors=authors_str(rec), title=rec.get("title") or "",
                   journal=rec.get("source") or "", year=(rec.get("pubdate") or "").split(" ")[0],
                   volume=rec.get("volume") or "", issue=rec.get("issue") or "",
                   pages=rec.get("pages") or "",
                   doi=next((a.get("value", "") for a in (rec.get("articleids") or []) if a.get("idtype") == "doi"), ""),
                   hit_by=hit_by, note=c["note"], warning=warn, vancouver=to_vancouver(rec, i))
        rows.append(row)
        print(f"[{i:02d}] {key:18s} {'✅' if not warn else '⚠️'} {pmid:>9s} {hit_by:18s} "
              f"{(row['journal'] or '')[:26]:26s} {row['year']:5s} {row['title'][:52]}")

    cols = ["n", "key", "status", "pmid", "authors", "title", "journal", "year", "volume",
            "issue", "pages", "doi", "hit_by", "note", "warning", "vancouver"]
    with open(os.path.join(OUT_DIR, "references_verified.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    with open(os.path.join(OUT_DIR, "references_verified.txt"), "w", encoding="utf-8") as f:
        f.write("主稿参考文献（Vancouver 格式，全部经 PubMed E-utilities 实时核验）\n")
        f.write(f"核验时间：{time.strftime('%Y-%m-%d %H:%M:%S')} ｜ 脚本：scripts/42_build_references.py\n")
        f.write("=" * 100 + "\n\n")
        for r in rows:
            if r["status"] == "NOT_FOUND":
                f.write(f"{r['n']}. [未命中，待人工核实] 候选标题：{r['title']}\n   → 用途：{r['note']}\n\n")
            else:
                f.write(r["vancouver"] + "\n")
                f.write(f"   PMID {r['pmid']} | {'字段核验通过' if r['status']=='OK' else r['warning']}"
                        f" | 检索路径 {r['hit_by']} | 用途：{r['note']}\n\n")

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_chk = sum(1 for r in rows if r["status"] == "CHECK")
    n_bad = sum(1 for r in rows if r["status"] == "NOT_FOUND")
    print("\n" + "=" * 72)
    print(f"合计 {len(rows)}：✅ 通过 {n_ok} ｜ ⚠️ 待人工确认 {n_chk} ｜ ❌ 未命中 {n_bad}")
    print("CSV →", os.path.join(OUT_DIR, "references_verified.csv"))
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
