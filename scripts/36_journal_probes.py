#!/usr/bin/env python3
"""36_journal_probes.py —— 候选期刊「三探针法」预例筛查（journal-due-diligence Step 4a）

对每本候选刊跑 5 条 PubMed esearch，用返回条数排序：
  P1 病种先例   : 胆石症/胆囊切除
  P2 设计先例   : Mendelian randomization
  P3 邻近病种   : 胆道癌/胆管癌
  P4 饮食×MR    : 饮食/营养标题 + MR
  P5 最严交集   : MR AND (胆石 OR 胆道癌)

输出 TSV：journal / p1 / p2 / p3 / p4 / p5
访问频率 ≤3 req/s（sleep 0.34）。
"""
import json
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

JOURNALS = [
    ("Hepatology International", "Hepatol Int"),
    ("Liver International", "Liver Int"),
    ("European Journal of Epidemiology", "Eur J Epidemiol"),
    ("Journal of Gastroenterology", "J Gastroenterol"),
    ("Alimentary Pharmacology & Therapeutics", "Aliment Pharmacol Ther"),
    ("Clinical Nutrition", "Clin Nutr"),
    ("Digestive and Liver Disease", "Dig Liver Dis"),
    ("United European Gastroenterology Journal", "United European Gastroenterol J"),
    ("Hepatology Communications", "Hepatol Commun"),
    ("Journal of Gastroenterology and Hepatology", "J Gastroenterol Hepatol"),
    ("European Journal of Nutrition", "Eur J Nutr"),
    ("European Journal of Clinical Nutrition", "Eur J Clin Nutr"),
    ("Annals of Hepatology", "Ann Hepatol"),
    ("Clinics and Research in Hepatology and Gastroenterology", "Clin Res Hepatol Gastroenterol"),
    ("Frontiers in Nutrition", "Front Nutr"),
    ("Nutrition (Elsevier)", "Nutrition"),
    ("Nutrition & Metabolism", "Nutr Metab (Lond)"),
    ("Lipids in Health and Disease", "Lipids Health Dis"),
    ("BMC Gastroenterology", "BMC Gastroenterol"),
    ("Cancer Medicine", "Cancer Med"),
    ("Scientific Reports", "Sci Rep"),
    ("BMC Medicine", "BMC Med"),
]

TA_FMT = '"{ta}"[ta]'

PROBES = {
    "p1_gallstone": '(cholelithiasis[tiab] OR gallstone*[tiab] OR "gall stones"[tiab] OR cholecystectomy[tiab])',
    "p2_mr": '("Mendelian randomization"[tiab] OR "Mendelian randomisation"[tiab] OR "mendelian randomization"[tiab])',
    "p3_btc": '("biliary tract cancer"[tiab] OR cholangiocarcinoma[tiab] OR "bile duct cancer"[tiab] OR "biliary tract carcinoma"[tiab])',
    "p4_diet_mr": '((diet[ti] OR dietary[ti] OR diet*[ti] OR nutrition*[ti] OR food[ti]) AND ("Mendelian randomi*"[tiab]))',
}


def count(term: str) -> int:
    url = BASE + "?db=pubmed&retmode=json&retmax=0&term=" + urllib.parse.quote(term)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as fh:
                return int(json.load(fh)["esearchresult"]["count"])
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return -1


def main() -> None:
    rows = []
    for name, ta in JOURNALS:
        t = TA_FMT.format(ta=ta)
        vals = {}
        for key, cond in PROBES.items():
            vals[key] = count(f"{t} AND {cond}")
            time.sleep(0.34)
        # P5：MR AND (胆石 OR 胆道癌)
        vals["p5_mr_gall_btc"] = count(
            f"{t} AND {PROBES['p2_mr']} AND (cholelithiasis[tiab] OR gallstone*[tiab] "
            f"OR cholangiocarcinoma[tiab] OR \"biliary tract cancer\"[tiab])"
        )
        time.sleep(0.34)
        vals["journal"], vals["ta"] = name, ta
        rows.append(vals)
        print(
            f"{name:58s} P1={vals['p1_gallstone']:5d} P2={vals['p2_mr']:5d} "
            f"P3={vals['p3_btc']:4d} P4={vals['p4_diet_mr']:4d} P5={vals['p5_mr_gall_btc']:3d}",
            flush=True,
        )

    cols = ["journal", "ta", "p1_gallstone", "p2_mr", "p3_btc", "p4_diet_mr", "p5_mr_gall_btc"]
    out = "results/pathA_journal_probes.tsv"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\n已写出 {out}")


if __name__ == "__main__":
    main()
