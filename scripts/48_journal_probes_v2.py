#!/usr/bin/env python3
"""48_journal_probes_v2.py —— 候选期刊「多探针法」预例筛查 v2（journal-due-diligence Step 4a）

相对 36_journal_probes.py 的变化：
  1. 候选池从 22 刊扩到 34 刊，**剔除 Liver International**（李医生 2026-09-17 决定不投）；
     新增肝病专科（Liver Cancer / J Clin Transl Hepatol / World J Gastroenterol /
     Clin Gastroenterol Hepatol / Therap Adv Gastroenterol）、营养（AJCN / J Nutr / Nutrients）、
     流行病学与遗传学方法（IJE / AJE / Genet Epidemiol / Hum Mol Genet / Clin Transl Med）、
     综合与转化（EBioMedicine / J Intern Med / J Transl Med）、肿瘤（Int J Cancer / Br J Cancer /
     Cancer Epidemiol Biomarkers Prev）等方向。
  2. 新增两条探针，对准本稿的两个方法学卖点：
     P6 共定位  : colocali*（trait–trait 共定位是本文核心方法之一）
     P7 介导-MR : mediation AND Mendelian randomi*（两步介导 MR，PubMed 中罕见）

探针口径：
  P1 病种先例   : 胆石症/胆囊切除
  P2 设计先例   : Mendelian randomization
  P3 邻近病种   : 胆道癌/胆管癌
  P4 饮食×MR    : 饮食/营养标题 + MR
  P5 最严交集   : MR AND (胆石 OR 胆道癌)
  P6 共定位     : colocali*
  P7 介导-MR    : mediation AND MR

输出 TSV：results/journal_probes_v2.tsv
访问频率 ≤3 req/s（sleep 0.34），每次请求 timeout=30 并重试 4 次。

用法：
  PY=python3
  $PY scripts/48_journal_probes_v2.py
"""
import json
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
OUT = "results/journal_probes_v2.tsv"

# (显示名, PubMed 刊名缩写)
JOURNALS = [
    # —— 肝病 / 消化专科 ——
    ("Hepatology International", "Hepatol Int"),
    ("Annals of Hepatology", "Ann Hepatol"),
    ("Journal of Clinical and Translational Hepatology", "J Clin Transl Hepatol"),
    ("World Journal of Gastroenterology", "World J Gastroenterol"),
    ("Hepatology Communications", "Hepatol Commun"),
    ("Clinical Gastroenterology and Hepatology", "Clin Gastroenterol Hepatol"),
    ("Digestive and Liver Disease", "Dig Liver Dis"),
    ("Clinics and Research in Hepatology and Gastroenterology", "Clin Res Hepatol Gastroenterol"),
    ("United European Gastroenterology Journal", "United European Gastroenterol J"),
    ("Therapeutic Advances in Gastroenterology", "Ther Adv Gastroenterol"),  # 已校验：缩写正确形式
    ("Liver Cancer", "Liver Cancer"),
    ("Journal of Gastroenterology", "J Gastroenterol"),
    ("Journal of Gastroenterology and Hepatology", "J Gastroenterol Hepatol"),
    ("Alimentary Pharmacology & Therapeutics", "Aliment Pharmacol Ther"),
    # —— 营养 / 膳食 ——
    ("Clinical Nutrition", "Clin Nutr"),
    ("American Journal of Clinical Nutrition", "Am J Clin Nutr"),
    ("European Journal of Nutrition", "Eur J Nutr"),
    ("European Journal of Clinical Nutrition", "Eur J Clin Nutr"),
    ("Frontiers in Nutrition", "Front Nutr"),
    ("Journal of Nutrition", "J Nutr"),
    ("Nutrients", "Nutrients"),
    ("Nutrition & Metabolism", "Nutr Metab (Lond)"),
    ("Lipids in Health and Disease", "Lipids Health Dis"),
    # —— 流行病学 / 遗传学方法 ——
    ("European Journal of Epidemiology", "Eur J Epidemiol"),
    ("International Journal of Epidemiology", "Int J Epidemiol"),
    ("American Journal of Epidemiology", "Am J Epidemiol"),
    ("Genetic Epidemiology", "Genet Epidemiol"),
    ("Human Molecular Genetics", "Hum Mol Genet"),
    ("Clinical and Translational Medicine", "Clin Transl Med"),
    # —— 综合 / 转化 / 肿瘤 ——
    ("Journal of Internal Medicine", "J Intern Med"),
    ("BMC Medicine", "BMC Med"),
    ("EBioMedicine", "EBioMedicine"),
    ("Journal of Translational Medicine", "J Transl Med"),
    ("BMC Gastroenterology", "BMC Gastroenterol"),
    ("Cancer Medicine", "Cancer Med"),
    ("International Journal of Cancer", "Int J Cancer"),
    ("British Journal of Cancer", "Br J Cancer"),
    ("Cancer Epidemiology Biomarkers & Prevention", "Cancer Epidemiol Biomarkers Prev"),
    ("Scientific Reports", "Sci Rep"),
]

TA_FMT = '"{ta}"[ta]'

PROBES = {
    "p1_gallstone": '(cholelithiasis[tiab] OR gallstone*[tiab] OR "gall stones"[tiab] OR cholecystectomy[tiab])',
    "p2_mr": '("Mendelian randomization"[tiab] OR "Mendelian randomisation"[tiab])',
    "p3_btc": '("biliary tract cancer"[tiab] OR cholangiocarcinoma[tiab] OR "bile duct cancer"[tiab] OR "biliary tract carcinoma"[tiab])',
    "p4_diet_mr": '((diet[ti] OR dietary[ti] OR nutrition*[ti] OR food[ti]) AND ("Mendelian randomi*"[tiab]))',
    "p6_coloc": '(colocalization[tiab] OR colocalisation[tiab] OR colocalized[tiab] OR colocalised[tiab] OR colocalize[tiab] OR "co-localization"[tiab])',
    "p7_med_mr": '(mediation[tiab] AND ("Mendelian randomi*"[tiab]))',
}

P5_TERM = ('{ta} AND ("Mendelian randomization"[tiab] OR "Mendelian randomisation"[tiab]) '
           'AND (cholelithiasis[tiab] OR gallstone*[tiab] OR cholangiocarcinoma[tiab] '
           'OR "biliary tract cancer"[tiab])')


def count(term: str) -> int:
    url = BASE + "?db=pubmed&retmode=json&retmax=0&term=" + urllib.parse.quote(term)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as fh:
                return int(json.load(fh)["esearchresult"]["count"])
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                print(f"    ! 失败: {exc}", file=sys.stderr)
            time.sleep(1.5 * (attempt + 1))
    return -1


def main() -> None:
    keys = list(PROBES) + ["p5_mr_gall_btc"]
    rows = []
    for idx, (name, ta) in enumerate(JOURNALS, 1):
        t = TA_FMT.format(ta=ta)
        vals = {}
        for key, cond in PROBES.items():
            vals[key] = count(f"{t} AND {cond}")
            time.sleep(0.34)
        vals["p5_mr_gall_btc"] = count(P5_TERM.format(ta=t))
        time.sleep(0.34)
        rows.append((name, ta, vals))
        print(f"  [{idx:>2}/{len(JOURNALS)}] {name:<52} " +
              " ".join(f"{k.split('_')[0].upper()}={vals[k]:>4}" for k in keys), flush=True)

    hdr = ["journal", "ta"] + keys
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\t".join(hdr) + "\n")
        for name, ta, vals in rows:
            fh.write("\t".join([name, ta] + [str(vals[k]) for k in keys]) + "\n")
    print(f"\n已写出 {OUT}（{len(rows)} 刊 × {len(keys)} 探针）")


if __name__ == "__main__":
    main()
