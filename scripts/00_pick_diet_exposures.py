#!/usr/bin/env python3
"""
从 OpenGWAS 全量元数据中筛出 UK Biobank 饮食相关暴露（即文献所称 "226 dietary factors"）。

用法:
  python3 00_pick_diet_exposures.py /tmp/all_meta.json out/diet_traits.csv

规则:
  1. 仅保留 ukb-a-* (Neale 2017 round1) 与 ukb-b-* (Neale 2018 round2)
  2. 用精心设计的正则白名单命中"饮食"语义
  3. 命中黑名单的排除（体成分、家族史、用药等被宽松关键词误伤的表型）
  4. 同名表型在 round1/round2 都出现时，保留样本量更大的那个（通常 round2）
"""
import csv
import json
import re
import sys
from collections import OrderedDict

INCLUDE = [
    r"\bintake\b", r"\bdiet\b", r"\bdietary\b", r"\bfood\b",
    r"coffee", r"\btea\b", r"alcohol", r"\bbeer\b", r"cider", r"\bwine\b",
    r"spirits", r"champagne", r"fortified wine", r"drinking",
    r"fruit", r"vegetable", r"salad", r"\bmeat\b", r"\bfish\b", r"poultry",
    r"\bbeef\b", r"\bpork\b", r"\blamb\b", r"mutton",
    r"cheese", r"\bmilk\b", r"yogurt", r"yoghurt", r"bread", r"cereal",
    r"\brice\b", r"pasta", r"salt", r"sugar", r"chocolate", r"cake",
    r"biscuit", r"\begg", r"\bnut\b", r"nuts\b", r"legume", r"\bbean",
    r"potato", r"butter", r"margarine", r"spread", r"\boil\b", r"cooking",
    r"supplement", r"vitamin", r"mineral", r"fizzy", r"juice", r"smoothie",
    r"\bwater\b", r"\bdrink", r"hot drink", r"processed meat",
    r"spoon", r"tablespoon", r"mouthful", r"heaped", r"level tablespoon",
    r"never eat", r"food intolerance", r"roast", r"fried", r"grilled",
    r"raw vegetable", r"cooked vegetable", r"dried fruit", r"fresh fruit",
    r"main type of", r"variation in diet", r"major dietary changes",
    r"meat consumption", r"fish consumption", r"number of",
]
INCLUDE_RE = re.compile("|".join(INCLUDE), re.I)

EXCLUDE = [
    r"fat-free mass", r"fat mass", r"fat percentage", r"body fat",
    r"illnesses of", r"non-cancer illness", r"medication", r"treatment/medication",
    r"family history", r"father", r"mother", r"sibling",
    r"cholesterol", r"blood", r"height", r"\bweight\b", r"\bbmi\b",
    r"whole body", r"\barm\b", r"\bleg\b", r"trunk", r"bone",
    r"walking", r"sleep", r"smoking", r"vascular", r"diastolic", r"systolic",
    # 关键词误伤（非饮食表型）
    r"number of", r"matches in round", r"diagnoses", r"self-reported cancers",
    r"gas or solid-fuel", r"cooking/heating", r"operations",
    r"alcohol use\b", r"icd10", r"cancer code",
    r"heating type", r"invitation to complete", r"reason for reducing",
]
EXCLUDE_RE = re.compile("|".join(EXCLUDE), re.I)


def norm(trait: str) -> str:
    """用于跨 round 去重的归一化 trait 名。"""
    t = trait.lower()
    t = re.sub(r"^(average weekly|average |mean |number of )", "", t)
    t = re.sub(r"\s+", " ", t)
    return t.strip()


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    src = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "diet_traits.csv"

    with open(src) as f:
        meta = json.load(f)
    print(f"元数据总条目: {len(meta)}")

    kept = OrderedDict()
    for gid, d in meta.items():
        if not (gid.startswith("ukb-a-") or gid.startswith("ukb-b-")):
            continue
        trait = d.get("trait") or ""
        if not INCLUDE_RE.search(trait) or EXCLUDE_RE.search(trait):
            continue
        key = norm(trait)
        n = d.get("sample_size") or 0
        try:
            n = int(n)
        except (TypeError, ValueError):
            n = 0
        row = {
            "id": gid,
            "trait": trait,
            "sample_size": n,
            "ncase": d.get("ncase") or "",
            "ncontrol": d.get("ncontrol") or "",
            "population": d.get("population") or "European",
            "author": d.get("author") or "",
            "year": d.get("year") or "",
            "sex": d.get("sex") or "both",
            "category": d.get("category") or "",
        }
        prev = kept.get(key)
        if prev is None or row["sample_size"] > prev["sample_size"]:
            kept[key] = row

    rows = sorted(kept.values(), key=lambda r: r["trait"].lower())
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["id"])
        w.writeheader()
        w.writerows(rows)

    r1 = sum(1 for r in rows if r["id"].startswith("ukb-a-"))
    r2 = sum(1 for r in rows if r["id"].startswith("ukb-b-"))
    print(f"去重后饮食暴露: {len(rows)}  (ukb-a round1: {r1}, ukb-b round2: {r2})")
    print(f"已写出: {out}")
    for r in rows[:15]:
        print(f"  {r['id']} | n={r['sample_size']} | {r['trait'][:60]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
