#!/usr/bin/env bash
# 32_watch_extract.sh —— 在 MVMR 并集观察名单上，从 5 份 FinnGen 全量 GWAS 取效应
#
# 为什么要这一步：MVMR 要求所有暴露共用同一套工具变量。每个 SNP 必须在
# **全部暴露与结局**上都有效应估计，缺任何一侧该 SNP 都要丢弃；只用自己的
# 显著位点集去结局文件取数，会让胆囊炎/其他胆囊疾病的条件 F 直接塌掉。
#
# 前置：Rscript scripts/31_clump_mvmr_ivs.R 已生成 results/mvmr/watch_all_rsids.txt
set -euo pipefail

cd "${MR_BASE:?set MR_BASE}"
PY=python3
OUT=results/mvmr
W=$OUT/watch_all_rsids.txt

[ -f "$W" ] || { echo "缺观察名单 $W，请先跑 31_clump_mvmr_ivs.R"; exit 1; }
echo "观察名单: $(wc -l < "$W" | tr -d ' ') 个 rsID"

declare -a PAIRS=(
  "K11_CHOLELITH:watch_cholelith.tsv"
  "K11_CHOLECYST:watch_cholecyst.tsv"
  "K11_GALLBLAOTH:watch_gallblaoth.tsv"
  "K11_CHOLECYSTECTOMY:watch_cholecstec.tsv"
  "C3_BILIARY_GALLBLADDER_EXALLC:watch_btc.tsv"
  "BMI_IRN:watch_bmi.tsv"
)

for pair in "${PAIRS[@]}"; do
  ep=${pair%%:*}; out=${pair##*:}
  f="data/gwas/finngen_R12_${ep}.gz"
  if [ ! -f "$f" ]; then echo "!! 缺文件 $f，跳过"; continue; fi
  $PY scripts/30_mvmr_extract.py watch "$f" "$W" "$OUT/$out"
  echo
done

echo "完成。产物："
ls -la "$OUT"/watch_*.tsv
