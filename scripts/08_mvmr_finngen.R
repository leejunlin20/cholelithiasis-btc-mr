#!/usr/bin/env Rscript
# 08_mvmr_finngen.R —— MVMR（FinnGen 同源数据版）：校正 BMI 后 胆囊结石病 → 胆道癌
#
# 与 07_mvmr.R 的关键区别：
#   1. BMI 改用 FinnGen 自有端点 BMI_IRN（n=500,348），与胆石/BTC 同源
#      （同 GRCh38、同 dbSNP 版本）→ rsID 可完美匹配，不再依赖跨库 proxy
#   2. 工具变量取「胆石 IV ∪ BMI IV」并集 —— 这是 MVMR 识别的前提：
#      若只用胆石 IV，BMI 侧条件 F 会极低（实测 0.14），估计不可信
#
# 用法: Rscript scripts/08_mvmr_finngen.R <chole_snps.tsv> <bmi_snps.tsv> <btc_snps.tsv> <outdir>

suppressMessages({
  library(data.table)
  library(MVMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("用法: Rscript 08_mvmr_finngen.R <chole.tsv> <bmi.tsv> <btc.tsv> <outdir>")
f_ch <- args[1]; f_bm <- args[2]; f_bt <- args[3]; outdir <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) {
  d <- fread(f)
  d[, rsid := sub(",.*$", "", rsids)]
  d[, `:=`(b = as.numeric(beta), s = as.numeric(sebeta), p = as.numeric(pval),
           a1 = alt, a2 = ref)]
  d[, .(rsid, chr = as.character(`#chrom`), pos, a1, a2, b, s, p)]
}

ch <- rd(f_ch); bm <- rd(f_bm); bt <- rd(f_bt)
cat(sprintf("胆石位点 %d | BMI 位点 %d | BTC 位点 %d\n", nrow(ch), nrow(bm), nrow(bt)))

# 三方共有
m <- merge(ch, bm, by = "rsid", suffixes = c("_ch", "_bm"))
m <- merge(m, bt[, .(rsid, a1_bt = a1, a2_bt = a2, b_bt = b, s_bt = s)], by = "rsid")
cat(sprintf("三方共有 rsID: %d\n", nrow(m)))

# 以胆石等位基因为参照对齐（flip 符号；回文且频率接近的丢弃）
flip_if <- function(ae, ao, a1, a2, b) {
  ifelse(a1 == ae & a2 == ao, b,
  ifelse(a1 == ao & a2 == ae, -b, NA_real_))
}
m[, b_bt := flip_if(a1_ch, a2_ch, a1_bt, a2_bt, b_bt)]
m[, b_bm := flip_if(a1_ch, a2_ch, a1_bm, a2_bm, b_bm)]
m <- m[!is.na(b_bt) & !is.na(b_bm)]
cat(sprintf("等位基因对齐后可用 SNP: %d\n", nrow(m)))
if (nrow(m) < 10) stop("可用 SNP 太少（<10），MVMR 不可靠")

BXGs   <- as.matrix(m[, .(b_ch, b_bm)])
BYG    <- m$b_bt
seBXGs <- as.matrix(m[, .(s_ch, s_bm)])
seBYG  <- m$s_bt
mv <- format_mvmr(BXGs, BYG, seBXGs, seBYG, m$rsid)

cat("\n===== MVMR-IVW（校正 BMI）=====\n")
res <- ivw_mvmr(mv)
res_dt <- data.table(exposure = c("Cholelithiasis", "BMI (FinnGen BMI_IRN)"),
                     b = res[, 1], se = res[, 2])
res_dt[, `:=`(or = exp(b), or_lci = exp(b - 1.96 * se), or_uci = exp(b + 1.96 * se),
              pval = 2 * pnorm(-abs(b / se)))]
print(res_dt, digits = 3)

cat("\n===== 条件 F 统计量（弱工具检验，>10 为佳）=====\n")
sf <- tryCatch(strength_mvmr(mv, gencov = 0), error = function(e) NULL)
if (!is.null(sf)) print(as.data.frame(sf))

cat("\n===== 多效性（Q）检验 =====\n")
pm <- tryCatch(pleiotropy_mvmr(mv, gencov = 0), error = function(e) NULL)
if (!is.null(pm)) print(as.data.frame(pm))

fwrite(res_dt, file.path(outdir, "mvmr_finngen_chole_bmi_to_btc.csv"))
cat(sprintf("\n已写出: %s\n", file.path(outdir, "mvmr_finngen_chole_bmi_to_btc.csv")))
