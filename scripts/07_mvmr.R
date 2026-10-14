#!/usr/bin/env Rscript
# 07_mvmr.R —— 多变量 MR（MVMR）：校正 BMI 后，胆囊结石病 → 胆道癌 的直接效应
#
# 动机：胆石与 BTC 都与肥胖强相关（BMI 是二者共同上游），必须做 MVMR 排除 BMI 混杂。
#
# 用法: Rscript scripts/07_mvmr.R <chole_iv_clumped.tsv> <btc_for_chole.tsv> <outdir> [bmi_id]

suppressMessages({
  library(data.table)
  library(ieugwasr)
  library(MVMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("用法: Rscript 07_mvmr.R <chole_iv.tsv> <btc_tsv> <outdir> [bmi_id=ieu-b-40]")
chole_f <- args[1]; btc_f <- args[2]; outdir <- args[3]
bmi_id  <- if (length(args) >= 4) args[4] else "ieu-b-40"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

chole <- fread(chole_f)
btc   <- fread(btc_f)
cat(sprintf("胆石 IV: %d | BTC 位点: %d\n", nrow(chole), nrow(btc)))

## ---- 统一 rsid 与等位基因列名 ----
setnames(chole, c("rsid", "alt", "ref", "beta", "se"),
                c("rsid", "ae", "ao", "b_ch", "s_ch"))
btc[, rsid := sub(",.*$", "", rsids)]
btc[, `:=`(ae_b = alt, ao_b = ref, b_bt = beta, s_bt = sebeta)]
btc_s <- btc[, .(rsid, ae_b, ao_b, b_bt, s_bt)]

## ---- 拉取 BMI 效应（允许 proxy，但只保留精确匹配）----
bmi <- associations(variants = chole$rsid, id = bmi_id)
bmi <- as.data.table(bmi)
cat(sprintf("BMI(%s) 返回: %d 行，其中 proxy 命中 %d\n",
            bmi_id, nrow(bmi), sum(bmi$proxy %in% TRUE)))
bmi <- bmi[is.na(proxy) | proxy == FALSE | proxy_snp == rsid]
bmi_s <- bmi[, .(rsid, ae_m = ea, ao_m = nea, b_bm = beta, s_bm = se)]
cat(sprintf("精确匹配保留: %d\n", nrow(bmi_s)))

## ---- 合并 ----
m <- merge(chole[, .(rsid, ae, ao, b_ch, s_ch)], btc_s, by = "rsid")
m <- merge(m, bmi_s, by = "rsid")
cat(sprintf("三方共有 SNP: %d\n", nrow(m)))

## ---- 以胆石等位基因为参照对齐（flip 符号）----
align <- function(ae, ao, ae_x, ao_x, b_x) {
  ok    <- ae == ae_x & ao == ao_x
  flip  <- ae == ao_x & ao == ae_x
  out   <- ifelse(ok, b_x, ifelse(flip, -b_x, NA_real_))
  out
}
m[, b_bt := align(ae, ao, ae_b, ao_b, b_bt)]
m[, b_bm := align(ae, ao, ae_m, ao_m, b_bm)]
m <- m[!is.na(b_bt) & !is.na(b_bm)]
cat(sprintf("等位基因对齐后可用 SNP: %d\n", nrow(m)))

if (nrow(m) < 5) stop("可用 SNP 太少，MVMR 不可靠")

## ---- MVMR ----
BXGs   <- as.matrix(m[, .(b_ch, b_bm)])
BYG    <- m$b_bt
seBXGs <- as.matrix(m[, .(s_ch, s_bm)])
seBYG  <- m$s_bt
mv <- format_mvmr(BXGs, BYG, seBXGs, seBYG, m$rsid)

cat("\n===== MVMR-IVW =====\n")
res <- ivw_mvmr(mv)
print(res)
res_dt <- data.table(exposure = c("Cholelithiasis", "BMI"),
                     b = res[, 1], se = res[, 2])
res_dt[, `:=`(or = exp(b), or_lci = exp(b - 1.96 * se), or_uci = exp(b + 1.96 * se),
              pval = 2 * pnorm(-abs(b / se)))]
print(res_dt, digits = 3)

cat("\n===== 条件 F 统计量（弱工具检验，>10 为佳）=====\n")
sf <- tryCatch(strength_mvmr(mv, gencov = 0), error = function(e) NULL)
if (!is.null(sf)) print(as.data.frame(sf))

cat("\n===== 多效性检验 =====\n")
pm <- tryCatch(pleiotropy_mvmr(mv, gencov = 0), error = function(e) NULL)
if (!is.null(pm)) print(as.data.frame(pm))

fwrite(res_dt, file.path(outdir, "mvmr_chole_bmi_to_btc.csv"))
fwrite(m, file.path(outdir, "mvmr_input_snps.csv"))
cat(sprintf("\n已写出: %s\n", file.path(outdir, "mvmr_chole_bmi_to_btc.csv")))
