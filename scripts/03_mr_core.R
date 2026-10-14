#!/usr/bin/env Rscript
# 03_mr_core.R —— 两样本 MR 主分析 + 敏感性分析
#
# 用法:
#   Rscript scripts/03_mr_core.R <exposures.rds> <outcome.tsv> <outcome_label> <outdir>
#
# <outcome.tsv> 由 02_extract_outcome_snps.py 从 FinnGen 汇总统计抽取而来
# 产物: <outdir>/mr_results_<label>.csv, mr_detail_<label>.rds

suppressMessages({
  library(TwoSampleMR)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("用法: Rscript 03_mr_core.R <exposures.rds> <outcome.tsv> <label> <outdir>")
expo_file <- args[1]; out_file <- args[2]; label <- args[3]; outdir <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

expo <- readRDS(expo_file)
cat(sprintf("暴露: %s 个数据集, %d 条 SNP-暴露记录\n",
            length(unique(expo$id.exposure)), nrow(expo)))

## ---- 读结局（FinnGen 格式）----
ff <- fread(out_file)
setnames(ff, gsub("^#", "", names(ff)))
cat(sprintf("结局文件变异数: %d\n", nrow(ff)))

# FinnGen release 文件没有 rsid 缺失时用 rsids 首值
ff[, SNP := tstrsplit(rsids, ",", keep = 1L)]
ff <- ff[!is.na(beta) & !is.na(sebeta) & !is.na(pval)]
ff[, N := fifelse(is.na(af_alt_cases), NA_real_, NA_real_)]  # 占位，后面用 manifest 的 N

outcome_dat <- data.frame(
  SNP                  = ff$SNP,
  beta.outcome         = as.numeric(ff$beta),
  se.outcome           = as.numeric(ff$sebeta),
  effect_allele.outcome = ff$alt,
  other_allele.outcome  = ff$ref,
  eaf.outcome          = as.numeric(ff$af_alt),
  pval.outcome         = as.numeric(ff$pval),
  outcome              = label,
  mr_keep.outcome      = TRUE,
  id.outcome           = label,
  stringsAsFactors     = FALSE
)
outcome_dat <- outcome_dat[!duplicated(outcome_dat$SNP), ]
cat(sprintf("与暴露可匹配的结局 SNP: %d\n", sum(outcome_dat$SNP %in% expo$SNP)))

## ----  harmonise ----
dat <- tryCatch(
  harmonise_data(expo, outcome_dat, action = 2),
  error = function(e) { cat("harmonise 失败:", conditionMessage(e), "\n"); NULL }
)
if (is.null(dat)) quit(status = 1)

dat <- dat[dat$mr_keep, ]
cat(sprintf("harmonise 后可用于 MR 的记录: %d\n", nrow(dat)))

# 每个暴露的 IV 数与 F 统计量
DT <- as.data.table(dat)
iv_stat <- DT[, .(nsnp = .N,
                  F_mean = mean((beta.exposure / se.exposure)^2, na.rm = TRUE)),
              by = .(id.exposure, exposure)]
print(iv_stat[order(-nsnp)][1:min(10, .N)])

## ----  MR 主分析 ----
# ⚠️ 不用 mr_ivw_mre：TwoSampleMR 的实现按 se*sqrt(Q/df) 缩放，当 Q < df 时
#    会把 SE 缩小到远低于固定效应（实测 se 0.438 → 0.028，P 从 0.27 变成 5e-68），
#    属于反保守。这里以固定效应 IVW 为基准，再按 max(1, sqrt(Q/df)) 手工做随机效应校正。
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode"))
res <- generate_odds_ratios(res)

## 敏感性分析
het  <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
loo  <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
stei <- tryCatch(directionality_test(dat), error = function(e) NULL)

saveRDS(list(dat = dat, res = res, het = het, plei = plei, loo = loo, steiger = stei),
        file.path(outdir, paste0("mr_detail_", label, ".rds")))

out <- as.data.table(res)
out <- merge(out, as.data.table(iv_stat), by = c("id.exposure", "exposure"), all.x = TRUE)
if (!is.null(plei)) {
  pl <- as.data.table(plei)[, .(id.exposure, egger_intercept, se_egger = se, pval_pleiotropy = pval)]
  out <- merge(out, pl, by = "id.exposure", all.x = TRUE)
}
if (!is.null(het)) {
  h <- as.data.table(het)[method == "Inverse variance weighted",
                          .(id.exposure, Q = Q, Q_df = Q_df, Q_pval)]
  out <- merge(out, h, by = "id.exposure", all.x = TRUE)
}

# 手工随机效应校正（Q>df 放大 SE；Q<df 不缩小）
out[, se_ivw_fe := se]
out[method == "Inverse variance weighted" & !is.na(Q_df) & Q_df > 0,
    se := se * pmax(1, sqrt(Q / Q_df))]
out[method == "Inverse variance weighted",
    `:=`(pval     = 2 * pnorm(-abs(b / se)),
         lo_ci    = b - 1.96 * se,
         up_ci    = b + 1.96 * se,
         or       = exp(b),
         or_lci95 = exp(b - 1.96 * se),
         or_uci95 = exp(b + 1.96 * se))]

PRIMARY <- "Inverse variance weighted"
out[, fdr_primary := NA_real_]
out[method == PRIMARY, fdr_primary := p.adjust(pval, method = "BH")]
out <- out[order(pval)]

fwrite(out, file.path(outdir, paste0("mr_results_", label, ".csv")))
cat(sprintf("\n已写出: %s\n", file.path(outdir, paste0("mr_results_", label, ".csv"))))

sig <- as.data.frame(out)[out$method == PRIMARY & !is.na(out$fdr_primary) & out$fdr_primary < 0.05, , drop = FALSE]
cat(sprintf("FDR<0.05 的暴露: %d 个（主方法: %s，随机效应校正后）\n", nrow(sig), PRIMARY))
if (nrow(sig)) {
  keep <- intersect(c("exposure", "nsnp.x", "nsnp.y", "b", "se", "se_ivw_fe", "pval",
                      "or", "or_lci95", "or_uci95", "Q_pval", "pval_pleiotropy", "fdr_primary"), names(sig))
  print(sig[, keep, drop = FALSE], digits = 3)
}
