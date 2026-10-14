#!/usr/bin/env Rscript
# 06_medi_mr.R —— 两步介导 MR：饮食 → 胆囊结石病 → 胆道癌
#
# 设计（本项目核心创新点）：
#   β1  = 饮食 → 胆石（源自 03_mr_core.R 的 Cholelithiasis 结果）
#   β2  = 胆石 → 胆道癌（本脚本 Step 2 计算）
#   βtot= 饮食 → 胆道癌（总效应）
#   间接效应 = β1 × β2 ；介导比例 = β1×β2 / βtot （delta 法求 SE）
#
# 用法: Rscript scripts/06_medi_mr.R <chole_iv_raw.tsv> <btc_outcome_tsv> <outdir>

suppressMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("用法: Rscript 06_medi_mr.R <chole_iv_raw.tsv> <btc_outcome.tsv> <outdir>")
iv_file <- args[1]; btc_file <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------- Step 0: 读入胆石 IV 并 clump ----------
iv <- fread(iv_file)
cat(sprintf("原始显著位点: %d\n", nrow(iv)))

cl <- tryCatch(
  ld_clump(data.table(rsid = iv$rsid, pval = iv$pval, id = "chole"),
           clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8, pop = "EUR"),
  error = function(e) { cat("clump 失败:", conditionMessage(e), "\n"); NULL }
)
if (is.null(cl) || !nrow(cl)) stop("clump 未产生结果，无法继续")
cat(sprintf("clump 后独立 IV: %d\n", nrow(cl)))
iv <- iv[rsid %in% cl$rsid]
fwrite(iv, file.path(outdir, "chole_iv_clumped.tsv"), sep = "\t")
writeLines(iv$rsid, file.path(outdir, "chole_iv_rsids.txt"))

## ---------- Step 2: 胆石 → 胆道癌 ----------
exp_dat <- data.table(
  SNP = iv$rsid, beta.exposure = iv$beta, se.exposure = iv$se,
  effect_allele.exposure = iv$alt, other_allele.exposure = iv$ref,
  eaf.exposure = iv$eaf, pval.exposure = iv$pval,
  id.exposure = "Cholelithiasis", exposure = "Cholelithiasis (49,834 cases)",
  mr_keep.exposure = TRUE, samplesize.exposure = 487252
)

out_dat <- fread(btc_file)
cat(sprintf("结局位点文件: %d 行, 列: %s\n", nrow(out_dat), paste(names(out_dat), collapse = ",")))

# FinnGen 原始列名 → TwoSampleMR 所需列名（rsids 列可能是 "rs1,rs2" 形式）
if ("rsids" %in% names(out_dat)) {
  out_dat[, SNP := sub(",.*$", "", rsids)]
  setnames(out_dat,
           c("alt", "ref", "beta", "sebeta", "pval", "af_alt"),
           c("effect_allele.outcome", "other_allele.outcome",
             "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome"),
           skip_absent = TRUE)
}
need <- c("SNP","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome")
missing <- setdiff(need, names(out_dat))
if (length(missing)) stop("结局文件缺列: ", paste(missing, collapse = ", "))
out_dat[, outcome := "Biliary tract cancer (2,298 cases)"]
out_dat[, id.outcome := "BTC"]

dat <- harmonise_data(exp_dat, out_dat, action = 2)
dat <- dat[dat$mr_keep, ]
cat(sprintf("harmonise 后可用 SNP: %d\n", nrow(dat)))
if (nrow(dat) < 3) stop("可用 SNP 太少，无法做稳健 MR")

res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode"))
res <- generate_odds_ratios(as.data.frame(res))
fwrite(as.data.table(res), file.path(outdir, "step2_chole_to_btc.csv"))

het  <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
if (!is.null(het))  fwrite(as.data.table(het),  file.path(outdir, "step2_heterogeneity.csv"))
if (!is.null(plei)) fwrite(as.data.table(plei), file.path(outdir, "step2_pleiotropy.csv"))
saveRDS(list(dat = dat, res = res), file.path(outdir, "step2_detail.rds"))

cat("\n===== Step 2: 胆囊结石病 → 胆道癌 =====\n")
print(res[, c("method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")], digits = 3)
if (!is.null(plei)) cat(sprintf("\nMR-Egger 截距 P = %.3f\n", plei$pval[1]))
if (!is.null(het))  cat(sprintf("Cochran Q P = %.3g\n", het$Q_pval[1]))

cat(sprintf("\n已写出: %s\n", file.path(outdir, "step2_chole_to_btc.csv")))
