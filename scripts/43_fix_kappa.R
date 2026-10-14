#!/usr/bin/env Rscript
# 43_fix_kappa.R —— 修正 IVW 行的「双重随机效应膨胀」并输出 v2 权威版本
#
# 【缺陷】scripts/03_mr_core.R 第 96–99 行：
#     out[, se_ivw_fe := se]                       # ← mr_ivw 的 se 本身已经是随机效应 se
#     out[method == "Inverse variance weighted",   # ← 于是这里再乘一次 κ = 乘了两次
#         se := se * pmax(1, sqrt(Q / Q_df))]
# 后果：列 `se_ivw_fe` 实际是【单次 κ 的随机效应 SE】（标签错误），
#       列 `se` 及其派生的 pval / lo_ci / up_ci / or_lci95 / or_uci95 是【κ²】，过度保守。
#
# 验证（酒精 → 胆石症，93 SNP）：
#   手算 FE se = 0.035245；Q = 453.9 (df 92)；κ = 2.2213
#   单次 κ（正确） = 0.078287  ← 与 CSV 的 `se_ivw_fe` 列一致
#   二次 κ（错误） = 0.173896  ← 与 CSV 的 `se` 列一致
#
# 【本脚本】从 mr_detail_*.rds 的 SNP 级数据独立复算，输出 *_v2.csv：
#   se_fe      固定效应 IVW SE（1/√Σw，w = β_exp²/se_out²）
#   kappa      max(1, √(Q/df))
#   se         随机效应 IVW SE = se_fe × kappa（**单次**）
# 其余方法（Egger / 加权中位数 / 加权众数）的 se 不受影响，原样保留。
# 同时按修正后的 p 值重算 fdr_primary。

suppressMessages({library(data.table)})

ROOT <- Sys.getenv("MR_BASE", unset = ".")
FULL <- file.path(ROOT, "results", "full")

run_one <- function(label) {
  rds <- file.path(FULL, paste0("mr_detail_", label, ".rds"))
  csv <- file.path(FULL, paste0("mr_results_", label, ".csv"))
  if (!file.exists(rds) || !file.exists(csv)) {
    cat("[跳过]", label, "：缺文件\n"); return(invisible(NULL))
  }
  d <- readRDS(rds)
  dat <- as.data.table(d$dat)
  old <- fread(csv)
  if ("mr_keep" %in% names(dat)) dat <- dat[mr_keep == TRUE]

  rec <- rbindlist(lapply(split(dat, dat$id.exposure), function(z) {
    if (nrow(z) < 2) return(NULL)
    w  <- z$beta.exposure^2 / z$se.outcome^2
    rt <- z$beta.outcome / z$beta.exposure
    b  <- sum(w * rt) / sum(w)
    se_fe <- 1 / sqrt(sum(w))
    Q  <- sum(w * (rt - b)^2); df <- nrow(z) - 1
    kap <- max(1, sqrt(Q / df))
    data.table(id.exposure = z$id.exposure[1], exposure = z$exposure[1],
               n_recalc = nrow(z), b_recalc = b, se_fe = se_fe,
               Q_recalc = Q, Q_df_recalc = df, kappa = kap, se_re = se_fe * kap)
  }), fill = TRUE)

  m <- merge(old, rec, by = c("id.exposure", "exposure"), all.x = TRUE)

  # 只改 IVW 行的 se 与派生量；se_ivw_fe 更名为 se_singlekappa 以消除标签歧义
  isw <- m$method == "Inverse variance weighted"
  m[, se_singlekappa := fifelse(isw, se_re, NA_real_)]
  m[, se_old := fifelse(isw, se, NA_real_)]
  m[isw & !is.na(se_re), `:=`(
    se       = se_re,
    pval     = 2 * pnorm(-abs(b / se_re)),
    lo_ci    = b - 1.96 * se_re,
    up_ci    = b + 1.96 * se_re,
    or       = exp(b),
    or_lci95 = exp(b - 1.96 * se_re),
    or_uci95 = exp(b + 1.96 * se_re))]
  m[isw, fdr_primary := p.adjust(pval, method = "BH")]
  m[, se_ivw_fe := NULL]           # 标签错误，删除以免误用
  m[, ratio_old_new := se_old / se]

  n_change <- sum(isw & !is.na(m$ratio_old_new) & abs(m$ratio_old_new - 1) > 1e-6)
  out <- file.path(FULL, paste0("mr_results_", label, "_v2.csv"))
  fwrite(m, out)

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("IVW 行 %d 条，其中被修正 %d 条（κ>1 者）\n",
              sum(isw), n_change))
  if (n_change) {
    z <- m[isw & abs(ratio_old_new - 1) > 1e-6][order(-ratio_old_new)]
    cat("修正幅度最大的 5 条（se 缩小倍数 = 原 κ）：\n")
    print(z[1:min(5, nrow(z)), .(exposure = substr(exposure, 1, 42), n = nsnp.x,
                                 kappa = round(ratio_old_new, 3),
                                 b = round(b, 4), se_old = round(se_old, 5),
                                 se_new = round(se, 5),
                                 p_old = signif(pval, 3), or_old = round(or_old <- or, 3))],
          digits = 4)
    cat("\nFDR<0.05（修正后）:", sum(m$method == "Inverse variance weighted" &
                                     !is.na(m$fdr_primary) & m$fdr_primary < 0.05), "个\n")
  }
  cat("已写出:", out, "\n")
  invisible(m)
}

for (lab in c("Cholelithiasis", "BTC")) run_one(lab)
