#!/usr/bin/env Rscript
# 04_mr_presso.R v2 —— MR-PRESSO 离群值检测与校正 IVW
#
# 修正（v1 的 bug）：
#   1. Main MR results 的 MR Analysis 取值是 "Raw" / "Outlier-corrected"，
#      不是 "Causal Estimate"（v1 因此取到 0 行，列被丢弃 → 汇总时报 object not found）
#   2. global test 在 r$`MR-PRESSO results`$`Global Test`$Pvalue，不在顶层
#   3. 离群值索引取自 r$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
#   4. 增量落盘：每处理若干个暴露就重写 CSV，避免长任务全损
#
# 用法: Rscript scripts/04_mr_presso.R <mr_detail.rds> <outdir> [min_snp=5] [nbdist=1000]

suppressMessages({
  library(MRPRESSO)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法: Rscript 04_mr_presso.R <mr_detail.rds> <outdir> [min_snp=5] [nbdist=1000]")
detail  <- args[1]
outdir  <- args[2]
min_snp <- if (length(args) >= 3) as.integer(args[3]) else 5L
nbdist  <- if (length(args) >= 4) as.integer(args[4]) else 1000L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "mr_presso_results.csv")

d <- as.data.frame(readRDS(detail)$dat)
ids <- names(which(table(d$id.exposure) >= min_snp))
cat(sprintf("待做 MR-PRESSO 的暴露: %d 个（IV >= %d, NbDistribution=%d）\n",
            length(ids), min_snp, nbdist))

num <- function(v) {
  if (is.null(v)) return(NA_real_)
  v <- as.character(v)[1]
  v <- gsub("[<>]", "", v)          # "<0.002" -> "0.002"
  suppressWarnings(as.numeric(v))
}

res   <- vector("list", length(ids))
t0    <- Sys.time()

flush_out <- function(res, path) {
  R <- rbindlist(res[!vapply(res, is.null, logical(1))], fill = TRUE)
  if (!nrow(R)) return(invisible(NULL))
  for (col in c("b_raw","se_raw","p_raw","b_corr","se_corr","p_corr","global_p")) {
    if (!col %in% names(R)) R[, (col) := NA_real_]
  }
  R[, or_corr     := exp(b_corr)]
  R[, or_lci_corr := exp(b_corr - 1.96 * se_corr)]
  R[, or_uci_corr := exp(b_corr + 1.96 * se_corr)]
  R[, fdr_corr    := p.adjust(p_corr, method = "BH")]
  setorder(R, p_corr, na.last = TRUE)
  fwrite(R, path)
  invisible(R)
}

for (i in seq_along(ids)) {
  id <- ids[i]
  x  <- d[d$id.exposure == id, ]
  nm <- sub(" \\|\\|.*$", "", x$exposure[1])

  out <- tryCatch(
    suppressWarnings(
      mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = x, NbDistribution = nbdist, SignifThreshold = 0.05)
    ),
    error = function(e) structure(conditionMessage(e), class = "presso_err")
  )

  if (inherits(out, "presso_err")) {
    res[[i]] <- data.table(id = id, exposure = nm, nsnp = nrow(x),
                           status = "failed", err = as.character(out))
    next
  }

  main <- as.data.table(out$`Main MR results`, keep.rownames = FALSE)
  pr   <- out$`MR-PRESSO results`
  raw  <- main[`MR Analysis` == "Raw", ]
  corr <- main[`MR Analysis` == "Outlier-corrected", ]

  outliers <- pr$`Distortion Test`$`Outliers Indices`
  nout <- if (is.null(outliers)) 0L else length(unlist(outliers))
  distort_p <- num(pr$`Distortion Test`$Pvalue)
  distort_b <- num(pr$`Distortion Test`$`Distortion Coefficient`)

  res[[i]] <- data.table(
    id = id, exposure = nm, nsnp = nrow(x),
    nsnp_outlier = nout,
    global_p = num(pr$`Global Test`$Pvalue),
    # 注意：点估计列名是带空格的 `Causal Estimate`，用 $Estimate 会静默取到 NULL
    b_raw  = if (nrow(raw))  num(raw[["Causal Estimate"]])  else NA_real_,
    se_raw = if (nrow(raw))  num(raw[["Sd"]])                else NA_real_,
    p_raw  = if (nrow(raw))  num(raw[["P-value"]])           else NA_real_,
    b_corr  = if (nrow(corr)) num(corr[["Causal Estimate"]]) else NA_real_,
    se_corr = if (nrow(corr)) num(corr[["Sd"]])               else NA_real_,
    p_corr  = if (nrow(corr)) num(corr[["P-value"]])          else NA_real_,
    distort_p = distort_p, distort_b = distort_b,
    status = "ok"
  )

  if (i %% 5 == 0 || i == length(ids)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  进度 %d/%d  (%.1f 分钟, 约 %.1f 秒/暴露)\n",
                i, length(ids), el, el * 60 / i))
    flush_out(res, outfile)
  }
}

R <- flush_out(res, outfile)
cat(sprintf("\n完成。已写出 %s\n", outfile))
if (!is.null(R)) {
  cat(sprintf("成功 %d 个，失败 %d 个\n",
              sum(R$status == "ok", na.rm = TRUE), sum(R$status == "failed", na.rm = TRUE)))
  sig <- R[!is.na(fdr_corr) & fdr_corr < 0.05]
  cat(sprintf("校正后 FDR<0.05 的暴露: %d 个\n", nrow(sig)))
  if (nrow(sig)) {
    print(sig[order(p_corr),
              .(exposure, nsnp, nsnp_outlier, or_corr, or_lci_corr, or_uci_corr, p_corr, fdr_corr)],
          digits = 3)
  }
}
