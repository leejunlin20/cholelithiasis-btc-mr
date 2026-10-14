#!/usr/bin/env Rscript
# 01_extract_instruments.R —— 从 OpenGWAS 抽取饮食暴露的工具变量（IV）
#
# 用法: Rscript scripts/01_extract_instruments.R <diet_traits.csv> <输出目录>
# 产物: <输出目录>/exposures.rds   (TwoSampleMR 格式的 exposure data frame)
#       <输出目录>/exposure_log.csv (每个数据集的 IV 数量 / 失败原因)

suppressMessages({
  library(TwoSampleMR)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法: Rscript 01_extract_instruments.R <diet_traits.csv> <outdir>")
trait_csv <- args[1]
outdir    <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 0)

tr <- fread(trait_csv)
cat(sprintf("待抽取暴露: %d 个数据集\n", nrow(tr)))

P1  <- 5e-8     # 基因组显著
R2  <- 0.001    # LD clumping
KB  <- 10000    # 窗口

log_rows <- list()
expo_list <- list()

for (i in seq_len(nrow(tr))) {
  id <- tr$id[i]
  cand <- file.path(outdir, "cache", paste0(id, ".rds"))
  dir.create(dirname(cand), showWarnings = FALSE, recursive = TRUE)

  x <- NULL
  if (file.exists(cand)) {
    x <- readRDS(cand)
  } else {
    x <- tryCatch(
      extract_instruments(outcomes = id, p1 = P1, clump = TRUE, r2 = R2, kb = KB),
      error = function(e) { cat(sprintf("  [ERR] %s: %s\n", id, conditionMessage(e))); NULL }
    )
    if (!is.null(x) && nrow(x) > 0) saveRDS(x, cand)
    Sys.sleep(0.3)
  }

  nsnp <- if (is.null(x)) 0L else nrow(x)
  log_rows[[length(log_rows) + 1]] <- data.table(
    id = id, trait = tr$trait[i], sample_size = tr$sample_size[i], nsnp = nsnp,
    status = if (nsnp == 0) "no_iv" else "ok"
  )
  if (nsnp > 0) expo_list[[id]] <- x
  if (i %% 10 == 0) cat(sprintf("  进度 %d/%d（当前: %s，累计有效 %d）\n", i, nrow(tr), id, length(expo_list)))
}

if (length(expo_list) == 0) stop("没有抽到任何 IV，检查 token / 数据集 ID")

expo <- rbindlist(expo_list, fill = TRUE)
expo <- as.data.frame(expo)
saveRDS(expo, file.path(outdir, "exposures.rds"))
fwrite(rbindlist(log_rows), file.path(outdir, "exposure_log.csv"))

cat(sprintf("\n完成: %d/%d 个数据集抽到 IV，合计 %d 条 SNP-暴露记录\n",
            length(expo_list), nrow(tr), nrow(expo)))
cat("IV 数量分布:\n"); print(summary(rbindlist(log_rows)$nsnp))
