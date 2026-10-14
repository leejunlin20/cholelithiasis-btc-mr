#!/usr/bin/env Rscript
# ============================================================
# 12_coloc.R  ——  共定位分析（coloc.abf）
# ============================================================
# 目的：检验胆石病、BMI 与胆道癌（BTC）是否由同一因果变异驱动，
#       为 MR 因果链排除 LD 混杂 / 独立信号。
#
# 数据：FinnGen R12（同一 Finnish 人群、同一 GRCh38 build、同一 imputation 面板）
#   CHOLE  K11_CHOLELITH                  49,834 / 437,418  (cc)
#   BMI    BMI_IRN                        500,348           (quant, IRN → sdY≈1)
#   BTC    C3_BILIARY_GALLBLADDER_EXALLC  2,298 / 378,749   (cc)
#
# 区域：胆石病 84 个独立先导 SNP 为中心 ±500kb 合并 → 69 个区域
# 先验：p12 敏感性 1e-6 / 1e-5 / 1e-4（主分析 p12=1e-5）
# ============================================================

suppressMessages({
  library(coloc)
  library(data.table)
})

BASE <- Sys.getenv("MR_BASE", unset = ".")
DIR_COLOC <- file.path(BASE, "results/coloc")
DIR_DATA  <- file.path(BASE, "data/coloc")

## ---------- 样本量（FinnGen R12 manifest，2026-09-17 核对） ----------
SPECS <- list(
  CHOLE = list(file = "CHOLE_regions.tsv", type = "cc",
               N = 487252, s = 49834 / 487252, label = "Cholelithiasis (49,834/437,418)"),
  BTC   = list(file = "BTC_regions.tsv",   type = "cc",
               N = 381047, s = 2298 / 381047,  label = "Biliary tract cancer (2,298/378,749)"),
  BMI   = list(file = "BMI_regions.tsv",   type = "quant",
               N = 500348, s = NULL, sdY = 1,   label = "BMI (IRN, N=500,348)")
)

regions <- fread(file.path(DIR_COLOC, "regions.tsv"), sep = "\t", header = TRUE,
                 data.table = FALSE, check.names = FALSE)
names(regions) <- c("region_id","chrom","start","end","n_lead","lead_snps","top_pval","nearest_genes")

read_trait <- function(key) {
  sp <- SPECS[[key]]
  d <- fread(file.path(DIR_DATA, sp$file), sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE,
             colClasses = list(character = "rsids"))
  names(d)[1] <- "region_id"
  names(d)[2] <- "chrom"
  d$rsids[is.na(d$rsids) | d$rsids == ""] <- "."
  d$snp <- ifelse(d$rsids == ".", paste0(d$chrom, ":", d$pos, "_", d$ref, "_", d$alt), d$rsids)
  # 每个区域内按位置去重（保留 p 最小者）
  d <- d[order(d$pval), ]
  d <- d[!duplicated(paste(d$region_id, d$pos)), ]
  d
}

cat("读入区域数据 ...\n")
DATA <- list(CHOLE = read_trait("CHOLE"), BTC = read_trait("BTC"), BMI = read_trait("BMI"))
for (k in names(DATA)) cat(sprintf("  %-6s %d 变异\n", k, nrow(DATA[[k]])))

## ---------- 通用单区域 coloc ----------
build_dataset <- function(d, sp, pos, idx) {
  out <- list(
    snp      = d$snp[idx],
    position = pos,
    beta     = d$beta[idx],
    varbeta  = d$sebeta[idx]^2,
    type     = sp$type
  )
  if (sp$type == "cc") {
    out$s <- sp$s
    out$N <- sp$N
  } else {
    out$N <- sp$N
    out$sdY <- sp$sdY
  }
  if ("af_alt" %in% names(d)) {
    maf <- pmin(d$af_alt[idx], 1 - d$af_alt[idx])
    if (all(is.finite(maf))) out$MAF <- maf
  }
  out
}

run_coloc <- function(rid, k1, k2, p12) {
  d1 <- DATA[[k1]]; d2 <- DATA[[k2]]
  a <- d1[d1$region_id == rid, ]
  b <- d2[d2$region_id == rid, ]
  if (nrow(a) < 50 || nrow(b) < 50) return(NULL)

  m <- merge(a[, c("chrom","pos","ref","alt","snp","beta","sebeta","pval","af_alt")],
             b[, c("chrom","pos","ref","alt","snp","beta","sebeta","pval","af_alt")],
             by = c("chrom","pos","ref","alt"), suffixes = c(".1", ".2"))
  m <- m[!is.na(m$beta.1) & !is.na(m$beta.2) &
         !is.na(m$sebeta.1) & !is.na(m$sebeta.2) &
         m$sebeta.1 > 0 & m$sebeta.2 > 0, ]
  m <- m[order(m$pos), ]
  if (nrow(m) < 50) return(NULL)

  s1 <- build_dataset(data.frame(snp = paste0(m$chrom,":",m$pos), beta = m$beta.1,
                                 sebeta = m$sebeta.1, af_alt = m$af_alt.1),
                      SPECS[[k1]], m$pos, seq_len(nrow(m)))
  s2 <- build_dataset(data.frame(snp = paste0(m$chrom,":",m$pos), beta = m$beta.2,
                                 sebeta = m$sebeta.2, af_alt = m$af_alt.2),
                      SPECS[[k2]], m$pos, seq_len(nrow(m)))

  res <- tryCatch(suppressWarnings(coloc.abf(dataset1 = s1, dataset2 = s2, p12 = p12)),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  sm <- res$summary
  ranked <- res$results
  top_shared <- ranked$snp[which.max(ranked$SNP.PP.H4)]

  data.frame(
    region_id = rid, chrom = m$chrom[1], start = min(m$pos), end = max(m$pos),
    p12 = p12, nsnps = as.integer(sm[["nsnps"]]),
    PP.H0 = sm[["PP.H0.abf"]], PP.H1 = sm[["PP.H1.abf"]],
    PP.H2 = sm[["PP.H2.abf"]], PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]],
    top_shared_snp = top_shared,
    min_p_1 = min(m$pval.1), top_snp_1 = m$snp.1[which.min(m$pval.1)],
    min_p_2 = min(m$pval.2), top_snp_2 = m$snp.2[which.min(m$pval.2)],
    stringsAsFactors = FALSE
  )
}

## ---------- 批量 ----------
PAIRS <- list(
  c("CHOLE", "BTC"),   # 核心：胆石 ↔ 胆道癌
  c("CHOLE", "BMI"),   # 阳性对照：BMI ↔ 胆石
  c("BMI",   "BTC")    # BMI ↔ 胆道癌
)
P12S <- c(1e-6, 1e-5, 1e-4)

all_rows <- list()
for (pr in PAIRS) {
  nm <- paste0(pr[1], "_vs_", pr[2])
  cat(sprintf("\n>>> %s (%s vs %s)\n", nm, SPECS[[pr[1]]]$label, SPECS[[pr[2]]]$label))
  rows <- list()
  for (i in seq_len(nrow(regions))) {
    rid <- regions$region_id[i]
    for (p12 in P12S) {
      r <- run_coloc(rid, pr[1], pr[2], p12)
      if (!is.null(r)) { r$pair <- nm; rows[[length(rows) + 1]] <- r }
    }
    if (i %% 20 == 0) cat(sprintf("   ... %d/%d 区域\n", i, nrow(regions)))
  }
  out <- do.call(rbind, rows)
  fwrite(out, file.path(DIR_COLOC, paste0("coloc_", nm, ".csv")))
  cat(sprintf("   -> coloc_%s.csv (%d 行)\n", nm, nrow(out)))
  all_rows[[nm]] <- out
}

final <- do.call(rbind, all_rows)
final <- merge(regions[, c("region_id","n_lead","lead_snps","nearest_genes")], final,
               by = "region_id", all.y = TRUE)
final <- final[order(final$pair, final$p12, -final$PP.H4), ]
fwrite(final, file.path(DIR_COLOC, "coloc_all_summary.csv"))
cat(sprintf("\n汇总: coloc_all_summary.csv (%d 行)\n", nrow(final)))

main <- final[final$p12 == 1e-5, ]
fwrite(main, file.path(DIR_COLOC, "coloc_main_p12_1e-5.csv"))

cat("\n================ 主分析 p12 = 1e-5 ================\n")
for (p in unique(main$pair)) {
  sub <- main[main$pair == p, ]
  cat(sprintf("\n[%s]  H4>0.8: %d 区域 | H4>0.5: %d | H3>0.8: %d | 共 %d 区域\n",
              p, sum(sub$PP.H4 > 0.8), sum(sub$PP.H4 > 0.5), sum(sub$PP.H3 > 0.8), nrow(sub)))
  top <- head(sub[order(-sub$PP.H4), c("region_id","chrom","start","end","nearest_genes",
                                       "nsnps","PP.H3","PP.H4","top_shared_snp",
                                       "min_p_1","min_p_2")], 10)
  print(top, row.names = FALSE, digits = 3)
}

cat("\n完成。\n")
