#!/usr/bin/env Rscript
# ============================================================
# 59_coloc_eqtl_gwas.R —— eQTL × GWAS 共定位（QTL 层因果检验 · 最强形式）
# ============================================================
# 目的：把现稿四个共定位区（R07/R51/R67/R69）从"区域级 GWAS-GWAS 共定位"
#       推进到"基因表达层":该基因的 cis-eQTL 与胆石症/胆道癌 GWAS 信号
#       是否共享同一因果变异。
#
# 数据
#   dataset1：eQTL Catalogue · GTEx v8 分组织 cis-eQTL（Liver / Whole_Blood）
#             results/qtl/raw/<region>_<tissue>.tsv.gz
#   dataset2：FinnGen R12 区域切片
#             data/coloc/CHOLE_regions.tsv (K11_CHOLELITH 49,834/437,418)
#             data/coloc/BTC_regions.tsv   (C3_BILIARY_GALLBLADDER_EXALLC 2,298/378,749)
#
# 说明：两侧同属 GRCh38；coloc.abf 为单因果变异假设，用 MAF 近似，
#       与现稿 12_coloc.R 的口径一致（现稿亦未使用 LD 矩阵）。
# 先验：p12 敏感性 1e-6 / 1e-5 / 1e-4（R07 的 PP.H4 对先验高度敏感，必须同时报）
# ============================================================

suppressMessages({
  library(coloc)
  library(data.table)
})

BASE     <- Sys.getenv("MR_BASE", unset = ".")
DIR_Q    <- file.path(BASE, "results/qtl")
DIR_RAW  <- file.path(DIR_Q, "raw")
DIR_DATA <- file.path(BASE, "data/coloc")

GENES <- list(
  ABCG5    = list(ensg = "ENSG00000138075", chrom = "2",  tss = 43806213, end = 43838865),
  ABCG8    = list(ensg = "ENSG00000143921", chrom = "2",  tss = 43831942, end = 43882988),
  LRPPRC   = list(ensg = "ENSG00000138095", chrom = "2",  tss = 43886224, end = 43996226),
  SERPINA1 = list(ensg = "ENSG00000197249", chrom = "14", tss = 94376745, end = 94391024),
  HNF4A    = list(ensg = "ENSG00000101076", chrom = "20", tss = 44355439, end = 44432845),
  PNPLA3   = list(ensg = "ENSG00000100344", chrom = "22", tss = 43923761, end = 43964488)
)

REGIONS <- list(
  R07 = list(chrom = "2",  genes = c("ABCG5", "ABCG8", "LRPPRC")),
  R51 = list(chrom = "14", genes = "SERPINA1"),
  R67 = list(chrom = "20", genes = "HNF4A"),
  R69 = list(chrom = "22", genes = "PNPLA3")
)

ORGANS <- c("Liver", "Whole_Blood")

SPECS <- list(
  CHOLE = list(file = "CHOLE_regions.tsv", type = "cc", N = 487252,
               s = 49834 / 487252, label = "Cholelithiasis"),
  BTC   = list(file = "BTC_regions.tsv",   type = "cc", N = 381047,
               s = 2298 / 381047,  label = "Biliary tract cancer")
)

WINDOW <- 1000000L     # cis 窗口 = TSS ±1 Mb（与 58 脚本一致）
P12S   <- c(1e-6, 1e-5, 1e-4)

## ---------------- 读 GWAS 区域切片 ----------------
read_gwas <- function(key) {
  sp <- SPECS[[key]]
  d <- fread(file.path(DIR_DATA, sp$file), sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE,
             colClasses = list(character = "rsids"))
  names(d)[1] <- "region_id"; names(d)[2] <- "chrom"
  d <- d[d$region_id %in% names(REGIONS), ]
  for (cc in c("pos","beta","sebeta","pval","af_alt"))
    d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  d$chrom <- as.character(d$chrom)
  d$maf <- pmin(d$af_alt, 1 - d$af_alt)
  d
}
G <- list(CHOLE = read_gwas("CHOLE"), BTC = read_gwas("BTC"))
for (k in names(G)) cat(sprintf("GWAS %-6s %d 变异\n", k, nrow(G[[k]])))

## ---------------- 读 eQTL ----------------
read_eqtl <- function(region, tissue, gene) {
  f <- file.path(DIR_RAW, sprintf("%s_%s.tsv.gz", region, tissue))
  if (!file.exists(f)) return(NULL)
  d <- read.delim(gzfile(f), sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, check.names = FALSE)
  names(d)[names(d) == "chromosome"] <- "chrom"
  for (cc in c("beta","se","maf","an","pvalue","position"))
    d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  d$chrom <- as.character(d$chrom)
  gi <- GENES[[gene]]
  d <- d[d$gene_id == gi$ensg, ]
  d <- d[d$position >= gi$tss - WINDOW & d$position <= gi$end + WINDOW, ]
  if (nrow(d) == 0) return(NULL)
  d <- d[order(d$pvalue), ]
  d <- d[!duplicated(d$position), ]
  d$chrom <- as.character(d$chrom)
  d
}

## ---------------- 单次 coloc ----------------
run_one <- function(region, tissue, gene, outcome, p12) {
  e <- read_eqtl(region, tissue, gene)
  if (is.null(e) || nrow(e) < 50) return(NULL)
  g <- G[[outcome]][G[[outcome]]$region_id == region, ]
  g <- g[!duplicated(g$pos), ]
  m <- merge(e[, c("chrom","position","ref","alt","beta","se","maf","an","pvalue","variant","rsid")],
             g[, c("chrom","pos","ref","alt","beta","sebeta","pval","maf")],
             by.x = c("chrom","position"), by.y = c("chrom","pos"),
             suffixes = c(".e", ".g"))
  same <- m$ref.e == m$ref.g & m$alt.e == m$alt.g
  flip <- m$ref.e == m$alt.g & m$alt.e == m$ref.g
  m <- m[same | flip, ]
  if (nrow(m) < 50) return(NULL)
  m$beta.g[m$ref.e == m$alt.g & m$alt.e == m$ref.g] <-
    -m$beta.g[m$ref.e == m$alt.g & m$alt.e == m$ref.g]
  m <- m[order(m$position), ]
  m <- m[!duplicated(m$position), ]
  # 严格清洗：MAF ∈ (0,1)、beta/se 有限且 se>0
  ok <- is.finite(m$maf.e) & is.finite(m$maf.g) &
        m$maf.e > 0 & m$maf.e < 1 & m$maf.g > 0 & m$maf.g < 1 &
        is.finite(m$beta.e) & is.finite(m$se) & m$se > 0 &
        is.finite(m$beta.g) & is.finite(m$sebeta) & m$sebeta > 0
  m <- m[ok, ]
  if (nrow(m) < 50) return(NULL)

  n_e <- median(m$an, na.rm = TRUE) / 2
  d1 <- list(snp = paste0(m$chrom, ":", m$position),
             position = m$position, beta = m$beta.e,
             varbeta = m$se^2, MAF = m$maf.e, N = n_e, type = "quant")
  sp <- SPECS[[outcome]]
  d2 <- list(snp = paste0(m$chrom, ":", m$position),
             position = m$position, beta = m$beta.g,
             varbeta = m$sebeta^2, MAF = m$maf.g, N = sp$N, s = sp$s, type = "cc")

  r <- tryCatch(suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2, p12 = p12)),
                error = function(x) { cat(sprintf("    !! coloc 失败 %s/%s/%s/%s/%s: %s\n",
                                                  region, tissue, gene, outcome, p12,
                                                  conditionMessage(x))); NULL })
  if (is.null(r)) return(NULL)
  sm <- r$summary
  top <- r$results$snp[which.max(r$results$SNP.PP.H4)]
  data.frame(region = region, gene = gene, tissue = tissue, outcome = outcome,
             p12 = p12, n_snp = as.integer(sm[["nsnps"]]), n_eqtl = n_e,
             PP.H0 = sm[["PP.H0.abf"]], PP.H1 = sm[["PP.H1.abf"]],
             PP.H2 = sm[["PP.H2.abf"]], PP.H3 = sm[["PP.H3.abf"]],
             PP.H4 = sm[["PP.H4.abf"]], top_shared_snp = top,
             min_p_eqtl = min(m$pvalue), lead_eqtl_snp = m$variant[which.min(m$pvalue)],
             min_p_gwas = min(m$pval), stringsAsFactors = FALSE)
}

## ---------------- 批量 ----------------
rows <- list()
for (region in names(REGIONS)) {
  for (gene in REGIONS[[region]]$genes) {
    for (tissue in ORGANS) {
      for (outcome in names(SPECS)) {
        for (p12 in P12S) {
          r <- run_one(region, tissue, gene, outcome, p12)
          if (!is.null(r)) rows[[length(rows) + 1]] <- r
        }
      }
    }
  }
  cat(sprintf("   ... %s 完成 (%d 行)\n", region, length(rows)))
}

out <- do.call(rbind, rows)
if (!is.null(out) && nrow(out) > 0) {
  out <- out[order(out$region, out$gene, out$tissue, out$outcome, out$p12), ]
  fwrite(out, file.path(DIR_Q, "coloc_eqtl_gwas.csv"))
} else {
  cat("!! 无结果行\n")
}

cat("\n================ eQTL × GWAS 共定位（全部） ================\n")
print(out[, c("region","gene","tissue","outcome","p12","n_snp",
              "PP.H3","PP.H4","top_shared_snp")],
      row.names = FALSE, digits = 3)

cat("\n================ 主口径 p12 = 1e-5 ================\n")
main <- out[out$p12 == 1e-5, ]
print(main[, c("region","gene","tissue","outcome","n_snp","PP.H3","PP.H4")],
      row.names = FALSE, digits = 3)

cat("\n完成 -> results/qtl/coloc_eqtl_gwas.csv\n")
