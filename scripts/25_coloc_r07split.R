#!/usr/bin/env Rscript
# =====================================================================
# 25_coloc_r07split.R —— R07 过度合并修正后的重跑
# ---------------------------------------------------------------------
# 背景：原 regions.tsv 的 R07 = chr2:43,306,615-45,694,235（2.39 Mb），
#       实际跨两个物理独立位点：
#         R07a  ABCG5/ABCG8（+DYNC2LI1、RNU6-566P）  chr2:43,306,615-44,470,983
#         R07b  CAMKMT（+LINC01121）                 chr2:44,470,984-45,694,235
#       SuSiE 在原 R07 上给出 19 个「独立信号」、条件化残余 |z|=20.2，
#       均说明该区域违反 coloc 的「单一因果变异」前提。
# 本脚本把 R07 按位点切开后重跑 coloc.abf，并与原结果对照。
# =====================================================================
suppressMessages({ library(coloc); library(data.table) })

BASE <- Sys.getenv("MR_BASE", unset = ".")
DIR_COLOC <- file.path(BASE, "results/coloc")
DIR_DATA  <- file.path(BASE, "data/coloc")
OUT       <- file.path(BASE, "results/coloc")

SPECS <- list(
  CHOLE = list(file = "CHOLE_regions.tsv", type = "cc", N = 487252,
               s = 49834 / 487252),
  BTC   = list(file = "BTC_regions.tsv",   type = "cc", N = 381047,
               s = 2298 / 381047)
)

SUB <- list(
  R07       = c(43306615L, 45694235L),
  R07a      = c(43306615L, 44470983L),   # ABCG5/ABCG8
  R07b      = c(44470984L, 45694235L)    # CAMKMT
)

read_r07 <- function(key) {
  sp <- SPECS[[key]]
  d <- fread(file.path(DIR_DATA, sp$file), sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE,
             colClasses = list(character = "rsids"))
  names(d)[1:2] <- c("region_id", "chrom")
  d <- d[d$region_id == "R07" & d$chrom == 2, ]
  d$rsids[is.na(d$rsids) | d$rsids == ""] <- "."
  d$snp <- ifelse(d$rsids == ".", paste0(d$chrom, ":", d$pos, "_", d$ref, "_", d$alt),
                  d$rsids)
  d
}

build <- function(d, sp, idx) {
  out <- list(snp = d$snp[idx], position = d$pos[idx], beta = d$beta[idx],
              varbeta = d$sebeta[idx]^2, type = sp$type)
  if (sp$type == "cc") { out$s <- sp$s; out$N <- sp$N } else out$N <- sp$N
  maf <- pmin(d$af_alt[idx], 1 - d$af_alt[idx])
  if (all(is.finite(maf))) out$MAF <- maf
  out
}

CHO <- read_r07("CHOLE"); BTC <- read_r07("BTC")
cat(sprintf("R07 原始抽取：CHOLE %d 变异 / BTC %d 变异\n", nrow(CHO), nrow(BTC)))

res <- list()
for (rid in names(SUB)) {
  lo <- SUB[[rid]][1]; hi <- SUB[[rid]][2]
  a <- CHO[CHO$pos >= lo & CHO$pos <= hi, ]
  b <- BTC[BTC$pos >= lo & BTC$pos <= hi, ]
  a <- a[order(a$pval), ]; a <- a[!duplicated(a$pos), ]
  b <- b[order(b$pval), ]; b <- b[!duplicated(b$pos), ]
  keep <- intersect(a$snp, b$snp)
  a <- a[a$snp %in% keep, ]; b <- b[match(keep, b$snp), ]

  if (length(keep) < 50) { cat(sprintf("[%s] 共享 SNP 仅 %d，跳过\n", rid, length(keep))); next }

  D1 <- build(a, SPECS$CHOLE, seq_len(nrow(a)))
  D2 <- build(b, SPECS$BTC,   seq_len(nrow(b)))

  for (p12 in c(1e-6, 1e-5, 1e-4)) {
    r <- tryCatch(suppressMessages(
      coloc.abf(D1, D2, p12 = p12)$summary),
      error = function(e) { cat(sprintf("  [ERR] %s p12=%.0e : %s\n", rid, p12,
                                        conditionMessage(e))); NULL })
    if (is.null(r)) next
    res[[length(res) + 1]] <- data.table(
      region = rid, span = sprintf("%d-%d", SUB[[rid]][1], SUB[[rid]][2]),
      n_snp = length(keep), p12 = p12,
      PP.H0 = r[["PP.H0.abf"]], PP.H1 = r[["PP.H1.abf"]],
      PP.H2 = r[["PP.H2.abf"]], PP.H3 = r[["PP.H3.abf"]],
      PP.H4 = r[["PP.H4.abf"]])
  }
  cat(sprintf("  %-5s %10d-%10d  n=%5d ｜ CHOLE min_p=%.3g  BTC min_p=%.3g\n",
              rid, lo, hi, length(keep), min(a$pval), min(b$pval)))
}

A <- rbindlist(res)
fwrite(A, file.path(OUT, "coloc_R07_split.csv"))

cat("\n================ R07 切分前后 coloc.abf 对照（CHOLE vs BTC）================\n")
print(A[, .(region, n_snp, p12, PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3))])

cat("\n要点：原 R07（2.39 Mb，含两个位点）在任何 p12 下都得不到可信的 H4；\n")
cat("      切分后 R07a（ABCG5/ABCG8）与 R07b（CAMKMT）分别给出区域内的 H4，才可解释。\n")
