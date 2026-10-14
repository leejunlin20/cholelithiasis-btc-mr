#!/usr/bin/env Rscript
# ============================================================
# 12b_coloc_window_sensitivity.R
# ============================================================
# 质疑点：1 Mb 区域内可能含多个独立信号，coloc.abf 的 H4 是否被"平均掉"？
# 对策  ：把窗口收窄到 ±100 kb / ±250 kb（以区域内胆石病最强 SNP 为中心）
#          重跑 coloc，比较 PP.H4 稳定性。
# 输出  : results/coloc/coloc_window_sensitivity.csv
# ============================================================

suppressMessages({ library(coloc); library(data.table) })

BASE <- Sys.getenv("MR_BASE", unset = ".")
DIR_COLOC <- file.path(BASE, "results/coloc")
DIR_DATA  <- file.path(BASE, "data/coloc")

SPECS <- list(
  CHOLE = list(file = "CHOLE_regions.tsv", type = "cc", N = 487252, s = 49834 / 487252),
  BTC   = list(file = "BTC_regions.tsv",   type = "cc", N = 381047, s = 2298 / 381047),
  BMI   = list(file = "BMI_regions.tsv",   type = "quant", N = 500348, sdY = 1)
)

read_trait <- function(key) {
  sp <- SPECS[[key]]
  d <- fread(file.path(DIR_DATA, sp$file), sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE,
             colClasses = list(character = "rsids"))
  names(d)[1] <- "region_id"; names(d)[2] <- "chrom"
  d$rsids[is.na(d$rsids) | d$rsids == ""] <- "."
  # 注意：必须分两步。写成 d[order(p),][!duplicated(...),] 会让逻辑索引错位，
  # 多等位位点不会被去掉，coloc 会报 "duplicated snps found"。
  d <- d[order(d$pval), ]
  d <- d[!duplicated(paste(d$region_id, d$pos)), ]
  d
}

DATA <- lapply(names(SPECS), read_trait); names(DATA) <- names(SPECS)

run <- function(rid, k1, k2, halfwin, p12 = 1e-5) {
  a <- DATA[[k1]][DATA[[k1]]$region_id == rid, ]
  if (!nrow(a)) return(NULL)
  ctr <- a$pos[which.min(a$pval)]
  a <- a[abs(a$pos - ctr) <= halfwin, ]
  b <- DATA[[k2]][DATA[[k2]]$region_id == rid, ]
  b <- b[abs(b$pos - ctr) <= halfwin, ]
  if (nrow(a) < 30 || nrow(b) < 30) return(NULL)

  m <- merge(a[, c("chrom","pos","ref","alt","beta","sebeta","pval")],
             b[, c("chrom","pos","ref","alt","beta","sebeta","pval")],
             by = c("chrom","pos","ref","alt"), suffixes = c(".1", ".2"))
  m <- m[order(m$pos), ]
  if (nrow(m) < 30) return(NULL)

  mk <- function(beta, se, sp) {
    l <- list(snp = paste0(m$chrom, ":", m$pos), position = m$pos,
              beta = beta, varbeta = se^2, type = sp$type, N = sp$N)
    if (sp$type == "cc") l$s <- sp$s else l$sdY <- sp$sdY
    l
  }
  D1 <- mk(m$beta.1, m$sebeta.1, SPECS[[k1]])
  D2 <- mk(m$beta.2, m$sebeta.2, SPECS[[k2]])
  res <- tryCatch(suppressWarnings(coloc.abf(D1, D2, p12 = p12)),
                  error = function(e) { message(sprintf("  [ERR] %s %gkb: %s", rid, halfwin/1000, conditionMessage(e))); NULL })
  if (is.null(res)) return(NULL)
  sm <- res$summary
  data.frame(region_id = rid, pair = paste0(k1, "_vs_", k2),
             halfwin_kb = halfwin / 1000, center = ctr, nsnps = as.integer(sm[["nsnps"]]),
             PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]],
             min_p_1 = min(m$pval.1), min_p_2 = min(m$pval.2), stringsAsFactors = FALSE)
}

regions <- fread(file.path(DIR_COLOC, "regions.tsv"), sep = "\t", header = TRUE,
                 data.table = FALSE, check.names = FALSE)
regions <- regions[[1]]
PAIRS <- list(c("CHOLE","BTC"), c("CHOLE","BMI"), c("BMI","BTC"))
WINS <- c(100e3, 250e3, 500e3)

rows <- list()
for (pr in PAIRS) {
  cat(sprintf(">>> %s_vs_%s\n", pr[1], pr[2]))
  for (rid in regions) for (w in WINS) {
    r <- run(rid, pr[1], pr[2], w)
    if (!is.null(r)) rows <- c(rows, list(r))
  }
}
out <- do.call(rbind, rows)
fwrite(out, file.path(DIR_COLOC, "coloc_window_sensitivity.csv"))
cat(sprintf("\n输出: coloc_window_sensitivity.csv (%d 行)\n", nrow(out)))

wide <- reshape(out[, c("region_id","pair","halfwin_kb","PP.H4")],
                idvar = c("region_id","pair"), timevar = "halfwin_kb", direction = "wide")
names(wide) <- gsub("PP.H4\\.", "H4_w", names(wide))
cat("\n==== CHOLE_vs_BTC：窄窗 vs 宽窗 PP.H4（前 12，按 ±100kb 排序）====\n")
s <- wide[wide$pair == "CHOLE_vs_BTC", ]
s <- s[order(-s$H4_w100), ]
print(head(s, 12), row.names = FALSE, digits = 3)

cat("\n==== 一致性统计（H4>0.8 的区域在各窗口下是否保持）====\n")
for (p in unique(wide$pair)) {
  sb <- wide[wide$pair == p, ]
  a <- sum(sb$H4_w100 > 0.8); b <- sum(sb$H4_w250 > 0.8); c3 <- sum(sb$H4_w500 > 0.8)
  cat(sprintf("  %-14s ±100kb: %d | ±250kb: %d | ±500kb: %d\n", p, a, b, c3))
}
