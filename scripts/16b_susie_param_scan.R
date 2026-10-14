# 16b_susie_param_scan.R —— susie 参数扫描
#
# 目的：coloc.susie 要求**两侧都有 credible set**。R51(SERPINA1)/R67(HNF4A) 的
#       胆道癌侧 min p 仅 4.9e-6 / 2.0e-5（未达全基因组显著），在 susie 默认
#       min_abs_corr=0.5 下建不出 CS，导致只有 R69 能算。
#
#       本脚本扫描 min_abs_corr / coverage / L / estimate_residual_variance，
#       判定弱信号侧能否在不放松「同一 CS 内纯度」这一实质约束的前提下恢复。
#       min_abs_corr 是 CS 内部两两 |r| 的下限（纯度），文献中处理弱信号时常降至
#       0.1~0.2；它不是「让结果变好看的旋钮」，而是显式声明 CS 允许的最大异质性。
#
# 产出 results/coloc_susie/susie_param_scan.tsv

suppressMessages({
  library(coloc)
  library(data.table)
})

BASE <- Sys.getenv("MR_BASE", unset = ".")
OUT  <- file.path(BASE, "results/coloc_susie")

REGIONS <- list(R69 = "PNPLA3", R51 = "SERPINA1", R67 = "HNF4A",
                R07 = "ABCG5/ABCG8", R54 = "FTO")
TRAITS <- list(
  CHOLE = list(N = 49834 + 437418, s = 49834 / (49834 + 437418)),
  BTC   = list(N = 2298 + 378749,  s = 2298 / (2298 + 378749))
)

rd_gz <- function(f) as.data.table(read.delim(gzfile(f), header = TRUE,
                                              check.names = FALSE, stringsAsFactors = FALSE))
ld_long_path <- function(rid) {
  p <- file.path("/tmp/ld", sprintf("%s_ld_long.tsv", rid))
  if (file.exists(p)) p else file.path(BASE, "data/ld", sprintf("%s_ld_long.tsv.gz", rid))
}
build_ld <- function(rid, snp_names) {
  f <- ld_long_path(rid)
  d <- fread(f, header = FALSE, col.names = c("i1", "i2", "r", "r2"))
  d <- d[i1 != i2]
  n <- length(snp_names); R <- diag(1, n)
  ii <- cbind(d$i1 + 1L, d$i2 + 1L)
  R[ii] <- d$r; R[ii[, c(2, 1)]] <- d$r
  R <- (R + t(R)) / 2; diag(R) <- 1
  dimnames(R) <- list(snp_names, snp_names)
  R
}
ensure_psd <- function(R) {
  dn <- dimnames(R); e <- eigen(R, symmetric = TRUE); ev <- e$values
  neg <- ev < 1e-6
  if (!any(neg)) return(R)
  R <- R - e$vectors[, neg, drop = FALSE] %*% (ev[neg] * t(e$vectors[, neg, drop = FALSE]))
  dd <- sqrt(diag(R)); R <- R / outer(dd, dd); diag(R) <- 1; dimnames(R) <- dn
  R
}

rows <- list()
for (rid in names(REGIONS)) {
  snps <- rd_gz(file.path(BASE, "data/ld", sprintf("%s_snps.tsv.gz", rid)))
  n <- nrow(snps)
  R <- build_ld(rid, snps$snp)
  flip <- ifelse(snps$af < 0.5, 1, -1)
  R <- R * outer(flip, flip); diag(R) <- 1
  dimnames(R) <- list(snps$snp, snps$snp)
  R <- ensure_psd(R)
  cat(sprintf("\n=== %s %s (n=%d, 翻转 %d) ===\n", rid, REGIONS[[rid]], n, sum(flip < 0)))

  for (tw in names(TRAITS)) {
    tr <- TRAITS[[tw]]
    which <- if (tw == "CHOLE") 1 else 2
    d <- list(snp = snps$snp, beta = snps[[paste0("beta_", which)]],
              varbeta = snps[[paste0("se_", which)]]^2, type = "cc",
              s = tr$s, N = tr$N, MAF = snps$maf, LD = R, position = snps$pos)

    for (L in c(1L, 5L, 10L)) {
      for (mac in c(0.5, 0.2, 0.1)) {
        for (erv in c(FALSE, TRUE)) {
          S <- tryCatch(suppressWarnings(
                 runsusie(d, L = L, maxit = 200, estimate_residual_variance = erv,
                          min_abs_corr = mac, coverage = 0.95)),
                 error = function(e) NULL)
          if (is.null(S)) { rows[[length(rows) + 1]] <- data.table(
            region = rid, trait = tw, L = L, min_abs_corr = mac,
            est_resid_var = erv, n_cs = NA_integer_, cs_max_pip = NA_real_,
            cs_size = NA_integer_, err = "error"); next }
          ncs <- length(S$sets$cs)
          mp  <- if (ncs > 0) max(S$pip) else NA_real_
          sz  <- if (ncs > 0) length(S$sets$cs[[1]]) else NA_integer_
          rows[[length(rows) + 1]] <- data.table(
            region = rid, trait = tw, L = L, min_abs_corr = mac,
            est_resid_var = erv, n_cs = ncs, cs_max_pip = mp,
            cs_size = sz, err = "")
        }
      }
    }
    cat(sprintf("  %s: min p = %.2e\n", tw, min(snps[[paste0("p_", which)]])))
  }
}

DF <- rbindlist(rows)
fwrite(DF, file.path(OUT, "susie_param_scan.tsv"), sep = "\t")

cat("\n\n========== 汇总：各区域首个能建出 CS 的最小 min_abs_corr ==========\n")
for (tw in c("CHOLE", "BTC")) {
  cat(sprintf("\n---- %s ----\n", tw))
  z <- DF[trait == tw & est_resid_var == FALSE & n_cs > 0,
          .(n_cs_min = min(n_cs), n_cs_max = max(n_cs)),
          by = .(region, min_abs_corr)]
  print(dcast(z, region ~ min_abs_corr, value.var = "n_cs_max"))
}
cat("\n已写出", file.path(OUT, "susie_param_scan.tsv"), "\n")
