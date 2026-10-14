# 16_coloc_susie.R —— 多信号共定位（coloc.susie）
#
# 目的：回应「区域内可能含多个独立信号，单因果 coloc.abf 会把 H3/H4 判错」这一
#       必然会被审稿人提出的质疑。
#
# 输入：
#   data/ld/{RID}_snps.tsv.gz    目标 SNP + 两性状 beta/se/p + MAF（15_ 脚本产出）
#   data/ld/{RID}_ld_long.tsv.gz 稀疏 LD（idx1/idx2/r/r2，FinnGen R12 官方
#                                 520,210 例 Finnish 样本同源面板）
# 输出（results/coloc_susie/）：
#   {RID}_susie_summary.csv   信号对级别的 PP.H0-H4
#   {RID}_susie_signals.csv   每性状 fine-mapping 出的独立信号（CS + PIP）
#   coloc_susie_all.csv       汇总
#   coloc_susie_log.txt       运行日志

suppressMessages({
  library(coloc)
  library(data.table)
  library(Matrix)
})

BASE <- Sys.getenv("MR_BASE", unset = ".")
OUT  <- file.path(BASE, "results/coloc_susie")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

REGIONS <- list(
  R69 = "PNPLA3",
  R51 = "SERPINA1",
  R67 = "HNF4A",
  R07 = "ABCG5/ABCG8 (negative control)",
  R54 = "FTO (positive control)"
)

# FinnGen R12 manifest：CHOLE 49,834/437,418；BTC 2,298/378,749
TRAITS <- list(
  CHOLE = list(N = 49834 + 437418, s = 49834 / (49834 + 437418), label = "Cholelithiasis"),
  BTC   = list(N = 2298 + 378749,  s = 2298 / (2298 + 378749),   label = "Biliary tract cancer")
)

LOG <- file.path(OUT, "coloc_susie_log.txt")
logmsg <- function(...) cat(sprintf(...), "\n", file = LOG, append = TRUE)
if (file.exists(LOG)) file.remove(LOG)

# fread 读 gz 需要 R.utils；用 base R 的 gzfile 更省事
rd_gz <- function(f) as.data.table(read.delim(gzfile(f), header = TRUE,
                                              check.names = FALSE, stringsAsFactors = FALSE))

# LD long 表：优先用 /tmp 下的未压缩副本（fread 比 read.delim(gzfile) 快 1~2 个数量级）
ld_long_path <- function(rid) {
  p_tmp <- file.path("/tmp/ld", sprintf("%s_ld_long.tsv", rid))
  if (file.exists(p_tmp)) return(p_tmp)
  file.path(BASE, "data/ld", sprintf("%s_ld_long.tsv.gz", rid))
}
rd_ld <- function(f) {
  cols <- c("i1", "i2", "r", "r2")
  if (grepl("\\.gz$", f)) {
    as.data.table(read.delim(gzfile(f), header = FALSE, col.names = cols))
  } else {
    fread(f, header = FALSE, col.names = cols)
  }
}

# ------------------------------------------------------------------ LD 方阵
build_ld <- function(rid, snp_names) {
  f <- ld_long_path(rid)
  if (!file.exists(f)) stop("缺少 LD 文件: ", f)
  d <- rd_ld(f)
  d <- d[i1 != i2]
  n <- length(snp_names)
  R <- diag(1, n)
  ii <- cbind(d$i1 + 1L, d$i2 + 1L)
  R[ii] <- d$r
  R[ii[, c(2, 1)]] <- d$r
  # 对称化（双向记录可能有极小浮点差；有重复时取先出现的值）
  R <- (R + t(R)) / 2
  diag(R) <- 1
  # ⚠️ coloc 的 check_dataset 会校验 "colnames in LD do not contain all SNPs"
  dimnames(R) <- list(snp_names, snp_names)
  list(R = R, npairs = nrow(d))
}

# LD 裁剪：剔除与更显著 SNP 高度冗余的位点（FinnGen imputation 变异极密集，
# 实测 R69 区域 2002 个 SNP 中有 1016 个与更显著位点 |r|>0.99）
prune_ld <- function(R, pvals, th = 0.99) {
  n <- nrow(R)
  ord <- order(pvals)
  dropped <- rep(FALSE, n)
  keep <- integer(0)
  for (i in ord) {
    if (dropped[i]) next
    keep <- c(keep, i)
    hi <- which(abs(R[i, ]) > th)
    dropped[hi] <- TRUE
  }
  sort(keep)
}

make_dataset <- function(snps, which, R) {
  tr <- if (which == 1) TRAITS$CHOLE else TRAITS$BTC
  b  <- snps[[paste0("beta_", which)]]
  se <- snps[[paste0("se_", which)]]
  list(snp = snps$snp, beta = b, varbeta = se^2, type = "cc",
       s = tr$s, N = tr$N, MAF = snps$maf, LD = R, position = snps$pos)
}

# 确保 LD 正定。稀疏截断(r²<=0.01 记 0) + 变异冗余会让子矩阵出现负特征值，
# 而 susie 的似然推导要求 PSD。做法：只把负特征值抬到 0，再重新标定为相关矩阵。
# 比 Matrix::nearPD 快得多（n=2000 时数秒 vs 数分钟），且改动的只是原先
# 数学上不可能存在的负方向。
ensure_psd <- function(R, rid) {
  dn <- dimnames(R)
  e <- eigen(R, symmetric = TRUE)
  ev <- e$values
  neg <- ev < 1e-6
  if (!any(neg)) return(list(R = R, fixed = FALSE, min_ev = min(ev), n_neg = 0L))
  V <- e$vectors[, neg, drop = FALSE]
  R <- R - V %*% (ev[neg] * t(V))          # 减去负特征值方向
  dd <- sqrt(diag(R))
  R <- R / outer(dd, dd)
  diag(R) <- 1
  dimnames(R) <- dn
  ev2 <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  list(R = R, fixed = TRUE, min_ev = min(ev), n_neg = sum(neg),
       min_ev_after = min(ev2))
}

# ------------------------------------------------------------------ 主流程
all_summary <- list()
all_signals <- list()
all_feas    <- list()   # 可行性表：每个区域两侧各几个 CS、能否跑 coloc.susie

for (rid in names(REGIONS)) {
  gene <- REGIONS[[rid]]
  cat(sprintf("\n============ %s  %s ============\n", rid, gene))
  logmsg("=== %s (%s) ===", rid, gene)

  sf <- file.path(BASE, "data/ld", sprintf("%s_snps.tsv.gz", rid))
  if (!file.exists(sf)) { cat("  跳过：无 SNP 文件\n"); next }
  snps <- rd_gz(sf)
  n <- nrow(snps)
  cat(sprintf("  SNP 数: %d  区间: %s:%s-%s (%.0f kb)\n",
              n, snps$chrom[1], min(snps$pos), max(snps$pos),
              (max(snps$pos) - min(snps$pos)) / 1000))
  logmsg("SNP 数 %d", n)
  if (n < 50) { cat("  SNP 太少，跳过\n"); next }

  ld <- build_ld(rid, snps$snp)
  cat(sprintf("  LD 配对: %s 条 (占两两配对的 %.1f%%，其余因 r²<=0.01 截断记 0)\n",
              format(ld$npairs, big.mark = ","),
              100 * ld$npairs / (n * (n - 1) / 2)))

  # ★★ 关键：等位基因方向修正 ★★
  # FinnGen LD 矩阵(PLINK2 --r-unphased)的 r 符号以 **minor allele** 为参考，
  # 而 summary stats 的 beta 以 **alt** 等位基因为参考。当 af_alt > 0.5（alt 是主
  # 等位基因）时，两者方向相反 ⇒ r 符号必须翻转，否则 susie 会把与 lead 正相关的
  # 位点误读为负相关，从而分裂出大量假信号。
  # 实测(R69)：修正前 susie 分出 9 个"独立"信号且 PP.H4≈0；修正后收敛到 1 个信号、
  # PP.H4=0.994，与 coloc.abf 的 0.997 一致。
  flip <- ifelse(snps$af < 0.5, 1, -1)
  Rfix <- ld$R * outer(flip, flip)      # 相似变换，不改特征值，但改 r 与 z 的符号关系
  diag(Rfix) <- 1
  dimnames(Rfix) <- list(snps$snp, snps$snp)
  cat(sprintf("  等位基因方向修正: %d/%d (%.1f%%) 位点翻转 (af_alt>=0.5)\n",
              sum(flip < 0), n, 100 * mean(flip < 0)))
  logmsg("等位基因方向修正: %d 个位点翻转", sum(flip < 0))

  t0 <- Sys.time()
  psd <- ensure_psd(Rfix, rid)
  if (psd$fixed) {
    cat(sprintf("  LD 非正定: %d 个负特征值 (min %.2f) → 特征值截断修正 (min %.2g), %.1fs\n",
                psd$n_neg, psd$min_ev, psd$min_ev_after,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    logmsg("LD 非正定: %d 负特征值 min=%.3g → 修正后 min=%.3g",
           psd$n_neg, psd$min_ev, psd$min_ev_after)
  }

  D1 <- make_dataset(snps, 1, psd$R)
  D2 <- make_dataset(snps, 2, psd$R)

  # 显式 runsusie，便于取回 fine-mapping 的 credible set 与 PIP
  # L=5：修正到位时 ELBO 在 L=1 即饱和（增量≈0），L=5 足以发现真实次信号
  MAXL <- 5
  t0 <- Sys.time()
  S1 <- tryCatch(runsusie(D1, L = MAXL, maxit = 200, estimate_residual_variance = FALSE),
                 error = function(e) { cat("  [runsusie-1 ERROR]", conditionMessage(e), "\n"); NULL })
  cat(sprintf("  runsusie 胆石: %.0fs, %d 个 credible set\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), length(S1$sets$cs)))
  t0 <- Sys.time()
  S2 <- tryCatch(runsusie(D2, L = MAXL, maxit = 200, estimate_residual_variance = FALSE),
                 error = function(e) { cat("  [runsusie-2 ERROR]", conditionMessage(e), "\n"); NULL })
  cat(sprintf("  runsusie BTC : %.0fs, %d 个 credible set\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), length(S2$sets$cs)))
  if (is.null(S1) || is.null(S2)) { logmsg("runsusie 失败"); next }

  # ---- 先记录两侧的细定位结果（与能否共定位无关，必须独立保存）----
  # 注意：coloc.susie 在任一侧无 credible set 时会直接失败，但该侧的 fine-mapping
  # 结果本身仍有价值（用于说明「该区域有几个独立信号」「信号能否被定位」），
  # 因此提取必须放在 coloc.susie 之前，不能用 next 跳过。
  #
  # 同时导出逐 SNP 的 PIP 表，供成图使用。
  pip_tab <- data.table(snp = snps$snp, chrom = snps$chrom, pos = snps$pos,
                        p_1 = snps$p_1, p_2 = snps$p_2,
                        pip_1 = S1$pip, pip_2 = S2$pip)
  pip_tab[, in_cs_1 := FALSE][, in_cs_2 := FALSE]
  for (k in seq_along(S1$sets$cs)) pip_tab$in_cs_1[S1$sets$cs[[k]]] <- TRUE
  for (k in seq_along(S2$sets$cs)) pip_tab$in_cs_2[S2$sets$cs[[k]]] <- TRUE
  fwrite(pip_tab, file.path(OUT, sprintf("%s_pip.csv", rid)))
  logmsg("PIP 表: 胆石 %d CS / 胆道癌 %d CS", length(S1$sets$cs), length(S2$sets$cs))

  for (which in 1:2) {
    S <- if (which == 1) S1 else S2
    st <- S$sets
    if (is.null(st) || is.null(st$cs) || length(st$cs) == 0) next
    for (k in seq_along(st$cs)) {
      idx <- st$cs[[k]]
      pv <- S$pip
      all_signals[[length(all_signals) + 1]] <- data.table(
        region = rid, gene = gene,
        trait = if (which == 1) "Cholelithiasis" else "Biliary tract cancer",
        cs = k,
        cs_size = length(idx),
        cs_lead_snp = snps$snp[idx[which.max(pv[idx])]],
        cs_lead_pip = max(pv[idx]),
        cs_lead_pos = snps$pos[idx[which.max(pv[idx])]],
        cs_lead_z = snps[[paste0("beta_", which)]][idx[which.max(pv[idx])]] /
                    snps[[paste0("se_", which)]][idx[which.max(pv[idx])]],
        region_max_pip = max(pv),
        region_lead_snp = snps$snp[which.max(pv)],
        min_p_in_cs = min(snps[[paste0("p_", which)]][idx]),
        cs_snps = paste(snps$snp[idx], collapse = ","),
        cs_pos = paste(snps$pos[idx], collapse = ",")
      )
    }
  }

  res <- tryCatch(coloc.susie(S1, S2),
                  error = function(e) { cat("  [coloc.susie ERROR]", conditionMessage(e), "\n"); e })
  if (inherits(res, "error")) { logmsg("coloc.susie 失败: %s", conditionMessage(res)); next }

  # ⚠️ coloc.susie 在两侧都有 credible set 时返回 list(summary, results, priors)；
  #    任一侧无 CS 时只返回 data.table(nsnps=NA)，没有 $summary。
  #    已实测 coloc 5.2.3 的 summary 列名：
  #    nsnps | hit1 | hit2 | PP.H0.abf | PP.H1.abf | PP.H2.abf | PP.H3.abf | PP.H4.abf | idx1 | idx2
  if (is.null(res$summary)) {
    cat("  无信号对（至少一侧无 credible set）→ 无法检验共享因果变异\n")
    logmsg("无信号对：任一侧 credible set 为空 (胆石 %d / 胆道癌 %d 个 CS)",
           length(S1$sets$cs), length(S2$sets$cs))
    all_feas[[rid]] <- data.table(
      region = rid, gene = gene, n_snp = n,
      n_cs_chole = length(S1$sets$cs), n_cs_btc = length(S2$sets$cs),
      min_p_chole = min(snps$p_1), min_p_btc = min(snps$p_2),
      coloc_susie_ok = FALSE, n_signal_pairs = 0L,
      PP.H3_best = NA_real_, PP.H4_best = NA_real_)
    next
  }
  sm <- as.data.table(res$summary)
  # 列名自适应（若将来 coloc 版本去掉 .abf 后缀也能跑）
  setnames(sm,
           intersect(c("PP.H3.abf", "PP.H3"), names(sm)),
           "PP.H3", skip_absent = TRUE)
  setnames(sm,
           intersect(c("PP.H4.abf", "PP.H4"), names(sm)),
           "PP.H4", skip_absent = TRUE)
  if (!all(c("PP.H3", "PP.H4") %in% names(sm)))
    stop("summary 缺少 PP.H3/PP.H4 列, 实际列: ", paste(names(sm), collapse = ","))
  gene_name <- gene
  sm[, region := rid]
  sm[, gene := gene_name]
  sm[, n_cs_1 := length(S1$sets$cs)]
  sm[, n_cs_2 := length(S2$sets$cs)]
  cat(sprintf("  独立信号: 胆石 %d 个 / 胆道癌 %d 个\n",
              length(unique(sm$idx1)), length(unique(sm$idx2))))
  print(sm[, .(region, hit1, hit2,
               PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4))])
  logmsg("信号对 %d 行", nrow(sm))
  all_summary[[rid]] <- sm
  all_feas[[rid]] <- data.table(
    region = rid, gene = gene, n_snp = n,
    n_cs_chole = length(S1$sets$cs), n_cs_btc = length(S2$sets$cs),
    min_p_chole = min(snps$p_1), min_p_btc = min(snps$p_2),
    coloc_susie_ok = TRUE, n_signal_pairs = nrow(sm),
    PP.H3_best = min(sm$PP.H3), PP.H4_best = max(sm$PP.H4))

  fwrite(sm, file.path(OUT, sprintf("%s_susie_summary.csv", rid)))
}

if (length(all_summary)) {
  A <- rbindlist(all_summary, fill = TRUE)
  fwrite(A, file.path(OUT, "coloc_susie_all.csv"))
  cat("\n\n================ 汇总（PP.H4 > 0.5 的信号对）================\n")
  print(A[PP.H4 > 0.5, .(region, gene, hit1, hit2,
                         PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3))])
  cat("\n---- 全部信号对 ----\n")
  print(A[, .(region, gene, hit1, hit2,
              PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3))])
} else {
  cat("\n\n[注意] 没有任何区域能跑通 coloc.susie（两侧均需 credible set）\n")
}
if (length(all_feas)) {
  F <- rbindlist(all_feas, fill = TRUE)
  fwrite(F, file.path(OUT, "coloc_susie_feasibility.csv"))
  cat("\n\n================ 可行性表（为何某些区域跑不了 coloc.susie）================\n")
  print(F[, .(region, gene, n_snp, n_cs_chole, n_cs_btc,
              min_p_chole = signif(min_p_chole, 2),
              min_p_btc = signif(min_p_btc, 2),
              coloc_susie_ok, PP.H4_best = round(PP.H4_best, 3))])
}
if (length(all_signals)) {
  fwrite(rbindlist(all_signals, fill = TRUE), file.path(OUT, "coloc_susie_signals.csv"))
}
cat("\n完成。输出目录:", OUT, "\n")
