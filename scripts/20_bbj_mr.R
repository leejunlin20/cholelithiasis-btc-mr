#!/usr/bin/env Rscript
# =====================================================================
# BBJ（东亚）跨队列复制：两样本 MR + 阴性/阳性对照
# ---------------------------------------------------------------------
# 设计（全部为「跨队列、无样本重叠」）
#   暴露 = FinnGen R12 胆石症（49,834 例 / 437,418 对照, 芬兰）
#   结局 = BBJ（Sakaue & Kanai 2021 Nat Genet, 日本）
#     · BtC            418 例  ← 主结局（肝内胆管癌）
#     · Urolithiasis 11,699 例  ← 阴性对照（另一种结石病）
#     · PaC            499 例  ← 阴性对照（邻近器官肿瘤）
#     · Cholecystitis 1,860 例  ← 阳性对照（胆囊炎）
#     · Cholelithiasis 9,305 例 ← 阳性对照（同表型、不同血统）
#
# 结局数据来源分两级
#   ① BtC：从官方全量文件 hum0197.v3.BBJ.BtC.v1.zip 直接抽取（Allele1=ref,
#      Allele2=alt, BETA 以 Allele2 计），按 rsID 精确匹配，零猜测 —— 权威。
#   ② 其余表型：PheWeb 单点 API（已用响应体 rsids 字段回验位点身份，
#      并按 bbj_alt 对效应等位基因做方向对齐）。
# =====================================================================
suppressMessages({library(data.table); library(TwoSampleMR)})

BASE <- Sys.getenv("MR_BASE", unset = ".")
OUT  <- file.path(BASE, "results/bbj_replication")
LOG  <- file.path(OUT, "bbj_mr_log.txt")
logmsg <- function(...) cat(sprintf(...), "\n", file = LOG, append = TRUE)

BBJ <- list(
  BtC            = list(ncase = 418,   nctrl = 159201, label = "肝内胆管癌(主结局)"),
  Cholelithiasis = list(ncase = 9305,  nctrl = 168253, label = "胆石症(阳性对照)"),
  Cholecystitis  = list(ncase = 1860,  nctrl = NA,     label = "胆囊炎(阳性对照)"),
  Urolithiasis   = list(ncase = 11699, nctrl = 167027, label = "尿路结石(阴性对照)"),
  PaC            = list(ncase = 499,   nctrl = 159201, label = "胰腺癌(阴性对照)"),
  GB_polyp       = list(ncase = 843,   nctrl = 167459, label = "胆囊息肉")
)

# ---------------- 暴露侧工具变量 ----------------
iv <- fread(file.path(BASE, "results/medi/chole_iv_clumped.tsv"))
setnames(iv, c("rsid","chr","pos","ref","alt","beta","se","pval","eaf","nearest_gene"))
iv[, `:=`(effect_allele = alt, other_allele = ref)]
iv <- iv[!is.na(beta) & !is.na(se) & se > 0 & !is.na(eaf)]
iv[, fstat := (beta / se)^2]
iv[, r2 := 2 * eaf * (1 - eaf) * beta^2]
cat(sprintf("暴露工具变量 %d 个 | 平均 F=%.1f 最小 F=%.1f | R2=%.4f (%.2f%%)\n",
            nrow(iv), mean(iv$fstat), min(iv$fstat), sum(iv$r2), 100 * sum(iv$r2)))
logmsg("暴露: %d SNP, R2=%.4f", nrow(iv), sum(iv$r2))
R2 <- sum(iv$r2)

# ---------------- 结局侧 ----------------
btc <- fread(file.path(OUT, "bbj_btc_instruments.tsv"))
btcOut <- btc[, .(rsid, beta_out = beta_aligned, se_out = bbj_se,
                  p_out = bbj_p, af_out = bbj_af)]

api <- fread(file.path(OUT, "bbj_lookup_instr.tsv"))
api[, `:=`(beta_out = beta_aligned, se_out = se_bbj, p_out = p_bbj, af_out = af_bbj)]
apiU <- unique(api, by = c("rsid", "pheno"))
cat(sprintf("结局侧：BtC 全量文件命中 %d SNP | API 表覆盖 %d 个 SNP/%d 表型\n",
            nrow(btcOut), uniqueN(apiU$rsid), uniqueN(apiU$pheno)))

results <- list(); detail <- list()
for (ph in names(BBJ)) {
  cfg <- BBJ[[ph]]
  o <- if (ph == "BtC") copy(btcOut)
       else apiU[pheno == ph, .(rsid, beta_out, se_out, p_out, af_out)]
  o <- unique(o, by = "rsid")
  m <- merge(iv[, .(rsid, effect_allele, other_allele, beta, se, pval, eaf, nearest_gene, fstat)],
             o, by = "rsid")
  m <- m[!is.na(beta_out) & !is.na(se_out) & se_out > 0]
  if (nrow(m) < 5) { cat(sprintf("[%s] 可用 SNP 仅 %d 个，跳过\n", ph, nrow(m))); next }

  # 方向已在两端统一到 effect_allele，此处直接构造，不依赖 harmonise 的链向猜测
  d <- data.table(
    SNP = m$rsid, beta.exposure = m$beta, se.exposure = m$se,
    pval.exposure = m$pval, eaf.exposure = m$eaf,
    effect_allele.exposure = m$effect_allele, other_allele.exposure = m$other_allele,
    id.exposure = "Cholelithiasis", exposure = "Cholelithiasis (FinnGen R12)",
    beta.outcome = m$beta_out, se.outcome = m$se_out, pval.outcome = m$p_out,
    eaf.outcome = m$af_out,
    effect_allele.outcome = m$effect_allele, other_allele.outcome = m$other_allele,
    id.outcome = ph, outcome = cfg$label, mr_keep = TRUE)

  res <- tryCatch(as.data.table(mr(
                    d, method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression",
                                       "mr_weighted_median", "mr_weighted_mode"))),
                  error = function(e) { cat("[MR ERR] ", ph, conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  # TwoSampleMR 0.7.9 的 mr() 只返回 9 列（无 lo_ci/hi_ci），此处显式重建
  if (!all(c("lo_ci", "hi_ci") %in% names(res)))
    res[, `:=`(lo_ci = b - 1.96 * se, hi_ci = b + 1.96 * se)]
  het <- tryCatch(mr_heterogeneity(d), error = function(e) NULL)
  ple <- tryCatch(mr_pleiotropy_test(d), error = function(e) NULL)
  qp  <- if (!is.null(het) && nrow(het))
           as.data.table(het)[method == "Inverse variance weighted", Q_pval][1] else NA_real_
  res[, `:=`(pheno = ph, pheno_label = cfg$label, n_snp = nrow(d),
             ncase = cfg$ncase, nctrl = cfg$nctrl, R2_exposure = R2,
             Q_pval = qp,
             egger_intercept = if (!is.null(ple)) ple$egger_intercept[1] else NA_real_,
             egger_p = if (!is.null(ple)) ple$pval[1] else NA_real_)]
  results[[ph]] <- res; detail[[ph]] <- d

  ivw <- res[method == "Inverse variance weighted"]
  cat(sprintf("\n== %-16s %-18s %5s 例 | SNP=%2d | IVW OR=%.3f (%.3f-%.3f) p=%.4g | Q_p=%.3g | Egger_p=%.3g\n",
              ph, cfg$label, ifelse(is.na(cfg$ncase), "NA", cfg$ncase), nrow(d),
              exp(ivw$b), exp(ivw$lo_ci), exp(ivw$hi_ci), ivw$pval, qp,
              res$egger_p[1]))
  logmsg("%-16s n=%d IVW b=%.5f se=%.5f p=%.5g", ph, nrow(d), ivw$b, ivw$se, ivw$pval)
}

A <- rbindlist(results, fill = TRUE)
fwrite(A, file.path(OUT, "bbj_mr_results.csv"))
fwrite(rbindlist(detail, fill = TRUE, idcol = "pheno"), file.path(OUT, "bbj_mr_input_snps.csv"))

cat("\n\n================ MR 汇总（IVW 主结果）================\n")
print(A[method == "Inverse variance weighted",
        .(pheno_label, n_snp, ncase, OR = round(exp(b), 3),
          lci = round(exp(lo_ci), 3), uci = round(exp(hi_ci), 3),
          p = signif(pval, 3), Q_pval = signif(Q_pval, 2),
          Egger_int_p = signif(egger_p, 2))])

cat("\n================ 全部方法（供敏感性分析）================\n")
print(A[, .(pheno = pheno, method, n_snp, OR = round(exp(b), 3),
            lci = round(exp(lo_ci), 3), uci = round(exp(hi_ci), 3), p = signif(pval, 2))])

# ---------------- 效能评估（用实际观测到的 MR 标准误，不做解析近似）----------------
cat("\n================ 效能评估 ================\n")
cat("暴露侧观测尺度 R2 = ", sprintf("%.4f", R2),
    "（二分类暴露的观测尺度 R2 会高估责任尺度，故效能用实测 SE 而非该 R2 推算）\n", sep = "")
ivwA <- A[method == "Inverse variance weighted"]
for (i in seq_len(nrow(ivwA))) {
  se_obs <- ivwA$se[i]
  b_eu <- 0.198                                # 欧裔基线 logOR (OR=1.22)
  pw <- pnorm(qnorm(0.025) + b_eu / se_obs)
  cat(sprintf("  %-16s 实测 SE(MR)=%.3f  ⇒  可检出(80%% 把握) OR≥%.2f  |  对欧裔效应 OR=1.22 的检出力=%.0f%%\n",
              ivwA$pheno[i], se_obs,
              exp((qnorm(0.975) + qnorm(0.8)) * se_obs), 100 * pw))
}
logmsg("done")
