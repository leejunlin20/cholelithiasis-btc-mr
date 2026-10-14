#!/usr/bin/env Rscript
# 40_finngen_outcome_panel.R —— 同一暴露（胆囊结石病 84 个 IV）× FinnGen 多终点面板
#
# 目的：为 FinnGen 侧补上与 BBJ 报告（§四 效能表 / §五 管线有效性）对称的证据：
#   · 阳性对照：胆囊切除、胆囊炎 —— 应与胆石症强相关（同一疾病谱）
#   · 主结局    ：胆道癌 —— 效应量小、病例少
#   · 阴性对照  ：BMI（连续）—— 不应被胆石 IV 强驱动（或应显示方向差异）
#   同时给出每个终点的**实测 MR 标准误 → 可检出最小 OR → 对 OR 1.22 的检出力**。
#
# 数据来源：results/mvmr/watch_*.tsv（此前为 MVMR 生成的 505 SNP 观察名单命中结果，
#          已含全部 84 个胆石 clumped IV，因此**无需重新下载或扫描任何 GWAS 文件**）
# 病例/对照数：FinnGen R12 官方 finngen_R12_pheno_n.tsv
#
# 用法: Rscript scripts/40_finngen_outcome_panel.R
# 输出: results/power/

suppressMessages({
  library(data.table)
  library(TwoSampleMR)
})

ROOT <- normalizePath(".")
OUT  <- file.path(ROOT, "results", "power")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(OUT, "outcome_panel_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG) }

Z_A <- qnorm(0.975); SUM <- qnorm(0.975) + qnorm(0.80)
or_min <- function(se) exp(SUM * se)
pw     <- function(beta, se) pnorm(abs(beta) / se - Z_A) + pnorm(-abs(beta) / se - Z_A)

## ── 暴露：胆囊结石病 84 个 clumped IV ────────────────────────────────
iv <- fread(file.path(ROOT, "results/mvmr/iv_cholelith_clumped.tsv"))
say("================================================================")
say("FinnGen 多终点面板：胆囊结石病（", nrow(iv), " 个 clumped IV）")
say("运行时间: ", format(Sys.time()))
say("================================================================")

## ── 终点定义（观察名单文件名 → 标签 / 病例数来源端点）────────────────
PANEL <- list(
  # ⚠️ K11_CHOLELITH 与暴露是同一份 GWAS → β 必然 = 1、Q 必然 = 0。
  #    这不是阳性对照，而是**等位基因对齐与合并流程的技术校验**，不进主表。
  list(file = "watch_cholelith.tsv",  label = "【技术校验】自身回归（同一 GWAS）", pheno = "K11_CHOLELITH", type = "selfchk"),
  list(file = "watch_cholecstec.tsv", label = "胆囊切除（阳性对照）",   pheno = "K11_CHOLECYSTECTOMY", type = "bin"),
  list(file = "watch_cholecyst.tsv",  label = "胆囊炎（阳性对照）",     pheno = "K11_CHOLECYST",      type = "bin"),
  list(file = "watch_btc.tsv",        label = "胆道癌（主结局）",       pheno = "C3_BILIARY_GALLBLADDER_EXALLC", type = "bin"),
  list(file = "watch_gallblaoth.tsv", label = "胆囊其他疾病 K82",       pheno = "K11_GALLBLAOTH",     type = "bin"),
  list(file = "watch_bmi.tsv",        label = "BMI（阴性对照，连续）",  pheno = "BMI_IRN",            type = "cont")
)

man <- fread(file.path(ROOT, "data/finngen/finngen_R12_pheno_n.tsv"))
setnames(man, gsub("^﻿", "", names(man)))
man[, `:=`(num_cases = as.numeric(num_cases), num_controls = as.numeric(num_controls))]

run_one <- function(spec) {
  ob <- fread(file.path(ROOT, "results/mvmr", spec$file))
  m <- merge(iv[, .(rsid, beta_exp = beta, se_exp = se, ref_exp = ref, alt_exp = alt, eaf_exp = eaf)],
             ob[, .(rsid, ref_o = ref, alt_o = alt, beta_o = beta, se_o = se, pval_o = pval, eaf_o = eaf)],
             by = "rsid")
  m <- m[!is.na(beta_o) & !is.na(se_o)]
  m[, ok := (alt_exp == alt_o & ref_exp == ref_o)]
  m[alt_exp == ref_o & ref_exp == alt_o, beta_o := -beta_o]
  m[, flip := (alt_exp == ref_o & ref_exp == alt_o)]
  n_flip <- sum(m$flip, na.rm = TRUE)
  m <- m[ok == TRUE]
  if (nrow(m) < 3) return(NULL)

  m[, r := beta_o / beta_exp][, w := beta_exp^2 / se_o^2]
  b <- sum(m$w * m$r) / sum(m$w)
  se_fe <- 1 / sqrt(sum(m$w))
  Q <- sum(m$w * (m$r - b)^2); df <- nrow(m) - 1
  kap <- max(1, sqrt(Q / df))
  se_re <- se_fe * kap

  # Egger / weighted median（用 TwoSampleMR 保证实现一致）
  d <- data.table(SNP = m$rsid, beta.exposure = m$beta_exp, se.exposure = m$se_exp,
                  beta.outcome = m$beta_o, se.outcome = m$se_o,
                  effect_allele.exposure = m$alt_exp, other_allele.exposure = m$ref_exp,
                  effect_allele.outcome = m$alt_o, other_allele.outcome = m$ref_o,
                  eaf.exposure = m$eaf_exp, eaf.outcome = m$eaf_o,
                  exposure = "Cholelithiasis", outcome = spec$label,
                  id.exposure = "chole", id.outcome = "out",
                  mr_keep.exposure = TRUE, mr_keep.outcome = TRUE, mr_keep = TRUE,
                  remove = FALSE, palindromic = NA, ambiguous = NA,
                  action = 2, SNP_index = seq_len(nrow(m)),
                  samplesize.exposure = NA_real_, samplesize.outcome = NA_real_)
  eg  <- tryCatch(mr_egger_regression(m$beta_exp, m$beta_o, m$se_exp, m$se_o), error = function(e) NULL)
  wm  <- tryCatch(mr_weighted_median(m$beta_exp, m$beta_o, m$se_exp, m$se_o), error = function(e) NULL)

  nc <- man[phenocode == spec$pheno]
  # ⚠️ Egger 的 P 值必须直接用包内返回值（$pval / $pval_i，基于 t 分布），
  #    不要用 2*pnorm(...) 重算 —— 正态近似会系统性偏小（实测 0.00243 vs t 分布 0.00325）。
  list(b = b, se_fe = se_fe, se_re = se_re, kap = kap, Q = Q, Q_df = df,
       Q_p = pchisq(Q, df, lower.tail = FALSE), n = nrow(m), n_flip = n_flip,
       eg_b = if (is.null(eg)) NA_real_ else eg$b,
       eg_se = if (is.null(eg)) NA_real_ else eg$se,
       eg_p = if (is.null(eg)) NA_real_ else eg$pval,
       eg_int = if (is.null(eg)) NA_real_ else eg$b_i,
       eg_int_p = if (is.null(eg)) NA_real_ else eg$pval_i,
       wm_b = if (is.null(wm)) NA_real_ else wm$b,
       wm_se = if (is.null(wm)) NA_real_ else wm$se,
       wm_p = if (is.null(wm)) NA_real_ else wm$pval,
       case = nc$num_cases[1], ctrl = nc$num_controls[1], type = spec$type,
       label = spec$label)
}

res <- rbindlist(lapply(PANEL, function(s) {
  z <- tryCatch(run_one(s), error = function(e) { say("  [失败] ", s$label, ": ", conditionMessage(e)); NULL })
  if (is.null(z)) return(NULL)
  data.table(label = z$label, type = z$type,
             ncase = z$case, nctrl = z$ctrl, nsnp = z$n, n_flip = z$n_flip,
             or_ivw = exp(z$b), lci = exp(z$b - 1.96 * z$se_re), uci = exp(z$b + 1.96 * z$se_re),
             p_ivw = 2 * pnorm(-abs(z$b / z$se_re)),
             se_fe = z$se_fe, kappa = z$kap, se_re = z$se_re,
             Q = z$Q, Q_df = z$Q_df, Q_p = z$Q_p,
             or_egger = exp(z$eg_b), p_egger = z$eg_p,
             egger_int = z$eg_int, egger_int_p = z$eg_int_p,
             or_wm = exp(z$wm_b), p_wm = z$wm_p)
}), fill = TRUE)

res[, `:=`(or_min80 = or_min(se_re),
           power_or120 = pw(log(1.20), se_re),
           power_or122 = pw(log(1.22), se_re),
           power_or150 = pw(log(1.50), se_re))]
fwrite(res, file.path(OUT, "finngen_outcome_panel.csv"))

say("\n================================================================")
say("胆囊结石病 IV → FinnGen 各终点（84 clumped IV，等位基因 0 不一致）")
say("================================================================")
say(sprintf("%-26s %9s %6s %23s %10s %11s %12s %11s %10s", "终点", "病例数", "nSNP",
            "IVW OR (95% CI)", "IVW P", "Egger 截距", "截距 P（多效性）", "Egger 斜率P", "WMed P"))
for (i in seq_len(nrow(res))) {
  r <- res[i]
  say(sprintf("%-26s %9s %6d %23s %10.3g %11.6f %12.3g %11.3g %10.3g",
              r$label, format(r$ncase, big.mark = ","), r$nsnp,
              sprintf("%.3f (%.3f–%.3f)", r$or_ivw, r$lci, r$uci),
              r$p_ivw, r$egger_int, r$egger_int_p, r$p_egger, r$p_wm))
}
say("\n注：判读多效性看『截距 P』；『Egger 斜率P』是 MR-Egger 因果估计本身的 P 值，与前者的含义不同。")

say("\n================================================================")
say("量表 P4　各终点的统计效能（用实测 MR 标准误反推）")
say("================================================================")
say(sprintf("%-26s %9s %9s %10s %9s %9s %9s", "终点", "病例数", "SE(RE)", "最小OR80", "P(1.20)", "P(1.22)", "P(1.50)"))
for (i in seq_len(nrow(res))) {
  r <- res[i]
  say(sprintf("%-26s %9s %9.4f %10.3f %8.0f%% %8.0f%% %8.0f%%",
              r$label, format(r$ncase, big.mark = ","), r$se_re, r$or_min80,
              100 * r$power_or120, 100 * r$power_or122, 100 * r$power_or150))
}

say("\n--- Cochran Q 异质性（RE 校正因子的来源）---")
for (i in seq_len(nrow(res))) {
  r <- res[i]
  say(sprintf("  %-26s Q = %8.2f on %3d df, P = %-10.3g, sqrt(Q/df) = %.3f",
              r$label, r$Q, r$Q_df, r$Q_p, r$kappa))
}

say("\n完成 → ", OUT)
close(LOG)
cat("\n完成。日志:", file.path(OUT, "outcome_panel_log.txt"), "\n")
