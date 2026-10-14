#!/usr/bin/env Rscript
# 37_mediation.R —— 两步介导 MR：饮食 → 胆囊结石病 → 胆道癌
#
# ── 设计 ──────────────────────────────────────────────────────────────
#   Step 1（路径 a）：饮食暴露 → 胆囊结石病   —— 暴露侧 IV（UKB 饮食），结局为 FinnGen 胆石症
#   Step 2（路径 b）：胆囊结石病 → 胆道癌     —— 中介侧 IV（FinnGen 胆石 84 个 clumped IV）
#   总效应（路径 c）：饮食暴露 → 胆道癌       —— 暴露侧 IV（与 a 同一组 IV）
#
#   间接效应 = a × b         直接效应 = c − a × b         中介比例 = a × b / c
#
# ── 标准误口径（★ 关键，已与主稿对齐）─────────────────────────────────
#   全部路径一律采用 **随机效应校正 IVW**：
#       se_RE = se_FE × max(1, sqrt(Q/df)),  se_FE = 1/sqrt(Σ w_k)
#   这与 TwoSampleMR::mr_ivw() 默认行为及 03_mr_core.R 的手工校正一致。
#   本数据下 a 步的 Q 极大（饮酒频率 Q=453.9/92 df, P=3e-49），
#   FE 会把 se 低估 2.2 倍，因此**必须**用 RE；报告同时给出 FE 作对照。
#
# ── 三个估计量（互为敏感分析）─────────────────────────────────────────
#   IVW-RE        ：全体 IV，随机效应校正
#   MR-PRESSO     ：剔除离群 SNP 后的 IVW（主稿 §三 主结果口径）
#   Weighted med. ：对多效性稳健、不需数据依赖的离群剔除
#
# ── 相关性处理（★ 本脚本相对常规做法的改进）───────────────────────────
#   a 与 c 用**同一组饮食 IV**、且两步都用到同一个 FinnGen 胆石症 GWAS，
#   因此三个估计量彼此不独立。汇总数据层面无法直接求出抽样相关性，处理方式：
#     · 把「两 FinnGen GWAS 因样本重叠导致的每 SNP 效应相关」参数化为 rho_g；
#     · cov(a, c) 由 rho_g **解析**求出（cov_ac()，无需个体数据）：
#         cov(a,c) = rho_g · κ1·κ2 · Σ w1k·w2k·(se_chol,k·se_btc,k / β_diet,k²) / (Σw1·Σw2)
#       其中 κ = sqrt(Q/df) 为各自 RE 膨胀因子（floor 1）。已证该式给出的是真相关系数（Cauchy–Schwarz ≤ 1）。
#     · rho_g 未知 → 主结果取 0（标准/保守假设），并给 rho_g ∈ {0,.25,.5,.75,.95} 网格。
#     · a 与 b 的相关系数 rho_ab 同样网格（默认 0：两步的 IV 集合与 GWAS 角色均不同）。
#
# 用法: Rscript scripts/37_mediation.R
# 输出: results/medi2/

suppressMessages({
  library(data.table)
  library(TwoSampleMR)
})

ROOT <- normalizePath(".")
OUT  <- file.path(ROOT, "results", "medi2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LOG <- file(file.path(OUT, "mediation_log.txt"), open = "wt")
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); writeLines(msg, LOG) }

say("================================================================")
say("两步介导 MR：饮食 → 胆囊结石病 → 胆道癌")
say("运行时间: ", format(Sys.time()))
say("标准误口径: 随机效应校正 IVW（与主稿 03_mr_core.R 一致）")
say("================================================================")

## ───────────────────────────── 1. 读入 ─────────────────────────────
d_chol <- readRDS(file.path(ROOT, "results/full/mr_detail_Cholelithiasis.rds"))
d_btc  <- readRDS(file.path(ROOT, "results/full/mr_detail_BTC.rds"))
dt1 <- as.data.table(d_chol$dat)[mr_keep == TRUE]   # 饮食 IV × 胆石症
dt3 <- as.data.table(d_btc$dat)[mr_keep == TRUE]    # 饮食 IV × 胆道癌

chole_iv <- fread(file.path(ROOT, "results/medi/chole_iv_clumped.tsv"))
btc_ch   <- fread(file.path(ROOT, "results/medi/btc_for_chole.tsv"))
btc_ch[, rsid := tstrsplit(rsids, ",", keep = 1L)]

presso_a <- fread(file.path(ROOT, "results/full/mr_presso_results.csv"))
presso_c <- fread(file.path(ROOT, "results/presso_btc/mr_presso_results.csv"))
num <- function(x) suppressWarnings(as.numeric(x))

say("\n[输入] 饮食×胆石症: ", nrow(dt1), " 条；饮食×胆道癌: ", nrow(dt3), " 条")
say("[输入] 胆石症 clumped IV: ", nrow(chole_iv), "；btc_for_chole: ", nrow(btc_ch), " 行")

## ───────────────────── 2. 路径 b：胆囊结石病 → 胆道癌 ─────────────────────
b_dat <- merge(
  chole_iv[, .(rsid, beta_chol = beta, se_chol = se, ref_chol = ref, alt_chol = alt)],
  btc_ch[, .(rsid, btc_ref = ref, btc_alt = alt, beta_btc = beta, se_btc = sebeta)],
  by = "rsid"
)
b_dat[, `:=`(ok   = (alt_chol == btc_alt & ref_chol == btc_ref),
             flip = (alt_chol == btc_ref & ref_chol == btc_alt))]
b_dat[flip == TRUE, beta_btc := -beta_btc]
b_dat <- b_dat[ok == TRUE]
say("\n[Step2] 匹配 ", nrow(b_dat), " 个 IV；等位基因直接匹配 ", sum(b_dat$ok),
    "，翻转 ", sum(b_dat$flip, na.rm = TRUE), "，丢弃 0")

b_dat[, r := beta_btc / beta_chol][, w := beta_chol^2 / se_btc^2]
b       <- sum(b_dat$w * b_dat$r) / sum(b_dat$w)
se_b_fe <- 1 / sqrt(sum(b_dat$w))
Qb      <- sum(b_dat$w * (b_dat$r - b)^2); dfb <- nrow(b_dat) - 1
kapp_b  <- max(1, sqrt(Qb / dfb))
se_b    <- se_b_fe * kapp_b
say(sprintf("[Step2] b = %.6f | FE SE %.6f | Q = %.2f on %d df (P=%.3g, κ=%.4f) | RE SE %.6f",
            b, se_b_fe, Qb, dfb, pchisq(Qb, dfb, lower.tail = FALSE), kapp_b, se_b))
say(sprintf("[Step2] OR = %.4f (%.4f–%.4f), P = %.3g, n = %d",
            exp(b), exp(b - 1.96 * se_b), exp(b + 1.96 * se_b), 2 * pnorm(-abs(b / se_b)), nrow(b_dat)))

fwrite(data.table(term = "b", b = b, se_fe = se_b_fe, se_re = se_b, kappa = kapp_b,
                  Q = Qb, Q_df = dfb, Q_pval = pchisq(Qb, dfb, lower.tail = FALSE),
                  nsnp = nrow(b_dat), or = exp(b),
                  lci = exp(b - 1.96 * se_b), uci = exp(b + 1.96 * se_b),
                  p = 2 * pnorm(-abs(b / se_b))),
       file.path(OUT, "mediation_path_b.csv"))

## ───────────────── 3. 逐暴露计算 a、c（IVW-RE，含 Σ） ─────────────────
ids <- intersect(unique(dt1$id.exposure), unique(dt3$id.exposure))
say("\n[暴露] 可同时估计 a 与 c: ", length(ids), " 个")

calc_one <- function(id) {
  x1 <- dt1[id.exposure == id]; x3 <- dt3[id.exposure == id]
  K <- intersect(x1$SNP, x3$SNP)
  x1 <- x1[SNP %in% K][order(SNP)]; x3 <- x3[SNP %in% K][order(SNP)]
  if (nrow(x1) < 2 || nrow(x3) < 2) return(NULL)
  if (!all(x1$effect_allele.exposure == x3$effect_allele.exposure)) return(NULL)
  bd <- x1$beta.exposure; sd_ <- x1$se.exposure
  bc <- x1$beta.outcome;  sc <- x1$se.outcome
  bt <- x3$beta.outcome;  st <- x3$se.outcome
  if (any(sd_ == 0) || any(!is.finite(c(bc, sc, bt, st)))) return(NULL)

  r1 <- bc / bd; w1 <- bd^2 / sc^2
  r2 <- bt / bd; w2 <- bd^2 / st^2
  a <- sum(w1 * r1) / sum(w1); se_a_fe <- 1 / sqrt(sum(w1))
  c_ <- sum(w2 * r2) / sum(w2); se_c_fe <- 1 / sqrt(sum(w2))
  Q1 <- sum(w1 * (r1 - a)^2);  df1 <- length(r1) - 1
  Q2 <- sum(w2 * (r2 - c_)^2); df2 <- length(r2) - 1
  k1 <- if (df1 > 0) max(1, sqrt(Q1 / df1)) else 1
  k2 <- if (df2 > 0) max(1, sqrt(Q2 / df2)) else 1
  se_a <- se_a_fe * k1; se_c <- se_c_fe * k2

  list(id = id, exposure = x1$exposure[1], nsnp = nrow(x1),
       a = a, se_a_fe = se_a_fe, se_a = se_a, k1 = k1, Q1 = Q1, Q1_df = df1,
       c = c_, se_c_fe = se_c_fe, se_c = se_c, k2 = k2, Q2 = Q2, Q2_df = df2,
       W = data.table(w1 = w1, w2 = w2, se_chol = sc, se_btc = st, b_diet = bd))
}

res_list <- Filter(Negate(is.null), lapply(ids, function(id) tryCatch(calc_one(id), error = function(e) NULL)))
say("[暴露] 成功估计: ", length(res_list), " 个")

## cov(a,c) 解析式（含 RE 膨胀因子 κ）
cov_ac <- function(W, k1, k2, rho_g) {
  if (rho_g == 0) return(0)
  rho_g * k1 * k2 * sum(W$w1 * W$w2 * (W$se_chol * W$se_btc / W$b_diet^2)) / (sum(W$w1) * sum(W$w2))
}

base <- rbindlist(lapply(res_list, function(z)
  data.table(id = z$id, exposure = z$exposure, nsnp = z$nsnp,
             a = z$a, se_a_fe = z$se_a_fe, se_a = z$se_a, kappa_a = z$k1,
             Q_a = z$Q1, Q_a_df = z$Q1_df,
             p_a = 2 * pnorm(-abs(z$a / z$se_a)), p_a_fe = 2 * pnorm(-abs(z$a / z$se_a_fe)),
             c = z$c, se_c_fe = z$se_c_fe, se_c = z$se_c, kappa_c = z$k2,
             Q_c = z$Q2, Q_c_df = z$Q2_df,
             p_c = 2 * pnorm(-abs(z$c / z$se_c)))))

# 两性状 GWAS 在饮食 IV 上的实际共线性（rho_g = 1 时的 cor(a,c)，用于判断重叠上限）
base[, corr_ac_max := mapply(function(i) {
  W <- res_list[[i]]$W
  s <- sum(W$w1 * W$w2 * (W$se_chol * W$se_btc / W$b_diet^2)) /
       (sqrt(sum(W$w1)) * sqrt(sum(W$w2)))
  s
}, seq_len(nrow(base)))]

fwrite(base, file.path(OUT, "mediation_paths_all.csv"))
say("已写出 mediation_paths_all.csv（", nrow(base), " 个暴露 × 路径 a/c）")

## ───────────── 4. 间接效应 + 中介比例（三种估计量） ─────────────
add_ind <- function(a, se_a, c_, se_c, rho_ab = 0, rho_g = 0, W = NULL,
                    k1 = 1, k2 = 1, se_a_ref = NULL, se_c_ref = NULL) {
  ind    <- a * b
  se_ind <- sqrt(a^2 * se_b^2 + b^2 * se_a^2 + 2 * a * b * rho_ab * se_a * se_b)
  cav    <- if (is.null(W)) 0 else cov_ac(W, k1, k2, rho_g)
  cadj   <- if (is.null(se_a_ref)) se_a else se_a_ref   # cov 用未校正 SE 求定 → 按 rho_g=1 的最大值线性缩放
  cbc    <- if (is.null(se_c_ref)) se_c else se_c_ref
  covic  <- b * cav * (se_a / cadj) * (se_c / cbc)
  se_dir <- sqrt(pmax(0, se_c^2 + se_ind^2 - 2 * covic))
  list(ind = ind, se_ind = se_ind, dir = c_ - ind, se_dir = se_dir,
       prop = ind / c_, cov_ic = covic)
}

## 组装三套估计量
est_rows <- list()

for (i in seq_len(nrow(base))) {
  r <- base[i]
  W <- res_list[[i]]$W
  # (1) IVW-RE
  z1 <- add_ind(r$a, r$se_a, r$c, r$se_c, 0, 0, W, r$kappa_a, r$kappa_c, r$se_a, r$se_c)
  est_rows[[length(est_rows) + 1]] <- data.table(
    id = r$id, exposure = r$exposure, nsnp = r$nsnp, estimator = "IVW-RE",
    a = r$a, se_a = r$se_a, p_a = r$p_a, c = r$c, se_c = r$se_c, p_c = r$p_c,
    indirect = z1$ind, se_indirect = z1$se_ind, p_indirect = 2 * pnorm(-abs(z1$ind / z1$se_ind)),
    direct = z1$dir, se_direct = z1$se_dir, prop = z1$prop, corr_ac_max = r$corr_ac_max)
}

## MR-PRESSO 校正值（仅部分暴露可运行）
for (i in seq_len(nrow(presso_a))) {
  pa <- presso_a[i]
  if (is.na(num(pa$b_corr))) next
  pc <- presso_c[id == pa$id]
  cval <- if (nrow(pc) && !is.na(num(pc$b_corr[1]))) num(pc$b_corr[1]) else NA_real_
  cse  <- if (nrow(pc) && !is.na(num(pc$se_corr[1]))) num(pc$se_corr[1]) else NA_real_
  a_c  <- num(pa$b_corr); se_ac <- num(pa$se_corr)
  lv <- add_ind(a_c, se_ac, NA_real_, NA_real_, 0, 0)     # 只用其间接效应
  est_rows[[length(est_rows) + 1]] <- data.table(
    id = pa$id, exposure = pa$exposure, nsnp = as.integer(pa$nsnp), estimator = "MR-PRESSO",
    a = a_c, se_a = se_ac, p_a = num(pa$p_corr),
    c = cval, se_c = cse,
    p_c = if (is.na(cval)) NA_real_ else 2 * pnorm(-abs(cval / cse)),
    indirect = lv$ind, se_indirect = lv$se_ind,
    p_indirect = 2 * pnorm(-abs(lv$ind / lv$se_ind)),
    direct = if (is.na(cval)) NA_real_ else cval - lv$ind,
    se_direct = NA_real_,
    prop = if (is.na(cval)) NA_real_ else lv$ind / cval,
    corr_ac_max = NA_real_)
}

## 加权中位数（稳健估计量，a 与 c 同时可得时才算比例）
wm_a <- fread(file.path(ROOT, "results/full/mr_results_Cholelithiasis.csv"))[method == "Weighted median"]
wm_c <- fread(file.path(ROOT, "results/full/mr_results_BTC.csv"))[method == "Weighted median"]
wm <- merge(wm_a[, .(rsid_key = id.exposure, exposure, a = as.numeric(b), se_a = as.numeric(se), p_a = as.numeric(pval))],
            wm_c[, .(rsid_key = id.exposure, c = as.numeric(b), se_c = as.numeric(se), p_c = as.numeric(pval))],
            by = "rsid_key")
for (i in seq_len(nrow(wm))) {
  z <- add_ind(wm$a[i], wm$se_a[i], wm$c[i], wm$se_c[i], 0, 0)
  est_rows[[length(est_rows) + 1]] <- data.table(
    id = wm$rsid_key[i], exposure = wm$exposure[i], nsnp = NA_integer_, estimator = "Weighted median",
    a = wm$a[i], se_a = wm$se_a[i], p_a = wm$p_a[i],
    c = wm$c[i], se_c = wm$se_c[i], p_c = wm$p_c[i],
    indirect = z$ind, se_indirect = z$se_ind, p_indirect = 2 * pnorm(-abs(z$ind / z$se_ind)),
    direct = z$dir, se_direct = z$se_dir, prop = z$prop, corr_ac_max = NA_real_)
}

est <- rbindlist(est_rows, fill = TRUE)
est[, `:=`(or_indirect = exp(indirect), or_ind_lci = exp(indirect - 1.96 * se_indirect),
           or_ind_uci = exp(indirect + 1.96 * se_indirect))]
fwrite(est, file.path(OUT, "mediation_estimators.csv"))

## ───────────── 5. 参数 bootstrap：中介比例 CI（rho 网格） ─────────────
boot_prop <- function(a, se_a, c_, se_c, cov_ac_v, rho_ab, nboot = 10000, seed = 20260917) {
  cov_ab <- rho_ab * se_a * se_b
  cov_bc <- rho_ab * se_b * se_c
  S <- matrix(c(se_a^2, cov_ab, cov_ac_v,
                cov_ab,  se_b^2, cov_bc,
                cov_ac_v, cov_bc, se_c^2), 3, 3, byrow = TRUE)
  ee <- eigen(S, symmetric = TRUE)
  if (any(ee$values < 0)) { ee$values[ee$values < 0] <- 0
    S <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors) }
  set.seed(seed)
  L <- t(chol(S)); Z <- matrix(rnorm(nboot * 3), nboot, 3) %*% t(L)
  A <- a + Z[, 1]; B <- b + Z[, 2]; C <- c_ + Z[, 3]
  R <- (A * B) / C
  R <- R[is.finite(R) & abs(R) < 1e3]
  c(quantile(R, 0.025, names = FALSE), quantile(R, 0.975, names = FALSE), length(R))
}

GRID <- c(0, 0.25, 0.5, 0.75, 0.95)
sens <- rbindlist(lapply(seq_len(nrow(base)), function(i) {
  r <- base[i]; W <- res_list[[i]]$W
  rbindlist(lapply(GRID, function(rg) {
    cav <- cov_ac(W, r$kappa_a, r$kappa_c, rg)
    rbindlist(lapply(c(0, 0.25, 0.5), function(rab) {
      bs <- boot_prop(r$a, r$se_a, r$c, r$se_c, cav, rab, nboot = 4000)
      data.table(id = r$id, exposure = r$exposure, rho_g = rg, rho_ab = rab,
                 prop = (r$a * b) / r$c, b_lci = bs[1], b_uci = bs[2], n_ok = bs[3])
    }))
  }))
}))
fwrite(sens, file.path(OUT, "mediation_sensitivity_grid.csv"))

## ───────────── 6. 报告 ─────────────
PRIMARY_IDS <- c("ukb-b-5779", "ukb-b-1489", "ukb-b-16576")
NICE <- c("ukb-b-5779" = "Alcohol intake frequency", "ukb-b-1489" = "Cheese intake",
          "ukb-b-16576" = "Dried fruit intake")

say("\n================================================================")
say("路径 b：胆囊结石病 → 胆道癌")
say("================================================================")
say(sprintf("  IVW-RE  b = %.4f (SE %.4f)  OR %.3f (%.3f–%.3f)  P = %.3g  n = %d",
            b, se_b, exp(b), exp(b - 1.96 * se_b), exp(b + 1.96 * se_b),
            2 * pnorm(-abs(b / se_b)), nrow(b_dat)))

say("\n================================================================")
say("主暴露（MR-PRESSO 校正 FDR<0.05）：间接效应与中介比例")
say("================================================================")
for (e in c("MR-PRESSO", "IVW-RE", "Weighted median")) {
  say(sprintf("\n【%s】", e))
  for (eid in PRIMARY_IDS) {
    z <- est[id == eid & estimator == e]
    if (!nrow(z)) { say(sprintf("  %-26s —— 该估计量不可用", NICE[[eid]])); next }
    say(sprintf("  %-26s  a=%+.4f  c=%+.4f  indirect=%+.4f (SE %.4f, P=%.3g)  OR=%.3f (%.3f–%.3f)",
                NICE[[eid]], z$a, z$c, z$indirect, z$se_indirect, z$p_indirect,
                z$or_indirect, z$or_ind_lci, z$or_ind_uci))
    if (!is.na(z$prop)) {
      s <- sens[id == eid]
      say(sprintf("  %-26s  proportion = %+.3f | bootstrap CI: rho_g=0 [%+.2f,%+.2f] | rho_g=0.5 [%+.2f,%+.2f] | rho_g=0.95 [%+.2f,%+.2f]",
                  "", z$prop,
                  s[rho_g == 0 & rho_ab == 0]$b_lci, s[rho_g == 0 & rho_ab == 0]$b_uci,
                  s[rho_g == 0.5 & rho_ab == 0]$b_lci, s[rho_g == 0.5 & rho_ab == 0]$b_uci,
                  s[rho_g == 0.95 & rho_ab == 0]$b_lci, s[rho_g == 0.95 & rho_ab == 0]$b_uci))
    }
  }
}

say("\n================================================================")
say("直接效应（差值法：c − a×b），IVW-RE")
say("================================================================")
for (eid in PRIMARY_IDS) {
  z <- est[id == eid & estimator == "IVW-RE"]
  say(sprintf("  %-26s  total c=%+.4f (SE %.4f) | direct=%+.4f (SE %.4f) | indirect=%+.4f",
              NICE[[eid]], z$c, z$se_c, z$direct, z$se_direct, z$indirect))
}

say("\n================================================================")
say("全部 a 的 P<0.05 暴露（含 IV 数不足、无法做 MR-PRESSO 者）")
say("================================================================")
sub <- base[p_a < 0.05][order(p_a)]
for (i in seq_len(nrow(sub))) {
  s <- sub[i]
  e <- est[id == s$id & estimator == "IVW-RE"]
  say(sprintf("  %-40s nIV=%-3d a=%+.4f (P=%.3g) c=%+.4f (P=%.3g) ind=%+.4f 比例=%+.3f 同向=%s",
              substr(s$exposure, 1, 40), s$nsnp, s$a, s$p_a, s$c, s$p_c,
              e$indirect, e$prop, sign(s$a) == sign(s$c)))
}

say("\n================================================================")
say("完成 → ", OUT)
say("================================================================")
close(LOG)
cat("\n完成。日志:", file.path(OUT, "mediation_log.txt"), "\n")
