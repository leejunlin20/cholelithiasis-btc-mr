#!/usr/bin/env Rscript
# 33_mvmr_run.R —— 多变量 MR：胆石症 → 胆道癌，校正胆囊炎 / 胆囊疾病
#
# 动机（回应竞品）
# --------------
# Bai Y et al. J Adv Res 2026 (PMID 40064441) 的网络 MR 显示：MVMR 校正胆囊炎+胆囊息肉后，
# 胆石症→胆道癌 衰减至 OR 1.45 (0.60–3.52)、不再显著。我方此前只校正了 BMI 与 E4_OBESITY，
# 必须补做同一套校正，否则审稿人必问。
#
# ★ 端点定义的关键事实（Risteys 官方定义，2026-09-17 核实）
#   K11_CHOLELITH      = ICD-10 K80（胆石症）——K80.0/K80.1 本身即「结石性胆囊炎」
#   K11_CHOLECYST      = ICD-10 K81（胆囊炎）——结石性者归 K80，故与 K80 部分互斥
#   K11_GALLBLAOTH     = ICD-10 K82（其他胆囊疾病）——K82.8 即胆囊息肉所在码
#   K11_CHOLECYSTECTOMY= 胆囊切除术（36,556 例）——高功效替代校正项
#   三端点共同对照排除 K11_GALLBILPANC
#   ⇒ FinnGen **无**独立「胆囊息肉」端点；K82 含息肉码但全基因组仅 2 个显著位点
#     且均为 AF≈0.027 的罕见变异、不在 1000G EUR LD 面板中 → **clump 后 0 个可用 IV**
#
# 用法: Rscript scripts/33_mvmr_run.R

suppressMessages({
  library(data.table)
  library(MVMR)
  library(MendelianRandomization)
  library(TwoSampleMR)
})

BASE <- Sys.getenv("MR_BASE", unset = ".")
OUT  <- file.path(BASE, "results/mvmr")
LOG  <- file.path(OUT, "mvmr_run_log.txt")
sink(LOG, split = TRUE)

TT <- list(
  cholelith  = list(f = "watch_cholelith.tsv",  lab = "Cholelithiasis (K80)",
                    ncase = 49834,  nctrl = 437418),
  cholecyst  = list(f = "watch_cholecyst.tsv",  lab = "Cholecystitis (K81)",
                    ncase = 5851,   nctrl = 437418),
  gallblaoth = list(f = "watch_gallblaoth.tsv", lab = "Other gallbladder dis. (K82)",
                    ncase = 1655,   nctrl = 437418),
  cholecstec = list(f = "watch_cholecstec.tsv", lab = "Cholecystectomy",
                    ncase = 36556,  nctrl = 463792),
  bmi        = list(f = "watch_bmi.tsv",        lab = "BMI (BMI_IRN)",
                    ncase = NA,     nctrl = NA),
  btc        = list(f = "watch_btc.tsv",        lab = "Biliary tract cancer (EXALLC)",
                    ncase = 2298,   nctrl = 378749)
)

## 二分类性状的有效样本量（log-OR 方差尺度）；连续性状用总 N
eff_n <- function(ncase, nctrl) if (is.na(ncase)) NA_real_ else
  4 / (1 / ncase + 1 / nctrl)

## ---------------------------------------------------------------- 读入
cat("==================== 观察名单命中 ====================\n")
W <- list()
for (k in names(TT)) {
  f <- file.path(OUT, TT[[k]]$f)
  if (!file.exists(f)) { cat(sprintf("  -- %-10s 缺文件，跳过\n", k)); next }
  d <- fread(f)
  setnames(d, c("beta", "se", "pval", "eaf"), c("b", "s", "p", "af"),
           skip_absent = TRUE)
  W[[k]] <- d
  cat(sprintf("  %-10s %4d SNP 命中\n", k, nrow(d)))
}
stopifnot(all(c("cholelith", "cholecyst", "btc") %in% names(W)))

## MVMR 工具变量集：必须是「各暴露 IV 的并集再 LD 剪枝」后的独立集。
## 直接拿 505 个观察名单 SNP 进去会把条件 F 稀释（未与暴露关联的位点拉低均值），
## 且名单内不同性状的位点可能互相 LD → 信息矩阵病态。
IV_MAIN <- readLines(file.path(OUT, "mvmr_union_rsids.txt"))
ivbmi_f <- file.path(OUT, "mvmr_union_with_bmi_rsids.txt")
IV_BMI  <- if (file.exists(ivbmi_f)) readLines(ivbmi_f) else IV_MAIN
cat(sprintf("MVMR 工具变量集：主模型 %d 个 | 含 BMI 模型 %d 个\n",
            length(IV_MAIN), length(IV_BMI)))

## ---------------------------------------------------------------- 合并 + 等位基因对齐
## FinnGen 同一数据冻结、同一插补面板，ref/alt 以**字符串**给出，
## 故按等位基因字符串比对即可定号，不存在链方向歧义（无需按回文/MAF 丢弃）。
base <- W$cholelith[, .(rsid, chr, pos, gene = nearest_gene, ref0 = ref, alt0 = alt,
                        b_cholelith = b, s_cholelith = s,
                        p_cholelith = p, af_cholelith = af)]
m <- base
for (k in setdiff(names(W), "cholelith")) {
  ## 逐性状用唯一列名，避免 merge 时同名列被拆成 refx.x / refx.y
  x <- W[[k]][, .(rsid, refx = ref, altx = alt, b = b, s = s, p = p, af = af)]
  setnames(x, c("refx", "altx", "b", "s", "p", "af"),
           paste0(c("refx_", "altx_", "bb_", "ss_", "pp_", "aa_"), k))
  m <- merge(m, x, by = "rsid", all = FALSE)
  rx <- paste0("refx_", k); ax <- paste0("altx_", k)
  ok   <- m[[rx]] == m$ref0 & m[[ax]] == m$alt0
  flip <- m[[rx]] == m$alt0 & m[[ax]] == m$ref0
  m[[paste0("b_", k)]]  <- ifelse(ok, m[[paste0("bb_", k)]],
                                  ifelse(flip, -m[[paste0("bb_", k)]], NA_real_))
  m[[paste0("s_", k)]]  <- m[[paste0("ss_", k)]]
  m[[paste0("p_", k)]]  <- m[[paste0("pp_", k)]]
  m[[paste0("af_", k)]] <- m[[paste0("aa_", k)]]
  cat(sprintf("  %-10s 合并 %4d SNP | 等位基因不一致 %d | 需翻转 %d\n",
              k, nrow(m), sum(!(ok | flip)), sum(flip)))
  m[, paste0(c("refx_", "altx_", "bb_", "ss_", "pp_", "aa_"), k) := NULL]
}

palin <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")
m[, palindromic := palin(ref0, alt0)]
cat(sprintf("\n回文位点 %d 个（alt 等位基因以字符串给出 → 定号无歧义，不丢弃）\n",
            sum(m$palindromic, na.rm = TRUE)))

est_row <- function(b, s) data.table(
  b = b, se = s, or = exp(b), lci = exp(b - 1.96 * s), uci = exp(b + 1.96 * s),
  pval = 2 * pnorm(-abs(b / s)))

## ---------------------------------------------------------------- 单变量参照
cat("\n==================== 单变量 MR（各暴露 → 胆道癌）====================\n")
uni_out <- list()
EXPS_ALL <- c("cholelith", "cholecyst", "gallblaoth", "cholecstec")
for (k in EXPS_ALL) {
  bk <- paste0("b_", k); sk <- paste0("s_", k); pk <- paste0("p_", k)
  if (!(bk %in% names(m))) next
  ## ★ 必须用「已 LD 剪枝的独立 IV」，不能拿观察名单里所有 p<5e-8 的位点：
  ##   观察名单含其他性状与 BMI 的位点，其中不少在 LD 上互相冗余，
  ##   直接喂进 IVW 会低估 SE、并让条件 F 被稀释。
  ivf <- file.path(OUT, sprintf("iv_%s_clumped.tsv", k))
  if (!file.exists(ivf)) { cat(sprintf("\n-- %s 无 clumped IV 文件，跳过\n", k)); next }
  ivk <- fread(ivf)
  d <- m[rsid %in% ivk$rsid]
  if (nrow(d) < 2) {
    cat(sprintf("\n-- %-10s 自身独立 IV 仅 %d 个 → 单变量 MR 不做\n", k, nrow(d)))
    next
  }
  ## 用 TwoSampleMR 的 mr()，与主报告/two-step MR 保持同一估计量栈
  ## （β 已在上面按 alt0 统一对齐，故两侧效应等位基因一致）
  meth <- if (nrow(d) >= 3)
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode") else
    c("mr_ivw", "mr_weighted_median")
  dat <- data.frame(
    SNP = d$rsid,
    beta.exposure = d[[bk]], se.exposure = d[[sk]],
    effect_allele.exposure = d$alt0, other_allele.exposure = d$ref0,
    eaf.exposure = d[[paste0("af_", k)]], pval.exposure = d[[pk]],
    beta.outcome = d$b_btc, se.outcome = d$s_btc,
    effect_allele.outcome = d$alt0, other_allele.outcome = d$ref0,
    eaf.outcome = d$af_btc, pval.outcome = d$p_btc,
    id.exposure = k, id.outcome = "btc", exposure = k, outcome = "btc",
    mr_keep = TRUE, samplesize.outcome = 381047)
  res <- tryCatch(as.data.table(mr(dat, method_list = meth)),
                  error = function(e) { cat("  mr() 失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(res) || !nrow(res)) next
  res <- cbind(exposure = k, res[, .(method, b, se, pval)],
               est_row(res$b, res$se)[, .(or, lci, uci)])
  res[, nsnp := nrow(d)]
  setcolorder(res, c("exposure", "method", "nsnp", "or", "lci", "uci", "pval", "b", "se"))
  uni_out[[k]] <- res
  cat(sprintf("\n-- %s（%s）: %d 个自身 IV\n", k, TT[[k]]$lab, nrow(d)))
  print(res[, .(method, nsnp, or = round(or, 3),
                CI = sprintf("%.3f–%.3f", lci, uci), pval = signif(pval, 3))])
  ## IV 很少（≤3）时，逐位点 Wald 比更透明，必须列出
  if (nrow(d) <= 3) {
    cat("   逐位点 Wald 比（暴露→胆道癌）：\n")
    for (i in seq_len(nrow(d))) {
      wr <- d$b_btc[i] / d[[bk]][i]
      sewr <- abs(d$s_btc[i] / d[[bk]][i])
      cat(sprintf("     %-13s %-10s β_exp=%+.3f (P_exp=%.2g) β_out=%+.3f → OR %.2f (%.2f–%.2f), P=%.3g\n",
                  d$rsid[i], d$gene[i], d[[bk]][i], d[[pk]][i], d$b_btc[i],
                  exp(wr), exp(wr - 1.96 * sewr), exp(wr + 1.96 * sewr),
                  2 * pnorm(-abs(wr / sewr))))
    }
  }
}
if (length(uni_out)) fwrite(rbindlist(uni_out), file.path(OUT, "mvmr_univariable.csv"))

## ---------------------------------------------------------------- MVMR 封装
run_mvmr <- function(dt, exps, tag) {
  if (is.null(dt) || nrow(dt) < 3) {
    cat(sprintf("\n\n===== MVMR [%s] 跳过：可用 SNP %d < 3 =====\n", tag, nrow(dt)))
    return(NULL)
  }
  cat(sprintf("\n\n==================== MVMR [%s] 暴露 = %s | SNP = %d ====================\n",
              tag, paste(exps, collapse = " + "), nrow(dt)))
  BX   <- as.matrix(dt[, paste0("b_", exps), with = FALSE])
  seBX <- as.matrix(dt[, paste0("s_", exps), with = FALSE])
  mv   <- format_mvmr(BX, dt$b_btc, seBX, dt$s_btc, dt$rsid)

  ivw <- tryCatch(ivw_mvmr(mv), error = function(e) {
    cat("ivw_mvmr 失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(ivw)) return(NULL)
  dt_ivw <- data.table(exposure = exps, b = ivw[, 1], se = ivw[, 2])
  dt_ivw <- cbind(dt_ivw, est_row(dt_ivw$b, dt_ivw$se)[, .(or, lci, uci, pval)])
  cat("\n-- MVMR-IVW --\n")
  print(dt_ivw[, .(exposure, or = round(or, 3),
                   CI = sprintf("%.3f–%.3f", lci, uci), pval = signif(pval, 3))])

  cf <- tryCatch(strength_mvmr(mv, gencov = 0), error = function(e) NULL)
  if (!is.null(cf)) {
    cat("\n-- 条件 F（Sanderson；>10 为佳）--\n")
    print(as.data.frame(cf)[, 1:2])
  }
  cfr <- tryCatch(strhet_mvmr(mv, gencov = 0), error = function(e) NULL)
  if (!is.null(cfr)) { cat("\n-- 稳健条件 F（strhet_mvmr）--\n"); print(as.data.frame(cfr)) }
  q <- tryCatch(pleiotropy_mvmr(mv, gencov = 0), error = function(e) NULL)
  if (!is.null(q)) { cat("\n-- 多效性 Q 检验 --\n"); print(as.data.frame(q)) }

  ## --- MendelianRandomization 稳健估计量 ---
  ## 注：mr_mvpcgmm 需要 SNP 间 LD（相关）矩阵，本项目未拉取 86 个 SNP 的全对 LD，
  ##     故不运行；改用 MVMR-Lasso / MVMR-median 承担多效性稳健性的角色。
  minp <- mr_mvinput(bx = BX, bxse = seBX, by = dt$b_btc, byse = dt$s_btc, snps = dt$rsid)
  se_of <- function(o) {
    sn <- slotNames(o)
    for (s in c("StdError", "StdError.Est", "StdError.Int"))
      if (s %in% sn) return(o@Estimate * 0 + slot(o, s))
    rep(NA_real_, length(o@Estimate))
  }
  robust <- list()
  add <- function(nm, o) if (!is.null(o)) {
    ## Lasso 会把被保留的 SNP 全列出来（上百个），这里只记录个数，避免日志爆炸
    if ("ValidSNPs" %in% slotNames(o))
      nm <- sprintf("%s (保留 %s 个 SNP)", nm, paste(o@ValidSNPs, collapse = "+"))
    robust[[nm]] <<- data.table(method = nm, exposure = exps,
                                b = o@Estimate, se = se_of(o))
  }
  add("MVMR-IVW",    tryCatch(mr_mvivw(minp), error = function(e) NULL))
  add("MVMR-Egger",  tryCatch(mr_mvegger(minp, orientate = 1), error = function(e) NULL))
  add("MVMR-median", tryCatch(mr_mvmedian(minp, iterations = 10000), error = function(e) NULL))
  add("MVMR-Lasso",  tryCatch(mr_mvlasso(minp, orientate = 1), error = function(e) NULL))
  rb <- as.data.table(rbindlist(robust))
  if (nrow(rb)) {
    rb <- cbind(rb, est_row(rb$b, rb$se)[, .(or, lci, uci, pval)])
    cat("\n-- 稳健估计量汇总 --\n")
    print(rb[, .(method, exposure, or = round(or, 3),
                 CI = sprintf("%.3f–%.3f", lci, uci), pval = signif(pval, 3))])
  }

  fwrite(dt_ivw, file.path(OUT, sprintf("mvmr_ivw_%s.csv", tag)))
  if (nrow(rb)) fwrite(rb, file.path(OUT, sprintf("mvmr_robust_%s.csv", tag)))
  if (!is.null(cf)) fwrite(as.data.table(as.data.frame(cf)),
                           file.path(OUT, sprintf("mvmr_condF_%s.csv", tag)))
  if (!is.null(q))  fwrite(as.data.table(as.data.frame(q)),
                           file.path(OUT, sprintf("mvmr_Q_%s.csv", tag)))
  fwrite(dt, file.path(OUT, sprintf("mvmr_input_%s.csv", tag)))
  invisible(list(ivw = dt_ivw, robust = rb, condF = cf, Q = q))
}

## ---- 模型 1：胆石 + 胆囊炎（= Bai 2026 的核心校正）
r1 <- run_mvmr(m[rsid %in% IV_MAIN & !is.na(b_cholecyst)],
               c("cholelith", "cholecyst"), "chole_cholecyst")
## ---- 模型 2：胆石 + 胆囊炎 + 胆囊切除术
if ("b_cholecstec" %in% names(m)) {
  r2 <- run_mvmr(m[rsid %in% IV_MAIN & !is.na(b_cholecyst) & !is.na(b_cholecstec)],
                 c("cholelith", "cholecyst", "cholecstec"), "chole_cholecyst_ec")
}
## ---- 模型 3：胆石 + 胆囊炎 + 胆囊切除 + BMI
if (all(c("b_bmi", "b_cholecstec") %in% names(m))) {
  r3 <- run_mvmr(m[rsid %in% IV_BMI & !is.na(b_bmi) & !is.na(b_cholecyst) &
                     !is.na(b_cholecstec)],
                 c("cholelith", "cholecyst", "cholecstec", "bmi"), "with_bmi")
}
## ---- 模型 4（对照）：胆石 + BMI（此前已报过的模型，检查一致性）
if ("b_bmi" %in% names(m)) {
  r4 <- run_mvmr(m[rsid %in% IV_BMI & !is.na(b_bmi)],
                 c("cholelith", "bmi"), "chole_bmi")
}

## ---------------------------------------------------------------- 遗传共享度
cat("\n\n==================== 胆石症 vs 胆囊炎 的遗传共享度 ====================\n")
for (k in c("cholecyst", "cholecstec")) {
  bk <- paste0("b_", k); sk <- paste0("s_", k)
  if (!(bk %in% names(m))) next
  sh <- m[!is.na(get(bk))]
  w  <- 1 / (sh[[sk]]^2 + sh$s_cholelith^2)
  mx <- sum(w * sh$b_cholelith) / sum(w); my <- sum(w * sh[[bk]]) / sum(w)
  slope <- sum(w * (sh$b_cholelith - mx) * (sh[[bk]] - my)) /
           sum(w * (sh$b_cholelith - mx)^2)
  se_sl <- sqrt(1 / sum(w * (sh$b_cholelith - mx)^2))
  cat(sprintf("\n-- %s（n = %d）--\n", k, nrow(sh)))
  cat(sprintf("   加权回归斜率（%s β 对 胆石 β）: %.3f (SE %.3f, z = %.2f, P = %.3g)\n",
              k, slope, se_sl, slope / se_sl, 2 * pnorm(-abs(slope / se_sl))))
  cat(sprintf("   Pearson r（未加权）: %.3f\n", cor(sh$b_cholelith, sh[[bk]])))
}

## ---------------------------------------------------------------- 恰好识别 MVMR
## 胆囊炎全基因组只有 2 个显著位点。用「自己的 2 个位点」做两暴露 MVMR 属恰好识别
## （2 个工具变量、2 个待估参数），无过度识别检验可用，但可给出完全不依赖
## 胆石位点的、独立的胆石直接效应估计，用于与主模型（86 SNP）互证。
cat("\n\n===== 敏感性：仅用胆囊炎的 2 个全基因组显著位点（恰好识别）=====\n")
cc_iv <- tryCatch(fread(file.path(OUT, "iv_cholecyst_clumped.tsv")),
                  error = function(e) NULL)
if (!is.null(cc_iv)) {
  dd <- m[rsid %in% cc_iv$rsid]
  cat(sprintf("胆囊炎先导位点 %d 个，其中在本表可用的 %d 个：%s\n",
              nrow(cc_iv), nrow(dd), paste(dd$rsid, collapse = ", ")))
  if (nrow(dd) >= 2) {
    print(dd[, .(rsid, gene, chr, pos,
                 b_cholelith = round(b_cholelith, 4), p_cholelith = signif(p_cholelith, 3),
                 b_cholecyst = round(b_cholecyst, 4), p_cholecyst = signif(p_cholecyst, 3),
                 b_btc = round(b_btc, 4), p_btc = signif(p_btc, 3))])
    BX <- as.matrix(dd[, .(b_cholelith, b_cholecyst)])
    mv <- format_mvmr(BX, dd$b_btc, as.matrix(dd[, .(s_cholelith, s_cholecyst)]),
                      dd$s_btc, dd$rsid)
    ivw <- tryCatch(ivw_mvmr(mv), error = function(e) {
      cat("恰好识别 MVMR 失败:", conditionMessage(e), "\n"); NULL })
    if (!is.null(ivw)) {
      rr <- est_row(ivw[, 1], ivw[, 2])
      print(data.table(exposure = c("cholelith", "cholecyst"), or = round(rr$or, 3),
                       CI = sprintf("%.3f–%.3f", rr$lci, rr$uci),
                       pval = signif(rr$pval, 3)))
      cat("（恰好识别：自由度 0，无 Q 检验；SE 完全由这 2 个位点决定）\n")
    }
  }
}

## ---------------------------------------------------------------- 留一区域法
cat("\n\n==================== 敏感性：剔除关键位点区域（模型 1）====================\n")
REGIONS <- list(ABCG8 = c(2, 43350000, 43950000),
                UGT1A = c(2, 233300000, 234200000),
                FUT2  = c(19, 48400000, 49000000))
for (rg in names(REGIONS)) {
  v <- REGIONS[[rg]]
  base_d <- m[rsid %in% IV_MAIN & !is.na(b_cholecyst)]
  d <- base_d[!(chr == v[1] & pos >= v[2] & pos <= v[3])]
  cat(sprintf("\n-- 剔除 %s（chr%d:%.2f–%.2f Mb；去掉 %d 个 SNP，剩 %d）--\n",
              rg, v[1], v[2] / 1e6, v[3] / 1e6,
              nrow(base_d) - nrow(d), nrow(d)))
  BX <- as.matrix(d[, .(b_cholelith, b_cholecyst)])
  mv <- format_mvmr(BX, d$b_btc, as.matrix(d[, .(s_cholelith, s_cholecyst)]),
                    d$s_btc, d$rsid)
  ivw <- tryCatch(ivw_mvmr(mv), error = function(e) NULL)
  if (!is.null(ivw)) {
    rr <- est_row(ivw[, 1], ivw[, 2])
    print(data.table(exposure = c("cholelith", "cholecyst"), or = round(rr$or, 3),
                     CI = sprintf("%.3f–%.3f", rr$lci, rr$uci),
                     pval = signif(rr$pval, 3)))
  }
}
sink()
cat(sprintf("\n完成。日志: %s\n", LOG))
