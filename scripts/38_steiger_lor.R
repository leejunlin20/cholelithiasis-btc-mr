#!/usr/bin/env Rscript
# 38_steiger_lor.R —— 用 get_r_from_lor() 重算 Steiger 方向性检验
#
# ── 为什么要重算 ──────────────────────────────────────────────────────
#   TwoSampleMR::directionality_test() 在未提供 r.exposure / r.outcome 时，
#   一律用 get_r_from_pn() 由 p 值与样本量反推相关系数 ——
#   该做法**默认所有性状都是定量性状**。本项目里
#     暴露侧：FinnGen 胆囊结石病（二分类）、FinnGen 胆道癌（二分类）
#     结局侧：同上
#   都是二分类，直接套用会把责任尺度（liability scale）的 r² 严重低估，
#   从而使方向判断系统性偏向 FALSE（即"反向因果"）。
#
# ── 修正做法 ─────────────────────────────────────────────────────────
#   · 二分类性状：get_r_from_lor(lor = β, af = eaf, ncase, ncontrol, prevalence)
#   · 连续性状  ：get_r_from_bsen(β, se, n)
#   把逐 SNP 的 r.exposure / r.outcome 事先算好塞进数据框，
#   再调用 directionality_test() —— 它会直接采用这些值，不再自行反推。
#
# ── 病例/对照数来源（官方）────────────────────────────────────────────
#   FinnGen R12 官方表型清单 finngen_R12_pheno_n.tsv
#   （https://storage.googleapis.com/finngen-public-data-r12/finngen_R12_pheno_n.tsv）
#     K11_CHOLELITH               : 49,834 / 437,418
#     C3_BILIARY_GALLBLADDER_EXALLC: 2,298 / 378,749
#
# 用法: Rscript scripts/38_steiger_lor.R
# 输出: results/steiger/

suppressMessages({
  library(data.table)
  library(TwoSampleMR)
})

ROOT <- normalizePath(".")
OUT  <- file.path(ROOT, "results", "steiger")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(OUT, "steiger_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG) }

## FinnGen R12 官方病例/对照数
NC <- list(cholelith = list(case = 49834, ctrl = 437418),
           btc       = list(case =  2298, ctrl = 378749))
## 人群患病率（用于责任尺度换算；主值 + 敏感性）
PREV <- list(cholelith = c(0.05, 0.10, 0.15, 0.20),
             btc       = c(0.0005, 0.001, 0.002, 0.005))
PREV_MAIN <- list(cholelith = 0.10, btc = 0.002)

say("================================================================")
say("Steiger 方向性检验 —— 用 get_r_from_lor() 重算二分类性状")
say("运行时间: ", format(Sys.time()))
say("================================================================")
say("\n[病例数来源] FinnGen R12 finngen_R12_pheno_n.tsv（官方）")
say("  K11_CHOLELITH                : 49,834 / 437,418")
say("  C3_BILIARY_GALLBLADDER_EXALLC:  2,298 / 378,749")

## r 计算辅助
r_binary  <- function(lor, af, ncase, nctrl, prev)
  get_r_from_lor(lor = lor, af = af, ncase = ncase, ncontrol = nctrl,
                 prevalence = prev, model = "logit")
r_cont    <- function(b, se, n) get_r_from_bsen(b, se, n)

## 逐 (暴露, 结局) 的 r 计算规格
build_r <- function(d, exp_type, out_type, prev_exp, prev_out,
                    ncase_exp = NA, nctrl_exp = NA, ncase_out = NA, nctrl_out = NA) {
  d <- copy(d)
  d[, r.exposure := if (exp_type == "binary")
       r_binary(beta.exposure, eaf.exposure, ncase_exp, nctrl_exp, prev_exp)
     else r_cont(beta.exposure, se.exposure, samplesize.exposure)]
  d[, r.outcome := if (out_type == "binary")
       r_binary(beta.outcome, eaf.outcome, ncase_out, nctrl_out, prev_out)
     else r_cont(beta.outcome, se.outcome, samplesize.outcome)]
  d
}

run_dt <- function(d, label) {
  d <- as.data.frame(d)
  dt <- suppressMessages(directionality_test(d))
  if (is.null(dt)) stop("directionality_test 返回 NULL（缺少 p/samplesize/r 列）")
  dt <- as.data.table(dt)
  dt[, pair := label]
  dt
}

results <- list()

## ════════════ 1. 胆囊结石病 → 胆道癌（Limitations 第 6 条所指的那一个）════════════
s2 <- readRDS(file.path(ROOT, "results/medi/step2_detail.rds"))$dat
s2 <- as.data.table(s2)
s2 <- s2[mr_keep == TRUE]
s2[, samplesize.outcome := NULL]
s2[, samplesize.exposure := as.numeric(NC$cholelith$case + NC$cholelith$ctrl)]
s2[, samplesize.outcome  := as.numeric(NC$btc$case + NC$btc$ctrl)]
s2[, pval.exposure := as.numeric(pval.exposure)]
s2[, pval.outcome  := as.numeric(pval.outcome)]
say("\n[1] 胆囊结石病 → 胆道癌：n = ", nrow(s2), " 个 IV")

for (p in PREV$cholelith) {
  for (q in PREV$btc) {
    dd <- build_r(s2, "binary", "binary", p, q,
                  NC$cholelith$case, NC$cholelith$ctrl, NC$btc$case, NC$btc$ctrl)
    dt <- run_dt(dd, "Cholelithiasis → BTC")
    dt[, `:=`(prev.exp = p, prev.out = q)]
    results[[length(results) + 1]] <- dt
  }
}
# 对照：原方法（get_r_from_pn）—— 直接读取既有输出，因该方法对 p=0 的位点会直接报错
old <- fread(file.path(ROOT, "results/medi/step2_steiger.csv"))
dt0 <- data.table(pair = "Cholelithiasis → BTC",
                  method = "原方法 get_r_from_pn()",
                  prev.exp = NA_real_, prev.out = NA_real_,
                  snp_r2.exposure = as.numeric(old$snp_r2.exposure[1]),
                  snp_r2.outcome  = as.numeric(old$snp_r2.outcome[1]),
                  correct_causal_direction = as.logical(old$correct_causal_direction[1]),
                  steiger_pval = as.numeric(old$steiger_pval[1]))
say(sprintf("  [原方法] r²_exp = %.5f, r²_out = %.5f, 方向正确 = %s, P = %.4f",
            dt0$snp_r2.exposure[1], dt0$snp_r2.outcome[1],
            dt0$correct_causal_direction[1], dt0$steiger_pval[1]))
results[[length(results) + 1]] <- dt0

## ════════════ 2. 饮食 → 胆囊结石病 / 胆道癌 ════════════
d1 <- as.data.table(readRDS(file.path(ROOT, "results/full/mr_detail_Cholelithiasis.rds"))$dat)
d3 <- as.data.table(readRDS(file.path(ROOT, "results/full/mr_detail_BTC.rds"))$dat)
d1 <- d1[mr_keep == TRUE]; d3 <- d3[mr_keep == TRUE]
# ⚠️ data.table 无法用 := 就地改变列类型（logical→double 会被静默强转为 TRUE/FALSE），
#    必须先删除该列再新建。这一步踩过坑：samplesize.outcome 原为逻辑型 NA，
#    直接赋值导致 psych::r.test 收到 n2=1 → sqrt 负数 → P 全变 NaN。
d1[, samplesize.outcome := NULL]; d3[, samplesize.outcome := NULL]
d1[, samplesize.outcome := as.numeric(NC$cholelith$case + NC$cholelith$ctrl)]
d3[, samplesize.outcome := as.numeric(NC$btc$case + NC$btc$ctrl)]
d1[, samplesize.exposure := as.numeric(samplesize.exposure)]
d3[, samplesize.exposure := as.numeric(samplesize.exposure)]
d1[, pval.exposure := as.numeric(pval.exposure)]; d3[, pval.exposure := as.numeric(pval.exposure)]
d1[, pval.outcome  := as.numeric(pval.outcome)];  d3[, pval.outcome  := as.numeric(pval.outcome)]

PRIMARY <- c("ukb-b-5779", "ukb-b-1489", "ukb-b-16576")
NICE <- c("ukb-b-5779" = "Alcohol intake frequency", "ukb-b-1489" = "Cheese intake",
          "ukb-b-16576" = "Dried fruit intake")

run_subset <- function(d, ids, out_type, nc, prev_out, label) {
  rbindlist(lapply(ids, function(id) {
    x <- d[id.exposure == id]
    if (!nrow(x)) return(NULL)
    dd <- build_r(x, "continuous", out_type, NA, prev_out,
                  NA, NA, nc$case, nc$ctrl)
    dt <- run_dt(dd, label)
    dt[, prev.out := prev_out]
    dt
  }), fill = TRUE)
}

say("\n[2] 饮食（连续）→ 胆囊结石病（二分类）")
for (id in PRIMARY) {
  x <- d1[id.exposure == id]
  dd <- build_r(x, "continuous", "binary", NA, PREV_MAIN$cholelith,
                NA, NA, NC$cholelith$case, NC$cholelith$ctrl)
  dt <- run_dt(dd, "Diet → Cholelithiasis")
  dt[, prev.out := PREV_MAIN$cholelith]
  results[[length(results) + 1]] <- dt
  say(sprintf("  %-26s r²_exp=%.5f r²_out=%.5f 方向正确=%s P=%.4f",
              NICE[[id]], dt$snp_r2.exposure[1], dt$snp_r2.outcome[1],
              dt$correct_causal_direction[1], dt$steiger_pval[1]))
}

say("\n[3] 饮食（连续）→ 胆道癌（二分类）")
for (id in PRIMARY) {
  x <- d3[id.exposure == id]
  dd <- build_r(x, "continuous", "binary", NA, PREV_MAIN$btc,
                NA, NA, NC$btc$case, NC$btc$ctrl)
  dt <- run_dt(dd, "Diet → BTC")
  dt[, prev.out := PREV_MAIN$btc]
  results[[length(results) + 1]] <- dt
  say(sprintf("  %-26s r²_exp=%.5f r²_out=%.5f 方向正确=%s P=%.4f",
              NICE[[id]], dt$snp_r2.exposure[1], dt$snp_r2.outcome[1],
              dt$correct_causal_direction[1], dt$steiger_pval[1]))
}

## ════════════ 3. 汇总输出 ════════════
res <- rbindlist(results, fill = TRUE)
res[, method := fifelse(!is.na(method), method, "get_r_from_lor() 修正")]
fwrite(res[, .(pair, method, prev_exp = prev.exp, prev_out = prev.out,
               r2_exp = snp_r2.exposure, r2_out = snp_r2.outcome,
               ratio = snp_r2.exposure / snp_r2.outcome,
               correct_direction = correct_causal_direction, steiger_p = steiger_pval)],
       file.path(OUT, "steiger_lor_results.csv"))

say("\n================================================================")
say("主对照：胆囊结石病 → 胆道癌")
say("================================================================")
rows0 <- res[method == "原方法 get_r_from_pn()"]
row <- res[pair == "Cholelithiasis → BTC" & !is.na(prev.exp) &
           abs(prev.exp - 0.10) < 1e-9 & abs(prev.out - 0.002) < 1e-9]
say(sprintf("  原方法 (get_r_from_pn) : r²暴露 = %.5f, r²结局 = %.5f, 比值 = %.3f, 方向正确 = %s, P = %.4f",
            rows0$snp_r2.exposure[1], rows0$snp_r2.outcome[1],
            rows0$snp_r2.exposure[1] / rows0$snp_r2.outcome[1],
            rows0$correct_causal_direction[1], rows0$steiger_pval[1]))
say(sprintf("  修正后 (get_r_from_lor): r²暴露 = %.5f, r²结局 = %.5f, 比值 = %.3f, 方向正确 = %s, P = %.4f",
            row$snp_r2.exposure[1], row$snp_r2.outcome[1],
            row$snp_r2.exposure[1] / row$snp_r2.outcome[1],
            row$correct_causal_direction[1], row$steiger_pval[1]))

say("\n--- prevalence 敏感性（胆囊结石病 → 胆道癌）---")
say(sprintf("  %-8s %-8s %-12s %-12s %-10s %-8s", "prev_exp", "prev_out", "r2_exp", "r2_out", "方向正确", "P"))
sub <- res[pair == "Cholelithiasis → BTC" & !is.na(res$prev.exp)]
for (i in seq_len(nrow(sub))) {
  r <- sub[i]
  say(sprintf("  %-8.3f %-8.4f %-12.5f %-12.5f %-10s %-8.4f",
              r$prev.exp, r$prev.out, r$snp_r2.exposure, r$snp_r2.outcome,
              r$correct_causal_direction, r$steiger_pval))
}

say("\n--- 交叉核对：手动 Fisher-z 差异检验（不依赖 psych 的样本量假设）---")
manual_z <- function(r2e, r2o, n_e, n_o) {
  ze <- atanh(sqrt(r2e)); zo <- atanh(sqrt(r2o))
  se <- sqrt(1 / (n_e - 3) + 1 / (n_o - 3))
  z <- (ze - zo) / se
  c(z = z, p = 2 * pnorm(-abs(z)))
}
mz <- manual_z(row$snp_r2.exposure[1], row$snp_r2.outcome[1],
               NC$cholelith$case + NC$cholelith$ctrl, NC$btc$case + NC$btc$ctrl)
say(sprintf("  胆囊结石病 → 胆道癌：z = %.2f, P = %.3g（注意：两侧 r 来自同一样本、同一批 SNP，\n  故该 P 值在独立性假设下**偏乐观**，应作为方向判断的辅助而非精确 P 值）", mz["z"], mz["p"]))

say("\n--- 饮食暴露 ---")
sub2 <- res[grepl("Diet", res$pair)]
for (i in seq_len(nrow(sub2))) {
  r <- sub2[i]
  say(sprintf("  %-22s %-28s r²_exp=%.5f r²_out=%.5f 方向正确=%s P=%.4f",
              r$pair, substr(r$exposure, 1, 28), r$snp_r2.exposure, r$snp_r2.outcome,
              r$correct_causal_direction, r$steiger_pval))
}

say("\n完成 → ", OUT)
close(LOG)
cat("\n完成。日志:", file.path(OUT, "steiger_log.txt"), "\n")
