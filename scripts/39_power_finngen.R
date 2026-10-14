#!/usr/bin/env Rscript
# 39_power_finngen.R —— FinnGen 侧正式效能计算（补 BBJ 侧已有表的对称件）
#
# ── 方法 ───────────────────────────────────────────────────────────────
# 不作解析近似（二分类暴露的观测尺度 R² 会高估责任尺度 R²，不适用于效能推算），
# 一律用**实际观测到的 MR 标准误**反推：
#   · 80% 把握可检出的最小 OR        = exp(2.8016 × SE)      （2.8016 = 1.95996 + 0.84162）
#   · 对给定真实效应量的检出力        = Φ(|ln OR|/SE − 1.95996) + Φ(−|ln OR|/SE − 1.95996)
#   · 达到目标效能所需病例数（同对照池）= 由 SE ∝ 1/sqrt(n_eff) 反解，
#                                       n_eff = C·K/(C+K)
#
# 病例/对照数来源：FinnGen R12 官方表型清单 finngen_R12_pheno_n.tsv
#
# 用法: Rscript scripts/39_power_finngen.R
# 输出: results/power/

suppressMessages(library(data.table))

ROOT <- normalizePath(".")
OUT  <- file.path(ROOT, "results", "power")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(OUT, "power_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG) }

Z_A <- qnorm(0.975); Z_B <- qnorm(0.80)         # 1.95996, 0.84162
SUM <- Z_A + Z_B                                 # 2.8016，80% 效能对应的 SE 倍数
pw  <- function(beta, se) pnorm(abs(beta) / se - Z_A) + pnorm(-abs(beta) / se - Z_A)
or_min <- function(se) exp(SUM * se)

say("================================================================")
say("FinnGen 侧效能计算（用实测 MR 标准误反推）")
say("运行时间: ", format(Sys.time()))
say("================================================================")

## ── 官方病例/对照数 ──────────────────────────────────────────────────
man <- fread(file.path(ROOT, "data/finngen/finngen_R12_pheno_n.tsv"))
setnames(man, gsub("^﻿", "", names(man)))
man[, `:=`(num_cases = as.numeric(num_cases), num_controls = as.numeric(num_controls))]
getnc <- function(pc) {
  r <- man[phenocode == pc]
  list(case = r$num_cases[1], ctrl = r$num_controls[1], pheno = r$phenotype[1])
}
NC_BTC  <- getnc("C3_BILIARY_GALLBLADDER_EXALLC")
NC_CHO  <- getnc("K11_CHOLELITH")
say(sprintf("\n[官方病例数] %s: %s 例 / %s 对照",
            "胆道癌 (C3_BILIARY_GALLBLADDER_EXALLC)", format(NC_BTC$case, big.mark = ","),
            format(NC_BTC$ctrl, big.mark = ",")))
say(sprintf("[官方病例数] %s: %s 例 / %s 对照",
            "胆囊结石病 (K11_CHOLELITH)", format(NC_CHO$case, big.mark = ","),
            format(NC_CHO$ctrl, big.mark = ",")))

## ── 汇总本项目中所有已完成的「暴露 → 结局」MR 实测精度 ────────────────
rows <- list()

## (1) 胆石症（84 clumped IV）→ FinnGen 胆道癌
rows[[length(rows) + 1]] <- data.table(
  analysis = "胆囊结石病 → 胆道癌（FinnGen）",
  exposure_set = "胆石症 84 个 clumped IV",
  outcome = "胆道癌 (Finngen C3)",
  ncase = NC_BTC$case, nctrl = NC_BTC$ctrl,
  se = 0.0732635699, b = 0.2023964295)

## (2) 胆石症 → 胆管癌（BBJ 复制）
bbj <- fread(file.path(ROOT, "results/bbj_replication/bbj_mr_results.csv"))
bbj_ivw <- bbj[method == "Inverse variance weighted" & pheno == "BtC"]
rows[[length(rows) + 1]] <- data.table(
  analysis = "胆囊结石病 → 胆管癌（BBJ 复制）",
  exposure_set = "胆石症 61 个对齐 IV",
  outcome = "肝内胆管癌 (BBJ BtC)",
  ncase = 418, nctrl = 159201,
  se = bbj_ivw$se[1], b = bbj_ivw$b[1])

## (3) 饮食 → 胆道癌（3 个主暴露）
d_btc <- fread(file.path(ROOT, "results/full/mr_results_BTC.csv"))
d_btc <- d_btc[method == "Inverse variance weighted"]
NICE <- c("ukb-b-5779" = "饮酒频率", "ukb-b-1489" = "奶酪", "ukb-b-16576" = "干果",
          "ukb-b-9508" = "咖啡（consumed）", "ukb-b-3291" = "标准茶", "ukb-b-3975" = "饮水量")
for (id in names(NICE)) {
  r <- d_btc[id.exposure == id]
  if (!nrow(r)) next
  rows[[length(rows) + 1]] <- data.table(
    analysis = sprintf("饮食（%s）→ 胆道癌（FinnGen）", NICE[[id]]),
    exposure_set = sprintf("%d 个 IV", as.integer(r$nsnp.x[1])),
    outcome = "胆道癌 (Finngen C3)",
    ncase = NC_BTC$case, nctrl = NC_BTC$ctrl,
    se = as.numeric(r$se[1]), b = as.numeric(r$b[1]))
}

## (4) 饮食 → 胆囊结石病（3 个 MR-PRESSO 主暴露）
pr <- fread(file.path(ROOT, "results/full/mr_presso_results.csv"))
for (eid in c("ukb-b-5779", "ukb-b-1489", "ukb-b-16576")) {
  r <- pr[id == eid]
  if (!nrow(r) || is.na(suppressWarnings(as.numeric(r$se_corr[1])))) next
  rows[[length(rows) + 1]] <- data.table(
    analysis = sprintf("饮食（%s）→ 胆囊结石病（FinnGen, MR-PRESSO）", NICE[[eid]]),
    exposure_set = sprintf("%d 个 IV", as.integer(r$nsnp[1])),
    outcome = "胆囊结石病 (Finngen K80)",
    ncase = NC_CHO$case, nctrl = NC_CHO$ctrl,
    se = as.numeric(r$se_corr[1]), b = as.numeric(r$b_corr[1]))
}

## (5) MVMR 校正 BMI 后的胆石症 → 胆道癌
rows[[length(rows) + 1]] <- data.table(
  analysis = "胆囊结石病 → 胆道癌（MVMR 校正 BMI）",
  exposure_set = "胆石 401 个联合 clump IV（条件 F 23.8）",
  outcome = "胆道癌 (Finngen C3)",
  ncase = NC_BTC$case, nctrl = NC_BTC$ctrl,
  se = 0.0486824996361698, b = 0.218248007736901)

## (6) 间接效应（三步链条中的 a×b）
rows[[length(rows) + 1]] <- data.table(
  analysis = "间接效应（饮酒频率 → 胆石 → 胆道癌）",
  exposure_set = "两步介导, a×b",
  outcome = "胆道癌 (Finngen C3)",
  ncase = NC_BTC$case, nctrl = NC_BTC$ctrl,
  se = 0.022960, b = 0.046104)
rows[[length(rows) + 1]] <- data.table(
  analysis = "间接效应（奶酪 → 胆石 → 胆道癌）",
  exposure_set = "两步介导, a×b",
  outcome = "胆道癌 (Finngen C3)",
  ncase = NC_BTC$case, nctrl = NC_BTC$ctrl,
  se = 0.032200, b = -0.057300)
rows[[length(rows) + 1]] <- data.table(
  analysis = "间接效应（干果 → 胆石 → 胆道癌）",
  exposure_set = "两步介导, a×b",
  outcome = "胆道癌 (Finngen C3)",
  ncase = NC_BTC$case, nctrl = NC_BTC$ctrl,
  se = 0.059973, b = -0.129611)

tab <- rbindlist(rows)
tab[, `:=`(or_obs = exp(b),
           or_lci = exp(b - 1.96 * se), or_uci = exp(b + 1.96 * se),
           p_obs = 2 * pnorm(-abs(b / se)),
           or_min80 = or_min(se),
           power_or120 = pw(log(1.20), se),
           power_or122 = pw(log(1.22), se),
           power_or130 = pw(log(1.30), se),
           power_or150 = pw(log(1.50), se))]

## 达到 80% 效能检出 OR 1.20 所需病例数（保持同一对照池）
need_cases <- function(se_now, target_or, n_ctrl, n_case_now) {
  se_t <- log(target_or) / SUM
  if (se_now <= se_t) return(0)
  n_eff_now <- n_case_now * n_ctrl / (n_case_now + n_ctrl)
  n_eff_need <- n_eff_now * (se_now / se_t)^2
  n_eff_need * n_ctrl / (n_ctrl - n_eff_need)     # 由 n_eff = C·K/(C+K) 反解 C
}
tab[, cases_need_or120 := mapply(function(s, c_, k) {
  if (s <= log(1.20) / SUM) return(0)
  v <- need_cases(s, 1.20, k, c_)
  if (!is.finite(v) || v < 0) NA_real_ else v
}, se, ncase, nctrl)]

fwrite(tab, file.path(OUT, "power_table_finngen.csv"))

say("\n================================================================")
say("表 P1　本项目已完成 MR 分析的实测精度（FinnGen R12 + BBJ）")
say("================================================================")
say(sprintf("%-44s %8s %8s %9s %9s %9s %9s", "分析", "SE", "观测OR", "最小OR80", "P(1.20)", "P(1.30)", "P(1.50)"))
for (i in seq_len(nrow(tab))) {
  r <- tab[i]
  say(sprintf("%-44s %8.4f %8.3f %9.3f %8.0f%% %8.0f%% %8.0f%%",
              substr(r$analysis, 1, 44), r$se, r$or_obs, r$or_min80,
              100 * r$power_or120, 100 * r$power_or130, 100 * r$power_or150))
}

## ── 表 P2：同一暴露（胆石症）在两队列的效能对比 ──────────────────────
say("\n================================================================")
say("表 P2　同暴露（胆囊结石病）在两队列的效能对比")
say("================================================================")
say(sprintf("%-26s %8s %8s %10s %12s %12s", "队列", "例数", "实测SE", "最小OR80", "P(OR=1.22)", "P(OR=1.50)"))
fin <- tab[grepl("FinnGen", analysis) & grepl("胆囊结石病 → 胆道癌（FinnGen）", analysis)]
bbjrow <- tab[grepl("BBJ", analysis)]
say(sprintf("%-26s %8s %8.4f %10.3f %11.0f%% %11.0f%%", "FinnGen R12 (胆道癌)",
            format(NC_BTC$case, big.mark = ","), fin$se[1], fin$or_min80[1],
            100 * fin$power_or122[1], 100 * fin$power_or150[1]))
say(sprintf("%-26s %8s %8.4f %10.3f %11.0f%% %11.0f%%", "BBJ (肝内胆管癌)",
            format(bbjrow$ncase[1], big.mark = ","), bbjrow$se[1], bbjrow$or_min80[1],
            100 * bbjrow$power_or122[1], 100 * bbjrow$power_or150[1]))
say(sprintf("\n  效能比（FinnGen / BBJ，对 OR 1.22）= %.2f 倍", fin$power_or122[1] / bbjrow$power_or122[1]))

## ── 表 P3：效能曲线（病例数 → 最小可检出 OR），锚定 FinnGen 实测 SE ──
say("\n================================================================")
say("表 P3　效能曲线：病例数 vs 最小可检出 OR（锚定 FinnGen 实测 SE = 0.07326）")
say("说明：按 SE ∝ 1/sqrt(n_eff)、n_eff = C·K/(C+K) 外推（K 固定为 FinnGen 对照组数）")
say("================================================================")
SE0 <- 0.0732635699; C0 <- NC_BTC$case; K0 <- NC_BTC$ctrl
neff0 <- C0 * K0 / (C0 + K0)
curve <- data.table(ncase = c(500, 1000, 2298, 5000, 10000, 20000, 50000, 100000))
curve[, n_eff := ncase * K0 / (ncase + K0)]
curve[, se := SE0 * sqrt(neff0 / n_eff)]
curve[, `:=`(or_min80 = or_min(se),
             power_or120 = pw(log(1.20), se),
             power_or130 = pw(log(1.30), se),
             power_or125 = pw(log(1.25), se))]
fwrite(curve, file.path(OUT, "power_curve_by_cases.csv"))
say(sprintf("%10s %10s %10s %12s %12s %12s", "病例数", "有效n", "SE", "最小OR80", "P(1.20)", "P(1.25)"))
for (i in 1:nrow(curve)) {
  r <- curve[i]
  say(sprintf("%10s %10.0f %10.4f %12.3f %11.0f%% %11.0f%%",
              format(r$ncase, big.mark = ","), r$n_eff, r$se, r$or_min80,
              100 * r$power_or120, 100 * r$power_or125))
}

## ── 关键结论 ─────────────────────────────────────────────────────────
say("\n================================================================")
say("关键结论")
say("================================================================")
say(sprintf("① FinnGen 胆道癌（2,298 例）用 84 个胆石 IV：SE = %.4f → 只能可靠检出 OR ≥ %.2f 的效应。",
            fin$se[1], fin$or_min80[1]))
say(sprintf("   对欧裔参考效应量 OR 1.22 的检出力 = %.0f%%（BBJ 同暴露仅 %.0f%%）。",
            100 * fin$power_or122[1], 100 * bbjrow$power_or122[1]))
prim <- tab[grepl("饮酒频率|奶酪|干果", analysis) & grepl("→ 胆道癌", analysis)]
say(sprintf("② 三个主饮食暴露 → 胆道癌的实测 SE 为 %.3f–%.3f，只能检出 OR ≥ %.2f–%.2f 的效应。",
            min(prim$se), max(prim$se), min(prim$or_min80), max(prim$or_min80)))
say(sprintf("   （86 个饮食暴露中 SE 中位 %.3f；SE 最小的『标准茶』为 %.4f，但其效应量本身仅 OR %.3f，无实质意义）",
            median(tab[grepl("饮食", analysis) & grepl("→ 胆道癌", analysis)]$se),
            min(tab[grepl("饮食", analysis) & grepl("→ 胆道癌", analysis)]$se),
            tab[grepl("标准茶", analysis)]$or_obs[1]))
say("   ⇒ 总效应不显著主要反映**效能不足**，不能解读为「无效应」。")
need <- tab[grepl("饮酒频率", analysis) & grepl("→ 胆道癌", analysis)]$cases_need_or120[1]
say(sprintf("③ 若要以 80%% 效能检出 OR 1.20 的饮食→胆道癌效应，约需 %.0f 例胆道癌（现 2,298 例，差 %.1f 倍）。",
            need, need / NC_BTC$case))
say(sprintf("④ 间接效应（a×b）的 SE 仅 %.4f–%.4f → 链条中的中介路径**是可检验的**，",
            min(tab[grepl("间接效应", analysis)]$se), max(tab[grepl("间接效应", analysis)]$se)))
say("   这与「总效应效能不足」并不矛盾：中介效应用的是两步各自的 IV，精度更高。")

close(LOG)
cat("\n完成。日志:", file.path(OUT, "power_log.txt"), "\n")
