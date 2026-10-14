#!/usr/bin/env Rscript
# 41_mvmr_mediation.R —— MVMR 直接效应：(饮食暴露 + 胆囊结石病) → 胆道癌
#
# 目的（37_mediation.R 的补角）
#   37 号脚本用**差值法**（直接效应 = 总效应 − 间接效应）估计直接效应；
#   该差值法的 SE 依赖 a 与 c 的抽样协方差，而二者共用同一组饮食 IV，
#   在汇总数据层面无法把该协方差钉死（只能做 ρ 网格敏感性）。
#   多变量 MR（MVMR）则**直接把两个暴露同时放进模型**，从而：
#     · 直接给出饮食的**条件（直接）效应**，无需做差；
#     · 给出每个暴露的**条件 F 统计量** —— 这正是本稿方法学讨论的核心诊断量。
#
# 工具变量集合
#   饮食 IV（UKB，P<5e-8）× 胆囊结石病 IV（FinnGen，84 个 clumped）的**并集**，
#   再做**联合** LD clumping（同一个 r²/kb 阈值，避免跨集合残留 LD）。
#   ⚠️ 两套 IV 在本数据中**完全不相交**（0/84），因此条件 F 有望较高 —— 与
#      校正胆囊炎那次（胆囊炎 IV 全被胆石 IV 覆盖 → 条件 F ≈ 1）形成对照。
#
# 数据获取
#   · 饮食效应（并集中属于饮食 IV 的部分）：本地 results/full/mr_detail_*.rds
#   · 饮食效应（并集中仅属于胆石 IV 的部分）：OpenGWAS associations()（ukb-b-* 就在 IEU）
#   · 胆石效应 / 胆道癌效应：本地 FinnGen 抽取结果
#
# 用法: Rscript scripts/41_mvmr_mediation.R
# 输出: results/medi2/

suppressMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MVMR)
  library(MendelianRandomization)   # mr_mvinput / mr_mvivw / mr_mvegger / mr_mvmedian / mr_mvlasso
  library(ieugwasr)
})

ROOT <- normalizePath(".")
OUT  <- file.path(ROOT, "results", "medi2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(OUT, "mvmr_mediation_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG) }

EXPS <- c("ukb-b-5779", "ukb-b-1489", "ukb-b-16576")
NICE <- c("ukb-b-5779" = "Alcohol intake frequency", "ukb-b-1489" = "Cheese intake",
          "ukb-b-16576" = "Dried fruit intake")

say("================================================================")
say("MVMR 直接效应：(饮食暴露 + 胆囊结石病) → 胆道癌")
say("运行时间: ", format(Sys.time()))
say("================================================================")

## ── 基础数据 ────────────────────────────────────────────────────────
dt1 <- as.data.table(readRDS(file.path(ROOT, "results/full/mr_detail_Cholelithiasis.rds"))$dat)[mr_keep == TRUE]
dt3 <- as.data.table(readRDS(file.path(ROOT, "results/full/mr_detail_BTC.rds"))$dat)[mr_keep == TRUE]
chole_iv <- fread(file.path(ROOT, "results/medi/chole_iv_clumped.tsv"))
btc_ch   <- fread(file.path(ROOT, "results/medi/btc_for_chole.tsv"))
btc_ch[, rsid := tstrsplit(rsids, ",", keep = 1L)]
btc_ch   <- btc_ch[!duplicated(rsid)]

# 胆石 IV 侧的 β_chole / β_btc
chole_side <- merge(
  chole_iv[, .(rsid, b_chol = beta, s_chol = se, ref = ref, alt = alt, eaf_chol = eaf)],
  btc_ch[, .(rsid, b_btc_iv = beta, s_btc_iv = sebeta, ref_btc = ref, alt_btc = alt)],
  by = "rsid"
)
say("\n[输入] 胆石 IV 表: ", nrow(chole_side), " 行（含 β_胆石 与 β_胆道癌）")

## ── 等位基因对齐辅助 ────────────────────────────────────────────────
# 把 (b, eaf) 统一到参考等位基因 (ref_ea, ref_oa) 方向
align <- function(b, eaf, ea, oa, ref_ea, ref_oa, maf_thr = 0.42) {
  n <- length(b)
  out_b <- rep(NA_real_, n); out_eaf <- rep(NA_real_, n); out_flag <- rep("drop", n)
  palind <- function(x, y) (x == "A" & y == "T") | (x == "T" & y == "A") |
                           (x == "C" & y == "G") | (x == "G" & y == "C")
  for (i in seq_len(n)) {
    if (is.na(ea[i]) | is.na(oa[i]) | is.na(ref_ea[i]) | is.na(ref_oa[i])) next
    if (ea[i] == ref_ea[i] && oa[i] == ref_oa[i]) {
      out_b[i] <- b[i]; out_eaf[i] <- eaf[i]; out_flag[i] <- "ok"
    } else if (ea[i] == ref_oa[i] && oa[i] == ref_ea[i]) {
      out_b[i] <- -b[i]; out_eaf[i] <- 1 - eaf[i]; out_flag[i] <- "flip"
    } else next
    # 回文位点：用 EAF 一致性再核一次
    if (palind(i_ea <- ea[i], oa[i])) {
      maf <- min(out_eaf[i], 1 - out_eaf[i])
      if (maf > maf_thr) out_flag[i] <- "palindromic_drop"
    }
  }
  list(b = out_b, eaf = out_eaf, flag = out_flag)
}

results <- list()

# 复用开关：REUSE_MVMR_INPUT=1 时，直接读取上一步已落盘的 mvmr_input_<eid>.csv，
# 跳过 OpenGWAS 查询与 LD clumping（单次网络查询约 2 分钟 × 3 暴露）。
# 用途：接口修好后只需重跑「拟合」环节，且保证拟合所用的工具变量集与已报告的条件 F 完全一致。
# ⚠️ 首次运行或换了暴露/结局数据后**不要**开启，否则会读到过期输入。
REUSE <- Sys.getenv("REUSE_MVMR_INPUT") == "1"
if (REUSE) say("\n[模式] REUSE_MVMR_INPUT=1 → 复用已落盘的 MVMR 输入矩阵")

for (eid in EXPS) {
  say("\n================================================================")
  say("暴露：", NICE[[eid]], " (", eid, ")")
  say("================================================================")

  inp_f <- file.path(OUT, sprintf("mvmr_input_%s.csv", eid))
  if (REUSE && file.exists(inp_f)) {
    u <- fread(inp_f)
    say("  [复用] ", inp_f, " → ", nrow(u), " 个位点（饮食 IV ", sum(u$src == "diet_iv"),
        " + 胆石 IV ", sum(u$src == "chole_iv"), "）")
  } else {

  ## 1) 饮食 IV 侧
  x1 <- dt1[id.exposure == eid]; x3 <- dt3[id.exposure == eid]
  K <- intersect(x1$SNP, x3$SNP)
  x1 <- x1[SNP %in% K][order(SNP)]; x3 <- x3[SNP %in% K][order(SNP)]
  diet_side <- data.table(
    rsid = x1$SNP,
    b_diet = x1$beta.exposure, s_diet = x1$se.exposure,
    ea_diet = x1$effect_allele.exposure, oa_diet = x1$other_allele.exposure,
    eaf_diet = x1$eaf.exposure,
    b_chol = x1$beta.outcome, s_chol = x1$se.outcome,        # 胆石症（来自 mr_detail_Cholelithiasis）
    b_btc  = x3$beta.outcome, s_btc  = x3$se.outcome,
    ref_c = x1$other_allele.outcome, alt_c = x1$effect_allele.outcome,
    ref_b = x3$other_allele.outcome, alt_b = x3$effect_allele.outcome,
    src  = "diet_iv")
  say("  饮食 IV: ", nrow(diet_side))

  ## 2) 胆石 IV 侧补上饮食效应（OpenGWAS，带磁盘缓存）
  only_chole <- setdiff(chole_side$rsid, diet_side$rsid)
  say("  仅属胆石 IV 的位点: ", length(only_chole), " → 从 OpenGWAS 取饮食效应")
  og_f <- file.path(OUT, sprintf("og_cache_%s.csv", eid))
  if (file.exists(og_f)) {
    og <- fread(og_f)
    say("  [缓存命中] ", og_f, " → ", nrow(og), " 个位点（跳过网络请求）")
  } else {
    og <- tryCatch(associations(variants = only_chole, id = eid),
                   error = function(e) { say("  [警告] OpenGWAS 查询失败: ", conditionMessage(e)); NULL })
    if (is.null(og) || !nrow(og)) {
      say("  [跳过] 无法取得饮食效应，本暴露不做 MVMR")
      next
    }
    setDT(og); og <- og[!is.na(beta) & !is.na(se)]
    fwrite(og, og_f)   # 落盘缓存：单次查询约 2 分钟，务必复用
    say("  OpenGWAS 返回: ", nrow(og), " 个位点（已缓存）")
  }
  setDT(og)
  og <- og[!is.na(beta) & !is.na(se)]

  chole_part <- merge(
    chole_side[rsid %in% og$rsid, .(rsid, b_chol, s_chol, ref_c = ref, alt_c = alt,
                                    b_btc = b_btc_iv, s_btc = s_btc_iv, ref_b = ref_btc, alt_b = alt_btc)],
    og[, .(rsid, b_diet = beta, s_diet = se, ea_diet = ea, oa_diet = nea, eaf_diet = eaf)],
    by = "rsid")
  chole_part[, src := "chole_iv"]

  u <- rbind(diet_side, chole_part, fill = TRUE)
  u <- u[!duplicated(rsid)]
  say("  并集: ", nrow(u), " 个位点（饮食 IV ", sum(u$src == "diet_iv"),
      " + 胆石 IV ", sum(u$src == "chole_iv"), "）")

  ## 3) 联合 LD clumping（关键：并集一起 clump，而非各自 clump 后拼接）
  u[, pval_clump := 5e-8]
  cl <- tryCatch(
    ld_clump(data.table(rsid = u$rsid, pval = u$pval_clump, id = "joint"),
             clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8, pop = "EUR"),
    error = function(e) { say("  [警告] 联合 clump 失败: ", conditionMessage(e)); NULL })
  if (is.null(cl) || !nrow(cl)) { say("  [跳过] 联合 clump 未产生结果"); next }
  say("  联合 clump 后: ", nrow(cl), " 个位点（移除 ", nrow(u) - nrow(cl), " 个 LD 冗余）")
  u <- u[rsid %in% cl$rsid]

  ## 4) 全部对齐到「饮食效应等位基因」方向
  a1 <- align(u$b_chol, u$eaf_diet, u$alt_c, u$ref_c, u$ea_diet, u$oa_diet)
  a2 <- align(u$b_btc,  u$eaf_diet, u$alt_b, u$ref_b, u$ea_diet, u$oa_diet)
  u[, `:=`(b_chol_h = a1$b, b_btc_h = a2$b,
           flag_chol = a1$flag, flag_btc = a2$flag)]
  say(sprintf("  对齐：胆石 直接匹配/flip/丢弃 = %d/%d/%d ；胆道癌 = %d/%d/%d",
              sum(u$flag_chol == "ok"), sum(u$flag_chol == "flip"), sum(u$flag_chol == "drop"),
              sum(u$flag_btc == "ok"),  sum(u$flag_btc == "flip"),  sum(u$flag_btc == "drop")))
  u <- u[flag_chol %in% c("ok", "flip") & flag_btc %in% c("ok", "flip")]
  say("  可用于 MVMR: ", nrow(u), " 个位点")

  fwrite(u, file.path(OUT, sprintf("mvmr_input_%s.csv", eid)))

  }  # ← 结束 else（非复用分支）

  BXM  <- as.matrix(u[, .(diet = b_diet,  chole = b_chol_h)])
  seBXM<- as.matrix(u[, .(diet = s_diet,  chole = s_chol)])
  BYv  <- as.numeric(u$b_btc_h)
  seBYv<- as.numeric(u$s_btc)

  ## 5) 条件 F（识别性诊断）—— 用 MVMR::format_mvmr 的 S3 输出
  mvmr_chk <- tryCatch(format_mvmr(as.data.frame(BXM), as.data.frame(BYv),
                                   as.data.frame(seBXM), as.data.frame(seBYv)),
                       error = function(e) { say("  [错误] format_mvmr: ", conditionMessage(e)); NULL })
  condF <- NULL
  if (!is.null(mvmr_chk)) {
    condF <- tryCatch(strength_mvmr(mvmr_chk, gencov = 0), error = function(e) NULL)
    if (!is.null(condF))
      say(sprintf("  ★ 条件 F：饮食 = %.3f | 胆囊结石病 = %.3f", condF[1, 1], condF[1, 2]))
  }

  ## 6) 四类 MVMR 估计量
  # ⚠️ 关键接口事实（2026-09-17 踩坑）：MVMR::format_mvmr() 返回的是 S3 `mvmr_format`
  #    data.frame，**只能喂给 MVMR::strength_mvmr() / pleiotropy_mvmr() / qhet_mvmr()**；
  #    而 mr_mvivw / mr_mvegger / mr_mvmedian / mr_mvlasso 是 **MendelianRandomization 的
  #    S4 泛型**，需要 mr_mvinput() 构造的 `MRMVInput` 对象。
  #    直接喂 format_mvmr 的输出会报 "unable to find an inherited method ... object = data.frame"，
  #    若用 tryCatch(error=function(e) NULL) 包住就会**静默全部失败**（本次真实发生，三个暴露全军覆没）。
  #    另：by / byse 必须是 **numeric 向量**，给 matrix 会报 invalid class "MRMVInput"。
  mvmr_in <- tryCatch(mr_mvinput(bx = BXM, bxse = seBXM, by = BYv, byse = seBYv),
                      error = function(e) { say("  [错误] mr_mvinput: ", conditionMessage(e)); NULL })
  if (is.null(mvmr_in)) next

  se_of <- function(o) {
    sn <- slotNames(o)
    for (s in c("StdError", "StdError.Est", "StdError.Int"))
      if (s %in% sn) return(slot(o, s))
    rep(NA_real_, length(o@Estimate))
  }
  fits <- list(); n_fail <- 0L
  run_est <- function(lbl, expr) {
    r <- tryCatch(expr, error = function(e) {
      n_fail <<- n_fail + 1L
      say("  [估计量失败] ", lbl, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) fits[[lbl]] <<- r
  }
  run_est("MVMR-IVW",    mr_mvivw(mvmr_in))
  run_est("MVMR-Egger",  mr_mvegger(mvmr_in, orientate = 1))
  run_est("MVMR-median", mr_mvmedian(mvmr_in, iterations = 10000))
  run_est("MVMR-Lasso",  mr_mvlasso(mvmr_in, orientate = 1))
  say("  成功估计量: ", length(fits), " / 4", if (n_fail) paste0("（失败 ", n_fail, " 个）") else "")
  if (!length(fits)) { say("  [跳过] 四个估计量全部失败"); next }

  for (nm in names(fits)) {
    o <- fits[[nm]]; se <- se_of(o)
    est <- as.numeric(o@Estimate)
    rows <- data.table(exposure_id = eid, exposure = NICE[[eid]], model = nm,
                       nSNP = nrow(u),
                       condF_diet = if (is.null(condF)) NA_real_ else condF[1, 1],
                       condF_chole = if (is.null(condF)) NA_real_ else condF[1, 2],
                       term = c("diet", "cholelithiasis"),
                       b = est, se = se,
                       p = 2 * pnorm(-abs(est / se)))
    results[[length(results) + 1]] <- rows
  }
}

## ── 汇总 ────────────────────────────────────────────────────────────
if (!length(results)) { say("\n[终止] 没有任何暴露完成 MVMR"); close(LOG); quit(status = 0) }
res <- rbindlist(results, fill = TRUE)
res[, `:=`(or = exp(b), lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se))]
fwrite(res, file.path(OUT, "mvmr_mediation_results.csv"))

say("\n================================================================")
say("MVMR 结果汇总：(饮食 + 胆囊结石病) → 胆道癌")
say("================================================================")
for (eid in unique(res$exposure_id)) {
  r <- res[exposure_id == eid]
  say(sprintf("\n【%s】nSNP = %d | 条件 F：饮食 %.2f，胆石 %.2f",
              NICE[[eid]], r$nSNP[1], r$condF_diet[1], r$condF_chole[1]))
  for (i in seq_len(nrow(r))) {
    z <- r[i]
    say(sprintf("   %-11s %-11s  OR %8.3f (%7.3f–%8.3f)  P = %.3g",
                z$model, z$term, z$or, z$lci, z$uci, z$p))
  }
}

## ── 与差值法对照 ────────────────────────────────────────────────────
say("\n================================================================")
say("直接效应：MVMR vs 差值法（37_mediation.R）")
say("================================================================")
diff_tab <- fread(file.path(OUT, "mediation_estimators.csv"))[estimator == "IVW-RE"]
for (eid in unique(res$exposure_id)) {
  mv <- res[exposure_id == eid & model == "MVMR-IVW" & term == "diet"]
  df <- diff_tab[id == eid]
  if (!nrow(mv) || !nrow(df)) next
  say(sprintf("  %-26s MVMR 直接 = %+.4f (SE %.4f) | 差值法直接 = %+.4f (SE %.4f) | 差 = %+.4f",
              NICE[[eid]], mv$b[1], mv$se[1], df$direct[1], df$se_direct[1],
              mv$b[1] - df$direct[1]))
}

close(LOG)
cat("\n完成。日志:", file.path(OUT, "mvmr_mediation_log.txt"), "\n")
