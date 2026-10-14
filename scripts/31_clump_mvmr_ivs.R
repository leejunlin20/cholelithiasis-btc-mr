#!/usr/bin/env Rscript
# 31_clump_mvmr_ivs.R —— MVMR 工具变量构造：各暴露独立 clump → 取并集 → 输出观察名单
#
# 为什么必须取并集
# ----------------
# MVMR 的识别条件要求「所有暴露共用同一套工具变量」。若只用胆石症的 IV，
# 胆囊炎侧的条件 F 统计量会极低（08_mvmr_finngen.R 里对 BMI 实测只有 0.14），
# 此时胆囊炎的直接效应估计不可信。
#
# 为什么对并集再 clump 一次
# ------------------------
# 三个暴露各自的 IV 之间可能存在 LD（同一单倍型块内两性状各有显著位点）。
# 对并集按「跨暴露最小 P」去冗余，可保证最终 SNP 集近似独立，避免共线性
# 把 MVMR 的信息矩阵拉成病态。
#
# 用法: Rscript scripts/31_clump_mvmr_ivs.R

suppressMessages({
  library(data.table)
  library(ieugwasr)
})

BASE <- Sys.getenv("MR_BASE", unset = ".")
OUT  <- file.path(BASE, "results/mvmr")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EXPOS <- list(
  cholelith   = list(f = file.path(BASE, "results/medi/chole_iv_raw.tsv"),
                     label = "Cholelithiasis (K11_CHOLELITH, K80)"),
  cholecyst   = list(f = file.path(OUT,  "iv_cholecyst_raw.tsv"),
                     label = "Cholecystitis (K11_CHOLECYST, K81)"),
  gallblaoth  = list(f = file.path(OUT,  "iv_gallblaoth_raw.tsv"),
                     label = "Other gallbladder disease (K11_GALLBLAOTH, K82; incl. K82.8 polyp)"),
  cholecstec  = list(f = file.path(OUT,  "iv_cholecstec_raw.tsv"),
                     label = "Cholecystectomy (K11_CHOLECYSTECTOMY; 替代校正项)")
)

rd <- function(f) {
  if (!file.exists(f)) { cat(sprintf("  （缺文件，跳过）%s\n", f)); return(NULL) }
  d <- fread(f)
  stopifnot(all(c("rsid", "pval", "beta", "se", "alt", "ref") %in% names(d)))
  d[, .(rsid, chr = as.character(chr), pos, ref, alt,
        beta = as.numeric(beta), se = as.numeric(se),
        pval = as.numeric(pval), eaf = as.numeric(eaf))]
}

clump1 <- function(d, tag) {
  cl <- tryCatch(
    ld_clump(data.table(rsid = d$rsid, pval = d$pval, id = tag),
             clump_r2 = 0.001, clump_kb = 10000, clump_p = 5e-8, pop = "EUR"),
    error = function(e) { cat("  clump 失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(cl)) return(NULL)
  if (!nrow(cl)) {
    cat(sprintf("  !! %s 在 LD 参考面板中无可用独立位点（候选 %d 个全部被剔除）\n",
                tag, nrow(d)))
    return(d[0])
  }
  d[rsid %in% cl$rsid]
}

cat("======================= 各暴露独立 clump =======================\n")
ivs <- list()
for (k in names(EXPOS)) {
  d <- rd(EXPOS[[k]]$f)
  if (is.null(d)) next
  if (nrow(d) == 0) { cat(sprintf("%-12s 候选 0 个，跳过\n", k)); next }
  iv <- clump1(d, k)
  if (is.null(iv) || !nrow(iv)) {
    cat(sprintf("%-12s 原始 %6d → clump 后 0 个可用 IV（该暴露无法进入 MVMR）\n",
                k, nrow(d)))
    next
  }
  ivs[[k]] <- iv
  cat(sprintf("%-12s 原始 %6d → clump 后 %4d 个 IV | 最强 P=%.3g (%s)\n",
              k, nrow(d), nrow(iv), min(iv$pval),
              iv[which.min(pval), paste0(alt, " @", chr, ":", pos)]))
  fwrite(iv, file.path(OUT, sprintf("iv_%s_clumped.tsv", k)), sep = "\t")
}

cat("\n======================= 并集构造 =======================\n")
all_iv <- rbindlist(lapply(names(ivs), function(k) {
  x <- copy(ivs[[k]]); x[, src := k]; x
}))
# 同一 rsid 可能在两个暴露里都显著 → 保留跨暴露最小 P（后续按此 P clump）
uni <- all_iv[, .(pval = min(pval), src = paste(sort(unique(src)), collapse = "+")),
              by = .(rsid)]
uni <- merge(uni, all_iv[, .SD[which.min(pval)][1],
                         by = rsid][, .(rsid, chr, pos, ref, alt, beta, se, eaf)],
             by = "rsid", all.x = TRUE)
cat(sprintf("三暴露 IV 并集（去重后）: %d 个 rsID\n", nrow(uni)))

cl_u <- tryCatch(
  ld_clump(data.table(rsid = uni$rsid, pval = uni$pval, id = "union"),
           clump_r2 = 0.001, clump_kb = 10000, clump_p = 5e-8, pop = "EUR"),
  error = function(e) { cat("并集 clump 失败:", conditionMessage(e), "\n"); NULL })
if (is.null(cl_u)) stop("并集 clump 无结果")
uni <- uni[rsid %in% cl_u$rsid]
cat(sprintf("并集再 clump（r2<0.001, 10Mb）后: %d 个 SNP ← MVMR 工具变量集\n", nrow(uni)))

cat("\n各暴露在该并集中的贡献：\n")
for (k in names(ivs)) {
  n <- sum(grepl(k, uni$src, fixed = TRUE))
  cat(sprintf("  %-12s 持有 IV %4d 个，其中进入并集 %4d 个\n",
              k, nrow(ivs[[k]]), n))
}

fwrite(uni, file.path(OUT, "mvmr_union_clumped.tsv"), sep = "\t")
writeLines(uni$rsid, file.path(OUT, "mvmr_union_rsids.txt"))

## ---------------------------------------------------------------------------
## BMI 侧 IV：08_mvmr_finngen.R 实测「只用胆石 IV 时 BMI 条件 F 仅 0.14」，
## 即 BMI 的效应在模型里不可识别。因此 BMI 敏感性模型必须把 BMI 自身 IV 并进来。
## 只抽一次观察名单（含 BMI IV），三个暴露的模型在 R 里按 rsid 取子集即可。
## ---------------------------------------------------------------------------
cat("\n======================= BMI 工具变量（并入观察名单）=======================\n")
bmi_f <- file.path(BASE, "results/medi/bmi_iv_raw.tsv")
bmi_cache <- file.path(OUT, "iv_bmi_clumped.tsv")
bmi_iv <- NULL
if (file.exists(bmi_cache)) {
  ## BMI 候选位点 5.5 万个，走 API clump 较慢 → 缓存复用
  bmi_iv <- fread(bmi_cache)
  cat(sprintf("BMI_IRN 复用缓存 %s（%d 个 IV）\n", basename(bmi_cache), nrow(bmi_iv)))
} else if (file.exists(bmi_f)) {
  b <- rd(bmi_f)
  bmi_iv <- clump1(b, "bmi")
  cat(sprintf("BMI_IRN 原始 %6d → clump 后 %4d 个 IV\n", nrow(b), nrow(bmi_iv)))
  fwrite(bmi_iv, bmi_cache, sep = "\t")
} else {
  cat("未找到 bmi_iv_raw.tsv，跳过 BMI 敏感性模型\n")
}

## ---------------------------------------------------------------------------
## 含 BMI 模型的工具变量集：把「全部四暴露的候选位点」放在一起联合剪枝。
## 不能简单地把「三性状并集」与「BMI 自己 clump 出的 IV」直接相加——两个子集
## 各自独立，但合起来可能存在互相 LD 的位点，会让 MVMR 信息矩阵病态。
## ---------------------------------------------------------------------------
cat("\n======================= 含 BMI 模型的联合剪枝 =======================\n")
if (!is.null(bmi_iv)) {
  cand <- rbindlist(list(all_iv[, .(rsid, pval)],
                         bmi_iv[, .(rsid, pval = as.numeric(pval))]))
  cand <- cand[, .(pval = min(pval)), by = rsid]
  cl_b <- tryCatch(
    ld_clump(data.table(rsid = cand$rsid, pval = cand$pval, id = "union_bmi"),
             clump_r2 = 0.001, clump_kb = 10000, clump_p = 5e-8, pop = "EUR"),
    error = function(e) NULL)
  if (!is.null(cl_b)) {
    writeLines(cl_b$rsid, file.path(OUT, "mvmr_union_with_bmi_rsids.txt"))
    cat(sprintf("四暴露候选 %d → 联合剪枝后 %d 个（含 BMI 模型用）\n",
                nrow(cand), nrow(cl_b)))
  }
}

## 观察名单必须包含：
##   ① 并集（MVMR 主模型用）
##   ② BMI 自身 IV（否则 BMI 条件 F 塌掉）
##   ③ **各暴露自己的先导位点**——它们可能在并集再 clump 时被同 LD 窗口的
##      另一性状位点代表而落选，但「恰好识别」等敏感性分析仍需要其效应值
##   ④ 含 BMI 模型的联合剪枝集（其代表位点可能与 ① 不同）
watch <- unique(c(uni$rsid,
                  if (!is.null(bmi_iv)) bmi_iv$rsid else character(0),
                  unlist(lapply(ivs, function(x) x$rsid)),
                  if (file.exists(file.path(OUT, "mvmr_union_with_bmi_rsids.txt")))
                    readLines(file.path(OUT, "mvmr_union_with_bmi_rsids.txt")) else
                    character(0)))
writeLines(watch, file.path(OUT, "watch_all_rsids.txt"))
cat(sprintf("\n观察名单（并集 ∪ BMI IV ∪ 各暴露先导位点）: %d 个\n", length(watch)))
cat(sprintf("  mvmr_union_rsids.txt        %d 个（MVMR 主模型）\n", nrow(uni)))
cat(sprintf("  watch_all_rsids.txt         %d 个（抽取用）\n", length(watch)))
for (k in names(ivs)) cat(sprintf("    · %-12s 先导位点 %d 个已并入名单\n", k, nrow(ivs[[k]])))
