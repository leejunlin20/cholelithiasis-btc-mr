#!/usr/bin/env Rscript
# 34_diagnose_identification.R —— 识别性诊断：胆囊炎 vs 胆石症的遗传共线性
#
# 问题：MVMR 校正胆囊炎后，两暴露的条件 F 均 ≈ 1（远低于 10 的经验阈值）。
# 必须判定这是「真实的不可识别」还是「调用/列序错误」。
#
# 判定方法：
#   ① 手算 Sanderson 条件 F 的分子分母，看能否复现包的输出
#   ② 检验 β_胆石 与 β_胆囊炎 在 86 个工具变量上的关系（斜率/相关/条件残差）
#   ③ 分区域（ABCG8 / UGT1A / 其余）分别拟合，看比例关系是否全基因组一致
#   ④ 对照实验：换成「校正 BMI」的模型，条件 F 应为大值

suppressMessages({library(data.table); library(MVMR)})

BASE <- Sys.getenv("MR_BASE", unset = ".")
OUT  <- file.path(BASE, "results/mvmr")
sink(file.path(OUT, "identifiability_log.txt"), split = TRUE)

IV_MAIN <- readLines(file.path(OUT, "mvmr_union_rsids.txt"))
IV_BMI  <- readLines(file.path(OUT, "mvmr_union_with_bmi_rsids.txt"))

rdw <- function(f, pre, withaf = TRUE) {
  d <- fread(file.path(OUT, f))
  setnames(d, c("beta", "se", "pval", "eaf"),
           paste0(pre, c("b", "s", "p", "af")), skip_absent = TRUE)
  d
}

cl <- rdw("watch_cholelith.tsv", "ch_")
cc <- rdw("watch_cholecyst.tsv", "cc_")
bt <- rdw("watch_btc.tsv", "bt_")
bm <- rdw("watch_bmi.tsv", "bm_")

i <- match(cl$rsid, cc$rsid)                       # 以胆石为基准对齐
ok <- cc$ref[i] == cl$ref & cc$alt[i] == cl$alt
fl <- cc$ref[i] == cl$alt & cc$alt[i] == cl$ref
stopifnot(all(ok | fl))
cc_b <- ifelse(ok, cc$cc_b[i], -cc$cc_b[i])

m <- data.table(rsid = cl$rsid, chr = cl$chr, pos = cl$pos,
                ch_b = cl$ch_b, ch_s = cl$ch_s, ch_p = cl$ch_p,
                cc_b = cc_b, cc_s = cc$cc_s[i], cc_p = cc$cc_p[i])
m <- merge(m, bt[, .(rsid, bt_b, bt_s, bt_p)], by = "rsid")
m <- merge(m, bm[, .(rsid, bm_b, bm_s, bm_p)], by = "rsid")
mM <- m[rsid %in% IV_MAIN]

cat("================== ① 手算对照包输出（86 个工具变量）==================\n")
delta <- coef(lm(ch_b ~ -1 + cc_b, data = mM))[1]
mM[, rz := (ch_b - delta * cc_b) / sqrt(ch_s^2 + delta^2 * cc_s^2)]
cat(sprintf("无截距 OLS 斜率 delta = %.4f\n", delta))
cat(sprintf("手算 Q_strength = %.6f\n", mean(mM$rz^2)))
mv <- format_mvmr(as.matrix(mM[, .(ch_b, cc_b)]), mM$bt_b,
                  as.matrix(mM[, .(ch_s, cc_s)]), mM$bt_s, mM$rsid)
cat("包输出 strength_mvmr：\n"); print(strength_mvmr(mv, gencov = 0)["F-statistic", ])
cat("\n→ 两者一致即证明调用与列序无误，条件 F ≈ 1 是真实结论。\n")

cat("\n================== ② 两性状效应的共线性 ===================\n")
w <- 1 / (mM$ch_s^2 + mM$cc_s^2)
mx <- sum(w * mM$ch_b) / sum(w); my <- sum(w * mM$cc_b) / sum(w)
cat(sprintf("Pearson r（未加权）        = %.4f\n", cor(mM$ch_b, mM$cc_b)))
cat(sprintf("加权 r                    = %.4f\n",
            sum(w * (mM$ch_b - mx) * (mM$cc_b - my)) /
              sqrt(sum(w * (mM$ch_b - mx)^2) * sum(w * (mM$cc_b - my)^2))))
f0 <- lm(ch_b ~ -1 + cc_b, data = mM); f1 <- lm(ch_b ~ cc_b, data = mM)
cat(sprintf("无截距: 斜率 %.4f (SE %.4f), R2 = %.4f\n",
            coef(f0)[1], summary(f0)$coefficients[2], summary(f0)$r.squared))
cat(sprintf("有截距: 斜率 %.4f (SE %.4f), 截距 %.4f (P = %.3g), R2 = %.4f\n",
            coef(f1)[2], summary(f1)$coefficients[4], coef(f1)[1],
            summary(f1)$coefficients[8], summary(f1)$r.squared))
cat(sprintf("条件残差 z（delta 固定为 %.2f）：中位 %.2f，|z|>3 的位点 %d / %d\n",
            delta, median(mM$rz), sum(abs(mM$rz) > 3), nrow(mM)))

cat("\n================== ③ 分区拟合 ===================\n")
inABC <- mM$chr == 2 & mM$pos >= 43300000 & mM$pos <= 44000000
inUGT <- mM$chr == 2 & mM$pos >= 233300000 & mM$pos <= 234200000
groups <- list("ABCG5/8 区（chr2:43.3-44.0Mb）" = inABC,
               "UGT1A 区（chr2:233.3-234.2Mb）" = inUGT,
               "其余位点"                       = !inABC & !inUGT,
               "全部 86 个"                      = rep(TRUE, nrow(mM)))
for (g in names(groups)) {
  s <- mM[groups[[g]]]
  if (nrow(s) < 3) { cat(sprintf("%-30s n=%d（过少，跳过）\n", g, nrow(s))); next }
  fs <- lm(ch_b ~ -1 + cc_b, data = s)
  cat(sprintf("%-30s n=%2d  斜率=%.3f (SE %.3f)  R2=%.3f  r=%.3f\n",
              g, nrow(s), coef(fs)[1], summary(fs)$coefficients[2],
              summary(fs)$r.squared, cor(s$ch_b, s$cc_b)))
}
cat(sprintf("\n86 个位点中胆囊炎达全基因组显著（P<5e-8）的：%d 个 → %s\n",
            sum(mM$cc_p < 5e-8), paste(mM[cc_p < 5e-8]$rsid, collapse = ", ")))
cat(sprintf("达 P<1e-3 的：%d 个\n", sum(mM$cc_p < 1e-3)))

cat("\n================== ④ 对照：换成校正 BMI 的模型 ===================\n")
mbT <- m[rsid %in% IV_BMI & !is.na(bm_b)]
cat(sprintf("工具变量 %d 个\n", nrow(mbT)))
mvB <- format_mvmr(as.matrix(mbT[, .(ch_b, bm_b)]), mbT$bt_b,
                   as.matrix(mbT[, .(ch_s, bm_s)]), mbT$bt_s, mbT$rsid)
cat("[对照 A] 胆石 + BMI：\n"); print(strength_mvmr(mvB, gencov = 0)["F-statistic", ])
cat("\n[检验 B] 胆石 + 胆囊炎（同规模工具变量集）：\n")
mvC <- format_mvmr(as.matrix(mbT[, .(ch_b, cc_b)]), mbT$bt_b,
                   as.matrix(mbT[, .(ch_s, cc_s)]), mbT$bt_s, mbT$rsid)
print(strength_mvmr(mvC, gencov = 0)["F-statistic", ])
cat("\n→ 若 A 为大值而 B 仍 ≈1，说明诊断量本身有效，失败的是胆囊炎模型的识别性。\n")

fwrite(mM, file.path(OUT, "identifiability_input86.tsv"), sep = "\t")
sink()
cat("完成\n")
