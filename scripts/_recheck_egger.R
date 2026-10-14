## 57_recheck_egger.R —— 用 TwoSampleMR 复算 84-IV 的 Egger，确认与 mvmr_univariable.csv 一致
suppressMessages(library(TwoSampleMR))
cat("TwoSampleMR", as.character(packageVersion("TwoSampleMR")), "\n")

## 注意：FinnGen 表头以 "#chrom" 开头，read.table 默认 comment.char="#" 会把表头行当注释吃掉，
## 导致列名错位（报 'by' must specify a uniquely valid column）。必须显式关掉。
rd <- function(p) read.table(p, header = TRUE, sep = "\t", comment.char = "",
                             check.names = FALSE)

iv <- rd("results/mvmr/iv_cholelith_clumped.tsv")
ou <- rd("/tmp/outcome_btc_84iv.tsv")
## 只取需要列并显式改名，避免 merge 后缀歧义（se 只在暴露侧，不会带 .e 后缀）
iv <- data.frame(rsid = iv$rsid, bx = iv$beta, sx = iv$se,
                 ae = iv$alt, oe = iv$ref, stringsAsFactors = FALSE)
ou <- data.frame(rsid = ou$rsids, by = ou$beta, sy = ou$sebeta,
                 ao = ou$alt, oo = ou$ref, stringsAsFactors = FALSE)
d <- merge(iv, ou, by = "rsid")
cat("n =", nrow(d), "\n")
cat("列名:", paste(names(d), collapse = ", "), "\n")
stopifnot(all(d$ae == d$ao), all(d$oe == d$oo))

dat <- data.frame(
  SNP = d$rsid,
  beta.exposure = d$bx, se.exposure = d$sx,
  beta.outcome = d$by, se.outcome = d$sy,
  effect_allele.exposure = d$ae, other_allele.exposure = d$oe,
  effect_allele.outcome = d$ao, other_allele.outcome = d$oo,
  exposure = "Cholelithiasis", outcome = "BTC",
  id.exposure = "Cholelithiasis", id.outcome = "BTC")

cat("\n--- mr_egger_regression 原始返回（逐个元素打印，该对象含 summary.lm 无法直接 as.data.frame）---\n")
eg <- mr_egger_regression(dat$beta.exposure, dat$beta.outcome,
                          dat$se.exposure, dat$se.outcome)
cat("class:", paste(class(eg), collapse = ", "), "\n")
cat("names:", paste(names(eg), collapse = ", "), "\n")
for (nm in names(eg)) {
  v <- eg[[nm]]
  if (is.numeric(v) && length(v) == 1) cat(sprintf("  %-28s %s\n", nm, format(v, digits = 8)))
}

cat("\n--- Egger 内部 lm 的真实权重（决定斜率的关键）---\n")
m <- eg$mod
cat("class(mod):", paste(class(m), collapse = ", "), "\n")
cat("系数:\n"); print(round(coef(m), 8))
w <- weights(m)
cat("weights(head):", paste(round(head(w, 5), 6), collapse = ", "), "\n")
cat("是否等于 1/se_out^2 :", isTRUE(all.equal(w, 1/dat$se.outcome^2)), "\n")
cat("是否等于 bx^2/se_out^2:", isTRUE(all.equal(w, dat$beta.exposure^2/dat$se.outcome^2)), "\n")
cat("是否等于 1/se_out    :", isTRUE(all.equal(w, 1/dat$se.outcome)), "\n")
cat("是否全为 1（无权）   :", isTRUE(all.equal(w, rep(1, nrow(dat)))), "\n")
cat("cor(w, 1/se_out^2) =", round(cor(w, 1/dat$se.outcome^2), 6), "\n")
cat("cor(w, bx^2/se_out^2) =", round(cor(w, dat$beta.exposure^2/dat$se.outcome^2), 6), "\n")

cat("\nDONE\n")
