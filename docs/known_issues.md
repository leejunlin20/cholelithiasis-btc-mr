# Known issues and version notes

This file records defects that were found and corrected during the analysis, and scripts that must not be used. It is
kept in the repository so that the code history is as transparent as the manuscript.

---

## 1. Double random-effects inflation of the IVW standard error in `03_mr_core.R`

**Status: corrected. All estimates reported in the manuscript use the corrected values.**

`03_mr_core.R` computed the inverse-variance weighted (IVW) standard error as follows (lines 96–99):

```r
out[, se_ivw_fe := se]                       # mr_ivw() already returns a random-effects SE
out[method == "Inverse variance weighted",
    se := se * pmax(1, sqrt(Q / Q_df))]      # so kappa was applied a second time
```

Because `TwoSampleMR::mr_ivw()` returns a standard error that has already been inflated by
κ = max(1, √(Q/df)), multiplying by κ again produced a standard error inflated by κ². The column named `se_ivw_fe`
was in fact the single-κ random-effects standard error, and the column `se`, together with everything derived from it
(`pval`, `lo_ci`, `up_ci`, `or_lci95`, `or_uci95`), was over-conservative.

Worked example (alcohol intake frequency → cholelithiasis, 93 variants):

| Quantity | Value |
|---|---|
| Fixed-effect SE, 1/√Σw | 0.035245 |
| Cochran Q (df 92) | 453.9 |
| κ = max(1, √(Q/df)) | 2.2213 |
| Correct single-κ random-effects SE (**used in the manuscript**) | 0.078287 |
| Double-κ SE (the defect) | 0.173896 |

**Correction.** `43_fix_kappa.R` recomputes, from the SNP-level data stored in `mr_detail_*.rds`, the fixed-effect SE,
Q, κ and the single-κ random-effects SE, then rewrites the IVW rows and recomputes the Benjamini–Hochberg FDR. It
writes `mr_results_<trait>_v2.csv`, which is the authoritative result set: every IVW estimate quoted in the manuscript
and its supplementary tables comes from these `*_v2.csv` files. The other estimators (MR-Egger, weighted median,
weighted mode) do not use this convention and were unaffected. The manuscript Methods describe the corrected
convention, and the convention was verified to reproduce `TwoSampleMR::mr_ivw()` to ten decimal places.

---

## 2. `09_make_figure.py` is deprecated

The script carries the header marker *"已作废，勿用于投稿"* (deprecated, do not use for submission). It is retained only
so that the file history is complete. The submitted figures are produced by `53_make_figure1.py`,
`54_make_figure2.py`, `56_make_figure3.py`, `55_make_figure4.py`, `61_make_figure5.py` and the corresponding
supplementary figure scripts.

---

## 3. Two mediation implementations

Both `06_medi_mr.R` and `37_mediation.R` implement the two-step mediation (diet → cholelithiasis → BTC).
`06_medi_mr.R` is the first implementation; `37_mediation.R` is the final specification whose output is reported, with
`41_mvmr_mediation.R` providing the MVMR direct-effect check. Retaining both is deliberate, because the agreement
between the differencing estimator and the MVMR estimator is itself reported in the manuscript.

---

## 4. Deliberate non-reporting of two analyses

These are documented here so that a reader does not interpret their absence as an omission:

- **Reverse-direction MR (BTC → cholelithiasis)** was attempted but is not reported because 219 genome-wide
  significant BTC variants reduced to 2 independent variants after clumping, which is insufficient to build an
  instrument set. This is stated in the manuscript.
- **Winner's-curse correction** of the cross-ancestry instrument attenuation was considered but not applied, because
  the estimated attenuation depended on the estimator chosen (zero-intercept IVW 0.667 vs intercept-inclusive weighted
  regression 0.915). This is stated in the manuscript.
