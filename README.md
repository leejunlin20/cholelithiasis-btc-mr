# Cholelithiasis → biliary tract cancer: cross-ancestry MR and region-definition-aware colocalisation

Analysis code for the manuscript:

> **Cholelithiasis and biliary tract cancer risk: a cross-ancestry Mendelian randomisation and colocalisation study with region-definition diagnostics**

This repository contains every script used to produce the results reported in the manuscript and its supplementary
material. It is intended to allow full re-execution of the pipeline from publicly available summary-level data.

---

## 1. What the pipeline does

| Stage | What is computed | Manuscript location |
|---|---|---|
| A | Exposure selection and instrument extraction | Table 1, Figure 1 |
| B | Univariable two-sample MR of diet on cholelithiasis (primary estimator: MR-PRESSO outlier-corrected, BH-FDR) | Table 2 panel A, Figure 2 |
| C | Univariable two-sample MR of cholelithiasis on biliary tract cancer (BTC), four estimators | Table 2 panel B, Figure 3A–C |
| D | Multivariable MR (BMI-adjusted; cholecystitis-adjusted), with conditional F-statistics | Table 3 panel A, Figure 4A–B, Supplementary Table S1 |
| E | Bayesian colocalisation (coloc.abf, coloc.susie), region-definition diagnostics and split-region sensitivity | Table 3 panel B, Figure 4C, Supplementary Table S3 |
| F | Liver cis-eQTL × GWAS colocalisation at ABCG5/ABCG8 | Results, Colocalisation, Supplementary Figure S2 |
| G | Two-step mediation (diet → cholelithiasis → BTC) and MVMR direct-effect check | Table 2 panel C, Figure 3D |
| H | Cross-ancestry replication in BioBank Japan with positive and negative controls | Table 3 panel C, Figure 5 |
| I | Steiger directionality, power calculations, FinnGen outcome panel | Table 3 panel D, Supplementary Figures S1 and S3, Supplementary Table S4–S5 |
| J | Figure and manuscript production tooling (not part of the analysis) | — |

---

## 2. Data sources (none are redistributed here)

All input data are publicly available summary-level statistics. **No individual-level data were used and no data files
are included in this repository** (the raw inputs total ~6.8 GB).

| Data | Source | Access |
|---|---|---|
| Dietary exposures (379 UK Biobank phenotypes) | IEU OpenGWAS, identifiers of the form `ukb-b-*`, European ancestry, n ≈ 450,000, GRCh37 | https://gwas.mrcieu.ac.uk/ (requires an OpenGWAS token; the API was used through `ieugwasr`) |
| Cholelithiasis (K11_CHOLELITH, ICD-10 K80; 49,834 cases / 437,418 controls) | FinnGen release 12 | https://r12.finngen.fi/ |
| Biliary tract cancer (C3_BILIARY_GALLBLADDER_EXALLC; 2,298 cases / 378,749 controls; controls exclude all cancers) | FinnGen release 12 | https://r12.finngen.fi/ |
| Covariates: BMI (`BMI_IRN`, inverse-rank-normalised), cholecystitis (K11_CHOLECYST), cholecystectomy (K11_CHOLECYSTECTOMY), other gallbladder disease (K11_GALLBLAOTH) | FinnGen release 12 (same freeze as the outcome) | https://r12.finngen.fi/ |
| Cross-ancestry replication and control outcomes (BBJ.BtC 418 cases; cholelithiasis 9,305; cholecystitis 1,860; gallbladder polyp 843; urolithiasis 11,699; pancreatic cancer 499) | BioBank Japan, Japanese ancestry, GRCh37 | https://pheweb.jp/downloads (dataset `hum0197`) |
| Liver cis-eQTLs for the ABCG5/ABCG8 region | GTEx v8 liver, eQTL Catalogue re-processing | https://www.ebiqtl.org/ |
| LD reference panel | 1000 Genomes phase 3 European samples; FinnGen R12 LD panel | see `15_fetch_ld_regions.py` |
| External BMI instruments (pipeline validation) | GIANT consortium | as referenced in the manuscript |

Cross-build matching is performed by rsID (token-level exact match, not substring) because FinnGen R12 is GRCh38 while
UK Biobank and BioBank Japan are GRCh37; coordinates are verified against the Ensembl REST service.

---

## 3. Software and versions

- **R 4.6.1** with `TwoSampleMR` 0.7.9, `ieugwasr` 1.1.0.9000, `MRPRESSO` 1.0, `MVMR` 0.4.8,
  `MendelianRandomization` 0.10.0, `coloc` 5.2.3, `susieR` 0.14.2, `data.table`
- **Python 3** with `pandas`, `numpy`, `scipy`, `matplotlib`, `requests`, `pysam`/`tabix`, `pymupdf`, `python-docx`
- **PLINK** 1.9 / 2.x for LD clumping and data management

Install the R packages with `install.packages()` or `remotes::install_github("MRCIEU/TwoSampleMR")`; the OpenGWAS token
is read from the environment by `ieugwasr` (`OPENGWAS_JWT`).

---

## 4. Script index

Scripts are kept flat and retain their original numeric prefixes, which record execution order.

### Project root: the `MR_BASE` environment variable

No absolute paths are hard-coded anywhere in this repository. Every script that reads or writes project data resolves
the project root from the environment variable `MR_BASE`:

```r
BASE <- Sys.getenv("MR_BASE", unset = ".")   # R
```
```python
BASE = os.environ.get("MR_BASE", ".")        # Python
```
```sh
cd "${MR_BASE:?set MR_BASE}"                 # shell
```

Before running the pipeline, point `MR_BASE` at a directory containing `data/` and `results/` — that is, at the
directory in which the GWAS summary statistics were staged:

```sh
export MR_BASE=/path/to/project      # e.g. ~/WorkBuddy/mr-biliary
python3 scripts/23_finngen_leads.py  # run from the repository root
```

If `MR_BASE` is unset the scripts fall back to the current working directory, so running them from the repository root
also works provided `data/` and `results/` are staged there. Raw and intermediate data are deliberately excluded from
this repository by `.gitignore`.

> **Note on credentials.** No API tokens, keys or passwords are stored in this repository. The OpenGWAS token is read
> at run time by `ieugwasr` from the `OPENGWAS_JWT` environment variable (see `01_extract_instruments.R`).

### A. Data acquisition
| Script | Purpose |
|---|---|
| `download_finngen_parallel.py` | Parallel download of FinnGen R12 GWAS summary statistics |
| `14_tabix_fetch.py` | Region-wise retrieval from tabix-indexed GWAS files |
| `15_fetch_ld_regions.py` | Retrieval of LD reference regions |
| `18_bbj_fetch.py` | Retrieval of BioBank Japan summary statistics |

### B. Exposure selection and instruments
| Script | Purpose |
|---|---|
| `00_pick_diet_exposures.py` | Pre-specified keyword whitelist over UK Biobank diet phenotypes, followed by cross-batch de-duplication (yields the 379 screened phenotypes) |
| `01_extract_instruments.R` | Instrument extraction and LD clumping from OpenGWAS for each dietary trait |
| `05_extract_iv_from_gwas.py` | Direct instrument extraction from GWAS files |
| `21_bbj_extract_instruments.py` | BBJ-side instrument effects for replication |
| `22_verify_leads.py` | Lead-variant verification: Ensembl GRCh37 coordinates and token-level exact rsID matching |
| `23_finngen_leads.py` | FinnGen R12 lead-variant coordinates and effects (hg38) with Ensembl cross-check |

### C. Univariable MR
| Script | Purpose |
|---|---|
| `02_extract_outcome_snps.py` | Outcome-side SNP extraction and matching |
| `03_mr_core.R` | Main two-sample MR plus pleiotropy-robust sensitivity analyses |
| `04_mr_presso.R` | MR-PRESSO outlier detection and outlier-corrected IVW |
| `43_fix_kappa.R` | **Correction of a random-effects double-inflation defect in `03_mr_core.R`.** Produces the authoritative `*_v2.csv` results used in the manuscript (single κ = max(1, √(Q/df))). See `docs/known_issues.md` |
| `_check_84iv.py` | Independent recomputation of IVW / MR-Egger / weighted median from the authoritative 84-variant set |
| `_recheck_egger.R` | Cross-check of the 84-variant MR-Egger fit with `TwoSampleMR` |
| `05_make_report.py` | Assembly of all MR results (including MR-PRESSO and mediation) into a report |

### D. Mediation
| Script | Purpose |
|---|---|
| `06_medi_mr.R` | Two-step mediation, diet → cholelithiasis → BTC (first implementation) |
| `37_mediation.R` | Two-step mediation, final specification |
| `41_mvmr_mediation.R` | MVMR direct effect with the dietary exposure and cholelithiasis entered simultaneously |

### E. Multivariable MR
| Script | Purpose |
|---|---|
| `07_mvmr.R` | MVMR of cholelithiasis on BTC conditional on BMI |
| `08_mvmr_finngen.R` | MVMR using the same-source FinnGen freeze |
| `30_mvmr_extract.py` | Extraction of MVMR input effects |
| `31_clump_mvmr_ivs.R` | Joint clumping across exposures and construction of the union instrument set |
| `32_watch_extract.sh` | Extraction of union-set effects from the five FinnGen GWASs |
| `33_mvmr_run.R` | MVMR with cholecystitis / gallbladder-disease adjustment |
| `34_diagnose_identification.R` | Identification diagnostics: genetic collinearity between cholecystitis and cholelithiasis (conditional F-statistics) |

### F. Colocalisation and region definition
| Script | Purpose |
|---|---|
| `10_define_coloc_regions.py` | Region definition around independent lead variants |
| `11b_extract_coloc_regions_stream.py` | Streaming (memory-safe) extraction of the coloc input regions |
| `12_coloc.R` | `coloc.abf` across all independent cholelithiasis loci |
| `12b_coloc_window_sensitivity.R` | Window-width sensitivity ( ±100, ±250, ±500 kb) |
| `16_coloc_susie.R` | `coloc.susie` multi-signal colocalisation |
| `16b_susie_param_scan.R` | Parameter scan because credible sets require signals on both traits |
| `17_count_independent_signals.py` | SuSiE independent-signal counts (region-definition diagnostic 1) |
| `24_split_regions.py` | Splitting of over-merged regions at the largest inter-lead gap |
| `25_coloc_r07split.R` | Re-run of the over-merged chromosome 2 region after splitting (the 0.31 → 0.19/0.12 result) |
| `26_locus_replication.py` | Locus-level replication checks |
| `58_qtl_layer_check.py`, `59_coloc_eqtl_gwas.R` | Liver cis-eQTL × GWAS colocalisation for ABCG5/ABCG8 |

### G. Cross-ancestry replication
| Script | Purpose |
|---|---|
| `19_bbj_lookup.py` | BBJ effect lookup for the 84 Finnish-derived instruments |
| `20_bbj_mr.R` | Two-sample MR in BioBank Japan with pre-specified positive and negative control outcomes |
| `27_direction_concordance_v2.py` | Error-in-variables-corrected weighted regression of BBJ on FinnGen effects |
| `28_sign_concordance.py` | Leave-one-region-out direction concordance |

### H. Directionality, power, outcome panel
| Script | Purpose |
|---|---|
| `38_steiger_lor.R` | Steiger testing with liability-scale variance explained (`get_r_from_lor()`), across 16 prevalence combinations |
| `39_power_finngen.R` | Power from observed standard errors (FinnGen side) |
| `40_finngen_outcome_panel.R` | Same 84 instruments against the FinnGen outcome panel |

### I. Figures
`53_make_figure1.py`, `54_make_figure2.py`, `56_make_figure3.py`, `55_make_figure4.py`, `61_make_figure5.py`,
`13_make_coloc_figure.py`, `16c_make_susie_figure.py`, `29_make_bbj_figure.py`, `35_make_mvmr_figure.py`,
`60_make_figureS_qtl.py`, `62_make_figureS3.py`, and `57_audit_figures.py` (pre-submission figure checks).
`09_make_figure.py` is **deprecated** and must not be used.

### J. Manuscript production tooling (not part of the analysis)
`36_journal_probes.py`, `48_journal_probes_v2.py`, `49_parse_author_guide.py`, `42_build_references.py`,
`44_finalize_references.py`, `45_number_citations.py`, `47_audit_references.py`, `50_word_count.py`,
`51_render_reference_list.py`, `52_render_docx.py`, `63_build_supplementary.py`, `64_make_submission_letters.py`,
`65_assemble_submission.py`, `66_fix_checklist.py`, `67_final_check_submission.py`.

---

## 5. Reproducing the main numbers

```bash
# 1. Select dietary exposures and extract instruments
python  scripts/00_pick_diet_exposures.py
Rscript scripts/01_extract_instruments.R

# 2. Univariable MR, then the IVW standard-error correction (authoritative results)
Rscript scripts/03_mr_core.R
Rscript scripts/04_mr_presso.R
Rscript scripts/43_fix_kappa.R          # writes mr_results_*_v2.csv  <- values reported in the paper

# 3. Multivariable MR
Rscript scripts/31_clump_mvmr_ivs.R
bash    scripts/32_watch_extract.sh
Rscript scripts/33_mvmr_run.R
Rscript scripts/34_diagnose_identification.R

# 4. Colocalisation and the region-definition diagnostic
Rscript scripts/12_coloc.R
python  scripts/17_count_independent_signals.py
python  scripts/24_split_regions.py
Rscript scripts/25_coloc_r07split.R

# 5. Mediation
Rscript scripts/37_mediation.R
Rscript scripts/41_mvmr_mediation.R

# 6. Cross-ancestry replication
Rscript scripts/20_bbj_mr.R
python  scripts/27_direction_concordance_v2.py
```

## 6. Key result files

The CSV/TSV files under `key_results/` are the small derived outputs that back the manuscript's tables:
univariable MR results (`mr_results_*_v2.csv`), the MVMR union inputs and conditional F-statistics (`mvmr_*`), the
BioBank Japan replication results, and the power tables. They are included so that the reported numbers can be checked
without re-downloading the 6.8 GB of raw summary statistics.

## 7. Known issues

See [`docs/known_issues.md`](docs/known_issues.md). In brief: `03_mr_core.R` inflated the IVW standard error by κ
twice; `43_fix_kappa.R` corrects this by recomputing κ once from SNP-level data, and all estimates reported in the
manuscript are taken from the corrected `*_v2.csv` files. `09_make_figure.py` is superseded and must not be used.

## 8. Licence

MIT — see [`LICENSE`](LICENSE).
