#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""63_build_supplementary.py —— 生成 AJCN 投稿用补充材料（独立文件）

依 ASN Guide for Authors（2026-03-24 版）"Supplementary Data" 节：
  * 补充材料须与主稿分开成文件
  * 文件名与标题须自解释（例："Supplementary Figure 1"）
  * **所有补充材料文件须加页眉：稿件标题 + 第一作者姓名**
  * 每个条目须在主稿正文中呼唤（本稿已用 "Supplementary Table S1" 等）
  * 单文件 ≤ 2 MB；排版与已发表稿件一致

产出（results/submission/supplementary/）
  Supplementary_Table_S1.docx    MVMR 四个校正集（条件 F 统计量 + 可识别性标注）
  Supplementary_Table_S2.docx    全部 379 个膳食表型（含 226 个无可用工具变量者）
  Supplementary_Table_S3.docx    全部共定位区（窗口宽度 + p12 敏感性）
  Supplementary_Table_S4.docx    效能表（观测 SE → 最小可检出 OR → 参照 OR 的效能）
  Supplementary_Table_S5.docx    84 个胆石症工具变量在六个 FinnGen 结局上的面板
  Supplementary_Methods.docx     可复现性说明（软件版本 / 文献核验 / IVW 随机效应约定）
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(os.environ.get("MR_BASE", "."))
RES = ROOT / "results"
OUT = RES / "submission" / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)

FONT = "Times New Roman"
TITLE = ("Dietary intake, cholelithiasis and biliary tract cancer: a two-sample "
         "Mendelian randomisation study with cross-ancestry replication and "
         "region-definition-aware colocalisation")
FIRST_AUTHOR = "Jun-Lin Li"

TBL_PT = 7.5
FOOT_PT = 7.5


# ------------------------------------------------------------------ 工具
def _set_font(run, size, bold=False, italic=False):
    run.font.name = FONT
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    rpr = run._element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts")
        rpr.append(rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rf.set(qn(a), FONT)


def new_doc(landscape=True):
    d = Document()
    st = d.styles["Normal"]
    st.font.name = FONT
    st.font.size = Pt(11)
    sec = d.sections[0]
    if landscape:
        sec.orientation = WD_ORIENT.LANDSCAPE
        w, h = sec.page_width, sec.page_height
        sec.page_width, sec.page_height = h, w
    sec.left_margin = sec.right_margin = Inches(0.7)
    sec.top_margin = sec.bottom_margin = Inches(0.7)
    return d


def add_header(doc):
    """ASN 强制：所有补充材料文件加「稿件标题 + 第一作者姓名」页眉。"""
    hdr = doc.sections[0].header
    p = hdr.paragraphs[0] if hdr.paragraphs else hdr.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _set_font(p.add_run(TITLE), 8, italic=True)
    p.add_run().add_break()
    _set_font(p.add_run(f"{FIRST_AUTHOR} et al. — Supplementary material"),
              8, italic=True)
    return p


def add_footer_pagenum(doc):
    sec = doc.sections[0]
    p = sec.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    _set_font(run, 9)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    run._element.addnext(fld)


def para(doc, text, size=11, bold=False, italic=False, align=None, space_after=6):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    p.paragraph_format.space_before = Pt(0)
    _set_font(p.add_run(text), size, bold, italic)
    return p


def set_cell(cell, text, size=TBL_PT, bold=False, align="left"):
    cell.text = ""
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(1)
    p.paragraph_format.space_before = Pt(1)
    p.alignment = {"left": WD_ALIGN_PARAGRAPH.LEFT,
                   "right": WD_ALIGN_PARAGRAPH.RIGHT,
                   "center": WD_ALIGN_PARAGRAPH.CENTER}[align]
    _set_font(p.add_run(text), size, bold)


def repeat_header(row):
    tr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    tr.append(el)


def make_table(doc, headers, rows, aligns=None, size=TBL_PT, widths=None):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    aligns = aligns or ["left"] + ["center"] * (len(headers) - 1)
    for i, h in enumerate(headers):
        set_cell(t.rows[0].cells[i], h, size=size, bold=True, align="center")
    repeat_header(t.rows[0])
    for r in rows:
        cells = t.add_row().cells
        for i, v in enumerate(r):
            set_cell(cells[i], str(v), size=size, align=aligns[i])
    if widths:
        for row in t.rows:
            for i, w in enumerate(widths):
                row.cells[i].width = Inches(w)
    return t


# ------------------------------------------------------------------ 数字格式化
def fnum(v, nd=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    return f"{v:.{nd}f}"


def pfmt(p, nd=2):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "—"
    if p < 1e-4:
        m, e = f"{p:.1e}".split("e")
        return f"{m} × 10^{int(e)}"
    return f"{p:.{nd}f}"


def intc(v):
    """整数化；NaN / inf → '—'。"""
    if v is None or not np.isfinite(v):
        return "—"
    return f"{int(np.ceil(v)):,}"


def orci(o, l, u, nd=3):
    return f"{o:.{nd}f} ({l:.{nd}f}–{u:.{nd}f})"


# ================================================================== Table S1
def table_s1():
    d = new_doc()
    add_header(d)
    add_footer_pagenum(d)
    para(d, "Supplementary Table S1. Multivariable Mendelian randomisation under "
            "all adjustment sets, with conditional F-statistics", 11, bold=True,
         space_after=4)
    para(d, "Conditional F-statistics below 10 indicate that the model is not "
            "identified; estimates from those models are reported for completeness "
            "only and are explicitly excluded from interpretation (main text "
            "Sections 3.5 and 4.4). All estimates use the same single-inflation "
            "random-effects IVW convention as the main analysis (main text "
            "Section 2.4).", FOOT_PT, italic=True, space_after=8)

    SETS = [
        ("chole_bmi", "Set 1: cholelithiasis + BMI (union instruments)"),
        ("chole_cholecyst", "Set 2: cholelithiasis + cholecystitis"),
        ("chole_cholecyst_ec", "Set 3: cholelithiasis + cholecystitis + "
                               "cholecystectomy"),
        ("with_bmi", "Set 4: cholelithiasis + cholecystitis + cholecystectomy "
                     "+ BMI"),
    ]
    EXPO = {"cholelith": "Cholelithiasis", "bmi": "BMI",
            "cholecyst": "Cholecystitis", "cholecstec": "Cholecystectomy"}

    for key, label in SETS:
        cond = pd.read_csv(RES / "mvmr" / f"mvmr_condF_{key}.csv").iloc[0].tolist()
        q = pd.read_csv(RES / "mvmr" / f"mvmr_Q_{key}.csv").iloc[0]
        rob = pd.read_csv(RES / "mvmr" / f"mvmr_robust_{key}.csv")
        rob["method"] = rob["method"].str.replace(
            r"^MVMR-Lasso.*$", "MVMR-Lasso", regex=True)
        exps = sorted(rob.exposure.unique())
        if len(cond) != len(exps):
            raise AssertionError(f"{key}: 条件 F 数 {len(cond)} ≠ 暴露数 {len(exps)}")

        para(d, label, 10, bold=True, space_after=2)
        para(d, "Conditional F-statistic — " + "; ".join(
            f"{EXPO.get(e, e)} = {c:.2f}" for e, c in zip(exps, cond)) +
            f";  Cochran Q = {q['Qstat']:.1f} (P = {pfmt(q['Qpval'])})",
             FOOT_PT, space_after=4)

        idflag = all(c >= 10 for c in cond)
        heads = ["Method", "Exposure", "OR", "95% CI", "P"]
        rows = []
        for _, r in rob.iterrows():
            rows.append([r.method, EXPO.get(r.exposure, r.exposure),
                         fnum(r["or"]), f"{r.lci:.3f}–{r.uci:.3f}", pfmt(r.pval)])
        make_table(d, heads, rows, widths=[1.5, 1.6, 0.9, 1.5, 1.1])
        para(d, ("Model identified (all conditional F ≥ 10)." if idflag else
                 "⚠ Model NOT identified: at least one conditional F < 10. "
                 "Point estimates and confidence intervals are not interpretable "
                 "as conditional effects."), FOOT_PT,
             italic=not idflag, space_after=12)

    para(d, "Abbreviations: BMI, body mass index; CI, confidence interval; "
            "IVW, inverse-variance weighted; MVMR, multivariable Mendelian "
            "randomisation; OR, odds ratio.", FOOT_PT, space_after=0)
    p = OUT / "Supplementary_Table_S1.docx"
    d.save(p)
    return p


# ================================================================== Table S2
def table_s2():
    d = new_doc()
    add_header(d)
    add_footer_pagenum(d)
    para(d, "Supplementary Table S2. All 379 UK Biobank diet-related phenotypes "
            "screened, including those with no usable instruments", 11, bold=True,
         space_after=4)
    para(d, "Phenotypes are listed in the screened order. Instrument counts are "
            "the number of independent genome-wide significant variants "
            "available for the phenotype. Cholelithiasis estimates are against "
            "FinnGen R12 cholelithiasis (49,834 cases); BTC estimates are against "
            "FinnGen R12 biliary tract cancer (2,298 cases). MR-PRESSO "
            "outlier-corrected estimates are the pre-specified primary result for "
            "the dietary screen and are undefined where fewer than four variants "
            "were available. FDR is the Benjamini–Hochberg adjusted P value across "
            "the analysable dietary exposures. The 83 phenotypes not discussed in "
            "the main text are reported here without interpretation, as stated in "
            "Section 3.3.", FOOT_PT, italic=True, space_after=8)

    traits = pd.read_csv(RES / "diet_traits.csv").fillna({"category": "—"})
    log = pd.read_csv(RES / "full" / "exposure_log.csv")
    chole = pd.read_csv(RES / "full" / "mr_results_Cholelithiasis_v2.csv")
    btc = pd.read_csv(RES / "full" / "mr_results_BTC_v2.csv")
    pres = pd.read_csv(RES / "full" / "mr_presso_results.csv")

    chole_ivw = (chole[chole.method == "Inverse variance weighted"]
                 .set_index("id.exposure"))
    btc_ivw = (btc[btc.method == "Inverse variance weighted"]
               .set_index("id.exposure"))
    pres_idx = pres.set_index("exposure")
    log_idx = log.set_index("id")

    # ---- 闭合断言 ----
    assert len(traits) == 379, len(traits)
    assert (log.status == "no_iv").sum() == 226
    assert chole_ivw.index.nunique() == 86 and btc_ivw.index.nunique() == 86
    assert set(pres_idx.index) >= {"Alcohol intake frequency.", "Cheese intake",
                                   "Dried fruit intake"}
    assert abs(float(pres_idx.loc["Alcohol intake frequency.", "or_corr"])
               - 1.307994) < 1e-5

    heads = ["#", "Phenotype ID", "Dietary phenotype", "Category", "N",
             "Instruments (n)", "Chain-analysable",
             "Cholelithiasis OR (95% CI), MR-PRESSO",
             "Cholelithiasis P, MR-PRESSO", "FDR P",
             "Cholelithiasis OR (95% CI), IVW", "BTC OR (95% CI), IVW", "BTC P"]
    rows = []
    for i, t in enumerate(traits.itertuples(), 1):
        pid, name = t.id, t.trait
        n_iv = int(log_idx.loc[pid, "nsnp"]) if pid in log_idx.index else 0
        ana = "Yes" if pid in chole_ivw.index else "No"
        po = pp = pf = "—"
        if name in pres_idx.index:
            r = pres_idx.loc[name]
            po = orci(r.or_corr, r.or_lci_corr, r.or_uci_corr)
            pp = pfmt(r.p_corr)
            pf = pfmt(r.fdr_corr) if np.isfinite(r.fdr_corr) else "—"
        co = "—"
        if pid in chole_ivw.index:
            r = chole_ivw.loc[pid]
            co = orci(r["or"], r["or_lci95"], r["or_uci95"])
        bo = bp = "—"
        if pid in btc_ivw.index:
            r = btc_ivw.loc[pid]
            bo = orci(r["or"], r["or_lci95"], r["or_uci95"])
            bp = pfmt(r.pval)
        rows.append([i, pid, name, t.category, f"{int(t.sample_size):,}"
                     if np.isfinite(t.sample_size) else "—",
                     n_iv, ana, po, pp, pf, co, bo, bp])

    make_table(d, heads, rows,
               aligns=["left", "left", "left", "left", "right", "center",
                       "center", "center", "center", "center", "center",
                       "center", "center"],
               widths=[0.35, 0.85, 2.5, 1.2, 0.75, 0.7, 0.8, 1.5, 0.75, 0.7,
                       1.5, 1.5, 0.7])
    para(d, "Abbreviations: BTC, biliary tract cancer; CI, confidence interval; "
            "FDR, false discovery rate; IVW, inverse-variance weighted; "
            "MR-PRESSO, Mendelian randomisation pleiotropy residual sum and "
            "outlier; N, sample size; OR, odds ratio.", FOOT_PT, space_after=0)
    p = OUT / "Supplementary_Table_S2.docx"
    d.save(p)
    return p


# ================================================================== Table S3
def table_s3():
    d = new_doc()
    add_header(d)
    add_footer_pagenum(d)
    para(d, "Supplementary Table S3. Colocalisation results for all regions, with "
            "window-width and prior-probability sensitivity analyses", 11,
         bold=True, space_after=4)
    para(d, "Posterior probability of hypothesis 4 (PP.H4; a single causal variant "
            "shared between the two traits) is reported under the main analysis "
            "(region-defined window, p12 = 1 × 10⁻⁵) and under two sensitivity "
            "families: the prior probability of a shared causal variant "
            "(p12 = 1 × 10⁻⁶, 1 × 10⁻⁵, 1 × 10⁻⁴) and the window half-width "
            "(100, 250, 500 kb around the lead variant). PP.H3 is the posterior "
            "probability of two distinct causal variants and is given under the "
            "main analysis settings. Regions were defined by merging lead variants "
            "within 1 Mb; because the merging rule itself changes the answer for "
            "some regions (main text Section 4.3), the split definition for the "
            "chromosome 2 locus is given in the accompanying text and in "
            "Supplementary Figure S2.", FOOT_PT, italic=True, space_after=8)

    allsum = pd.read_csv(RES / "coloc" / "coloc_all_summary.csv")
    win = pd.read_csv(RES / "coloc" / "coloc_window_sensitivity.csv")

    assert set(allsum.p12.unique()) == {1e-6, 1e-5, 1e-4}, allsum.p12.unique()
    assert set(win.halfwin_kb.unique()) == {100, 250, 500}, win.halfwin_kb.unique()

    main = allsum[allsum.p12 == 1e-5].copy()
    key = ["region_id", "pair"]
    p12map = {p: allsum[allsum.p12 == p].set_index(key)["PP.H4"] for p in (1e-6, 1e-4)}
    winmap = {w: win[win.halfwin_kb == w].set_index(key)["PP.H4"]
              for w in (100, 250)}
    assert len(main) == len(p12map[1e-6]) == len(winmap[100])

    heads = ["Region", "Trait pair", "Chr", "Lead variant(s)", "Nearest gene(s)",
             "n SNPs", "PP.H3", "PP.H4", "p12 = 1e-6", "p12 = 1e-4",
             "100 kb", "250 kb"]
    rows = []
    for _, r in main.sort_values(["pair", "region_id"]).iterrows():
        k = (r["region_id"], r["pair"])
        rows.append([
            r["region_id"], str(r["pair"]).replace("_", " vs "), str(r["chrom"]),
            r["lead_snps"], str(r["nearest_genes"]).replace(",", ", "),
            int(r["nsnps"]), fnum(r["PP.H3"]), fnum(r["PP.H4"]),
            fnum(p12map[1e-6].get(k, np.nan)), fnum(p12map[1e-4].get(k, np.nan)),
            fnum(winmap[100].get(k, np.nan)), fnum(winmap[250].get(k, np.nan))])

    make_table(d, heads, rows,
               aligns=["left", "left", "center", "left", "left", "right"] +
                      ["center"] * 6,
               widths=[0.6, 1.35, 0.4, 1.9, 1.5, 0.55, 0.6, 0.6, 0.75, 0.75,
                       0.65, 0.65])
    para(d, f"Rows: {len(rows)} region × trait-pair combinations. "
            "Abbreviations: Chr, chromosome; PP.H3, posterior probability of "
            "hypothesis 3 (two distinct causal variants); PP.H4, posterior "
            "probability of hypothesis 4 (one shared causal variant); "
            "p12, prior probability of a shared causal variant.",
         FOOT_PT, space_after=0)
    p = OUT / "Supplementary_Table_S3.docx"
    d.save(p)
    return p


# ================================================================== Table S4
def table_s4():
    d = new_doc()
    add_header(d)
    add_footer_pagenum(d)
    para(d, "Supplementary Table S4. Statistical power for each analysis, and the "
            "case count needed for 80% power", 11, bold=True, space_after=4)
    para(d, "The smallest detectable odds ratio at 80% power is "
            "exp(2.8016 × SE), computed from the observed standard error of each "
            "analysis rather than assumed. Power at reference odds ratios uses the "
            "same observed standard error. The bottom block gives the power curve "
            "obtained by scaling the effective sample size to the case counts "
            "shown, holding the instrument set and effect size fixed.",
         FOOT_PT, italic=True, space_after=8)

    pw = pd.read_csv(RES / "power" / "power_table_finngen.csv")
    cv = pd.read_csv(RES / "power" / "power_curve_by_cases.csv")
    assert abs(float(pw[pw.analysis.str.contains("胆囊结石病 → 胆道癌")].se.iloc[0])
               - 0.073264) < 1e-6
    assert abs(float(pw[pw.analysis.str.contains("BBJ")].se.iloc[0])
               - 0.156578) < 1e-6
    assert len(cv) == 8

    para(d, "Block A. Power of each reported analysis", 10, bold=True, space_after=2)
    heads = ["Analysis", "Exposure set", "Outcome", "Cases", "Controls", "SE",
             "OR (95% CI)", "P", "Smallest detectable OR", "Power at 1.20",
             "Power at 1.22", "Power at 1.30", "Power at 1.50",
             "Cases for 80% power"]
    rows = []
    for r in pw.itertuples():
        rows.append([r.analysis, r.exposure_set, r.outcome, f"{int(r.ncase):,}",
                     f"{int(r.nctrl):,}", f"{r.se:.4f}",
                     f"{r.or_obs:.3f} ({r.or_lci:.3f}–{r.or_uci:.3f})",
                     pfmt(r.p_obs), f"{r.or_min80:.3f}", f"{r.power_or120:.2f}",
                     f"{r.power_or122:.2f}", f"{r.power_or130:.2f}",
                     f"{r.power_or150:.2f}", intc(r.cases_need_or120)])
    make_table(d, heads, rows,
               aligns=["left", "left", "left", "right", "right", "center",
                       "center", "center", "center"] + ["center"] * 5,
               widths=[1.6, 1.35, 1.15, 0.65, 0.7, 0.55, 1.35, 0.65, 0.8,
                       0.65, 0.65, 0.65, 0.65, 0.8])
    para(d, "", 8, space_after=6)

    para(d, "Block B. Power curve by case count (biliary tract cancer)",
         10, bold=True, space_after=2)
    heads2 = ["Cases", "Effective sample size", "SE", "Smallest detectable OR",
              "Power at 1.20", "Power at 1.25", "Power at 1.30"]
    rows2 = []
    for r in cv.itertuples():
        rows2.append([f"{int(r.ncase):,}", f"{r.n_eff:,.0f}", f"{r.se:.4f}",
                      f"{r.or_min80:.3f}", f"{r.power_or120:.3f}",
                      f"{r.power_or125:.3f}", f"{r.power_or130:.3f}"])
    make_table(d, heads2, rows2, widths=[0.95, 1.5, 0.8, 1.3, 1.0, 1.0, 1.0])
    para(d, "Abbreviations: CI, confidence interval; OR, odds ratio; "
            "SE, standard error. Effective sample size is "
            "cases × controls / (cases + controls).", FOOT_PT, space_after=0)
    p = OUT / "Supplementary_Table_S4.docx"
    d.save(p)
    return p


# ================================================================== Table S5
def table_s5():
    d = new_doc()
    add_header(d)
    add_footer_pagenum(d)
    para(d, "Supplementary Table S5. Outcome panel: the 84 cholelithiasis "
            "instruments across six FinnGen endpoints", 11, bold=True,
         space_after=4)
    para(d, "Each endpoint was analysed with the same 84 instruments, the same "
            "harmonisation and the same single-inflation random-effects IVW "
            "convention. The self-regression row uses cholelithiasis as its own "
            "outcome and is a technical positive control for the pipeline: because "
            "exposure and outcome are the same GWAS, an odds ratio of exactly "
            "e^ln(2.718) and a Cochran Q of zero are expected. Cholecystectomy, "
            "cholecystitis and gallbladder polyp are positive controls; "
            "urolithiasis and pancreatic cancer are negative controls.",
         FOOT_PT, italic=True, space_after=8)

    op = pd.read_csv(RES / "power" / "finngen_outcome_panel.csv")
    assert len(op) == 6, len(op)
    self_chk = op[op.type == "selfchk"]
    assert abs(float(self_chk.or_ivw.iloc[0]) - 2.718282) < 1e-5
    assert float(self_chk.Q.iloc[0]) < 1e-6

    heads = ["Endpoint", "Role", "Cases", "Controls", "n IV",
             "OR (95% CI), IVW", "P, IVW", "Infl.", "Q (df)", "P, Q",
             "OR, MR-Egger", "P, MR-Egger", "Egger intercept (P)",
             "OR (95% CI), weighted median", "P, weighted median",
             "Smallest detectable OR"]
    rows = []
    for r in op.itertuples():
        rows.append([r.label,
                     "Technical control" if r.type == "selfchk" else
                     ("Positive control" if "阳性" in r.label else
                      ("Negative control" if "阴性" in r.label else "Main outcome")),
                     f"{int(r.ncase):,}", f"{int(r.nctrl):,}", int(r.nsnp),
                     f"{r.or_ivw:.3f} ({r.lci:.3f}–{r.uci:.3f})", pfmt(r.p_ivw),
                     f"{r.kappa:.2f}", f"{r.Q:.1f} ({int(r.Q_df)})", pfmt(r.Q_p),
                     f"{r.or_egger:.3f}", pfmt(r.p_egger),
                     f"{r.egger_int:+.4f} ({pfmt(r.egger_int_p)})",
                     f"{r.or_wm:.3f}", pfmt(r.p_wm), f"{r.or_min80:.3f}"])
    make_table(d, heads, rows,
               aligns=["left", "left"] + ["right"] * 2 + ["center"] * 12,
               widths=[1.5, 0.95, 0.55, 0.65, 0.4, 1.35, 0.7, 0.5, 0.8, 0.6,
                       0.7, 0.7, 1.1, 1.2, 0.8, 0.85])
    para(d, "Abbreviations: CI, confidence interval; IV, instrumental variable; "
            "IVW, inverse-variance weighted; OR, odds ratio; SE, standard error. "
            "\"Infl.\" is the random-effects inflation factor max(1, √(Q/df)).",
         FOOT_PT, space_after=0)
    p = OUT / "Supplementary_Table_S5.docx"
    d.save(p)
    return p


# ================================================================== Methods
def supp_methods():
    d = new_doc(landscape=False)
    add_header(d)
    add_footer_pagenum(d)
    para(d, "Supplementary Methods. Reproducibility notes", 12, bold=True,
         space_after=8)

    BLOCKS = [
        ("1. Software and versions",
         "All analyses were performed in R 4.6.1. The following packages were used: "
         "TwoSampleMR 0.7.9, ieugwasr 1.1.0.9000, MRPRESSO 1.0, MVMR 0.4.8, "
         "MendelianRandomization 0.10.0, coloc 5.2.3 and susieR 0.14.2. "
         "Figure assembly and supplementary tables were produced in Python 3.13 "
         "with matplotlib 3.11, pandas 2.3 and python-docx 1.2. "
         "Summary statistics were retrieved from FinnGen release 12, the IEU "
         "OpenGWAS database, BioBank Japan (PheWeb, dataset hum0197) and the "
         "eQTL Catalogue GTEx v8 release."),
        ("2. Inverse-variance weighted estimates and the random-effects convention",
         "Every inverse-variance weighted estimate in this paper — univariable, "
         "multivariable, mediation and cross-ancestry — uses one convention. "
         "For each variant k the causal estimate is the Wald ratio "
         "r_k = β_out,k / β_exp,k, and the weights are the delta-method weights "
         "w_k = β_exp,k² / se_out,k², so that SE_FE = 1/√Σw_k. Because "
         "between-variant heterogeneity is substantial for dietary exposures, the "
         "standard error is inflated by a single factor κ = max(1, √(Q/df)) where "
         "Q = Σw_k(r_k − b)², giving SE = SE_FE × κ. The inflation is never "
         "permitted to shrink the standard error below the fixed-effect value. "
         "For the reference analysis (cholelithiasis → biliary tract cancer, 84 "
         "instruments) this yields b = 0.2023964295 and SE = 0.0732635699, "
         "Q = 232.394 on 83 df and κ = 1.6733, reproducing the default of "
         "TwoSampleMR::mr_ivw() to ten decimal places. Multivariable standard "
         "errors were computed from the inverse of the weighted "
         "variance–covariance matrix of the genetic associations rather than by "
         "bootstrapping, which is why the leave-one-out panel in Supplementary "
         "Figure S3 and the multivariable estimates are directly comparable."),
        ("3. Reference verification procedure",
         "Every one of the 50 references was verified field by field against the "
         "PubMed record for the cited article (authors, title, journal, year, "
         "volume, pages and DOI) before the manuscript was assembled. No reference "
         "was generated, paraphrased or reconstructed from memory. Where a record "
         "in PubMed carried no author list, the citation was replaced with a "
         "verifiable source rather than completed by inference. Corrections, "
         "errata and retractions are indexed separately from the articles they "
         "concern and were excluded from the citation set."),
        ("4. Instrument selection and allele harmonisation",
         "Instruments were selected at P < 5 × 10⁻⁸ and clumped with an r² "
         "threshold of 0.001 within a 10,000 kb window, using the 1000 Genomes "
         "European reference panel; palindromic variants with intermediate allele "
         "frequencies were removed. Outcome-side alleles were harmonised to the "
         "exposure effect allele, with strand flips resolved through complementary "
         "allele matching. For the cross-ancestry replication, FinnGen variants "
         "were matched to BioBank Japan by rsID using token-level exact matching "
         "rather than substring matching, because substring matching silently "
         "maps rsIDs onto neighbouring variants and inflates agreement."),
        ("5. Steiger directionality with binary traits",
         "Because both cholelithiasis and biliary tract cancer are binary, the "
         "variance explained on the liability scale was computed with "
         "get_r_from_lor() using case and control counts from the FinnGen R12 "
         "phenotype manifest. Prevalence was varied over 16 combinations "
         "(exposure 0.05–0.20 × outcome 0.0005–0.005). The uncorrected "
         "get_r_from_pn() computation, which assumes a continuous trait, is "
         "reported alongside as a negative control: it returns a variance-explained "
         "ratio below 1 and therefore reverses the direction of the test, which is "
         "why it was not used."),
        ("6. Files provided",
         "Supplementary Table S1 (multivariable MR, all adjustment sets), "
         "Supplementary Table S2 (all 379 dietary phenotypes), Supplementary "
         "Table S3 (colocalisation, all regions, with sensitivity analyses), "
         "Supplementary Table S4 (power), Supplementary Table S5 (outcome panel), "
         "Supplementary Figure S1 (multivariable MR identifiability), "
         "Supplementary Figure S2 (regional association plots) and Supplementary "
         "Figure S3 (Steiger directionality with prevalence sensitivity, "
         "MR-PRESSO outlier diagnostics and leave-one-out analyses). Each figure "
         "and table is a separate file, as required."),
    ]
    for head, body in BLOCKS:
        para(d, head, 11, bold=True, space_after=2)
        para(d, body, 11, space_after=10)
    p = OUT / "Supplementary_Methods.docx"
    d.save(p)
    return p


# ================================================================== main
def main():
    import shutil
    produced = [table_s1(), table_s2(), table_s3(), table_s4(), table_s5(),
                supp_methods()]

    # 图 S1 / S2 / S3：均为带「稿件标题 + 第一作者」页眉的投稿版
    #   S1 由 scripts/35_make_mvmr_figure.py 生成（图内主标题已按 ASN 删除）
    #   S2 由 scripts/13_make_coloc_figure.py 生成
    #   S3 由 scripts/62_make_figureS3.py 生成
    # ⚠️ 这三个文件**由各自出图脚本直写**本目录，本脚本只校验存在性，不覆盖。
    for n in (1, 2, 3):
        for ext in ("png", "pdf"):
            f = OUT / f"Supplementary_Figure_S{n}.{ext}"
            if not f.exists():
                raise FileNotFoundError(
                    f"{f} 不存在 —— 请先跑对应的出图脚本（S1: 35, S2: 13, S3: 62）")
            produced.append(f)

    print("产出补充材料：")
    for p in produced:
        mb = p.stat().st_size / 1e6
        flag = "OK " if mb <= 2 else "!! >2MB"
        print(f"  {flag} {p.name:38s} {mb:5.2f} MB")
    return produced


if __name__ == "__main__":
    main()
