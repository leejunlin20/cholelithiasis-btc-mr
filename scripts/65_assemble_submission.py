#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""65_assemble_submission.py —— 把全部投稿件汇总到 ~/Desktop/投稿/AJCN/

产出结构（扁平，文件名自解释，便于 Editorial Manager 逐项上传）
  README_投稿上传清单.md
  AJCN_Manuscript.docx                     投稿主件（不含图片）
  AJCN_Manuscript_with_figures.docx        内部审阅用，勿投
  AJCN_Cover_Letter.docx
  AJCN_Declaration_of_Interests.docx
  AJCN_STROBE-MR_Checklist.docx
  AJCN_Figure1..5.pdf                      矢量投稿图（逐张单独文件）
  AJCN_Figure1..5.png                      400 dpi 位图备用
  AJCN_Supplementary_Table_S1..S5.docx
  AJCN_Supplementary_Figure_S1..S3.pdf
  AJCN_Supplementary_Methods.docx
"""
from __future__ import annotations
import os

import shutil
from pathlib import Path

ROOT = Path(os.environ.get("MR_BASE", "."))
RES = ROOT / "results"
MAIN = RES / "主稿"
SUB = RES / "submission"
SUPP = SUB / "supplementary"
LETTERS = SUB / "letters"
DEST = Path.home() / "Desktop" / "投稿" / "AJCN"

FIG_PREFIX = "AJCN_Figure"

PLAN = [
    (MAIN / "AJCN_Manuscript.docx", "AJCN_Manuscript.docx"),
    (MAIN / "AJCN_Manuscript_with_figures.docx", "AJCN_Manuscript_with_figures.docx"),
    (MAIN / "AJCN_STROBE-MR_Checklist.docx", "AJCN_STROBE-MR_Checklist.docx"),
    (LETTERS / "Cover_Letter_AJCN.docx", "AJCN_Cover_Letter.docx"),
    (LETTERS / "Declaration_of_Interests.docx", "AJCN_Declaration_of_Interests.docx"),
]
for i, stem in enumerate([
        "Figure1_study_design", "Figure2_diet_cholelithiasis",
        "Figure3_cholelithiasis_btc_chain", "Figure4_mvmr_and_failure_modes",
        "Figure5_cross_ancestry_replication"], 1):
    PLAN.append((SUB / f"{stem}.pdf", f"{FIG_PREFIX}{i}.pdf"))
    PLAN.append((SUB / f"{stem}.png", f"{FIG_PREFIX}{i}.png"))

for n in range(1, 6):
    PLAN.append((SUPP / f"Supplementary_Table_S{n}.docx",
                 f"AJCN_Supplementary_Table_S{n}.docx"))
for n in range(1, 4):
    PLAN.append((SUPP / f"Supplementary_Figure_S{n}.pdf",
                 f"AJCN_Supplementary_Figure_S{n}.pdf"))
PLAN.append((SUPP / "Supplementary_Methods.docx",
             "AJCN_Supplementary_Methods.docx"))

README = """# AJCN 投稿上传清单（The American Journal of Clinical Nutrition）

**稿件**： Dietary intake, cholelithiasis and biliary tract cancer: a two-sample
Mendelian randomisation study with cross-ancestry replication and
region-definition-aware colocalisation

**通讯作者**：Jun-Lin Li ｜ leejunlin@live.cn ｜ +86 746 6236441
**投稿系统**：https://www.editorialmanager.com/ajcn/
**稿约来源**：ASN Guide for Authors（官方 PDF，更新日期 **2026-03-24**）
<https://legacyfileshare.elsevier.com/promis_misc/asn-guide-for-authors.pdf>
**打包日期**：2026-09-27

---

## 一、逐项上传对照表

| # | 上传项 | 文件 | 系统内类型 | 状态 |
|---|---|---|---|---|
| 1 | Cover Letter | `AJCN_Cover_Letter.docx` | Cover Letter | ✅ |
| 2 | Manuscript | `AJCN_Manuscript.docx` | Manuscript | ✅ |
| 3 | Declaration of Interests（全体作者） | `AJCN_Declaration_of_Interests.docx` | Declaration of Interests | ⚠️ **须三位作者确认后上传** |
| 4 | Research report checklist（STROBE-MR） | `AJCN_STROBE-MR_Checklist.docx` | **Health Research Reporting Checklist** | ✅ |
| 5 | Figure 1 | `AJCN_Figure1.pdf` | Figure | ✅ |
| 6 | Figure 2 | `AJCN_Figure2.pdf` | Figure | ✅ |
| 7 | Figure 3 | `AJCN_Figure3.pdf` | Figure | ✅ |
| 8 | Figure 4 | `AJCN_Figure4.pdf` | Figure | ✅ |
| 9 | Figure 5 | `AJCN_Figure5.pdf` | Figure | ✅ |
| 10 | Supplementary Table S1 | `AJCN_Supplementary_Table_S1.docx` | Supplementary Material | ✅ |
| 11 | Supplementary Table S2 | `AJCN_Supplementary_Table_S2.docx` | Supplementary Material | ✅ |
| 12 | Supplementary Table S3 | `AJCN_Supplementary_Table_S3.docx` | Supplementary Material | ✅ |
| 13 | Supplementary Table S4 | `AJCN_Supplementary_Table_S4.docx` | Supplementary Material | ✅ |
| 14 | Supplementary Table S5 | `AJCN_Supplementary_Table_S5.docx` | Supplementary Material | ✅ |
| 15 | Supplementary Figure S1 | `AJCN_Supplementary_Figure_S1.pdf` | Supplementary Material | ✅ |
| 16 | Supplementary Figure S2 | `AJCN_Supplementary_Figure_S2.pdf` | Supplementary Material | ✅ |
| 17 | Supplementary Figure S3 | `AJCN_Supplementary_Figure_S3.pdf` | Supplementary Material | ✅ |
| 18 | Supplementary Methods | `AJCN_Supplementary_Methods.docx` | Supplementary Material | ✅ |

### 系统里几个必答项怎么选（2026-09-27 投稿实操）

| 系统提问 | 怎么填 |
|---|---|
| **Health Research Reporting Checklist** | 上传 **`AJCN_STROBE-MR_Checklist.docx`**（本稿适用的报告清单就是 STROBE-MR，20 项已逐条填写，含 8 条 Comment 与 2 项如实标注"未做"） |
| `The Acknowledgements section includes the contributions of each author` | **Yes**（稿约 §7 强制；稿件 `## Acknowledgments` 下已有 contributions 段）。Explain 框建议填指向性说明，而非只写 N/A |
| `If this paper involves a dietary supplement…` | **Not applicable**（本稿为膳食摄入/食物组的观察性 MR，无补充剂干预） |
| `Is author … an Editorial Board Member of AJCN?` | **No**（三位作者均非编委；下拉第一项是"有一位作者是编委"，**勿误选**） |
| `All systematic reviews must be pre-registered in PROSPERO … Please add the review registration URL and registration number` | **Not applicable** —— 本稿是 Original Research（两样本 MR），非系统综述。文本框填：`Not applicable - this manuscript is not a systematic review; it is an original two-sample Mendelian randomisation study.` |

> **participant flowchart 不必另做文件**：ASN 原文要求初投随交 research report checklist 与 participant flowchart；
> 本设计**无受试者**，流程图职责由 **Figure 1 面板 B**（379 → 153 → 86 → 38 筛选漏斗）承担，
> 这一说明**已写入 Cover Letter**（"there is no participant flowchart in the conventional sense; Panel B of Figure 1 is the screening flow diagram"），
> 故编辑部若问起，指向投稿信即可。Cover Letter 同时已说明本稿**未经前瞻性注册**、PROSPERO 不适用。

> **不要上传**这两个：
> - `AJCN_Manuscript_with_figures.docx` —— 供你通读的图文合并版（期刊要求图片单独上传）
> - `AJCN_Figure1.png` ～ `AJCN_Figure5.png` —— 400 dpi 位图备用件
>   （期刊接受 PDF/EPS 矢量与 TIFF/JPEG 位图，**未列 PNG**；若系统拒收 PDF 再用这些转 TIFF）

---

## 二、本稿对 ASN 稿约的逐条落点

> ASN 自 2026 年起改为**初次投稿免排版**（"hassle-free submission … Only at revision will
> authors be asked to format their manuscript in compliance with journal style"），
> 但下列项目**初次投稿即强制**。本稿全部满足。

| 稿约要求（官方原文要点） | 本稿落点 |
|---|---|
| 初投即须含**每位作者**的利益冲突声明 | `AJCN_Declaration_of_Interests.docx` |
| 初投即须含**连续行号** | Word 自 **Abstract** 起连续编号，**止于参考文献之前**；标题页不编号 |
| 初投即须含**页码** | 每页页脚居中 |
| MR 研究须随初投提交 **STROBE-MR 清单**与流程图 | `AJCN_STROBE-MR_Checklist.docx`；流程图职责由 **Figure 1 面板 B**（379→153→86 筛选漏斗）承担 |
| 章节表固定为 Title Page / Abstract / Introduction / Methods / Results / Discussion / Acknowledgments / References | 完全一致；**无独立 Conclusions**，已并入 Discussion 4.7 |
| 标题页须含 a–g 七项 | 标题、全体作者、单位、Disclaimers、Sources of support、通讯作者姓名+邮寄地址+电话+邮箱、20 条字母序缩写表 |
| 摘要 ≤300 词、结构式 Background / Objective / Methods / Results / Conclusions | **297 词**，五段式；缩写在首次出现处定义 |
| 摘要下 5–10 个关键词，须含研究对象 | **10 个**，含 `humans`、`gallstone disease` |
| 文内引用用括号阿拉伯数字、置于标点前 | 全稿 `(1,2)` 体例 |
| 参考文献 >6 作者取前 6 位 + et al. | 50 条全部符合；逐字段经 PubMed 核验 |
| 致谢须含作者贡献（官方七描述符 + 全体作者缩写 + 已读批准声明），**不得写 obtained funding** | 已按规范重写 |
| 数据共享须按官方三选一原句照录 | 采用**选项 2** 原句：`Data described in the manuscript, code book, and analytic code will be made available upon request pending application and approval.` |
| 生成式 AI 须在参考文献前设固定标题的声明章节 | 章节 `Declaration of Generative AI and AI-assisted technologies in the writing process`，用官方模板句式 |
| 表格须置于参考文献之后的正文文件中、每表另起一页、脚注用上标阿拉伯数字 | Table 1–3 位于 References 之后的横向节，逐表另页 |
| 图件须逐张单独成文件；标题与图注写在正文文件里而非图上 | 5 张图各自单独文件；图注在正文 `Figure legends` 节；图内无主标题 |
| 补充材料须单独成文件、加「稿件标题 + 第一作者姓名」页眉、单件 ≤2 MB | 9 件全部满足（最大 0.33 MB）；Figure S1–S3 的页眉经 **PDF 文本层核验**存在 |
| Electronic Artwork：`Embed the fonts used if the application provides that option`；`Make sure you use uniform lettering and sizing`；`Aim to use … Arial, Courier, Times New Roman, Symbol, or … fonts that look similar` | 8 张图**统一为单一无衬线族 DejaVu Sans**、一律 **TrueType（pdf.fonttype=42）嵌入**；无 Type 3 字体、无单图内跨族混排（`57`/`67` 已把这两项做成阻断级检查） |
| 图件缩小后字号约 8 pt、多为单栏 < 8.5 cm | ⚠️ **唯一未达标项**，属**返修阶段**要求（见第四节） |

---

## 三、投稿前您需要做的三件事

1. **确认利益冲突声明**（诚信项）
   打开 `AJCN_Declaration_of_Interests.docx`，请三位作者各自确认；如您使用 ASN 在线
   Declaration of Interests 工具生成的文件，可直接替换本件，内容等价。

2. **在系统中填写与标题页完全一致的元数据**
   作者姓名、单位、通讯作者邮箱与电话必须与 Word 标题页**逐字一致**（ASN 明确比对）。

3. **按需填写系统里的关键词（供选审稿人）**
   系统会另问一组"用于寻找审稿人"的关键词，可与稿件中的 10 个不同。

---

## 四、已知未办事项（不阻断初投）

**图件印幅等效字号**：ASN 要求图件缩小到发表尺寸后文字约 8 pt，多数图为单栏 < 8.5 cm。
本稿 5 张图为 11.4 in 宽画布，折到全宽 17.5 cm 后最小字号等效 3.4–3.8 pt。
这是**返修阶段**要求（ASN 明文：排版合规"仅在返修时"要求），不影响初投；
若编辑部在返修时提出，需按发表尺寸重制（减文字或增大面板占比，**单纯缩放画布无效**）。
体检脚本：`scripts/57_audit_figures.py`（当前阻断级 0 项，此项属卫生级）。

**图件字体族（已处理，留档备查）**：官方点名的首选字体是 Arial / Courier / Times New
Roman / Symbol。本稿 8 张图统一使用 **DejaVu Sans**（无衬线，属官方措辞
"or … fonts that look similar" 允许范围），换来的是**全库单一族 + 全部 TrueType 嵌入 +
无跨族混排**。若返修阶段编辑部点名要求 Arial，只需把各出图脚本 rcParams 里的
`"font.family"` 由 `"DejaVu Sans"` 改为 `"Arial"` 后重跑对应脚本即可（**逐图目视复核版式**，
Arial 与 DejaVu 字宽不同，个别面板的文字可能需重排）。
> 此前状态（2026-09-27 修复前）：`35`、`62` 因未设 `pdf.fonttype` 输出 **Type 3 字体**
> （文字不可检索、生产系统常拒收）；`13`、`62` 因 Arial 缺字形静默回退 DejaVu，
> 同一张图内 **Arial 与 DejaVu 两族并存**。两项均已消除。

---

## 五、本包的校验状态（2026-09-27 全绿）

```
三对齐           16/16 通过（图 5/5、表 3/3、补充 Table S1–S5 / Figure S1–S3 各自连续）
                正文引用 50 键 == 文献表 50 条；无僵尸文献、无缺号、无占位符残留
文献表体检        阻断级 0 项、卫生级 0 项（50 条均为 PubMed 逐字段核验）
字数             摘要 297（限 300）｜ 主文 7,479（ASN 不设正文字数上限）
图件体检          阻断级 0 项（400 dpi、无裁切、图内无仓库路径、图件清单 1–5 == 1–5）
图件字体          8 张图单一族 DejaVu Sans、TrueType 全嵌入、无 Type 3、无跨族混排
Word 排版         行号自摘要起、止于参考文献前；页码齐；Times New Roman 12 pt 双倍行距
补充材料          9 件，全部加页眉（Table/Methods 查 Word 页眉、Figure S1–S3 查 PDF 文本层），最大 0.33 MB（限 2 MB）
```

**复现命令**

```bash
PY=python3
cd $MR_BASE
$PY scripts/45_number_citations.py        # 引用编号 + 三对齐
$PY scripts/50_word_count.py              # 字数
$PY scripts/47_audit_references.py        # 文献表体检
$PY scripts/51_render_reference_list.py   # 同步 references_final.txt（不联网）
$PY scripts/35_make_mvmr_figure.py        # 补充图 S1（带 ASN 页眉）
$PY scripts/13_make_coloc_figure.py       # 补充图 S2（带 ASN 页眉）
$PY scripts/62_make_figureS3.py           # 补充图 S3（带 ASN 页眉）
$PY scripts/63_build_supplementary.py     # 补充材料 9 件
$PY scripts/64_make_submission_letters.py # 投稿信 + 利益冲突声明
$PY scripts/52_render_docx.py             # Word 三份
$PY scripts/57_audit_figures.py           # 图件体检（dpi / 贴边 / 印幅字号 / 字体与嵌入）
$PY scripts/65_assemble_submission.py     # 汇总到本目录
$PY scripts/67_final_check_submission.py  # 投稿包终检
```

> ⚠️ 不要重跑 `scripts/44_finalize_references.py`（联网抓 PubMed，可能覆盖已核验条目）。
"""


def main():
    DEST.mkdir(parents=True, exist_ok=True)

    missing = [(s, d) for s, d in PLAN if not s.exists()]
    if missing:
        for s, d in missing:
            print(f"[缺失] {s}")
        raise SystemExit("有源文件缺失，未做任何拷贝")

    total = 0
    for src, name in PLAN:
        shutil.copy2(src, DEST / name)
        total += 1
    (DEST / "README_投稿上传清单.md").write_text(README, encoding="utf-8")

    print(f"已汇总 {total} 个投稿件 + README 到 {DEST}\n")
    for p in sorted(DEST.iterdir()):
        if p.is_file():
            print(f"  {p.name:46s} {p.stat().st_size/1024:8.0f} KB")


if __name__ == "__main__":
    main()
