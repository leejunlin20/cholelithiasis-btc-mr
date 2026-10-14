# -*- coding: utf-8 -*-
"""修掉 STROBE-MR 清单里两处已过期的「卡在用户决定上」注释。"""
import io

P = "STROBE-MR_checklist.md"
s = io.open(P, encoding="utf-8").read()

REPS = [
    ("| \u23f3 See Comment G |", "| \u2705 |"),
    ("| \u23f3 See Comment H |", "| \u2705 |"),
    ("plus the ASN option-1 statement for code)",
     "using the ASN option-2 statement verbatim)"),
    ("**Comment G (item 18).** Blocked on a user decision: the "
     "`[[FUNDING_STATEMENT]]` placeholder. See the hand-off document.",
     "**Comment G (item 18).** Resolved. The study received no specific grant "
     "from any funding agency in the public, commercial or not-for-profit "
     "sectors; the same sentence appears on the title page under *Sources of "
     "support* and in the Acknowledgments under *Funding*. Because there is no "
     "funder, no funder role can be described. The contributing biobanks are "
     "publicly funded research resources, and their own funders are "
     "acknowledged in the source publications cited in \u00a72.2."),
    ("**Comment H (item 19).** Blocked on a user decision: the `[[REPO_URL]]` "
     "placeholder. Both the ASN data-sharing statement and \u00a72.10 "
     "(software) refer to it.",
     "**Comment H (item 19).** Resolved. All summary-level data used are "
     "public and are cited with URLs in the Acknowledgments (*Data sharing*), "
     "\u00a72.2 and **Table 1**. No repository was created for the analysis "
     "code, so the manuscript uses the ASN option-2 statement verbatim "
     "\u2014 \u201cData described in the manuscript, code book, and analytic "
     "code will be made available upon request pending application and "
     "approval\u201d \u2014 and \u00a72.10 states that the analysis scripts "
     "are available from the corresponding author upon reasonable request. "
     "No placeholder remains anywhere in the manuscript."),
]

for old, new in REPS:
    n = s.count(old)
    if n != 1:
        raise SystemExit("未命中 (%d 次): %r" % (n, old[:70]))
    s = s.replace(old, new)

io.open(P, "w", encoding="utf-8").write(s)

chk = io.open(P, encoding="utf-8").read()
print("hourglass residual:", chk.count("\u23f3"))
print("placeholder residual:", chk.count("[["))
print("option-1 residual:", chk.count("option-1"))
print("Blocked-on-user residual:", chk.count("Blocked on a user decision"))
