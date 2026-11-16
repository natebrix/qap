#!/usr/bin/env python3
"""Assembles README.md from its skeleton (README.md with @@ANSWER@@, @@STAGE3@@), answer_text.md and
stage3_text.md. Run once; kept for the record of how the README was put together."""
import sys
tai_row, tai_text, tai_note = sys.argv[1], sys.argv[2], sys.argv[3]
r = open("README.md").read()
a = open("answer_text.md").read().strip()
s = open("stage3_text.md").read().strip()
a = a.replace("2,248 certified values", "3,248 certified values")
a = a.replace("concentrated on one fixed facility\n   (88 to 98% of the mass on facility 8 at corner 0 in the tho40 subtree `26:10,8:0`)",
              "concentrated on one or two fixed facilities\n   (88 and 98% of the mass on facility 8 at corner 0 at two of the three tho40 nodes)")
a = a.replace("@@N10GAIN@@", "**+2.7** on the best round (-32, -12 in the two rounds after the triangles' own overshoot)")
a = a.replace("@@TAINOTE@@", tai_note)
s = s.replace("@@N10ROW@@", "| wil50 node 10 `39:27,30:48,25:49` | 3 | 47 | 6,627 | 355 (5.4%) | 68 | 79% on 30@48 | `logs/s3_wil50_n10_stats.log` |")
s = s.replace("@@TAIROW@@", tai_row)
s = s.replace("@@D4EROW@@\n", "")
s = s.replace("@@N10TEXT@@", """Node 10 (`39:27,30:48,25:49`; 3 rounds, no quiet round, `logs/s3_wil50_n10_jT.log`,
`..._jTXm.log`): 48,667.7 against 48,670.4 after round 1 (**+2.7**); in round 2 the library's
own triangle step overshoots in both arms (48,501.2 / 48,469.4) and round 3 has not recovered
(48,606.9 / 48,595.4), so the best certified values are round 1's. Node 10 shows 68 units of
violation, 4.5 times node 13's; at the tenth that tho40 and node 13 realise that is about +7
[e], still under the bar.""")
s = s.replace("@@TAITEXT@@", tai_text)
s = s.replace(" @@PATH8@@", "")
r = r.replace("@@ANSWER@@", a).replace("@@STAGE3@@", s)
assert "@@" not in r, [x for x in r.split() if "@@" in x]
open("README.md", "w").write(r)
