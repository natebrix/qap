# branch.xbar_lines: rule 9 with more candidate lines (Nov 9 2026)

Since 6335dea rule 9's `branch.xbar_lines` key admits the top-k rows and the
top-k columns by DNN marginal score as the lookahead's candidates (before, only
the single top row and column, whatever the key said), and their children's
QPB bounds pick the line. This sweep asks whether k > 1 helps. Motivation: the
line-truth studies (`../tho30_solve/line_truth/`, `../tho40_line_truth/`)
found rule 9's line 2nd to 5th of 8 at five of six nodes, with the best line
at 0.60x to 0.81x its subtree.

8 workers on the M2 Pro (Accelerate), binary at 6335dea, `bench/parallel.py`,
`sweep.sh`; logs in `logs/`. tho30 uses the rule 9 variant of the record
recipe with the CP cut phase through depth 2
(`../tho30_solve/line_truth/cutd2_cp_r9.par`); nug30 the Phase 1 recipe
(`infra/phase1/recipes/nug30.par`, the record with rule 9). The first tho30
run overlapped a video call, so its CPU is a little high; node counts are
unaffected.

| instance, lines per side | nodes | fathomed | FW iterations | CPU s | wall s |
|---|---|---|---|---|---|
| tho30, 1 (today's rule 9) | 14,256 | 13,673 | 1,572,943 | 4,980 | 746 |
| tho30, 2 | 11,912 | 11,417 | 1,768,538 | 4,283 | 625 |
| tho30, 4 | 11,915 | 11,420 | 2,539,257 | 4,223 | 605 |
| tho30, 8 | **11,205** | 10,737 | 3,646,040 | **3,942** | 544 |
| nug30, 1 (today's rule 9, the record) | **3,716** | 3,574 | 382,572 | **1,283** | 190 |
| nug30, 4 | 4,334 | 4,169 | 905,278 | 1,426 | 202 |
| nug30, 8 | 4,602 | 4,427 | 1,535,190 | 1,518 | 215 |

For reference, tho30's earlier best is rule 9 alone, Nov 8: 14,822 nodes,
4,733 CPU s; rule 8 with the same CP cut phase: 17,470.

## Reading

- **tho30: a new best.** 8 lines per side gives 11,205 nodes and 3,942 CPU s,
  21% fewer nodes and 21% less CPU than 1 line under the same par, and 24%
  fewer nodes than the Nov 8 best. Most of the gain is there at 2 lines
  (11,912); 4 adds nothing; 8 adds another 6%. The lookahead's FW iterations
  rise 2.3x (1.57M to 3.65M), but they are QPB iterations at m - 1 and cheap
  next to a DNN solve: total CPU still falls.
- **nug30: worse at every k > 1.** 4 lines gives 17% more nodes and 11% more
  CPU; 8 lines 24% and 18%. On nug30 the two marginal lines are better than
  anything the QPB children's bounds pick out of a wider shortlist.
- So the QPB child bound is not a reliable arbiter across instances: on the
  3 x 10 grid, where rows (facilities) and columns (locations) differ a lot,
  a wider shortlist finds better lines; on the 5 x 6 grid it is led astray.
  The line-truth studies scored children by a DNN first solve, which did
  correlate 0.7 to 0.9 with the subtree on both instances; the QPB bound is a
  much weaker proxy for the DNN bound that actually fathoms. The variant the
  settings plan reserves, scoring the shortlist's children by a short DNN
  solve from the node's warm state, is the version that could help both.
- Practical: keep `branch.xbar_lines = 1` as the default (the nug30 record
  stands). For tho30-like instances, 2 is most of the gain at a fifth of the
  extra lookahead; 8 is the best measured. For tho40 (5 x 8), which of the
  two regimes applies is not known; a depth-limited A/B at k = 1 and 2 on a
  few depth-2 subtrees is the cheap way to find out before a long run.
