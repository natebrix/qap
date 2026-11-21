# Rule 9's failures: what rule 9 saw at the c47 descent's labelled nodes (Nov 20-21 2026)

The coordinator's question for Nathan (Nov 20, 23:15 UTC): when do rule 9's bad failures happen? This part uses the
states this session holds. A census of the traces is in `experiments/branching_lab/failures/census/` on
claude/cool-babbage-ilnogj, and is not repeated here. Every number is measured [m] unless marked [e].

## Answer

**Every failure is a column whose children are near-interchangeable facilities of one cluster, and rule 9 takes it
for one large entry it shares with a row. The good line was never one of rule 9's two candidates.** [m]

- **The failures.** In this set there are four: three on c30 and one on c26. At 11:29 (depth 12, added from the
  sibling scan) c30 left 13 branching children against rule 2's c54 at 2. At 8:49 and 8:57 (depth 17) c30 left 12
  and 6, against r2's 0. At 22:31 (depth 17, on the path) c26 left 7, against c54's 3 and r13's 2. At the other 24
  nodes, rule 9's line is within 2 branching children of the best alternative bounded.
- **What rule 9 saw.** Its candidates are the top row and the top column by sum of x^2, and a QPB lookahead chooses
  between them. At c30's three nodes the two candidates are r51 and c30. They share one dominant cell, 51:30
  (x = 0.80-0.90), which alone makes both lines "most decided". The lookahead then splits them by 0.1-0.6%. Rule
  9's arbitration is that close at almost every node: the two values differ by 0.1% to 5%.
- **Shortlist or arbitration?** The failure is in the shortlist: the better line was never a candidate. r2 ranks
  9th and 15th by marginal score at 8:49 and 8:57, though 2nd by the full lookahead. Rule 2's c54 ranks 26th at
  11:29, and at 22:31 c54 ranks 18th and r13 8th. Rule 2's lines rank 1st to 4th by the full QPB lookahead almost
  everywhere (the table's "by LA"). Within the shortlist, the arbitration's pick at 8:49 (c30, 12) was worse than
  the other candidate (r51, 5); at 8:57 it was worse by 6 against 1.
- **Is c30 "decided" by one large entry or by several medium ones?** One large entry. At 8:49 c30 has 2 entries
  above 0.02, 2 above 0.05 and 1 above 0.1 (51:30 = 0.856). r2 has 2, 2 and 2 (2:23 = 0.546 and 2:57 = 0.455).
  The difference is in the tail. c30 has 22 entries above 0.001; r2 has 2, and its other 41 are exactly 0.
- **Clusters** (the census's question). The relaxation puts c30's mass on cluster 4 (facilities 48-51: 0.80-0.91,
  almost all of it on facility 51). But **12 of c30's 13, 11 of its 12 and all 6 of its branching children are
  cluster-1 facilities** (14-25). Those carry only 0.07-0.18 of the mass, spread thinly over 9-12 of cluster 1's 12
  free members, every one above 0.001 (`clusters_xbar.py`, `out/clusters_xbar.txt`). **So the mass does not identify
  the cluster whose members survive; the non-zero support does.** The cluster with many members carrying small but
  non-zero x is where the survivors are. At 22:31, c26's mass and its 14 non-zero entries are all in cluster 2,
  and its 7 survivors are all cluster 2.
- **Does anything separate the 4 failures from the 24 ordinary nodes?** Not cleanly, and four failures cannot
  support a fitted threshold. Shortfall, depth, path, the top score and its margin, and the row-or-column choice
  do not separate them (table). The nearest descriptive feature is that **all four are columns with at least 14
  entries above 0.001, 9 to 14 of them in one cluster.** Among rule 9's ordinary choices, columns have 3-10 such
  entries, with one exception: s17_9_29's c6, with the same shape as c30 (0.88 on cluster 4, 11 non-zero cluster-1
  members), left 2. Rule 2's habitual c54 has 9-10 non-zero cluster-1 entries everywhere and leaves 1-5 branching
  children. Rows never fail here: rule 9's rows leave 0 to 3.
- **"Support"** (fewest entries above 1/m, ties by rule 9's order; fixed before the data). **It fails its decision
  rule.** At 8:49 it picks r51 (5 branching children, against at most 2 allowed). At 8:57 it picks c30, rule 9's
  own failing line (6). Entries above 1/m do not see the tail: c30 has only 1-2 there, as r2 does. So it is not
  worth a wider test. Its lines at the other failures were bounded for completeness: c4 at 22:31 and c38 at 11:29
  (below). The six ordinary-node arms were not run, since the decision was already settled by 8:49 and 8:57.
- **Siblings.** 8:49 and 8:57 are near-copies: children of 4:29 under r8, differing only in facility 8's location
  (49 against 57). Their bounds are 608,024,831 and 607,991,775, and their c30 columns are the same shape (51:30 =
  0.856 and 0.896, 22 non-zero entries each, the same cluster split). So it is one event seen twice. The same thing
  happens at depth 12. Under 3:35 the siblings 11:29 and 11:59 both chose c30, leaving 13 and 1. In the whole
  descent (`siblings.py`, `out/siblings.txt`), these are the only sibling groups in which the same line leaves many
  branching children. In every other repeated-line group each sibling leaves 0-3. A c30 failure then makes the
  bulge one level down. Its branching children are near-interchangeable cluster-1 facilities at location 30, and
  they in turn choose one line together: r13 for 12 of 11:29's 13, c16 for 11 of 8:49's 12. Under 11:29 those 13
  then left no branching children at all (depth 13). Under 8:49 they were not bounded further.

## The 28 nodes

Rule 9's line and its branching children (B) are the descent's. The "alternative" lines are those bounded at the node
(rule 2's line, and at 8:49, 8:57 and 22:31 the signal's). Ranks are by sum of x^2 and by the full QPB lookahead (LA,
`QAP_LINE_SCORES_LA`). Entries are counts of x-bar above 0.02 / 0.05 / 0.1 / 1/m / 0.001. Shortfall = cutoff - the
node's bound; "sf/root slack" divides it by cutoff - the root bound (605,315,696). f12 is 11:29 at depth 12, the c30
failure found by the sibling scan (rebuilt as `runs/rebuild/b335`; its rule-2 arm is `runs/scan/f12_rule2`).

| node | depth | path | shortfall | sf/root slack | rule 9: line (B) | cands row/col (LA row, LA col) | top score, margin | rule 9 line: entries >0.02/0.05/0.1/>1/m/>0.001 | alternative: line (B) | alt rank by x^2 / by LA | alt entries >0.02/0.05/0.1/>1/m/>0.001 | support line (B if known) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| d11_a | 11 | no | 186k | 0.064 | r11 (3) | r11/c59 (1.56e+12, 1.53e+12) | 0.663, 0.022 | 3/3/2/3/5 | c54 5 | 40 / 1 | 6/4/3/6/12 | c59 |
| d11_b | 11 | no | 156k | 0.054 | r11 (3) | r11/c59 (1.56e+12, 1.53e+12) | 0.858, 0.002 | 3/2/1/3/3 | c54 3 | 31 / 1 | 6/3/3/6/10 | c59 |
| d11_c | 11 | yes | 134k | 0.046 | c40 (3) | r32/c40 (1.51e+12, 1.52e+12) | 0.604, 0.014 | 3/3/2/3/6 | c54 3 | 42 / 1 | 5/4/3/5/10 | c59 |
| f12 | 12 | no | 94k | 0.032 | c30 (13) | r51/c30 (1.52e+12, 1.53e+12) | 0.598, 0.009 | 6/1/1/6/21 | c54 2 | 26 / 7 | 5/3/3/5/11 | c38 |
| d13_a | 13 | yes | 137k | 0.047 | r11 (3) | r11/c59 (1.52e+12, 1.49e+12) | 0.717, 0.017 | 3/3/1/3/6 | c54 4 | 34 / 1 | 5/4/3/5/10 | c59 |
| d13_b | 13 | no | 116k | 0.040 | r11 (3) | r11/c59 (1.52e+12, 1.49e+12) | 0.858, 0.003 | 3/1/1/3/3 | c54 3 | 39 / 1 | 5/4/3/4/10 | c59 |
| d13_c | 13 | no | 92k | 0.032 | r12 (1) | r12/c0 (1.55e+12, 1.52e+12) | 0.725, 0.001 | 3/3/1/3/5 | r8 1 | 10 / 3 | 3/3/3/3/4 | c38 |
| d15_a | 15 | no | 61k | 0.021 | r4 (2) | r4/c23 (1.53e+12, 1.5e+12) | 0.724, 0.000 | 3/2/2/3/3 | r8 1 | 3 / 1 | 3/3/2/3/3 | c38 |
| d15_b | 15 | yes | 112k | 0.039 | c6 (2) | r48/c6 (1.44e+12, 1.44e+12) | 0.563, 0.008 | 4/2/2/4/6 | c54 3 | 28 / 1 | 5/3/2/5/10 | r22 |
| d15_c | 15 | no | 45k | 0.015 | r51 (1) | r51/c30 (1.43e+12, 1.43e+12) | 0.710, 0.008 | 3/2/2/3/9 | c54 1 | 29 / 1 | 4/3/3/4/9 | r51 (1) |
| s16_12_0 | 16 | no | 32k | 0.011 | r9 (2) | r9/c55 (1.49e+12, 1.43e+12) | 0.842, 0.000 | 2/2/1/2/4 | r8 1 | 12 / 3 | 3/3/3/3/3 | c55 |
| s16_48_6 | 16 | yes | 90k | 0.031 | c31 (2) | r22/c31 (1.42e+12, 1.42e+12) | 0.613, 0.000 | 3/2/2/3/7 | c54 3 | 26 / 1 | 5/3/2/3/10 | c56 |
| s16_4_23 | 16 | no | 27k | 0.009 | r8 (1) | r8/c29 (1.51e+12, 1.43e+12) | 0.882, 0.001 | 2/2/1/2/2 | r8 1 | 1 / 3 | 2/2/1/2/2 | r8 (1) |
| s16_4_29 | 16 | no | 318k | 0.110 | r8 (2) | r8/c23 (1.54e+12, 1.52e+12) | 0.698, 0.002 | 3/3/2/3/3 | c23 0 | 2 / 4 | 3/3/1/3/3 | c31 |
| s16_50_6 | 16 | no | 111k | 0.038 | r9 (2) | r9/c55 (1.47e+12, 1.42e+12) | 0.569, 0.007 | 3/2/2/2/6 | c54 2 | 19 / 1 | 5/3/2/5/10 | r9 (2) |
| s16_50_6a | 16 | no | 46k | 0.016 | r9 (1) | r9/c55 (1.47e+12, 1.42e+12) | 0.566, 0.004 | 2/2/2/2/6 | c54 1 | 13 / 1 | 5/3/2/5/10 | r9 (1) |
| s16_51_30c | 16 | no | 51k | 0.018 | c6 (1) | r48/c6 (1.41e+12, 1.42e+12) | 0.753, 0.001 | 3/2/1/2/6 | c54 1 | 36 / 1 | 3/3/2/3/10 | c6 (1) |
| p2231 | 17 | yes | 86k | 0.030 | c26 (7) | r45/c26 (1.39e+12, 1.39e+12) | 0.537, 0.021 | 6/3/1/6/14 | c54 3, r13 2 | 18 / 8 | 3/3/3/3/4 | c4 |
| p849 | 17 | no | 190k | 0.066 | c30 (12) | r51/c30 (1.46e+12, 1.46e+12) | 0.740, 0.002 | 2/2/1/2/22 | r2 0, r51 5 | 9 / 2 | 2/2/2/2/2 | r51 (5) |
| p857 | 17 | no | 223k | 0.077 | c30 (6) | r51/c30 (1.44e+12, 1.45e+12) | 0.804, 0.001 | 1/1/1/1/22 | r2 0, r51 1 | 15 / 2 | 2/2/2/2/2 | c30 (6) |
| s17_25_31 | 17 | no | 139k | 0.048 | c43 (1) | r18/c43 (1.38e+12, 1.43e+12) | 0.602, 0.004 | 5/2/1/4/10 | c54 2 | 10 / 1 | 5/4/1/4/9 | c9 |
| s17_48_6 | 17 | no | 38k | 0.013 | c11 (1) | r26/c11 (1.39e+12, 1.4e+12) | 0.712, 0.003 | 3/2/2/3/9 | c54 2 | 35 / 1 | 3/3/2/3/9 | c4 |
| s17_8_29 | 17 | no | 86k | 0.030 | r5 (0) | r5/c57 (1.52e+12, 1.44e+12) | 0.526, 0.000 | 2/2/2/2/2 | r2 1 | 4 / 2 | 2/2/2/2/2 | r5 (0) |
| s17_9_0 | 17 | no | 25k | 0.009 | c4 (0) | r49/c4 (1.41e+12, 1.45e+12) | 0.579, 0.000 | 2/2/2/2/5 | r8 1 | 13 / 1 | 4/4/2/4/4 | r49 |
| s17_9_29 | 17 | no | 201k | 0.069 | c6 (2) | r48/c6 (1.42e+12, 1.43e+12) | 0.661, 0.001 | 3/1/1/3/23 | c23 0 | 13 / 5 | 4/4/2/4/4 | c30 |
| s17_9_55a | 17 | no | 24k | 0.008 | r13 (0) | r13/c36 (1.46e+12, 1.43e+12) | 0.738, 0.000 | 2/2/2/2/3 | r8 0 | 7 / 2 | 3/2/2/3/4 | r13 (0) |
| s17_9_55p | 17 | no | 14k | 0.005 | r4 (1) | r4/c23 (1.48e+12, 1.45e+12) | 0.586, 0.000 | 3/2/2/3/3 | r8 1 | 5 / 1 | 3/3/3/3/3 | r4 (1) |
| s17_9_55q | 17 | no | 97k | 0.033 | r13 (1) | r13/c36 (1.46e+12, 1.43e+12) | 0.742, 0.002 | 2/2/2/2/3 | r8 1 | 8 / 2 | 3/2/2/3/4 | r13 (1) |

"Support" at the failures: 22:31's c4 and 11:29's c38 are in `runs/scan/sup_*` (results below when finished).

## How

`dumpnode.sh` re-bounds each node from its parent's state (the rebuilt states, every one checked against the run).
It runs under `tai60b.par` (rule 9 below depth 6, as in the descent) with `-d DEPTH`, which sets the children aside,
and dumps U, x-bar, the bound and every line's lookahead (`QAP_LINE_SCORES=DEPTH`, `QAP_LINE_SCORES_LA=1`), plus
rule 9's arbitration (`QAP_XBAR_DEBUG=2`'s `xbar9` line). `failures.py` builds the table, `clusters_xbar.py` the
cluster split, and `siblings.py` the sibling groups. The rule-2 lines are from the earlier arms (`runs/scan`,
`runs/arm17`). The clusters are the census's (`clusters.json` there: 0-13, 14-25, 26-42+45-47, 43-44, 48-51, 52-59).
