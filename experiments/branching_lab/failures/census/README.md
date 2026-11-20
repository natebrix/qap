# Where the branching rule's bad failures happen: a census from node traces (Nov 20 2026)

branching_lab agent "failures/census". Node traces only, no solver runs. Every number is measured **[m]** from
the traces by the script named (outputs in `out/`), unless marked **[e]**.

## Answer

**On the Taillard "b" instances a bad failure is a line whose children are the free members of one tight facility
cluster, at a node far from the cutoff, and it is paid for once per sibling.** In tai60b every column with two or
more branching children (1,265 of them) had all of those children in one facility cluster, the cluster the
incumbent puts at that location [m]. The bound cannot tell near-interchangeable cluster members apart, so a third
to a half of the free members survive. In tai60b the failing lines are columns into the locations of cluster 1
(facilities 14-25, its tightest large cluster). Rule 2 takes them 1,934 times below depth 6. Five such lines (c43,
c2, c13, c54, c18) carry 318 of the run's 412 nodes with b >= 4.

At equal depth band and shortfall, these columns leave **1.85x as many branching children as rows and 3.6x as many
as other columns** [m, `q2_strata.py`]. b grows as about the 0.55 power of the node's shortfall. Once shortfall,
depth and this line kind are in the model, b is close to Poisson (dispersion 1.17), so the heavy tail is mostly
explained by those three.

Siblings repeat the failure. 88% of sibling groups in tai60b branch on one line, and b is correlated within a group
(ICC 0.40). **71% of the high-b mass comes in same-line sibling sets** of two or more siblings with b >= 3 [m]. The
mirror chain then doubles it.

Where a counterfactual exists (63 nodes seen under two or more lines), the worse line is a column in 21 of the 24
pairs with a gap of 3 or more. In tai60b, a cluster-1 column did worse than the best alternative at 17 of 27 nodes,
tied at 9 and did better at 1 (sum b 116 against 67) [m, `q4_analyse.py`].

How sure: the descriptive facts on tai60b rest on thousands of nodes and are solid. The counterfactual support is
thin: 27 nodes, most of them in the one subtree under 10:39, with parent states that differ between arms. On the
other instances it is thinner still (tai50b has 13 high-b nodes; tai35b and tai40b have none).

## Which rule chose each line, and how that was established

The trace's `rule` field is the par file's level (`latest/bnb.C` `node_trace`: "the params level that chose the
line (pick_param)"). For each run I read the level-to-rule map from the par file **at the commit the run used**:

| run | how established | levels 0-2 (depths 0-6) | level 3 (bottom) |
|---|---|---|---|
| tai60b-261117 | `recipe.sh`: commit 316361fc; `git show 316361fc:infra/phase1/recipes/tai60b.par` | 9 | **2** |
| tai50b-261116 | `recipe.sh`: 749b36b2; `tai50b.par` there | 9 | **2** |
| tai50b-r9 | `recipe.sh`: e406a047, PARAMS `tai50b_r9.par` | 9 | 9 |
| tai35b_base, tai40b_base | `run.log`: `experiments/line_signals/tai35b.par` / `tai40b.par` (= `tai/*.par`), commit 8df95737 | 9 | 2 |
| tho40ab-32 | `launches.tsv`: f27dd94, recipe tho40ab-ship32 → `tho40-ship32.par` at f27dd94 | 9 | 2 |
| nug30-trace-a | `launches.tsv`: c79368d, recipe nug30trace → `nug30.par` at c79368d | 9 | 2 |
| tai60b-c47 (linetest descent) | linetest README "Correction": `tai60b.par` at 0e6c2ecf; depth 7 forced to c47 | 9 | 9 |
| tai50b-ls | its README: `recipes/tai50b.par` at 0e6c2ecf (rule 9 everywhere); the signal chose where `ls.log` says "chose" (`lsig` column) | 9 | 9 |

In every run above, level 3 appears only at depth 7 and below (`out/q1.txt`, depth x level counts in `build.py`'s
tables). So "rule 2 levels" in tai60b means depths 7 to 56.

**Fix convention** (`fixconv.py`, `out/fixconv.txt`): a fix `f:l` puts facility f, the row index of the FIRST
matrix in the .dat (A, symmetric), at location l, the index of the second (B, asymmetric). The best known
objective is reproduced as sum A[i][j] B[p(i)][p(j)] with p from the .sln (1-based). With that convention, 57
tai60b trace nodes lie on the best known path, down to depth 57; the inverse convention puts 1 node there. In
tai-b, A is the clustered matrix (Euclidean-like, block structure) and B the heavy-tailed one. **So a row r_i
(one facility, all its free locations) has location children, and a column c_j has facility children: the
interchangeable objects.** tho40ab's fixes do not match `tho40.sln` (the run used the transposed copy or a symmetry
file), so it has no path features.

**Clusters** (`clusters.py`, `out/clusters.txt`): single-linkage on A, cut at the largest gap between minimum
spanning tree edges. tai50b has 5 clusters (15, 4, 12, 3, 16; gap 63 → 356) and tai35b and tai40b have 3. For
tai60b the largest gap (186 → 572) leaves a 36-member group that row 0 of A shows is three, so I cut at 106 → 162:
six clusters. They are 0-13, **14-25 (cluster 1)**, 26-42+45-47 (cluster 2, 20), 43-44, 48-51 and 52-59. Tightness
(the mean relative L1 distance between two members' rows of A, small = interchangeable) is 0.027 for cluster 1,
0.065 for cluster 0, 0.074 for cluster 2 and 0.007 for the 8-member cluster 5. **Locations have no cluster
structure:** on 1 - correlation of their B profiles, the largest MST gap ratio is 1.05-1.09 in all four instances.
tai60b's mirror chain (44:58) is facilities 43 and 44 exchanged. They are the 2-member cluster.

## Q1. How the work is distributed over branched nodes

`q1_concentration.py` → `out/q1.txt`; `q1_sets.py` → `out/q1_sets.txt`. b = a branched node's children that
branch. Its mean is 1 in any complete tree (every branched node but the root is someone's branching child), so the
question is how b is spread. "Opened" = the node's branching children that are off the chain (best known path, and
in tai60b the mirror), with everything below them, counted as a union. The baseline draws random branched nodes in
the same number at the same depths (20 draws). Ties in b are broken by a hash of the node, never by its subtree.

| tree | branched | b >= 4: nodes (share of b) | top 1% by excess over the depth's median: opened nodes / CPU | random at same depths | top 5% by excess: opened nodes / CPU | random |
|---|---|---|---|---|---|---|
| tai60b-261117 | 7,496 | 412 (26%) | 78% / 80% (75 nodes) | 17% / 22% | 95% / 96% | 42% / 47% |
| tai50b-261116 | 214 | 13 (29%) | 11% / 12% (2) | 3% / 4% | 54% / 42% | 19% / 17% |
| tai50b-r9 | 141 | 8 (29%) | 8% / 8% (1) | 1% / 1% | 49% / 36% | 13% / 13% |
| tai35b_base | 46 | 0 | b <= 3 everywhere | | | |
| tai40b_base | 37 | 0 | a single chain, b <= 1 | | | |
| tho40ab-32 | 644 | 47 (62%) | 62% / 58% (6) | 14% / 13% | 85% / 82% | 35% / 34% |
| nug30-trace-a | 141 | 12 (59%) | 18% / 16% (1) | 9% / 8% | 56% / 50% | 32% / 30% |

tai60b by rule and depth band (`out/q1.txt`). Direct cost is bounding the band's children (additive over the tree).
var/mean is the dispersion of b (1 = Poisson):

| rule | depths | branched | mean b | var/mean | max b | direct nodes | direct CPU | top 10% of the band by b: share of the band's b | top 10% opened: nodes / CPU (random) |
|---|---|---|---|---|---|---|---|---|---|
| 9 | 0-2 | 10 | 2.10 | 2.8 | 6 | 0.2% | 0.5% | 29% | 1.7% / 3.4% (0.8% / 1.6%) |
| 9 | 3-6 | 216 | 1.39 | 3.0 | 12 | 3.6% | 9.1% | 46% | 29% / 51% (11% / 20%) |
| 2 | 7-10 | 1,880 | 1.22 | 2.2 | 10 | 28.9% | 52.3% | 43% | 81% / 81% (45% / 48%) |
| 2 | 11-15 | 1,965 | 0.85 | 1.0 | 5 | 28.1% | 21.3% | 30% | 45% / 24% (19% / 10%) |
| 2 | 16-22 | 2,861 | 0.96 | 1.7 | 6 | 34.6% | 16.0% | 38% | 26% / 12% (20% / 9%) |
| 2 | 23+ | 564 | 0.82 | 0.8 | 5 | 4.7% | 0.8% | 29% | 2.7% / 0.4% (1.6% / 0.2%) |

Reading: **tai60b's cost is a few hundred bad branchings, not an even spread.** b is overdispersed (var/mean 2-3)
at depths 3-10, where the cost is, and about Poisson below. The 75 nodes with the largest excess (1%) opened
subtrees holding 80% of the CPU, against 22% for 75 random branched nodes at the same depths. The 412 nodes with
b >= 4 (5.5%) carry 26% of all branching children. 323 of them are rule-2 columns into cluster-1 locations, and
those opened 83% of the nodes, against 50% at random [m, `out/q1_sets.txt`]. Caution: "opened" is nested. A node
high in the tree opens everything below it, and large shortfalls make large subtrees whatever the line. The
random-at-same-depth baseline controls for depth, not for shortfall, so these shares overstate what the line alone
decided. Under rule 9 alone (tai50b-r9) the tail is short: nothing below depth 6 has b > 2. At depths 7-10, rule 2
reaches 4 in tai50b-261116 and 10 in tai60b.

## Q2. What the high-b nodes have in common

`q2q3_features.py` → `out/q2.txt`; `q2_strata.py` → `out/q2_strata.txt`. "high" = b >= 4.

**Row or column, by rule (tai60b)**:

| rule | depths | side | nodes | mean b | max b | b >= 4 | rate |
|---|---|---|---|---|---|---|---|
| 9 | 3-6 | c / r | 107 / 109 | 1.29 / 1.50 | 12 / 11 | 9 / 10 | 8% / 9% |
| 2 | 7-10 | c / r | 1,464 / 416 | 1.23 / 1.18 | 10 / 4 | 196 / 12 | **13.4% / 2.9%** |
| 2 | 11-15 | c / r | 396 / 1,569 | 0.24 / 1.00 | 3 / 5 | 0 / 17 | 0% / 1.1% |
| 2 | 16-22 | c / r | 2,409 / 452 | 1.00 / 0.75 | 6 / 2 | 160 / 0 | **6.6% / 0%** |

In the c47 descent (rule 9 below depth 7), columns average 1.98 with a maximum of 13 and rows 0.81 with a maximum
of 3. In the signal run tai50b-ls, columns average 4.67 and rows 0.45 [m].

**The cluster, which matters more than the side.** Every tai60b column with b >= 2 (1,265) had all its branching
children in one facility cluster, the incumbent's cluster at that location. Every row with b >= 2 (663) had all its
branching children among the free locations that hold the row facility's cluster in the incumbent. The share of
those candidates that branch is 0.32 for columns (0.48 when b >= 4) and 0.49 for rows. In tai50b the same holds
for 96% of columns with b >= 2 (rule 2 run) and 86% (rule 9 run). So the line's candidate set is the free part of
one cluster, and the line fails when that set is large and its members interchangeable. By cluster, tai60b, rule 2
levels:

| line kind (rule 2) | nodes | mean b | max b | b >= 4 |
|---|---|---|---|---|
| column, incumbent's facility at the location in cluster 1 (14-25, tightness 0.027) | 1,934 | **1.77** | 10 | **323** |
| column, cluster 2 (20 members, 0.074) | 2,279 | 0.42 | 5 | 33 |
| column, cluster 0 / 4 / 5 | 267 / 201 / 98 | 0.22 / 0.88 / 1.04 | 3 / 3 / 5 | 0 / 0 / 5 |
| row, facility in cluster 0 (rule 2's rows are almost all cluster 0) | 2,437 | 0.99 | 5 | 29 |

Controlled for shortfall: within each depth band x shortfall quintile, cluster-1 columns have the largest mean b in 13 of the 18 cells that hold all three kinds
(`out/q2_strata.txt`). Depths 7-10, top quintile: cluster-1 columns 4.37, other columns 1.44, rows 1.96. A
quasi-Poisson regression of b on log(shortfall / root's), depth band and line kind (rule-2 levels, 7,270 nodes) gives
rate ratios **cluster-1 column 1.85, other column 0.51 (rows = 1)** and a shortfall elasticity of 0.55, with residual
dispersion 1.17 [m]. At rule 9 levels (226 nodes): cluster-1 column 1.76 (s.e. of the coefficient 0.20), other
column 1.15 (0.15), elasticity 1.3.

The incumbent-based label misses off the path. At the c47 descent's 8:49 and 8:57, c30's location holds cluster 4
in the incumbent, but 11 of c30's 12 branching children (and 6 of 6) are cluster-1 facilities. Off the chain the
subtree's own arrangement decides which cluster "owns" a location, and the trace does not show that before
branching.

**Which line recurs (tai60b)**: c43 (96 high-b nodes of 678), c2 (93 of 282), c13 (70 of 530), c54 (47 of 97), c18 (12)
are all cluster-1 locations in the incumbent (cluster-1 locations: 2, 9, 12, 13, 17, 18, 21, 25, 31, 43, 48, 54).
c50 (cluster 2) is used 1,515 times with mean b 0.47. In the c47 descent: c30 (3 high-b nodes, b up to 13), c26, c47.

**Shortfall** is the strongest single correlate. Within depth bands, Spearman(b, shortfall / root's) is 0.41-0.76 in
tai60b and 0.66-0.87 in tai50b. High-b nodes sit at 2-3x the band's median shortfall (depths 7-10: 0.065 against
0.022 of the root's). QPB's distance below the cutoff is weak and inconsistent (rho -0.16 to 0.43).

**Path and mirror**: chain nodes have the largest b at every band (mean 5-6.25 at depths 3-10) because they have
the largest shortfalls. Off the chain, nodes 1-2 fixes off the best known solution average 1.25-2.0, and 3-5 off
average 0.7-0.8 (depths 7-22).

**How the children went**: in tai60b the high-b nodes' survivors are real survivors, not near misses. Their
branching children's shortfall is a median 0.52 of the parent's (others 0.61), with 52% above half. The cut phase
fathoms 3.6% of the rest, against 3.8% elsewhere. 17-23% of tai60b's branching children have a bound below their
parent's (the child's first solve plus a 1,000-iteration phase does not regain the parent's strength). tai50b
differs: its high-b nodes have many near misses (10.8 children per node fathomed within 1% of the parent's
shortfall, against 1.9; 17% of fathoms by the cut phase against 3%; survivors a median 0.34 of the parent's). So in
tai50b part of the tail is an edge effect, and in tai60b it is not.

**Other instances**: tho40ab-32 (grid instance, no cluster analysis) has the heaviest tail (top 5% of branched nodes
hold 51% of b). Its high-b nodes are all rule-9 columns at depths 3-6 (c39, c8, c7, c24). In nug30, rows and columns
are alike.

## Q3. Sibling repetition

`q2q3_features.py` → `out/q3.txt`. A group = the branched children of one node (two or more).

| run | groups | all on one line | ICC of b | high-b mass in same-line sets of >= 2 siblings with b >= 3 | high-b nodes with no branched sibling |
|---|---|---|---|---|---|
| tai60b-261117 | 1,928 | **88%** | **0.40** | **71%** | 6% |
| tai60b-c47 | 29 | 24% | 0.36 | 42% (the 8:49 / 8:57 pair, c30: 12 + 6) | 12% |
| tai50b-261116 | 34 | 76% | 0.03 | 6% | 52% |
| tai50b-r9 | 22 | 27% | -0.19 | 0% | 57% |
| tho40ab-32 | 116 | 51% | 0.22 | 90% | 8% |
| nug30-trace-a | 30 | 7% | 0.54 | 60% | 11% |

In tai60b the rank correlation of b between two siblings is 0.52 at depths 7-10, -0.08 at 11-15 and 0.49 at
16-22. The largest repeated sets are the depth-8 children of 10:39 on both chains: c2 at all 10 siblings, b = 7,
4, 7, 7, 7, 5, 0, 7, 7, 5, on the path chain and the mirror chain alike. Then c54 at depth 7 under
`...,6:59,22:31,12:0`: 7 siblings, b 6-8. Rule 2's line is a property of the subtree (the U-sum tie-break lands on
the same location in every sibling), so a bad line is paid once per sibling, and the mirror doubles it. Rule 9's
lines vary more between siblings (24-27% one line in the rule-9 trees), and its failures come singly or in pairs.

## Q4. Against the counterfactuals

`q4_counterfactual.py` → `out/q4_obs.tsv` (only nodes seen more than once); `q4_analyse.py` → `out/q4.txt`. Sources:
production traces, the c47 descent, every linetest arm (scan, arm17, d9x, sc54 + lvl/c54r9, c54, c47, r31, kids),
tai60b_lines' depth-2 arms, the tai50b runs (rule 2, rule 9, signal), and tai35b/tai40b's line_signals arms. Nodes
are matched as sets of fixes.

- **Noise**: the same node under the same line in two sources, 283 pairs: |delta b| = 0 in 99%, at most 2 (tai60b
  10:39, cold against strong parent state: c47 5/7, c54 9-11).
- **63 nodes under two or more lines, 79 pairs. Of the 24 pairs with a gap of 3 or more, the worse line is a
  column in 21 (88%), and the better line in 7 (29%).** The worse line has more candidates (free members of its
  cluster) in 67% of them and fewer in 17%. In tai60b (23 pairs) the worse line is a cluster-1 column by the
  incumbent in 11 and has mostly cluster-1 children in 15.
- **Reverse test (tai60b)**: at the 27 nodes where a cluster-1 column met another line, the cluster-1 column was
  worse than the best alternative at 17, tied at 9 and better at 1 (sum b 116 against 67). The large gaps are at
  depths 7-8 under 10:39: c54 10 against c47 6; c2 4-7 against 1-4 at all 9 compared siblings; c54 8-9 against
  c10/c53 3. At depths 11-17 the gap is 0-2, as the shortfall elasticity predicts.
- **The exceptions are columns into the large, looser cluster 2.** At tai60b depth 2 (tai60b_lines), c26 gives 20
  and rows r1, r45, r8 give 6, 10 and 12. At tai60b depth 17 on the path, c26 gives 7 and c54 / r13 give 3 / 2. At
  tai50b's root, c26 gives 28 (clusters 2 and 4: 12 + 16) and r31 gives 10. Rows can also differ by 2x near the top
  (r1 6 against r8 12 at the same depth-2 node).
- **What does not separate the lines**: shortfall, QPB gap and position are properties of the node, the same for
  both lines up to the parent-state noise. They say when a failure is expensive, not which line fails.

The pairs are not independent. 9 of the 27 reverse-test nodes are siblings under one depth-7 node. The depth-8 pairs
compare production's parent state (rule 2's c2) with the strong arm's (rule 9's lines). linetest checked that 2 of
those parents, rebranched under rule 2 from the strong state, reproduce production's c2 counts exactly.

## Q5. What the traces say, and what they cannot

From the traces alone, the best description is this. On tai-b, a line's surviving children always come from one
facility cluster: the free members of the cluster that owns the location (a column), or the free locations of the
facility's own cluster (a row). How many survive scales with the node's shortfall (b ~ shortfall^0.55) and with how
interchangeable the cluster's members are. A column into the region of a tight, mostly unfixed cluster (tai60b's
14-25) asks which near-identical facility goes there, and the bound cannot answer: a third to a half survive, 7-12
at depths 7-10. Rule 2's U-sum tie-break lands on such columns (c54, c2, c43, c13) and repeats the choice at every
sibling, so the failure is multiplied by the sibling count and again by the mirror. Rule 9 falls into the same trap
more rarely and singly (c26 at depth 6, c30 off the chain at depths 12 and 17). Near the root, columns into the big
looser cluster (c26 at tai60b depth 2 and at tai50b's root) fail the same way.

I am confident of the cluster pattern and the shortfall scaling in tai60b (thousands of nodes, every case). I am
moderately confident that cluster-1 columns are worse than the alternative at the same node (27 nodes, mostly one
subtree). I am unsure how far this carries beyond tai60b: tai50b's tail is small and partly near misses, and tai35b
and tai40b barely branch.

To go further, the solver would have to dump, at each branched node:

- the DNN marginals of the chosen line and of rule 9's two candidates (the most decided row and column), to test
  whether a failing line's mass is spread over one cluster's members. That is the trace-free version of the
  cluster label, and it would catch the off-path cases where the incumbent's cluster is wrong (c30);
- rule 9's QPB lookahead values for both candidates;
- rule 2's U sums for its top lines and its count.

Then one can ask whether rule 9 saw a column with flat marginals over a cluster and took it anyway (lookahead noise),
or never saw the row. Per-child first-solve bounds for the line not taken would need a probe, not a dump
(line_signals: 30 warm iterations per child).

## What I could not determine

- **Causality beyond 27 nodes.** The census shows where high b is, and the counterfactuals that exist agree. But most
  of them are one subtree of tai60b (under 10:39), and several pair different parent states.
- **A predictive label off the chain.** The cluster that owns a location in an off-chain subtree is not in the trace
  (c30 at 8:49 holds cluster 4 in the incumbent but behaved as cluster 1). The marginals would show it.
- **The subtree cost of each failure as distinct from its shortfall.** The "opened" shares are nested and are not
  controlled for shortfall. No counterfactual subtree below depth 18 exists.
- **tai50b-ls and the c47 descent are incomplete trees.** They are used for b only, at nodes whose children are all
  traced (112 of 131 and 117 of 148). tai50b-r9 has 3 nodes with missing children, which are excluded.
- **tho40ab and nug30**: no path (tho40ab's fixes do not match `tho40.sln`) and no cluster analysis (grid
  instances).
- Whether rule 9 would have avoided the depth 7-10 cluster-1 columns in production tai60b is shown only where
  linetest bounded both rules (depths 8 and 9 under 10:39: rule 9 halves b, `../../linetest/README.md` 2x2).

## Files

`common.py` (trace loader; the run → files → level-to-rule map), `fixconv.py` (fix convention, `out/fixconv.txt`),
`clusters.py` (`out/clusters.txt`, `.json`), `build.py` (one row per branched node: `out/bn_<run>.tsv.gz`, with
column docs at its head), `q1_concentration.py` (`out/q1.txt`), `q1_sets.py` (`out/q1_sets.txt`),
`q2q3_features.py` (`out/q2.txt`, `out/q3.txt`), `q2_strata.py` (`out/q2_strata.txt`), `q4_counterfactual.py`
(`out/q4_obs.tsv`, nodes seen more than once; 50 kB), `q4_analyse.py` (`out/q4.txt`). Rerun in order: `clusters.py`,
`build.py`, then the rest (each under a minute; the tai60b trace loads in about 8 s).
