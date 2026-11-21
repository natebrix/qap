# The rule "tail": branch on the line with the fewest non-zero marginals

Nov 21 2026. Status: **not established on whole trees.** Ungated it helps on tai35b (13% fewer nodes) and
hurts on tai40b (at least 15% more: a weak root). Gated (the rule only where the count is at most m / 2) it
has no effect on either (+3%, -0.5%). Its support is still the one-level evidence at tai60b's tight nodes,
where the counts are 2 to 5 against 22; no whole tree here has such nodes going wrong under rule 9.
[m] = measured here, [e] = estimate or inference.

## The rule

At a node with DNN marginals x-bar: count, for each free row and column, the entries above 0.001. Branch on the
line with the fewest. Ties: the larger sum of squares (rule 9's score), then rows before columns, then the lower
index. It costs nothing: no lookahead, no QPB arbitration.

It came from rule 9's failures (`../linetest/failures.md`): rule 9 scores a line by its mass (sum of squares),
and a column with one large entry and a long tail of small non-zero entries scores well although the tail's
facilities (one tight cluster's members) all survive as children. The count sees the tail; the mass does not.

In the solver: `QAP_TAIL=THR[:DMIN:DMAX]` (`tail.patch`, on cognac 0e6c2ecf + `experiments/line_signals/
experiment_build.patch` + `../lsla/lsla.patch`). It replaces the params' rule at every node of depth DMIN..DMAX
that has marginals; `QAP_TAIL_LOG` gets one line per choice. It may choose a column, so it is for instances
without a symmetry file. Checked on tai25b's root and depth 1: the solver's five choices equal a separate
calculation from the dumped x-bar (`QAP_LINE_SCORES`) [m]. `QAP_TAIL_GATE=F` (added after tai40b): the rule
only where its fewest-entries count is at most F * m, the params' rule elsewhere. `lines.py` prints a node's
lines with their counts from such a dump.

## One level, tai60b

- **In the subtree where it was found** (`../linetest/failures.md`, "The tail rule") [m]: at rule 9's four
  failures the rule's line leaves 0, 0, 0 and 2 branching children (rule 9: 13, 12, 6, 7). Over the 18 nodes
  where both lines were bounded: 21 branching children against 63; worse than rule 9 at one node, by one.
  In-sample: the feature was noticed on these nodes.
- **Out of sample, the chain nodes at depths 7, 8, 9** (`../lsla/out/tail.md`) [m]: the rule picks c59 at all
  three, the best line bounded at each (3 branching children; the others 4 to 12). At depth 9 the landscape
  signal's c30 (12 branching children) has 23 entries and ranks 91st of 102.
- **Where it says nothing** [m]: at tai50b's root every line has 49 or 50 entries above 0.001 and at tai25b's
  root 23 to 25 of 25, so the choice falls to the tie-break and the rule is rule 9's first candidate. At a weak
  node the relaxation is spread and the count is flat.
- **Threshold** [m, linetest]: over 28 nodes 0.0001 changes the choice at 4 and 0.01 at 13. Ties at the minimum
  are common (up to 8 lines), so the tie-break does part of the work.

## Whole trees (`run.sh`)

The params of `../lsla` (rule 9 through depth 6, rule 2 below; DNN through depth 6), 2 slots, incumbent best
known + 1. The reference is `../lsla/t2/PROBLEM_prod`, the same params without the rule. Decision rule, fixed
before the runs: at most 95% of the reference's nodes, the rule helps; at least 105%, it hurts.

| tai35b | trace lines | solver's node count | CPU s | branched by depth 0-8 |
|---|---|---|---|---|
| production (rule 9, rule 2 below) | 1,017 | 987 | 6,935 | 1 3 3 4 5 2 2 2 1 |
| tail at every node with marginals | **884** | **822** | 6,527 | 1 2 3 4 2 2 2 1 1 |
| landscape signal, depths 0-2 (`../lsla`) | 983 | | 6,352 | |
| DNN lookahead (`line_signals` (c)) | 780 | | 8,900 | |

- **tai35b: 13% fewer nodes than production (884 against 1,017 trace lines), closed at 283,315,445** [m
  `runs/tai35b_tail`]. It passes the rule. 16 choices; the first three (c26, then r19 and r18) are the
  lookahead's and the signal's; it departs at depth 2 (r33).
- CPU is 6% lower, but the two runs were on different containers and this one shared its 4 cores with another
  run, so only the node count is a clean comparison [e].
- The saving is at depths 2, 5 and 8 (68 / 62 / 28 nodes against 102 / 155 / 56): one fewer branching node at
  depths 1, 4 (three fewer) and 7.
- tai35b's tree has about 20 branching nodes above depth 7, so this is one small tree, not a rate.

### tai40b: the ungated rule hurts

| tai40b | trace lines | CPU s | branching nodes at depths 1-7 |
|---|---|---|---|
| production | 731 (whole tree) | 6,925 | 1 1 1 1 1 1 1 |
| tail at every node with marginals | 743 through depth 7, **at least 843** | 15,293 so far | 3 3 4 4 3 3 3 |

- **Not closed: stopped twice by the two-hour limit on background commands and not restarted, because the
  partial tree decides it** [m `runs/tai40b_tail/nodes.partial.txt.gz`, the units' traces at the stop]. 743
  lines below the root, and the three branching nodes at depth 7 have 33 children each not yet bounded: at
  least 843 lines with the root, 115% of production. It hurts by the rule above.
- **The harm is the root's line** [m]. The rule chose c23, which leaves three branching children (28:23, 29:23,
  30:23); production's r23 leaves one (23:31). Each of the three runs its own chain (3 to 4 branching nodes per
  depth against production's 1). CPU by depth 1-7: 3,685 / 2,181 / 3,122 / 2,227 / 1,705 / 1,284 / 1,089 s
  against 867 / 598 / 503 / 466 / 385 / 335 / 277.
- **The root is a weak node and the count is flat there** [m `runs/tai40b_root`, `lines.py`]: of the 80
  lines 72 have 40 entries above 0.001 and 7 have 39; c23 has 35 (sum of squares 0.058, largest entry 0.105) and r23 39
  (0.084, 0.179). The minimum of a flat count is noise. Rule 9's two candidates are r23 and c2 (the QPB
  arbitration takes r23).
- Against tai50b's root (49-50 everywhere, where the tie-break gave rule 9's line): here one line is lower by
  four entries, enough to win and not enough to mean anything.

### The gate: no effect on either tree

Hypothesis, set before the runs: the count is informative only once it is small. `QAP_TAIL_GATE=0.5`: the rule
where the fewest-entries count is at most m / 2, the params' rule (rule 9, QPB-arbitrated) elsewhere. The 0.5
was not fitted.

| trace lines | production | tail, ungated | tail, gated at m / 2 |
|---|---|---|---|
| tai35b | 1,017 | 884 (87%) | 1,045 (103%) |
| tai40b | 731 | at least 843 (115%) | 727 (99.5%) |

Both gated trees closed at the best known value [m `runs/tai35b_tailg`, `runs/tai40b_tailg`]. (`run.sh`'s last
lines failed in both because the script was edited while they ran; the solves had finished and the traces were
copied by hand, as `cmd.txt` says.)

- **tai40b: the gate removes the harm** [m]. The root is gated (rule 9's r23), and below it the rule chooses
  at depths 1-6 (16, 7, 5, 4, 3, 2 entries) and the tree stays a single chain: the same node count at every
  depth through 23. The rule's lines differ from production's from depth 2 (c33 against r22); on a chain that
  changes nothing.
- **tai35b: the gate removes the gain** [m]. The gated nodes are the three at depth 2 (17, 17 and 19 entries
  of 33); rule 9's c5 there, then the rule's tight choices below (3 and 4 entries of 31-32), give 7 branching
  nodes at depth 4 against the ungated tree's 2 and production's 5.
- **So the hypothesis as stated is wrong**: a 17-of-33 line (r33) helped and the tight choices below rule 9's
  c5 did not.
- **And tai35b's gain hung on one entry** [m `runs/tai35b_d2`: the three depth-2 nodes bounded again in one
  process, `lines.py`]. There r33 has 18 entries and rule 9's c5 19, 22 and 18 (at the third node c5 ties r33
  and wins the tie-break); the whole-tree run saw r33 at 17. Both lines stand far below the field (quartiles
  28 / 30 / 32), so the count separates the pair from the rest and not one from the other. Marginals differ
  by an entry or two between runs of the same node (the ADMM's warm start), so which of the two the rule takes
  is not stable. r33 was the better line; the count did not know it.
- At tai40b's root nothing stands out (c23 35 against 40); at tai35b's depth 2 two lines do (18 against 30)
  and the rule picks between them by noise; at tai60b's tight nodes one line has 2 to 5 against rule 9's 22.
  Only the last is a margin the count can carry [e].

## Not done

- A tree where rule 9's failures occur (tight nodes deep in a big tree) with the gated rule in its place:
  tai50b (rule 9 at every depth, 6,629 nodes) on the laptop, `../laptop/tai50b_tail.sh`, now with the gate. On
  tai35b and tai40b the gated rule had nothing to repair.
- A margin in place of the count's minimum (the rule only where its line is below rule 9's candidates by more
  than the run-to-run noise): suggested by tai35b's depth 2, not tested.
- The rule as a shortlist for the lookahead at chain nodes, against the rule alone.
- Grids (symmetry: rows only where symmetry remains; not implemented).
