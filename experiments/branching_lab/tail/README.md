# The rule "tail": branch on the line with the fewest non-zero marginals

Nov 21 2026. Status: passes its one-level tests on tai60b and its first whole tree (tai35b); tai40b running.
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
calculation from the dumped x-bar (`QAP_LINE_SCORES`) [m].

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

tai40b: running (`runs/tai40b_tail`). Production's tree there is a single chain, one branching node per depth
(731 nodes), so it can show harm but little gain [e].

## Not done

- A tree where the rule replaces rule 9 deep (tai50b: rule 9 at every depth, 6,629 nodes). That is the test of
  the failures' cost; it needs the laptop.
- The rule as a shortlist for the lookahead at chain nodes, against the rule alone.
- Grids (symmetry: rows only where symmetry remains; not implemented).
