# The rule "tail": branch on the line with the fewest non-zero marginals

Nov 21 2026. Status: ungated, it helps on tai35b (13% fewer nodes) and hurts on tai40b (at least 15% more: a
weak root). A gated version (the rule only where the count is small) is running on both.
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

### The gate (running)

Hypothesis: the count is informative only once it is small. Where the rule worked the chosen line had 2 to 5
entries of about 50 (tai60b) or 9 of 35 (tai35b's root); where it failed, 35 of 40. `QAP_TAIL_GATE=0.5`: the
rule where the fewest-entries count is at most m / 2, the params' rule (rule 9, QPB-arbitrated) elsewhere. The
0.5 was set before any gated run and is not fitted: in tai35b's ungated tree four choices had 17 to 29 entries
of 32-33 and would be gated, and whether those helped is not known [e].

`runs/tai35b_tailg`, `runs/tai40b_tailg` (`QAP_TAIL=0.001 QAP_TAIL_GATE=0.5`): running, the same references and
the same 95% / 105% rule. On tai40b the root is gated (rule 9's r23), so the question there is whether the
rule keeps the single chain below it; on tai35b, whether the 13% survives the gate.

## Not done

- A tree where the rule replaces rule 9 deep (tai50b: rule 9 at every depth, 6,629 nodes). That is the test of
  the failures' cost; it needs the laptop (`../laptop/tai50b_tail.sh`, on hold until the gated runs are in).
- The rule as a shortlist for the lookahead at chain nodes, against the rule alone.
- Grids (symmetry: rows only where symmetry remains; not implemented).
