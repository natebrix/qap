# tho40 bound lab: the grid cut phase per node (Nov 15 2026)

Agent `tho40_nodes` of `experiments/bound_lab/` (see `../COMMON.md`). Question: is the
regional-triangle cut phase (CP-ADMM form, `latest/dnncut.h`) valid on tho40, and how much
does it add to the certified DNN bound, node by node? Tools as built in the main checkout:
`experiments/dnn_cuts/dnncuts`, `latest/dnntest`. Nothing under `latest/` was changed.

Marks: **[m]** measured here (the log is named), **[e]** estimated. The machine was shared
(load about 6 on 4 cores throughout), so iterations are the metric; CPU seconds are indicative.
This file is assembled by `scripts/make_readme.py` from `README.tmpl.md`; the tables are the
scripts' output on the logs.

## Summary

@@SUMMARY@@

## Files

| | |
|---|---|
| `scripts/path.py` | the best known solution (`QAPSOLN/tho40.sln`, inverted for `tho40t.dat`), its four grid images (all 240,516), node specs along it |
| `scripts/part1_jobs.py`, `scripts/p1_summary.py` | part 1's runs and table (`logs/p1_*.log`) |
| `scripts/sample_nodes.py`, `nodes.tsv` | part 2's sample from `infra/phase1/runs/tho40ab-32/nodes.txt.gz` |
| `scripts/p2_summary.py`, `part2.json` | part 2's table (`logs/p2pl_*.log` plain, `logs/p2cp_*.log` CP phase, `logs/p2u_*.log` and `logs/p2long_*.log` longer CP runs) |
| `scripts/cost_table.py`, `scripts/tree_model.py`, `scripts/level_gain.py` | the derived tables of part 2 |
| `scripts/p3_summary.py` | part 3's tables (`logs/p3r_*.log` rounds form, `logs/p3cpb_*.log` CP form with boxes) |
| `scripts/runq.sh`, `scripts/runseq.sh`, `scripts/*.jobs` | the job runners and lists; each log's first line is its command (the lists were edited as the session went: the log is the record) |
| `logs/validate_cp.out`, `logs/validate_cp/` | brute-force check of the CP path on small instances |
| `cp_d4.log`, `rounds_d4.log` | the gate test of `../COMMON.md` (not mine) |

## Part 1: validity

**Brute force** [m] (`logs/validate_cp.out`): the repository's CP-path validator
(`experiments/wil50_feasibility/roadmap/admm/cpadmm/validate_cp.py --sep --count 1`, 47 s) on
41 nodes of n = 7, 8, 9 random instances (prefix regions) and nug12 with a 3 x 4 grid,
optima by enumeration: 205 CP solves, **10,508 certified values, max(certified - optimum) =
-1e-10**; 3,547 of them from the arms that separate inside the solve (`-cpsep`, the form used
below), same maximum; no fathom with the cutoff at the optimum, 80 correct fathoms at
optimum - 1. (The script ends with "FAILURES (502)": every one is a per-sweep progress line of
self-test (e), which the script's "must end in ok" rule predates; all 123 x 5 self-test
verdict lines (a)-(e) end in `ok`, and no certified-value or fathom failure is listed.)
`tests/dnn_test.py` (`make dnncheck`) covers the plain bound only, so it was not rerun.

**The path to the best known solution** [m] (`logs/p1_*.log`).
`scripts/path.py` inverts `tho40.sln` for `tho40t.dat`, checks that it and its three grid
images evaluate to 240,516, and takes the image with facility 26 at location 2 (path_check's).
A node that fixes d assignments of that solution has optimum <= 240,516 (= 240,516 if the best
known value is optimal), so no certified value may pass 240,516. Each run is
`dnncuts -x SPEC -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 -cpsep 100:1000 -cpiters N -T 1 -v qaplib/tho40t.dat`
(cutoff far above, so nothing stops early): 300 plain iterations cold, then N CP iterations,
certified every 10. Orders: `warm` = the path_check warm dive's order
(26, 30, 2, 16, 35, 24, 34, 21, 33, 37, 39, 15), `cold` = the cold dive's (differs from depth 3
on), `index` = facilities 0, 1, 2, ... (not a branching order: other nodes of the same solution).
N = 2,000 at depths 8 to 12, where the bound comes within reach of the optimum and an error
would show; at depths 1 to 6, where the bound is thousands below, N = 500 (`warm`) and 300
(`index`), cut down from the planned 1,000 for time.

{{p1_summary.py}}

@@P1TEXT@@

## Part 2: what the CP phase gains, by depth and shortfall

**Sample** (`scripts/sample_nodes.py`, `nodes.tsv`): from the production trace of subtree
`26:10,8:0` (22,650 nodes, 644 branched), branched nodes at quantiles of the give-up shortfall:
five at depth 3, five at depth 4 plus the gate node of `../COMMON.md`, four each at depths 5 and
6; the depth-2 node itself, the depth-1 nodes `26:10` (its parent) and `26:2` (the best known
solution's), and the root.

**Runs.** Plain [m]: `dnntest -p experiments/estimates/tho40/tho40.par -x SPEC -u 240516 -i N -H 0 -T 1 -v`
(`logs/p2pl_*.log`), N = 1,500 for the first three, then 1,000 at depths 4 to 6, 1,200 at
depth 3, 1,000 at depths 1 and 2, 600 at the root (cut for time; the "last-300 gain" column
says how converged each is). CP [m]: the brief's command with the cutoff at 241,500 instead
of `-u 240516` and `-cpiters 800` (`logs/p2cp_*.log`); the fathom iteration is the first check
whose certified value passes 240,514.05. That is the same run: at the gate node the
`-u 240516` log (`cp_d4.log`) and the high-cutoff log pass the cutoff at the same CP iteration
(120) with the same value to all digits (240,514.98184145204); checks further from a run's own
cutoff read about 10 lower in the high-cutoff run (a cheaper certificate there), so the
high-cutoff numbers are if anything conservative. `-u 240516` confirmations:
the gate node (`cp_d4.log`), the four early-entry runs `logs/p2e_*.log` (all `fathomed`,
240,514.45 to 240,521.34) and the long run `logs/p2long_d5_d.log` (not fathomed in 2,000).

Columns: shortfalls are 240,514.05 less the certified bound (negative = past the cutoff).
"production shortfall" is the trace's give-up value (200 iterations, warm-started from the
parent), with the node's subtree size in the production run. "plain@300" is cold, by dnntest
with tho40.par and by dnncuts' own first solve (its DNN defaults; the CP solve starts from that
state). "CP gain" = CP value at its last check less the plain end value; where the CP run
stopped at its 241,500 cutoff or was still rising, it is a lower bound.

{{p2_summary.py}}

**Reading it** [m].

- *The plain bound is spent at the give-up.* At depths 3 to 6 the converged plain shortfall is
  3 to 10% below the production give-up shortfall (d6_d 1,538 -> 1,495; d4_e 3,554 -> 3,417;
  d3_e 5,039 -> 4,667) and gains 1 to 17 units over its last 300 iterations. More plain
  iterations fathom none of the 18 sampled give-up nodes (d3_a comes within 13). Only the
  depth-2 node, whose production give-up came 200 iterations from a far start, moves much
  (8,634 -> 6,494).
- *The CP phase fathoms every sampled node up to a production shortfall of 1,906*: 13 of the 18
  nodes at depths 3 to 6, in 70 to 740 CP iterations after the 300-iteration first solve:
  70 to 130 iterations up to a shortfall of 407, 190 to 280 for 650 to 940, 650 to 740 for
  1,540 to 1,910. Of the subtree's 644 branched nodes, 84% are short by at most 1,500 and 90% by
  at most 1,900 (depth 3: 16% / 26%; depth 4: 72% / 84%; depth 5: 95% / 97%; depth 6: 99% / 100%).
- *The gain is a number of units, not a share of the gap.* Over the converged plain value, at
  the far-miss nodes of every depth, the certified CP value is ahead by

  | CP iterations | 100 | 200 | 400 | 800 |
  |---|---|---|---|---|
  | d3_d (plain end 3,383 short) | +405 | +882 | +1,379 | +1,878 |
  | d4_e (3,417) | +341 | +668 | +1,078 | +1,544 |
  | d5_d (2,288) | +323 | +727 | +1,126 | +1,554 |
  | d6_d (1,495) | +313 | +742 | +1,134 | +1,539 |

  and at CP iteration 800: +1,706 at depth 1 (against a plain value still rising 87 per 300),
  +1,840 at depth 2, +1,660 to +1,880 at depth 3, +1,540 to +1,620 at depth 4, +1,550 to +1,690
  at depth 5, +1,540 at depth 6; +1,850 at the root (against 226,515.4, the 1,500-iteration
  plain value of `experiments/tho40_estimate/root_xbar.log`). As a share of the converged plain
  shortfall that is 1 or more (a fathom) up to about 1,650, then 0.75 to 0.77 at 2,050 to 2,170, 0.68
  at 2,288, 0.56 and 0.45 at 3,400, 0.36 at 4,667, 0.28 at 6,494 (depth 2), 0.15 at 11,400
  (depth 1), 0.13 at the root: all at 800 iterations and all still rising by 150 to 230 per 200.
- *The plateau.* `logs/p2long_d5_d.log` [m] runs d5_d (2,359 short in production, a
  371-node subtree) for 2,000 CP iterations with `-u 240516`: shortfall 734 at 800, 605 at
  1,000, 510 at 1,200, 433 at 1,400, 374 at 1,600, 324 at 1,800, **283 at 2,000, not
  fathomed**. That is +2,005 over the converged plain value, the increments shrinking by 18%
  per 200 iterations, so the plateau is about +2,150 to +2,200 [e]: this node sits just past
  what the triangles can fathom, and the practical reach of the phase at depths 4 to 6 is a
  shortfall of about 1,600 (800 CP iterations) to 2,000 (2,000 iterations). Shallower nodes
  have more room: the root is +2,797 at 2,000 CP iterations and +3,573 at 5,490 (below).
- *Where the CP solve starts.* It restarts below the first solve (the gate node reads 234,639
  at CP iteration 10 after 239,728 at the give-up) and is back at the give-up value by CP
  iteration 50 to 70, so those iterations re-converge rather than gain. Entering the
  CP phase earlier saves about 100 iterations in total [m, `logs/p2e_*.log`, `-u 240516`, all
  four end `fathomed`]: the gate node fathoms after 100 + 218 = 318 iterations (`-I 100`) or
  50 + 277 = 327 (`-I 50`) instead of 300 + 120 = 420; d6_c after 50 + 320 = 370 instead of
  490; d4_c after 50 + 830 = 880 instead of 980.

### What a fathom at the give-up node saves

{{cost_table.py}}

Iterations as the unit [m for the columns, e for the conversion: the CPU below a node also
holds its children's QPB and branching work]. At the give-up nodes whose children all fathomed
at once in production (subtree = 35 to 38 nodes: the class that is about 40% of tho40's CPU),
the CP fathom costs 2.2 to 6.4 times fewer iterations than the children it replaces; the same
at the nodes with larger subtrees (3.3 to 5.2). The saving is largest for the near misses
(100 CP iterations against 500 to 600) and smallest near a shortfall of 900 to 950.

### A threshold model of the subtree [e]

{{tree_model.py 1.2}}

This counts only what part 2 measured: a CP phase after the give-up, on nodes within the
sampled range, each CP iteration at 1.2 plain iterations (dnncuts' CPU per CP iteration is 1.19
times its plain one at every m, `logs/p2cp_*.log`: 1.12 to 1.30). With the threshold at 1,500 to
1,900 the subtree drops to **10 to 16% of its nodes** and to **57 to 62% of its CPU**: the
node count falls 6 to 10 fold, the CPU 1.6 to 1.75 fold. The CPU falls less because what is left
is the give-up nodes themselves: each still pays its 200-iteration first solve (the production
run's branched nodes cost 37 s each, 30% of its CPU, against 2.5 s for a fathomed one) plus its
CP iterations, which at a shortfall of 1,500 (740) are most of the children's cost. Not in the
model, and where the larger factor would have to come from: children that inherit the parent's
pool and start from its CP state (the cut bound from their first iteration), the CP solve
entered before the plain solve has run its 200 to 300 iterations, and fathoms above a shortfall
of 1,900 with longer CP solves. Those are the in-search questions (agent `tho40_search`).

### How much of a level is it worth?

{{level_gain.py}}

One level of branching raises the bound by about 6,300 from depth 2, 2,700 to 3,400 from
depth 3, 1,600 to 2,500 from depth 4 and 1,100 to 2,300 from depth 5 (the lower figures are
medians over all children, where a fathomed child's gain is cut off at the cutoff; the higher
ones are for far-miss parents, where less is cut off). The CP phase's +1,550 to +1,880 at 800
iterations is therefore **0.25 to 0.3 of a level at depth 2, about 0.5 to 0.65 at depth 3, 0.6 to
0.95 at depth 4, 0.7 to 1.4 at depth 5** [e], and its long-run gain (about +2,000 to +2,200, above)
a quarter more. In the "gap halves per level" picture: the median give-up node at depth 3 is
2,427 short; after 800 CP iterations it is 500 to 900 short, which is the median give-up
shortfall at depth 4 (930), so one level; from depth 4 down the median give-up node (930, 364,
342 short) is not moved a level but fathomed, in 280, 130 and 120 CP iterations. The gain does
not shrink with depth while the gap does, so the cut phase is worth more levels the deeper the
node, and it is worth most exactly where tho40's CPU is (depths 4 to 7).

### The root bound with cuts, run long

`logs/p2cp_root.log` [m]: `dnncuts -x - -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 -cpsep 100:1000 -cpiters 20000 -T 1 -v qaplib/tho40t.dat`
under `timeout 2400` (40 minutes of wall on the shared machine). Plain at 300 iterations:
225,681.44. The CP solve passes the plain DNN value (226,522) at CP iteration 120 and reads

| CP iteration | 500 | 1,000 | 2,000 | 3,000 | 4,000 | 5,000 | 5,490 (killed by the timeout) |
|---|---|---|---|---|---|---|---|
| certified | 227,870.2 | 228,603.7 | 229,312.3 | 229,668.5 | 229,884.6 | 230,031.7 | **230,088.31** |

3,568 triangles in the pool at the end (2,568 added by 54 in-solve separations), 1,260 CPU s
for the CP solve (0.23 s per iteration at m = 40), 5,790 iterations in all. Still rising: +57
over the last 490 iterations, the increments shrinking by about a third per 1,000, which
extrapolates to a plateau of 230,300 to 230,400 [e]. **230,088.31 is a certified lower bound for tho40**
(230,090 after rounding up to the next even value, the objective being even): 3,566 above the
DNN value of about 226,522 that the brief gives as the best published bound, and it takes the
root gap from 5.82% to 4.34% of 240,516. Plain DNN at the root for comparison [m]: 224,074.5 at
300 and 226,280.0 at 600 iterations (`logs/p2pl_root.log`, tho40.par), 226,515.4 at 1,500
(`experiments/tho40_estimate/root_xbar.log`). The claim rests on the validity of the regional
triangles (a counting argument on the grid's row and column prefix regions, `latest/dnncut.h`),
on the certificate (every value is certified from a dual point; 549 checks in this run), and on
part 1: I did not verify the published record independently.

## Part 3: the other families

@@P3INTRO@@

{{p3_summary.py}}

@@P3TEXT@@

## What it implies, and what I would do next

@@NEXT@@
