# Child warm start from the parent's Y as well as Z (Nov 3, cloud)

Idea 1 of `../roadmap/admm/README.md` (from `../roadmap/admm/splitting.md`).
Today a child starts its ADMM from the parent's final multiplier Z mapped to
the child (`dnn_warm_child`), with Y = e0 e0'. With Z near the parent's
multiplier the first R-step sees W = V'ZV/beta, nearly negative semidefinite,
returns R near 0, and the first few dozen iterations rebuild the positive
part the parent already had. The test also carries the parent's final Y.

Code: `dnn_warm_child(..., ymode)` in `latest/dnn.C` (0 = Z only, the default;
1 = the parent's Y restricted to the child's indices; 2 = that Y conditioned
on the child's new assignment x_ij = 1, row 0 from the parent's row (ij)
over x_ij and the off-diagonals rescaled so diag(Y) = row 0), used with
`DNNParams::warm_primal = 1`; `dnntest -K i:j -Q z,a,b` runs a child from each
start. The certificate is untouched: it comes from the R-step multiplier as
always. An earlier note in `dnn.C` said a conditioned Y "made the child
slower"; that was one nug30 child in the first port of the solver (Oct 28,
`../../dnn_bound/cpp_report.md`: "a hint, not a measurement").

`run.sh`: the 16 depth-4 children of `../stronger_bounds/depth4/nodes.txt`,
each bounded after its depth-3 parent (nodes 10, 13, 2, 7, all clear misses,
given up at 300 iterations as in the search), production settings (bundle D
of `../amx/README.md`: AMX, degree 3 with tight lambda, one-pass CholQR,
refresh 40, PRSM (1, 1), give-up after 300, at most 1,500 iterations), four
processes on the four cores. Per-child table in `summary.txt`, logs in `runs/`.

## Result

| start | outcome changes | fathoms (7 shared): iterations | CPU | misses (9 shared): CPU | misses' mean lb at 50 / 100 / 200 |
|---|---|---|---|---|---|
| Z only (today) | - | 1,800 | 604 s | 1,290 s | 48,498.0 / 48,604.5 / 48,672.0 |
| Z + Y restricted | none | 1,360 (-24%) | 465 s (1.30x) | 1,189 s | 48,610.1 / 48,664.0 / 48,700.1 |
| **Z + Y conditioned** | **none** | **1,360 (-24%)** | **467 s (1.29x)** | **1,145 s (1.13x)** | **48,606.9 / 48,671.3 / 48,705.6** |

- Every fathom comes earlier or at the same iteration with either Y start
  (node 8: 160 -> 90, node 5: 490 -> 340, node 3: 410 -> 340, node 13:
  530 -> 420), and no child changes between fathom and miss.
- On the misses the bound is +67 higher at iteration 100 and +34 at 200
  with the conditioned Y; the give-up fires sooner on most (1,145 s against
  1,290 s), one miss (node 9) runs 40 iterations longer.
- The kill rule (under 10% fewer iterations to fathom AND under +10 at
  iteration 100) is cleared on both counts. The conditioned variant is the
  better on the misses (node 10: 48,738 against 48,681 at iteration 100) and
  equal on the fathoms; it is the one to carry.

Over all 16 children the CPU is 1,612 s against 1,894 s (1.17x); on the
fathoms 1.29x. The earlier one-node nug30 result does not reproduce.

## What it would take in the search

The search keeps each missed node's `DNNWarm` for its children; `dnn_bound`
already fills `save->Y`, so the parent's Y is in memory today. Integration is
`bnb.C` calling `dnn_warm_child(..., 2)` with `dp.warm_primal = 1` behind an
environment switch, plus a check of the warm-start memory (Y doubles each
stored warm start: 39 MB more per missed node at m = 47) and of the state
files (`latest/state.C`), then a node-10 subtree A/B against the rule-9 arm.
