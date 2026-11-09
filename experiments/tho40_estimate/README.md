# tho40: how big is the tree? (November 8 2026)

tho40 (Thonemann and Bölte, 5 x 8 grid, QAPLIB best known 240,516, never
proven optimal) sized with the production solver, by tracing every node the
search bounds and fitting subtree size to the bound gap on the subtrees it
finished. Fujii, Ito, Kim, Kojima, Shinano and Toh (arXiv:2101.09629) estimate
3.7e7 nodes for tho40 with their rule P. The estimate here, after the depth-2
sweep below, is **6e5 to 3e6 nodes**, about 1e6 to 4e6 CPU s on this laptop's
cores: 1.3 to 6 days at 8 workers, 15 to 70 minutes on 1,000 cores. tho40 is
solvable with the resources the wil50 plan assumes, with a wide margin.

## Files

| | |
|---|---|
| `qaplib/tho40t.dat` | `tho40.dat` with the matrices swapped (flows first, the grid second), since symmetry acts on the second matrix |
| `latest/SYM/tho40t.sym` | the 5 x 8 grid's symmetry: 12 root orbits (Burnside: (40 + 8)/4), the middle row fixed by one reflection, the top three rows the other half-domain; `bnbdrv -d 1` gives 13 nodes |
| `root_xbar.log` | the root DNN bound, 1,500 iterations: 226,515.4, **gap 5.8%** (tho30 4.2%, nug30 2.8%); 95 ms per iteration at m = 40 |
| `traces/{nug30,tho30,tho40}.trace.gz` | `QAP_DNN_TRACE` output: one line per DNN-bounded node (pid, depth, fixed pairs, first-solve bound, cut-phase bound, cutoff, fathomed) |
| `{nug30,tho30,tho40}_rule9_j8_m2pro.log` | the three searches, the nug30 record recipe with rule 9 (`nug_deep_xbar9_all.par`), 8 workers |
| `estimate.py` | rebuilds the tree from the fixes, sums subtree sizes bottom-up, fits log10(size) to the relative gap per depth on complete subtrees, applies the fit to open nodes |

## Method

The record recipe with rule 9 is the production search, so a run of it on tho40
builds the real tree in the real order; the trace records every node's bound.
The run was stopped after 5 h (`timeout` did not stop the coordinator; a
SIGTERM did), with 29,254 nodes bounded, 5,440 open, 40,700 CPU s. The tree is
rebuilt from the fixes, each node's subtree size summed where the subtree is
finished, and for each depth d a line log10(size) = a_d + b_d (cutoff - lb)/cutoff
fitted on the finished, unfathomed subtrees (fathomed nodes are size 1). An
open node's subtree is the fit at its own gap and depth; a node with some
children open is 1 plus its children's estimates, missing children sharing the
parent's fitted size. The 3 root orbits the search had not reached count at
the mean of the 9 it had.

Cross-validation on the two solved instances, truncating the finished traces
at depth 2, 3 or 4 and estimating the rest:

| fit on | predict | actual | estimate (truncated at depth 2 / 3 / 4) |
|---|---|---|---|
| tho30 | nug30 | 3,705 | 2,900 / 3,800 / 4,200 |
| nug30 | tho30 | 14,822 | 38,000 / 40,000 / 41,000 |
| tho30 | tho30 | 14,822 | 15,600 / 18,900 / 20,800 |

So a model carried from a neighbouring instance is good to a factor of 3, and
an instance calibrated on its own finished subtrees to 1.5x. The slopes agree
across all three instances: 0.5 to 0.8 in log10(size) per 1% of relative gap
at depths 2-6 (tho40 pooled: log10 size = 1.42 + 0.69 gap% + 0.01 depth over
689 finished subtrees).

## The tho40 trace

| depth | bounded | fathomed | finished subtrees | open nodes used by the model, their gaps |
|---|---|---|---|---|
| 1 | 9 of 12 | 0 | 0 | 9, 4.8-6.4% |
| 2 | 32 of 468 | 0 | 3 | 29, 1.35-3.9% |
| 3 | 332 | 205 | 256 | 76, 0.07-3.1% |
| 4 | 3,239 | 2,913 | 3,146 | 93, up to 1.8% |
| 5 | 10,183 | 9,848 | 10,117 | 66, up to 1.3% |
| 6 | 10,627 | 10,503 | 10,610 | 17 |
| 7 | 3,936 | 3,906 | 3,932 | 4 |
| 8 | 896 | 896 | 896 | 0 |

Depth-1 bounds 225,070 to 228,881 (gaps 4.8-6.4%); the 32 depth-2 bounds
231,110 to 239,618 (median 233,837, gap 2.8%). The search is depth-first, so
after 5 hours it had been deep in 9 of the 12 root orbits and had opened 32 of
the 468 depth-2 nodes; 96.6% of bounded nodes fathom, like tho30 (96.1%) and
nug30 (96.3%), at a mean 1.39 CPU s per node (tho30 0.32 s).

## Estimates

| calibration | nodes |
|---|---|
| tho40's own finished subtrees | 4.0e6 |
| tho30's tree | 1.6e6 |
| nug30's tree (known to over-predict tho30 by 2.6x) | 1.8e8 |

Range: **1.6e6 to 4.0e6 nodes**, 2.2e6 to 5.6e6 CPU s at 1.39 s per node.
On this laptop at 8 workers fully used, 3.5 to 9 days; on 1,000 cores at
tho40's node cost, 40 minutes to 1.6 hours of wall. Fujii et al.'s 3.7e7 (rule P)
is 10x to 25x above this; the same ratio tho30 shows between this solver and
the older bound/branching, and in line with rule 9 plus the warm-started
bound. With a 3x model error either way, the 95% picture is 5e5 to 1e7 nodes:
tho40 is solvable, and not close to the margin.

Caveats, in order of weight:

- The top of the model is thin: the depth-2 line is fitted on 3 finished
  subtrees (gaps 0.4-0.7%) and applied to 29 open depth-2 nodes at gaps of
  1.3-3.9%, and 3 of 12 root orbits were never reached. The tho30-calibrated
  number does not have this problem (tho30's depth-2 fit covers 0-2%) and
  agrees within 2.5x. The fix is a cold sweep of all 468 depth-2 nodes with
  `dnntest` at 1,500 iterations (about a minute each, 1.1 h at 8 jobs), which
  removes the extrapolation at depth 2 entirely; not run.
- The run used about a quarter of the machine: `bench/parallel.py` with
  `--checkpoint-every 300` stops launching tasks while it waits for the running
  ones to finish before each checkpoint, and tho40 tasks are long, so the
  running count drained 7, 6, 5, ..., 0 before every checkpoint (40,700 CPU s in
  17,300 s of wall on 8 workers). A real tho40 run needs a long checkpoint
  interval or a checkpoint that does not wait. `timeout` sent SIGTERM to the
  coordinator at 3 h 40 min and it kept launching tasks; the TaskStop SIGTERM
  3 hours later stopped it with a checkpoint (3.0 GB, in the session's
  scratchpad, not kept).
- Nothing is tuned for n = 40: the DNN depths 1-7, the 200-iteration give-up,
  beta 0.21 m/3 and rule 9 through depth 6 are nug30's. The root line was the
  production rule's pick (facility 26 across the 12 J1 cells); the root-line
  studies in `../wil50_feasibility/tuning/lines/` suggest a shortlist-by-child-
  bounds choice could be worth 2-5x more.

## The depth-2 sweep (later on Nov 8): the extrapolation at the top removed

The production solver started at each of the 12 depth-1 nodes (`bnbdrv -F 26:c
-d 2`, the record environment, rule 9), so every depth-2 node gets the bound
the search would give it: the depth-1 node solved cold as in production, its
children warm from its shipped state, the line chosen by rule 9 at the node.
12 runs, 8 at a time, 17 min of wall (`depth2/d2_26_<c>.{trace,log}`,
`depth2/one.sh`; `depth2/tho40_depth2.trace` is the concatenation). On the 32
depth-2 nodes the 5 h search had also bounded, the sweep's bounds agree to
within 0.2 units, so the two traces are the same kind of number and are read
together (`estimate.py` takes `a.trace+b.trace`; `EXHAUSTIVE=1` tells it that
depth 1's children lists are complete).

404 depth-2 nodes, not 468: the four parents on the grid's middle row
(cells 16-19) keep one reflection, and the solver's symmetry reduction gives
them 23 children each, the other eight 39. Lines: rule 9 took cell 7 under
26:0, cell 0 under the other seven off-axis parents, and facility 15 or 5
across the cells under the middle-row parents. One depth-2 node fathoms.

| parent 26:c | children | median gap | max gap |
|---|---|---|---|
| c = 0 | 39 | 1.89% | 3.13% |
| 1, 8 | 39 | 2.33-2.53% | 3.55-3.67% |
| 2, 3, 9, 10, 11 | 39 | 2.86-3.24% | 4.02-4.38% |
| 16, 17, 18, 19 (middle row) | 23 | 3.13-3.68% | 4.03-4.63% |

Across all 404: median gap 2.87%, 10th-90th percentile 1.39-3.92%, max 4.63%.

With the sweep folded in, nothing above depth 2 is extrapolated, and the three
calibrations move together:

| calibration | estimate (search trace + sweep) | from the 404 depth-2 bounds alone |
|---|---|---|
| tho40's own finished subtrees | 6.1e5 | 5.4e5 |
| tho30's tree | 1.7e6 | 1.7e6 |
| nug30's tree | 2.9e6 | 2.9e6 |

The spread that remains is the slope of log10(size) against gap at depth 2
applied at gaps of 3-4%, above any finished depth-2 subtree on any of the
three instances (tho40's own depth-2 fit rests on 3 nodes at 0.4-0.7%; tho30's
covers 0-2%). The nug30 calibration no longer sits 100x above the others, as it
did when it had to extrapolate from the depth-1 gaps of 5-6%: that gap had
been the weakest step, and the sweep removed it. Range: **6e5 to 3e6 nodes**,
**0.9e6 to 4e6 CPU s** at 1.39 s per node; 1.3 to 6 days on this laptop with 8
workers fully used, 15 to 70 minutes on 1,000 cores. Fujii et al.'s 3.7e7 is
12x to 60x above.

What would tighten it further: finishing a handful of depth-2 subtrees at
gaps of 3-4% to pin the slope where it is applied. The sweep shows a typical
such node costs 10^3.5 to 10^4.3 nodes by the three fits, 1 to 7 CPU-hours
each, so three or four of them at 8 workers is an evening; the depth-first
search did this for the small ones only.

## Running it again

```sh
# the nug30 record environment (../nug30_solve/README.md), then
QAP_DNN_TRACE=$PWD/experiments/tho40_estimate/traces/tho40.trace \
python3 bench/parallel.py qaplib/tho40t.dat --jobs 8 --ub 240516 \
    --params PARAM/nug_deep_xbar9_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint tho40.ckpt --checkpoint-every 3600
# the depth-2 sweep: one depth-limited production run per depth-1 node
printf '%s\n' 0 1 2 3 8 9 10 11 16 17 18 19 | xargs -P 8 -n 1 experiments/tho40_estimate/depth2/one.sh
cd experiments/tho40_estimate; T=traces/tho40.trace+depth2/tho40_depth2.trace
EXHAUSTIVE=1 python3 -I estimate.py predict $T 40 12 $T 40 12          # self-calibrated
EXHAUSTIVE=1 python3 -I estimate.py predict traces/tho30.trace 30 10 $T 40 12   # tho30-calibrated
```
