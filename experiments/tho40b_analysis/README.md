# tho40-b analysis scripts (Nov 15)

The scripts the plan in `docs/tho40b_results_plan.md` calls for, written while
tho40-b runs and tried on the two finished subtree runs of the shipping A/B
(`infra/phase1/runs/tho40ab-32/` and `tho40ab-0/`: the depth-2 subtree
`-F 26:10,8:0` of `tho40t.dat`, about 22,500 nodes each, with the incumbent at the
best known, 240,516). They read the node trace (`QAP_NODE_TRACE`, format in
`latest/bnb.C`) and, for the give-ups, the DNN trace, as `qapctl export` writes
them (`coord/export.py`).

| script | question | input |
|---|---|---|
| `ntrace.py` | the trace as a tree: parents from the fixes, subtree nodes, CPU and depth, which subtrees are complete | |
| `bound_value.py` | how many nodes and how much CPU a bound d units better would save, by the depth of the new fathom; the near misses worth freezing for a bound lab | node trace |
| `estimator_score.py` | how well a run part-way predicts its own final tree: the first fraction of the trace (merge order, what the live tree showed then) through two estimators, scored against the real total | node trace |
| `tuning_set.py` | a tuning set: complete subtrees about one tuning run long on one instance, spread over the gaps, none inside another, as `-F` specs | node trace |
| `giveups.py` | what the DNN give-up costs and saves: give-ups per depth, and those whose children all fathomed at once | an export directory |

## Results on the A/B subtrees

**A better bound.** `python3 bound_value.py ../../infra/phase1/runs/tho40ab-32/nodes.txt.gz`
(shipped arm; the unshipped arm is within 1 point up to d = 1,000 and 7 points at 2,000):

| d (units) | d (% of cutoff) | nodes saved | CPU saved | the same, better bound only at depths <= 4 |
|---|---|---|---|---|
| 20 | 0.008% | 1.2% | 0.5% | 0% |
| 100 | 0.042% | 10.4% | 5.0% | 2.2% nodes, 1.3% CPU |
| 200 | 0.083% | 23.0% | 11.5% | 4.3%, 2.6% |
| 500 | 0.208% | 49.6% | 26.7% | 11.0%, 7.0% |
| 1,000 | 0.416% | 71.7% | 43.6% | 20.3%, 14.6% |
| 2,000 | 0.832% | 90.5% | 69.1% | 48.2%, 41.1% |

So a bound about 500 units (0.2%) better everywhere halves this subtree's nodes, but
saves only a quarter of its CPU: the nodes saved are the cheap deep ones, and
the expensive shallow nodes remain. A better bound used only down to depth 4
needs about 2,000 units to halve the tree. The bound recorded for a branched
node is what the recipe reached within its give-up (200 iterations), so part of
any gain may come from a later give-up rather than a new bound; `giveups.py` on
tho40-b's export separates the two.

**Estimating part-way.** `python3 estimator_score.py ../../infra/phase1/runs/tho40ab-32/nodes.txt.gz`,
the predicted total over the actual:

| prefix of the run | tho40ab-32 gap fit | rates | tho40ab-0 gap fit | rates |
|---|---|---|---|---|
| 5% | 0.06 | 0.18 | 0.07 | 0.20 |
| 10% | 0.18 | 0.26 | 0.16 | 0.36 |
| 20% | 0.28 | 0.34 | 0.36 | 0.40 |
| 30% | 0.59 | 0.42 | 0.48 | 0.46 |
| 50% | 0.73 | 0.59 | 0.71 | 0.64 |
| 75% | 0.82 | 0.79 | 0.84 | 0.79 |

Both estimators read low until the run is about half done, by 2 to 5 times
early on: the subtrees that finish first are the small ones, and the early
fathom rates are the easy nodes'. The Nov 8 estimate (6e5 to 3e6 nodes) was made
the same way, from a 5-hour prefix of a depth-first search, and a part-way
projection of tho40-b is low for the same reason: the per-depth rates projection
from its own published tables rose from 3.8e6 nodes at 4.9 h of wall (Nov 14,
21:44 UTC) to 6.3e6 at 13.3 h (Nov 15, 06:11 UTC). With these ratios it is more
likely 1e7 or more. At 3.7 CPU s per node and about 62 busy slots, 1e7 nodes
is about a week of wall: tho40-b will reach its 72-hour `--limit` and need a
resume launch (`infra/phase1/README.md`, Interrupted or stopped runs).

**A tuning set.** `python3 tuning_set.py ../../infra/phase1/runs/tho40ab-32/nodes.txt.gz --minutes 2 --n 6`
picks six depth-3 and depth-4 subtrees of 900 to 1,700 nodes (1 to 2 min on 64
slots). The A/B subtree itself is one 30-minute tuning run; in tho40-b's whole
tree the 30-minute subtrees are the larger depth-2 nodes.

**Give-ups.** No A/B export has DNN traces (they were not collected), so
`giveups.py` was checked on a small traced nug17 run only; it runs on
tho40-b's export.

## On tho40-b, when it ends

```
infra/phase1/collect.sh tho40-b          # or the phase1 workflow's collect: artifact traces-tho40-b
cd experiments/tho40b_analysis
python3 bound_value.py ../../infra/phase1/exports/tho40-b/nodes.txt.gz
python3 bound_value.py ../../infra/phase1/exports/tho40-b/nodes.txt.gz --max-depth 4
python3 estimator_score.py ../../infra/phase1/exports/tho40-b/nodes.txt.gz
python3 tuning_set.py ../../infra/phase1/exports/tho40-b/nodes.txt.gz --minutes 30 --n 8
python3 giveups.py ../../infra/phase1/exports/tho40-b
```
