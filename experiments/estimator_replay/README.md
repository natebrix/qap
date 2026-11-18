# Estimator designs replayed on finished trees (Nov 18 2026)

Step 1 of `docs/estimates_plan.md`. Pick the estimator's design offline, on trees whose totals are known,
before any of it runs on AWS. `replay.py` takes each design's dives on a finished run's tree instead of
solving nodes. It scores each estimate against the true totals, at the CPU a real run of the design would
spend: the recorded CPU of every node the design would bound, each node paid once.

    python3 -I replay.py TRACE --k 3 --budget 0.1,0.2,0.3 --reps 20 --designs y2000:1.5,survq:1,split-short-s90:1

TRACE is a node trace, or a directory of the website's tree shards. tai60b-261117's shards were
fetched with `fetch_tree.py tai60b-261117 DIR`; the collected trace,
`infra/phase1/runs/tai60b-261117/nodes.txt.gz` (committed Nov 18), reads the same: 331,149 nodes,
3,596.7 CPU h, no orphans.
`sanity.py TREE DIVES` checks that the dives are unbiased (from the root); `predict.py TREE` measures
how well the shortfall ranks a node's surviving children by their subtrees.

## Designs

Every design starts from the exact top: every node at depth < k is counted as it is, and every
depth-k node is bounded once.

- **y2000.** Today's `estimate.C` and the thesis's ch. 4. A dive starts at a depth-k node drawn by g^q
  and takes one child at a time, drawn by g^q. g is the relative gap of the child's own QPB bound
  (`QAP_EST_CHILD_QPB=1`), with 0.1 uniform mixed in.
- **strat.** The same dives, but every depth-k node that branches is its own stratum. Dives are
  allocated by Neyman at fixed cost, after a pilot of 2 per stratum.
- **surv, survq (survivor dives).** Stratified, and at every node of a dive all its children are
  bounded, as the search bounds them anyway:
  - the fathomed children are counted exactly;
  - the dive continues into a child that branched, drawn uniformly (surv) or by shortfall^q (survq).
    The shortfall is the cutoff minus the child's best bound.

  This is unbiased (Purdom's partial backtracking).
- **split-BY-sS.** Adaptive stratification with survivor dives. Until S% of the budget is spent, the
  stratum with the largest key is replaced by its children, which is one more exact breadth-first
  level under it. The key is one of:
  - `short`: the largest shortfall;
  - `bfs`: the shallowest, as a solve's frontier grows;
  - `var` or `mean`: from the pilot dives.

  The rest of the budget goes to Neyman dives on the strata left.

## Results

**tai60b-261117** (331,149 nodes, 3,597 CPU h, deepest 57). k = 3, 20 replays each
(`out/tai60b_grid.md`, `out/tai60b_split.md`). The columns are the estimate over the truth, for CPU:

| design | budget (share of the tree's CPU) | median | 10-90% | within 20% | within 2x | deepest level a dive reached |
|---|---|---|---|---|---|---|
| y2000, q = 1.5 | 20,000 dives (2.6%) | 0.09 | 0.06-0.12 | 0% | 0% | 7 |
| strat, q = 1.5 | 20,000 dives (2.9%) | 0.07 | 0.06-0.13 | 0% | 0% | 7 |
| surv | 10% | 0.32 | 0.23-0.69 | 5% | 15% | 12 |
| survq, q = 1 | 10% | 0.44 | 0.28-0.78 | 5% | 40% | 14 |
| split-short-s90, q = 1 | 10% | 0.68 | 0.37-1.92 | 5% | 60% | 18 |
| survq, q = 2 | 20% | 0.85 | 0.60-1.47 | 35% | 100% | 22 |
| split-bfs-s80, q = 1 | 20% | 0.87 | 0.75-1.42 | 60% | 95% | 18 |
| survq, q = 1 | 30% | 0.87 | 0.74-1.31 | 60% | 100% | 22 |
| split-short-s90, q = 1 | 30% | 0.90 | 0.79-1.03 | 75% | 100% | 22 |

- **Today's design cannot see this tree.** With 20,000 dives it gives 5-9% of the true CPU. The bound
  lab's estimate was 7.8% (279 of 3,597 h). The dives end within one level: 93-99% of a depth's nodes
  are fathomed, and the QPB gaps that weight the children barely differ. More dives do not help:
  20,000 cost 2.6% of the tree and the estimate does not move.
- **Survivor dives reach the bottom and are unbiased, but tai60b's tree is concentrated.** Two chains
  carry the cost (`../tai60b_postmortem/`). The sampling distribution is heavy-tailed: from the root,
  20,000 survivor dives have mean 0.43 to 0.72 and median 0.005 to 0.04, so the median estimate is
  low until a budget of about 20% of the tree.
- **Expanding the high-shortfall strata exactly helps**: split-short-s90 at 30% is within 20% three
  times in four. A partial solve's breadth-first frontier with survivor dives on its open nodes
  (split-bfs) does about as well at 20%.
- **The shortfall predicts survivors' subtrees weakly on tai60b**: Spearman 0.29; the largest
  shortfall has the largest subtree at 52% of nodes. On tai50b it is 0.87 and 82%.

**tai50b-261116** (8,216 nodes, 133 CPU h; `out/tai50b_grid.md`). k = 1 to 2. y2000 gives 0.26-0.51 of
the CPU at 10-30%. Survivor dives give 0.63-0.81 with 86-98% within 2x, but cost 42-55% of the tree:
the pilot dives' children are most of a tree this small.

**tai35b, tai40b** (the laptop's, about 2 CPU h each; `out/tai35b_grid.md`, `out/tai40b_grid.md`). Too
small to test a budget: survivor dives cover most of the tree. y2000 gives 0.23-0.44 of the CPU.

## What this means

1. **Today's estimator is not fit for these trees, at any number of dives.** That holds for the
   thesis's design (g^q from QPB) wherever the DNN bound fathoms nearly every child.
2. **A trustworthy estimate of a tai60b-like tree costs 20-30% of the tree** with these designs. That
   is about the cost of solving that share of the tree, and a partial solve does much of the same
   work. So the estimate run should **be the start of the solve**:
   - the coordinator solves breadth-first for a fixed wall-clock budget, and survivor dives estimate
     the open nodes' subtrees;
   - continuing loses only the dives' own extra nodes.
3. **The weights among the survivors are the lever** (next section). With the true subtree sizes as
   weights (the oracle), 10% of the tree is within 20% 95% of the time on tai60b. The best weights a
   node's own bound gives need 30% for 90%.
4. **Not tried yet:** the dives' own nodes kept for the solve. In the replay every node a survivor
   dive bounds is a node of the tree, so none of the estimate's CPU is outside the solve's own work.

## Round 2: weights, depth, another tree (Nov 18)

New in `replay.py`:
- survivor weights `xB`: exp(B x), with x the shortfall over the median shortfall of the branched
  nodes at the child's depth;
- `o`: the oracle, each survivor's true subtree CPU, the zero-variance ideal;
- `--bydepth`: the node estimate by depth band;
- subtree traces (k counts from the subtree's root).

**What predicts a survivor's subtree** (`size_model.py`, `size_features.py`). This is a within-sibling
regression of log subtree CPU on what a node records once it is bounded:

| tree | survivors | shortfall / depth's median: slope, R^2 | + own CPU, + QPB shortfall: R^2 |
|---|---|---|---|
| tai60b-261117 | 5,463 | 0.40, 0.19 | 0.24 |
| tai50b-261116 | 109 | 0.83, 0.63 | 0.63 |
| tho40ab-0 (subtree) | 560 | 0.36, 0.32 | 0.32 |
| nug30-c7a | 122 | 0.64, 0.68 | |

The shortfall is the only useful signal, and a weak one on tai60b. A node's own CPU has opposite
signs on two trees, and QPB's shortfall adds nothing.

**tai60b-261117, k = 3, 20 replays each** (`out/tai60b_weights.md`). CPU estimate over the truth:

| design, weights | budget | median | 10-90% | within 20% | within 2x |
|---|---|---|---|---|---|
| split-short-s90, exp(0.4 x) | 20% | 0.82 | 0.65-1.20 | 45% | 100% |
| split-short-s90, exp(0.4 x) | 30% | 0.92 | 0.85-0.99 | **90%** | 100% |
| split-short-s80, exp(0.8 x) | 30% | 0.85 | 0.74-0.97 | 70% | 100% |
| split-short-s80, oracle | 10% | 1.01 | 0.84-1.11 | 95% | 100% |
| split-short-s80, oracle | 20% | 0.98 | 0.95-1.07 | 100% | 100% |

- exp(0.4 x), tai60b's own slope, is better than shortfall^1 (90% against 75% within 20% at 30%).
- **By depth**, median estimate over the truth at 30% with exp(0.4 x):

  | depths | 0-3 | 4-8 | 9-13 | 14-23 | 24 and below |
  |---|---|---|---|---|---|
  | estimate / truth | 1.00 | 0.99 | 1.00 | 0.69 | 0 |

  The deep levels hold few nodes and little CPU. The oracle gets depths 14-23 right (0.99), but not
  24 and below.

**tho40ab-0** (a depth-2 subtree of tho40, 22,439 nodes, 43 CPU h, at most 8 levels deep;
`out/tho40ab_grid.md`):
- Even today's design is right there: 0.99 to 1.03 at 20-30% of the subtree, 77-100% within 20%.
- Survivor dives cost at least 26% of it: the pilot alone bounds most of a tree this small.
- A shallow, bushy tree is the case Knuth's dives were made for. tai60b's two deep chains are the
  case they miss.

## Where this leaves the estimate runs

- **The instances differ.** tho40's subtrees are shallow and bushy, and a Knuth estimate works there.
  tai-b's trees are chains with heavy subtrees off them. Of the estimates of record, wil50 and sko49
  are grids like tho40; tai60b is the chain case, and its truth is known.
- **For a chain-like tree, plan an estimate at 20-30% of the tree's CPU, and make it the first part
  of the solve:**
  - expand exactly the strata with the largest shortfall, in the order a best-first solve would;
  - give the rest to survivor dives weighted exp(0.4 x).
- **Every estimate should report its depth profile against the deepest level its dives reached.**
  When the deep bands are empty, the estimate is low, and the report says so.

## Round 3: what the as-built dives run costs after a pause (`dive_cost.py`, Nov 18)

The coordinator's flow is: solve, pause, then bound every open node and, where it branches, all
of its children (each stratum's setup), and at least 2 pilot dives per stratum with survivors.
None of that is capped by the per-unit budget. `dive_cost.py` prices it on the recorded trees,
with the solve approximated as breadth-first to a share of the CPU:

| tree | solve to | open nodes | setup | 2 pilot dives | solve + setup + pilot |
|---|---|---|---|---|---|
| tai50b-261116 (133 CPU h) | 25% | 773 | 32% | 23% | 79% of the tree |
| tai50b-261116 | 37% | 1,007 | 28% | 27% | 92% |
| tai60b-261117 (3,597 CPU h) | 10% | 7,742 | 26% | 35% | 71% |
| tai60b-261117 | 16% | 26,518 | 41% | 37% | 94% |
| tai60b-261117 | 25% | 42,843 | 37% | 28% | 90% |

**As built, an estimate of these trees costs nearly as much as the rest of the solve.** Most of
it is tree work the solve would do anyway: every open node bounded, and the walks down the heavy
subtrees. The replay's designs were cheaper for two reasons: they expanded only the
largest-shortfall strata, not the whole frontier, and their dives started from a few hundred
strata, not tens of thousands of open nodes. Before an estimate run on a tree like tai60b, the flow
needs three changes:
1. **Sample the open nodes** instead of bounding all of them: probabilities from the parent's
   shortfall, Horvitz-Thompson weights. This is still unbiased, at a cost proportional to the
   sample.
2. **A best-first solve before the pause**: the largest shortfalls first, as `split-short`, so the
   heavy part is solved exactly and the open nodes left are light.
3. **The budget caps the setup and the pilot too**, with the strata it never reached reported.

Each can be tried in the replay before any of it is built.

## Round 4: the flow with the fixes (`flow.py`, Nov 18)

Each replay runs the whole flow:
1. the solve, to a share of the tree's CPU, breadth-first (`bfs`) or best-first (`best`);
2. its open nodes Poisson-sampled with weights uniform, QPB (g^2 of each node's own QPB bound) or
   parent (exp(0.4 x) of its parent's shortfall), weighted 1/pi;
3. Knuth's or survivor dives below each sampled node, all under one total budget.

Knuth's dives go round the strata in turn. Neyman's allocation read low with their tails (the
estimator measured -6% on nug15 against -0.7% with equal turns). Grids: `out/flow_tai60b.md`,
`out/flow_tai50b.md`, `out/flow_tho40ab.md`.

| tree | best design | median est/true CPU | 10-90% | within 20% | cost / tree |
|---|---|---|---|---|---|
| tho40ab-0 (bushy) | solve 2%, sample 20% by QPB, Knuth, 10% budget | 0.98 | 0.85-1.22 | 85% | 0.10 |
| tho40ab-0 | the same at 20% | 1.05 | 0.93-1.25 | 90% | 0.20 |
| tho40ab-0 | solve 2%, every open node, Knuth | 0.95-1.01 | | 90-95% | 0.25-0.30 |
| tai60b-261117 (chains) | solve 10% bfs, sample 10% (parent or QPB), survivor, 20% budget | 0.84-0.88 | 0.74-1.31 | 70% | 0.19-0.20 |
| tai60b-261117 | the same, best-first, solve 5% | 0.87 | 0.60-1.42 | 50% | 0.20 |
| tai50b-261116 | (no design within 20% more than 25% of the time) | | | | |

- **Sampling the open nodes (fix 1) is what makes the flow affordable.** On tai60b it beats the
  earlier split designs at the same cost: 70% against 45-60% within 20% at 20% of the tree.
- **Best-first before the pause (fix 2) does not help.** Its best tai60b row is 50% within 20%,
  and the breadth-first solve is the coordinator's own order anyway.
- **Capping the setup (fix 3)** comes with the sampling: the setup is the sample's, not the frontier's.
- **The dives must fit the tree.** Knuth's dives on a bushy tree; survivor dives on a chain tree,
  where Knuth's never converge (tai50b: a median of 0.68 even after 200,000 dives).
- tai50b's tree is too small for any of this: its estimate costs about what solving it does.

**Built** (9f71ab2d): `QAP_EST_SURV=kQ` (Knuth's dives in strata), `QAP_EST_SAMPLE` /
`QAP_EST_SAMPLE_W` (open nodes sampled, each weighted 1/pi), `QAP_EST_ALLOC`, and
`QAP_ESTIMATE=SOLVE_S:DIVE_S:B:SAMPLE:W`. Two of the estimator's own early readings were biased low
and are fixed:
- nug15, Knuth's dives under Neyman's allocation: -5%;
- QPB-weighted sampling: -9%, because nodes their QPB bound fathoms were almost never drawn; they
  are now always kept.

After the fixes, nug15 reads 0.96-1.01 on the mean in every mode (`tests/estimate_test.py`), and
nug17 through the whole flow 0.998 with Knuth's dives and 20% sampling. **First runs**:
`tai50b-est` (survivor dives, every open node: the machinery's smoke test) and `tho40-est`
(`1:600:k1.5:0.2:qpb`: the first grid estimate).
