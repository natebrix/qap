# Pool inheritance in the search

Astra's study (`research/pool-inheritance`, `experiments/pool_inheritance/`)
measured the inherited cut pool at six fixed children of node 10: frozen
inherited pool (arm B) against the uncut solve (arm A) at 1,500 cold
iterations, median gain +54, gates passed; continuing the Polyak rounds from
the pool (arm C) adds +10 to +21 before the quiet round, which does the work.
This experiment puts it in the search and measures the end-to-end factor on
the node-10 subtree against the stacked arm (`stronger_bounds/search`, step
13: 322 nodes, 32,740 CPU-s).

## Implementation (`latest/dnncut.{h,C}`, `latest/bnb.C`, `QAP_DNN_POOL=1`)

- A node whose cut phase ends unfathomed hands its children the pool of the
  solve their DNN starts from: the cuts with a positive multiplier, at those
  multipliers, in original indices (`DNNPool`). That solve is the quiet
  round's when one ran, else the best round's; its state (Y, Z) replaces the
  uncut state the children would otherwise start from. A phase that did not
  inherit hands on a pool only if that solve beat its uncut bound (on nug12
  the rounds never did, and the frozen multipliers cost a child 90 units).
- A child with a pool solves with it frozen in place of the uncut solve:
  same iteration limit and give-up rule, warm-started from the parent's state
  as before (`QAP_DNN_WARM_Y=2`). Triangles with a newly fixed member are
  dropped; the rest stay valid on the child's part of their region (every
  regional triangle and box is valid for any region and members, so a remap
  can weaken a cut, never invalidate it).
- If the child's cut phase runs, it starts from that solve: the inherited
  cuts at their multipliers, the inherited solve as round 0 and the best
  round so far. Its quiet-round test uses the best round's primal instead of
  the last round's (in Astra's C runs the rounds after an inherited start
  overshoot; the quiet round restarts from the best round). Without a cut
  phase the child hands the inherited pool on.
- The pool travels in state files (version 5, one token on the warm line), so
  it survives the parallel driver's task boundaries and checkpoints; the DNN
  state does not (unchanged).
- Checks: `dnncuts -loadpool` on Astra's parent pool keeps 6,395 cuts with a
  positive multiplier on d1, d2 and s5, as her runs report; every assembled
  term is checked on 20 random permutations (`QAP_DNN_POOL_CHECK`); nug12 and
  nug20 smoke tests (traces here), state save/resume with pools,
  `tests/resume_test.py` and `tests/parallel_test.py` pass. With
  `QAP_DNN_POOL` unset the search is unchanged.

## Protocol (predeclared, Nov 4)

`run_pool_search.sh`: exactly `run_ab6.sh`'s stacked arm (rule 9 at every
node through depth 6, `-x 0:4:1:300:1`, EXT 20, cut zone 0.003, quiet 1,500,
stall 3, bundle D, `QAP_DNN_WARM_Y=2`, 4 workers, checkpoint every 30
minutes) plus `QAP_DNN_POOL=1`. Same incumbent (48,816), same subtree.

- **Primary:** total CPU seconds for the subtree, against 32,740; the factor
  is their ratio. Also nodes, DNN-fathomed nodes and wall time.
- **Host:** logged at every (re)start (`hosts.txt`). The stacked arm finished
  Nov 4 03:39 UTC, before the slower host appeared (about 05:19), so it ran
  on the 2.10 GHz AMX host this run starts on. If any part of this run lands
  on another host type, its CPU is not compared directly; the stacked arm is
  rerun on the same host first.
- **Per node (`trace.txt`):** first-solve bound with and without a pool,
  fathoms by the first solve, cut phases run, iterations and CPU.
- **Reading:** a factor of 1.25 or more: adopt into the production stack; 1.0
  to 1.25: within the noise of one subtree, keep opt-in; below 1.0: drop.

## Results (Nov 4, 15:36-17:39 UTC, one host, no restart)

**1.37x end to end: the same tree in 23,905 CPU-s against 32,740.** By the
predeclared reading (1.25 or more), inheritance is adopted into the
production stack.

| arm | nodes / fathomed / branched | fathomed by the first solve | fathomed by a cut phase | CPU s | DNN first solves s | cut phases s | wall s |
|---|---|---|---|---|---|---|---|
| stacked (`run_ab6.sh`, Nov 3-4) | 322 / 315 / 7 | 295 | 20 | 32,740 | 18,256 | about 14,440 | 5,873 for the last segment (two resumes) |
| **stacked + pool inheritance** | **322 / 315 / 7** | **307** (228 of them with a pool) | **8** | **23,905** | **16,854** | **7,016** | 7,342 (3 checkpoints) |

Both arms ran on a 2.10 GHz Xeon with AMX (`hosts.txt`; the stacked arm
finished before the slower host first appeared), with the same incumbent and
four workers. The node counts are identical: the tree is the same, and the
saving is all in the cost per node. The stacked arm's cut-phase time is its
total CPU less the DNN first solves and about 40 s of QPB work (this arm's
outside-DNN time is 36 s).

From the per-node trace (`trace_summary.txt`, `analyze_trace.py`):

- **All 228 nodes that inherited a pool fathomed on their first solve,** in
  101 iterations on average (59 at depth 6, 128 at depth 5), keeping a median
  1,795 of the parent's cuts. None needed a cut phase. In the stacked arm the
  same children start from the parent's uncut state, and 20 nodes of the
  subtree fathomed only in a cut phase.
- **Only 13 cut phases ran** (8 fathomed, 5 handed a pool on: +10.7, +26.1,
  +42.4, +50.6, +57.9 over their uncut bounds), all at nodes without a pool,
  for 7,016 CPU-s against about 14,440 in the stacked arm. The five pools
  reached all 228 pooled nodes, the children of those five misses (no pooled
  node branched, so no pool went further down).
- The first solves cost less too (16,854 against 18,256 s): the pooled
  children stop at 101 iterations on average.
- The validity guard (`QAP_DNN_POOL_CHECK=20`) checked every assembled term
  on 20 random permutations and never fired.

Caveats: one subtree (node 10), so the factor is measured where the stacked
arm was, not on a sample of the tree; the fraction of work in children of
cut-phase misses elsewhere decides how much it carries over. The stacked
arm's wall time is not comparable (two container restarts). The trace was
not run in the stacked arm, so its per-node split is inferred from totals.
