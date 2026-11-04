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

## Results

(running)
