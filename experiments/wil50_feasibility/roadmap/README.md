# wil50 in a week on 1,000 cores: the roadmap (Nov 3, 2026)

Three research memos written on Nov 2-3 from the measurements of Oct 28 to
Nov 2: `bound_strength.md`, `node_speed.md`, `tree_size.md`. This page puts
their factors together. Every number is an estimate unless it names a log.

## 1. The baseline is larger than we had been saying

The brief gave 2.2e9-7e9 core-seconds after the cut phase (26-83 days on
1,000 cores). `tree_size.md` re-derives the tree from the one end-to-end
subtree (node 10: 2,193 DNN nodes at 63 Xeon-s each, `stronger_bounds/
search/wil50n10ext_base.log`) and the `tuning/lines` child data, and finds
two things the earlier estimate missed:

- **The symmetry file was never used.** NOTES.md section 4 counts the
  depth-3 nodes modulo the 5x10 grid's order-4 symmetry (29,580 orbits), but
  every wil50 run printed "couldn't find [./SYM/wil50.sym]". As configured
  today the search is 3.98x the estimate. The file is half a day's work on
  the transposed instance (the solver applies symmetry to the second matrix;
  the sko49 session did the same for the 7x7 grid) and is a prerequisite,
  not a bonus.
- **The depth-6 tail sets the total.** Two of the 16 sampled depth-3 nodes,
  the ones 329-405 short, carry 85-91% of the work; the estimate's step
  function at the level boundaries swings the total 2x on one sample. With
  the symmetry file: 3e8-1.8e9 DNN nodes, 1.0e10-6e10 CPU-s with the cut
  phase, **120-700 days on 1,000 Xeon cores, center about 290**.

So a week (6e8 core-seconds) needs about 12x at the center and 3-45x across
the range, after the symmetry file.

## 2. What the three angles offer (factors on the total)

| angle | item | factor | effort | first test |
|---|---|---|---|---|
| tree | symmetry file for wil50 (transposed instance) | 3.98x against today's runs, already in the estimate | **done** (`sym/README.md`: `latest/SYM/wil50t.sym`, `qaplib/wil50t.dat`) | passed: 15 and 635 children at depths 1 and 2 (651 nodes against 2,501 plain); node 10's bound agrees to 0.02 after a cut round with `-side 1` |
| tree | branch on the DNN marginal line (rule 9) at every DNN depth, DNN from depth 1 | **measured on the node-10 subtree: 1.41x less CPU, 1.13x fewer nodes** (rule 8 gives the same tree at 1.23x; `stronger_bounds/README.md` step 12); the 2-5x estimate was too high because rule 2 already finds decent lines | done | the node-3 subtree would confirm; DNN from depth 1 untested |
| bound | children inherit the parent's cut pool and multipliers (and a phase at clear misses to equip their children) | 1.5-4x on the dominant class (bound memo); 1.1-1.3x (tree memo) | 2-4 days | one laptop-hour: node 10's pool reloaded frozen on d1, d2, d4 and three siblings; killed if the median gain over the child's plain-DNN value is under +40 |
| speed | AMX bf16 for the filter and Rayleigh-Ritz products, plus CPU trims | **measured Nov 3 on the 32-node sample: 1.71x on fathoms, 1.75x on misses, same 18 fathoms, same iterations** (bundle D: AMX kernel for the five W*Y products, degree 3 with tight lambda, one-pass CholQR, refresh every 40; `amx/README.md`). Closed as not worth it: guard K/40, fp32 R (certificate margin), the small products on AMX, the row fusion; neutral: the symmetrisation skip. Left: the kernel's last third (ceiling 2 TF against 1 TF reached) | done | production: `DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40` |
| speed | cut-phase trims: shorter inner solves with Y/Z restoration (Astra's Y12s: `-warmY`, 12 rounds x 200 iterations, quiet cap 1,000) | **measured on the phase: 1.56x** (64% of baseline CPU, bound and fathom gates kept on all six nodes; `same_node_continuation/budget_study/RESULTS.md`); 1.17-1.20x on the total at the 46% phase share | done; subtree A/B next | fewer rounds (Y6, Y8s) lose 4-16 units of bound on the short-of-cutoff nodes and are out |
| speed | coordinator instead of the checkpoint barrier and node-end steals | 1.15-1.35x of the fleet | in plan.md | 70-85% busy today against 93-97% |
| tree | confirm give-ups uncut before the phase (QAP_DNN_GU_ZONE was off in run_ab3), primal gate in the phase, zone 0.0035 | **replayed, about 1.1x; closed** (`trims/README.md`): the deferral has nothing to save (the extrapolation rule already lets slow uncut fathoms finish), the primal gate 4:20 is in the code (default off, no false stop on 88 logged phases) and worth about 6%, the zone a modelled 4%; the search A/B was abandoned as not worth four hours | done | switch `QAP_DNN_CUT_PGATE=4:20` on in the production configuration; no further test |
| bound | skew-symmetric certificate polishing at non-converged iterates | 1.05-1.15x on the phase | 2-3 days | offline ascent on a dumped residual; killed under 20% of the gap recovered |
| bound | regional M2 block as a Lagrangian decomposition (Astra's 10b) | +20 units (10-15% chance of +65) | 2-3 weeks | project node 13's capped iterate onto the M2 set offline; killed if the excess over the triangles is under 0.03 |
| speed | GPU port, 8-16 nodes per device | wall clock: 100 H100s = about 7x the 1,000 cores | 3-6 weeks | one-day CuPy prototype; dropped above 10 ms per node-iteration on an A100 |

Closed this week, do not revisit without new evidence: 5-point counts, triple
consistency, the dual method (bundle), the local-7 family as a pool cut type,
Anderson/Halpern, momentum, atomic triangles, eliminations from the
certificate. Another full level from new inequalities at depths 3-5 is rated
under 10% within a month (`bound_strength.md`).

## 3. Multiplying out

Centers, on top of the symmetry file, as estimated on Nov 2: rule-9 branching
3x, inheritance 2x, per-node speed 2x, phase trims and fleet efficiency 1.3x:
about **15x**. Three of these are now measured (Nov 3): rule 9 gives 1.41x, not
3x; the Y12s phase budget gives 1.17-1.20x on the total; per-node speed gives
1.71x on the 32-node sample (AMX kernel plus eigensolver trims), close to its
2x. With the measured factors the center is about **6-7x**, which puts a 1,000-core week out of reach
unless inheritance delivers its upper estimate or the tree is at the small end
of its range; GPUs or more cores remain the other way to the week.
From a center of 290 days that is about 20 days on 1,000 Xeon cores; the low
end of the tree range (120 days) lands at about 8 days, the high end (700) at
about 45. **A week on 1,000 Xeon-class cores is reachable only if the tree
is at the small end of its range or the branching lever delivers its upper
estimate.** With 100 H100-class GPUs in place of the cores (about 7x the
wall-clock throughput at similar dollars, `node_speed.md`), the center falls
to about 3 days and the whole range fits in a week.

The two measurements that decide between these outcomes are cheap:

1. **The sampling plan** (`tree_size.md` section 6): build the production
   frontier at depths 0-2 (1.6e5 CPU-s), take 400 frontier nodes to the
   give-up point, run stratified subtrees end to end, dive into the ten
   worst tail nodes, and run a cutoff ladder. About 1e7-2e7 CPU-s, 3-6 hours
   on 1,000 cores, 2-3% of the week, and every subtree finished is a piece of
   the proof. It pins the total within 1.5x instead of 15x.
2. **The two subtree A/Bs** on the laptop or the cloud box: rule 9 at the DNN
   depths, and inheritance of the pool. Together they are the 2-5x and
   1.5-4x of the table, and each is a day of work plus a day of compute.

## 4. Order of work

1. ~~Symmetry file for wil50~~ done Nov 3 (`sym/README.md`): every wil50 run
   from here uses `qaplib/wil50t.dat` with `QAP_DNN_CUT_SIDE=1` and `-F` specs
   as j:i; the node-10 subtree A/Bs are unaffected (no symmetry below node 10).
2. Inheritance test (laptop, one hour) and the rule-9 subtree A/B (cloud, one
   day); per-phase logging in the search so that failures are visible.
3. ~~Confirm give-ups uncut, the primal gate, zone 0.0035~~ closed Nov 3 by
   replay (`trims/README.md`, about 1.1x; the gate ships default off). Still
   open here: the Y12s phase budget (`-warmY`, 12 x 200, quiet 1,000) as a
   search option, priced on the node-10 subtree against the rule-9 arm.
4. AMX and the CPU trims; the coordinator from docs/distributed/plan.md.
5. The sampling plan on the first 1,000 cores, before the full run.
6. In parallel, gated: the M2 projection test and the GPU prototype.
