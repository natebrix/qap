# Convergence of the per-node DNN pipeline, second pass (Nov 4, 2026)

Question: beyond new cuts (closed), what else could make the DNN bound
converge faster on wil50? Running notes are in `NOTES.md`, tables in `logs/`,
and test logs in `runs/`.

**Answer: yes, three ideas worth building. None is in the ADMM step itself:
they are in the per-node pipeline, and all three are about warm starts.**

1. **Run a short cut phase at clear misses that will branch, to give their
   children a pool.** On node 10's children this turned all 5 of the
   search's phase fathoms into first-solve fathoms, at about 1/4 of the
   iterations.
2. **Keep the parent's DNN state when a node crosses a coordinator unit
   boundary.** Siblings that lost it took 2.8-6.9x more iterations to
   fathom. Storing that state in bf16 changed nothing in 4 children.
3. **Restore Y between cut rounds.** Today every round restarts from
   Y = e0e0'. That restart, not the certificate, causes most of the 25-90
   units of "certificate lag" in the phase. The switch is now wired;
   200-iteration rounds are the setting to carry, and 100-iteration rounds
   failed.

The ADMM step itself has nothing new beyond the Nov 3 memos' open items.
Every bound stays certified: these ideas change only start points,
multiplier schedules and which nodes run a phase.

[m] = measured (log named); [e] = estimate. Paths are relative to
`experiments/wil50_feasibility/` unless they start with `latest/`.

## 1. Where the DNN CPU goes (node-10 subtree, coordinator pool arm, `coord_search/trace.txt`)

Script: `roadmap/convergence2/trace_breakdown.py`. Total 18,465 CPU-s.

| class | nodes | mean iterations | CPU-s | share |
|---|---|---|---|---|
| cold first solves (no parent DNN state: lost at a unit boundary or at the initial split) | 93 | 162-1,007 | 7,523 | 41% |
| cut phases | 13 | 600-2,560 | 6,448 | 35% |
| warm first solves (228 of them with an inherited pool) | 228 | 48-68 | 4,379 | 24% |

In the `bench/parallel.py` run of the same tree, 146 of 322 nodes were cold.

## 2. What the code and logs show

**F1. The cut phase restarts every round from Y = e0e0'.** In the search,
`DNNCutParams::warm_primal` is 0, and `bnb.C` never set it. So every Polyak
round and the quiet round keep Z but reset Y.
- Budget-study logs, re-read (`round_lag.py`,
  `logs/round_lag_budget_study.txt`): the median of primal − certificate at
  the end of a 300-iteration round is 44-67 with Y reset (arm B, today's
  search) and 4-10 with Y restored (Y12). [m]
- The quiet round then adds 37-60 under B and 5.5-10 under Y12. [m]
- In the search, 25% of phase iterations are quiet rounds
  (`phase_split.py`). [m]
- Nodes 3.5-7.6 short after 990-1,190 uncut iterations still spend
  870-1,370 phase iterations to gain 8-9 units. [m]

So most of the "inner lag" in `research_bound_quality.md` is this restart,
not the certificate.

**F2. Cold nodes.** The coordinator ships a node's cut pool to other units,
but not its parent's DNN state (Y and Z). Compare siblings of the same
parent with the same pool status, differing only in whether that state
survived (`sibling_warm_cold.py`) [m]:
- fathoms take 49 vs 339 iterations (depth 5, no pool), 51 vs 174 (depth
  6, pooled), 70 vs 195 and 70 vs 268 (depth 5, pooled, parallel.py run);
- misses give up at 320-330 iterations warm vs 872-1,007 cold.

On the 16 depth-4 sample children, cold vs warm Z+Y: 2,440 vs 1,360
iterations on the 7 shared fathoms (`cold_vs_warm.py`). [m]

**F3. Node 10 ran no phase.** It is 232 short at its give-up, outside the
146-unit cut zone, so none of its 47 children inherited a pool. Its child
49:45 (186 short, also outside the zone) carries 137 nodes and 8,002 CPU-s,
43% of the subtree (`subtree_cost.py`). [m]

**F4. The "overshoot from inherited pools" never happened in the search.**
All 228 pooled nodes fathomed on their first solve, so none ran a phase. The
overshoot appears as soon as pooled nodes do run phases (test E2 below).

**F5. Every Polyak run overshoots the same way.** Every budget-study arm and
node 10's phase drop 18-75 units in round 4. In 100-iteration rounds the
drops reach −88 to −302. [m]

## 3. Ranked candidates

Factors are on wil50 DNN node CPU at depths 3-6.

| # | idea | factor | effort | validity | cheapest test; what kills it |
|---|---|---|---|---|---|
| 1 | **Equip phase:** a short cut phase (no quiet round, stall rule) at out-of-zone misses at depths 3-4 that will branch, so their children inherit a pool. | per child [m] (test E). On the subtree [e]: 1.2x if only node 10's direct children change, 1.5-1.8x if the 49:45 subtree shrinks as its children's pools suggest. | none for a first test (`QAP_DNN_CUT_ZONE=0.005`); 0.5 day for a separate equip zone | none: pool cuts are already certified, and `QAP_DNN_POOL_CHECK` guards their assembly | node-10 subtree through the coordinator with zone 0.005, same host as the 18,496 CPU-s pool arm; kill if under 1.2x |
| 2 | **Keep the parent's Y and Z across unit boundaries:** ship them in bf16 in `state.C` (upper triangles, about 10 MB per parent at m = 47, shared by its siblings), or re-bound the parent for 300 iterations when a unit starts. | about 1.3x on this subtree [e from m]: 41% of CPU is in cold nodes, 2.8-6.9x per sibling. At fleet scale it depends on how often units split. | 1-2 days (re-bounding: 0.5-1 day) | none: only the start point changes | test S passed [m]; then the node-10 coordinator rerun; kill if under 1.1x |
| 3 | **Phase with Y restored:** `QAP_DNN_CUT_WARMY=1` (wired now), 200-iteration rounds with the quiet round at most 1,000 (the budget study's Y12s). For phases that start from a pool, also damp the first step (θ about 0.3 or δ about 10-20; about 10 lines). | 1.56x on the phase [m, budget study]; about 1.15x on the total at today's 35% phase share [e]; more with #1 | hours | none | node-10 subtree A/B with the switch; kill if any fathom is lost or the gain is under 1.05x |
| 4 | The Nov 3 memos' open items, unchanged: restricted-PRSM pieces (Graham et al.), two-phase β, fp32 rank term off certificate iterations, complementary-subspace Rayleigh-Ritz. | 1.0-1.5x; 0-15%; 4% per iteration [m]; 5-9% per iteration [e] | 2-3 days; 2 days; 0.5 day; 3-5 days | as in the memos | as in the memos |

Ideas 1-3 stack, with some overlap: together about 1.8-2.5x on the node-10
subtree [e]. With #1, more pooled nodes run phases, so #3's first-step
damping becomes necessary.

**Closed or rejected:**
- **Parent investment** [m, test P]: with the parent at 1,000 iterations
  instead of 300, the 5 fathoming children fell only from 950 to 830
  iterations (−13%), and misses were unchanged. That roughly pays back the
  parent's extra 700 iterations over 47 children: about break-even.
- **100-iteration warm rounds** [m, test C]: misses end 13-23 below today's
  protocol.
- **Earlier hand-over from the uncut solve to the phase** [e]: the trims
  replay shows the extrapolation rule already lets slow uncut fathoms
  finish. Most of the cold-node overrun disappears with #2.
- **Adaptive certificate checks near the cutoff** [e]: 1-2%.
- **Sibling pool sharing:** dominated by #1. It reaches only siblings
  bounded after the first sibling's phase, and needs a new channel between
  units.
- **Cut multipliers updated inside the ADMM (CP-ADMM)** [e]: the waste it
  targets is the Y reset, which #3 removes for free. Below 200-iteration
  rounds the multiplier updates degrade (Y6, Y8s and Y12q all fail the bound
  gate).
- **A general adaptive-δ Polyak step:** worth about one 300-iteration round
  per phase; only the first-step damping from pools is worth doing.
- **Certifying an extrapolated Z:** a certificate-side question (see
  `../certificate/REPORT.md`).

## 4. Tests (cores 0-1, from 22:12 UTC; logs in `runs/`)

**E — equip phase at node 10** (`E_n10_phase.log`, `E_c*.log`,
`logs/E_P_S_summary.txt`)
- Node 10's phase reproduces the search exactly: 48,582.4671 at iteration
  300.
- It gains +63, +90 and +102 in rounds 1-3, drops −18 in round 4, and stops
  on the stall rule. Cost: 1,200 iterations, 446 CPU-s.
- 13 of its depth-4 children were then solved cold with that pool frozen,
  giving up after 300 iterations as in the search:

| child group | search (first solve + phase) | with node 10's pool (first solve only) |
|---|---|---|
| 5 children that fathomed in a phase | 3,030 / 2,600 / 1,860 / 1,340 / 2,020 iterations | fathom at 1,100 / 530 / 350 / 390 / 450 (3.8x fewer in total) |
| 3 phase misses | 106-130 short after the first solve; 64-79 short even after their phase | 42-44 short, inside the zone |
| 49:45 | 186 short, out of zone, branched | 82.6 short, inside the zone |
| 4 easy fathoms | 130 / 200 / 290 / 450 | 120 / 150 / 180 / 220 (−37%) |

**E2 / E2w — 49:45's own phase, started from the inherited pool**
(`E2_c49_phase.log`, partial; `E2w_c49_phase_warmY.log`)
- With Y reset (today): round 1 drops 48,731.5 → 48,615.6 (−116); the phase
  stalls after 3 rounds with no gain.
- With Y restored: round 1 still drops (−88), then reaches 48,746.5 (+15)
  before the stall.
- So the inherited-pool overshoot is the first Polyak step on near-optimal
  multipliers. Restoring Y helps the recovery, not the drop itself.

**C — warm Y, 12 rounds of 100, quiet round 500** (`C_*_Y12q.log`,
`logs/C_summary.txt`)
- Compared with the budget study's logs. Round 0 agrees within 0.013 (not
  bit-identical; the solver has changed since that study).

| node | today (B) | Y12s (200-iteration rounds) | Y12q (this test) |
|---|---|---|---|
| n14 fathom at | 4,750 | 2,110 | 1,680 |
| n6 fathom at | 1,260 | 780 | 780 |
| n9 fathom at | 990 | 720 | 700 |
| n13 final bound (miss) | 48,754.6 | +5.2 | **−13.2** |
| n10 final bound (miss) | 48,743.9 | +4.5 | **−23.1** |

Y12q is killed for its misses; carry Y12s instead.

**P — parent at 1,000 instead of 300 iterations** (`P_13_J1000.log`,
`P_7_J1000.log`)
- The parents' bounds at iteration 300 equal warm_y's exactly, so the
  comparison is like for like.
- Fathoming children: 950 → 830 iterations (−13%). Misses still give up at
  300-310, with +4 to +12 at iteration 200.

**S — parent state rounded to bf16** (`S_13_J300_bf16.log`)
- Parent 13's 4 children take 340 / 300 / 310 / 90 iterations, identical to
  the unrounded runs. Bounds at 50/100/200 agree within 0.8.

## 5. Code (opt-in only; defaults unchanged)

- `latest/bnb.C`: `QAP_DNN_CUT_WARMY=1` sets `cp.warm_primal`, so cut rounds
  and the quiet round restart from the previous Y as well as Z.
- `latest/dnntest.C`: `-J N` (the parent's own iteration budget), `-K
  i:j,i:j,...` (several children of one parent), `-Wr 16|32` (round the
  parent's Y and Z to bf16 or fp32 before the children start).
- `experiments/dnn_cuts/dnncuts.C`: `-H N` (give-up of the `-loadpool`
  first solve), `-vc` (certified bound at every check).
- Scripts here: `trace_breakdown.py`, `sibling_warm_cold.py`,
  `cold_vs_warm.py`, `round_lag.py`, `phase_split.py`, `subtree_cost.py`,
  `summarize_c.py`, `summarize_ep.py`, `run_tests.sh`, `run_tests4.sh`,
  `build.sh`.
- `E2_c49_phase.log` is partial: it was stopped in its quiet round, which
  the search would not run because the primal was below the cutoff.

Literature, checked only lightly: inexact-oracle and on-demand-accuracy
bundle methods (de Oliveira, Sagastizábal, Lemaréchal, Math. Prog. 148,
2014) support short inner solves plus one accurate final evaluation, which
is #3's structure. Whether BiqCrunch passes multipliers from parent to child
could not be confirmed.

## Measured vs estimated

**Measured [m]:**
- The Y reset causes the round lag: 44-67 vs 4-10; quiet-round gain 37-60
  vs 5.5-10.
- 25% of phase iterations are quiet rounds; near-near misses spend
  870-1,370 phase iterations.
- Cold vs warm siblings: 2.8-6.9x on fathoms; misses give up at 320-330 vs
  872-1,007; 41% of DNN CPU in cold nodes; 1.8x on the 7 depth-4 sample
  fathoms.
- Test E per child: 5 phase fathoms become first-solve fathoms (3.8x fewer
  iterations); misses move to 42-44 short; 49:45 moves from 186 to 83 short;
  easy fathoms −37%; phase cost 1,200 iterations.
- E2/E2w: first-step overshoot from a pool (−116 / −88); gain +0 / +15.
- Test C: Y12q 1.4-2.8x faster on near misses, but −13/−23 on misses
  (killed).
- Test P: −13%, break-even.
- Test S: bf16 makes no difference on 4 children.
- Budget study (re-read): Y12s 1.56x on the phase.

**Estimated [e]:**
- The equip phase's subtree factor (1.2x floor, 1.5-1.8x central).
- About 1.3x from keeping parent state across units, and its share at fleet
  scale.
- About 1.15x on the total from #3.
- 1.8-2.5x for #1-#3 combined.
- The rejected items' 1-2% (adaptive checks) and ~0 (earlier hand-over).
- The memo items' factors.

## Checks after the merge into cognac (Nov 4, 23:45 UTC)

- **F1 is in the code.** `dnn_cut_rounds` passes `par.warm_primal =
  cp.warm_primal && r > 0` (`latest/dnncut.C`). `cp.warm_primal` defaults
  to 0, and before this merge `bnb.C` set only the child warm start
  (`dp.warm_primal`), never the cut phase's. `dnn_bound` starts Y at
  e0e0' unless `warm_primal` is set. The certificate study's description
  ("PRSM restarts from the old (Y, Z)") is therefore wrong for the
  production protocol: only Z is kept. Its measurement still holds: 84% of
  the round-end lag at [9:45] is the certificate below its converged
  value. F1 explains why: the iterate is rebuilding Y.
- **Section 1 reproduces.** `trace_breakdown.py` on
  `coord_search/trace.txt` gives 322 nodes, 12,016 CPU-s of first solves
  and 6,448 CPU-s of cut phases. The cold rows (warm 0) sum to 94 nodes
  and 7,637 CPU-s with node 10 itself, 41%.
- **Test E agrees with the search.** The "search" column of
  `logs/E_P_S_summary.txt` matches the trace: 48:45 is 81.2 short after 470
  iterations and fathoms in a 2,560-iteration phase.
- **The 1.2x floor for #1 assumes the misses cost nothing extra.**
  - 49:45, 29:45, 19:45 and 9:45 start 42-83 short with node 10's pool,
    inside the zone, so each runs a phase.
  - E2 shows 49:45's phase from the pool stalling with no gain; it still
    branches.
  - Counting those phases (about 1,200 iterations each) and node 10's own,
    the measured part comes to about 1.15x on the subtree. The rest of the
    central estimate (1.5-1.8x) depends on 49:45's 137-node subtree doing
    better with the extra pool, and that is unmeasured.
- **Idea 3 is not new.** Astra's Y12s measured 1.56x on the phase on Nov
  3, and Nathan declined it then as too small (`../README.md`, section 2).
  What is new is the mechanism (F1), and that #1 raises the phase share,
  which raises #3's value.
- **Build and tests.** After the merge, `make solver dnntest`, `make -C
  latest dnncuts` and `make dnncheck` pass. `dnntest -K` takes both `,` and
  `;` (the certificate study added the same option with `;`).
