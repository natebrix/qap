# Pipeline A/B: the Nov 4 candidates on the node-10 subtree

The convergence and certificate studies (`../roadmap/convergence2/REPORT.md`,
`../roadmap/certificate/REPORT.md`) left four candidates. Each is tested
end to end on the node-10 subtree through the coordinator, against the
production configuration on the same machine at the same time.

| arm | what it adds | judged against | kill line (DNN iterations) |
|---|---|---|---|
| `base` | production: rule 9, bundle D, parent-Y warm start, pool inheritance (`../coord_search/run_coord_search.sh`) | - | - |
| `zone` | `QAP_DNN_CUT_ZONE=0.005` (from 0.003): misses up to 244 short run a cut phase, node 10 (232 short) among them, so their children inherit a pool (the "equip phase") | base | 1.15x |
| `warmy` | zone + Y kept between cut rounds (`QAP_DNN_CUT_WARMY=1`), 12 rounds of 200 iterations, quiet round at most 1,000 (Astra's Y12s), first Polyak step damped (0.3) for phases that start from an inherited pool | zone | 1.05x on the matched nodes and the phases both ran |
| `certx` | zone + the certificate polish within 10 units of the cutoff (`DNN_CERTX=32 DNN_CERTX_USE=1 ...`); on top of base until Nov 5 | zone | 1.03x on the matched nodes |
| `zone2`, `base2` | replicates: the run-to-run noise of one configuration on one machine | zone, base | - |

Iterations are the deciding measure: they do not depend on the machine.
Running two arms at once on one machine also makes their CPU seconds
comparable. The smoke runs on nug20 (`runs_smoke/`, not committed) checked
the plumbing only.

## Running it (laptop)

```
git pull origin cognac
make solver                      # latest/bnbdrv; on macOS with Accelerate
experiments/wil50_feasibility/pipeline_ab/run_pair.sh base zone
```

- Each arm gets half the performance cores (`SLOTS=N` overrides).
- The script keeps the Mac awake (`caffeinate`) while it runs. Keep the lid
  open, or the power and an external display connected.
- `kill PID` (the script prints it) or ^C stops it cleanly. The workers
  drain, and the same command resumes the arms where they stopped. A
  finished arm is skipped.
- At the end it writes `summary.txt` and prints the `git add/commit/push`
  line for the results: `runs/ARM/` trace, results, logs and hosts. The
  stores, work directories and binaries stay out of git.
- The second pair, `run_pair.sh warmy certx`, runs from the same build: the
  first-step damping it uses has been in since fdb845d. The script refuses
  `warmy` on a binary without it.

`summarize.py runs/` can be run at any time; unfinished arms are listed as
such.

**The Nov 6 levers** (`../speed/nov6_admm/`, `../speed/nov6_core/`), two
sessions of about 25 minutes each on the M2 Pro, each pair at the same time
so that its CPU ratio means something:

```
experiments/wil50_feasibility/pipeline_ab/run_pair.sh fast wlowc     # warm-call beta 0.25: kill line 1.15x in DNN iterations, matched nodes
experiments/wil50_feasibility/pipeline_ab/run_pair.sh best bestcore  # the three core switches on top: kill line 1.05x in DNN CPU, matched nodes
```

`best` against `wlowc` (the crossing certificate, kill line 1.02x in DNN
CPU) spans the two sessions, so its CPU ratio is weaker than the in-session
ones; its iterations are not. `core` (`fast` + the core switches) and `wlow`
(0.25 for the cut rounds too) are there for a third session if a result
needs splitting. `summarize.py` prints every pair and its kill line. The
gates for a binary with these switches are `make test`, `make dnncheck`,
`make shipcheck` and `make continuationcheck` with the switches unset; with
any `DNN_BETA_SCHED` exported, `tests/dnn_test.py --parts ACDE`, since part
B compares against the fixed-beta Python prototype and fails by design.

## The damping setting (Nov 5, cloud box, `damp_test.sh`, logs in `damp/`)

Child 49:45 of node 10 with node 10's equip-phase pool
(`../roadmap/convergence2/runs/E_n10.pool`), the warmy arm's phase: Y kept
between rounds, rounds of 200 iterations, stall rule 3, uncapped target,
quiet round of 1,000 run regardless (in the search it runs only when the best
round's primal is above the cutoff, which none of these is).

| damping (first step : later steps) | rounds 1-4 | best round | after the quiet round |
|---|---|---|---|
| none | −102.4, −14.3, +13.3 (stall) | 48,744.8 | 48,760.9 (+29.5) |
| 0.3 : 1 | +16.2, +31.8, −227.6, −68.9 | **48,763.3** | **48,776.2 (+44.8)** |
| 0.3 : 0.5 | +16.2, +30.8, −72.2, −10.9 | 48,762.2 | 48,775.1 (+43.7) |
| 0.3 : 0.3 | +16.2, +30.0, +2.8, +30.2 | 48,761.7 | 48,771.3 (+39.8) |

- The first step is what matters: damped, the inherited multipliers improve
  at once instead of losing 102 units. The phase ends 18 units higher with
  one more round (800 iterations against 600), and the multipliers it hands
  to 49:45's children are those of a better round.
- Damping the later steps removes the round-3 overshoot but does not improve
  the best round. The warmy arm uses `QAP_DNN_CUT_DAMP1=0.3` (first step
  only).
- 49:45 stays about 50 units short either way, so it branches in every
  variant. Whether the better pool helps its children is for the search A/B
  to show.

## Results

### base vs zone (laptop, Nov 4; `runs/base`, `runs/zone`, `summary.txt`)

Both arms at f6d448e, Apple M2 Pro, 3 slots each, run at the same time
(16:07-17:21 PST). **zone passes**: 1.635x fewer DNN iterations against a
kill line of 1.15x.

| arm | nodes / fathomed | cut phases (fathomed) | DNN iterations (first + phase) | DNN CPU | wall (3 slots) |
|---|---|---|---|---|---|
| base | 454 / 444 | 15 (6) | 83,640 (59,600 + 24,040) | 11,318 s | 4,431 s |
| zone | 366 / 358 | 8 (0) | 51,150 (40,650 + 10,500) | 7,985 s | 3,192 s |

- Iterations 1.635x, DNN CPU 1.417x, wall 1.388x. The CPU ratio is below
  the iteration ratio because the pooled iterations are dearer per step.
- The root's cut phase (zone only: 48,582 -> 48,685 in 1,200 iterations,
  still 129 short) fathoms nothing itself; the gain is in the children.
  Every zone node is a base node (0 nodes only in zone, 88 only in base).
  On the 366 shared nodes the iterations fall 74,390 -> 51,150 (1.454x):
  first solves are shorter from the inherited pool, and only 8 nodes need
  a phase of their own against 15.
- Two shared nodes change outcome, both to the good: a depth-5 node base
  branched after a 2,400-iteration phase (best 48,807.8) and zone fathomed
  in the first solve (48,814.1 at 1,570 iterations); a depth-7 node
  likewise (base phase best 48,813.4; zone first solve 48,814.8 at 190
  iterations). Those two account for the 88 extra base nodes.
### Checks of the base/zone result (cloud box, Nov 5)

**The base arm is the same configuration as the cloud's Nov 4 coordinator
run** (`../coord_search/`, 322 nodes, 58,870 DNN iterations, 18,496 CPU-s on
4 AMX Xeon slots). The tree differs not because of the configuration but
because of chance; an earlier bullet here said otherwise.
- On the 212 nodes the two base runs share with the same start (parent-Y warm
  start and pool status alike), the laptop needs 3% more iterations
  (35,520 against 34,310). The numerics (Accelerate against OpenBLAS and
  AMX) matter little.
- The difference comes from two places.
  - **Marginal nodes:** two nodes of the 49:45 subtree fathom in the cloud
    run and branch here: [0:8 49:45] misses by 6.3 after its phase, [0:25
    9:21 49:45] by 9.6, and its child [0:25 19:0 9:21 49:45] by 0.65. That
    adds depths 7-8: 132 nodes not in the cloud tree.
  - **Cold starts:** the laptop base split into 33 units (zone 14, the
    cloud 23), and 232 of its 454 nodes started without the parent's Y and
    Z (cloud 94 of 322).
- **Per iteration** on identical first solves (same node, same start, same
  iteration count; 168 nodes), the M2 Pro takes 2,235 CPU-s against 4,460
  for the AMX Xeon: **2.0x per core**.

**Zone, measured three ways:**

| comparison | DNN iterations | factor |
|---|---|---|
| zone against the laptop base, end to end (the same machine and hour) | 83,640 -> 51,150 | **1.64x** |
| the same, on the 366 nodes both trees have | 74,390 -> 51,150 | 1.45x |
| zone against the cloud base run, end to end (a luckier base: 322 nodes) | 58,870 -> 51,150 | 1.15x |
| node 10's 47 children alone (depth 4), against either base | 22,940 (laptop), 22,100 (cloud) -> 13,910 | **1.65x** |

- **The mechanism is the one test E predicted** (`../roadmap/convergence2/`).
  Node 10's phase (1,200 iterations, 48,582 -> 48,685) hands its 47
  children a pool. 43 of them fathom on the first solve against 38, and 4
  run a phase of their own against 9. The 184 depth-5 nodes fall from
  31,440 to 20,990 iterations.
- **What zone still wastes:** its phases from inherited pools gain little.
  On the 7 nodes with a phase in both arms the mean gain over the first
  solve is +7.7 against +30.8 in base; five zone phases of 900 iterations
  gain nothing at all (4,500 iterations, 9% of the arm). That is the
  first-step overshoot the warmy arm's damping targets (49:45: −102 -> +16,
  above).
- **Noise.** One configuration gives 322 or 454 nodes, depending on which
  marginal nodes fathom and which nodes start cold. ~~The end-to-end factor
  of a single pair is therefore uncertain by about 1.4x.~~ (Corrected after
  the replicates, below: the 1.4x was between machines and slot counts, the
  cloud's 4 AMX slots against the laptop's 3 Accelerate slots. On one
  machine a configuration repeats within 1%.) Zone passes its kill line
  even against the cloud base (1.15x). `summarize.py` also judges warmy and
  certx on the matched nodes (shared, same start) and, for warmy, on the
  phases both arms ran.
- **Adopted:** `QAP_DNN_CUT_ZONE=0.005` is the new reference. `certx` now
  runs on top of zone and is judged against it.

Next: `run_pair.sh warmy certx`, both judged against `runs/zone`. Then, if
the laptop is free, `run_pair.sh zone2 base2` for the noise.

### warmy and certx vs zone (laptop, Nov 4 evening; `runs/warmy`, `runs/certx`)

Both arms at ff7d517 (the base/zone pair was at f6d448e; the arms differ
only in the Nov 4 coordinator and damping commits), 3 slots each, run at
the same time (21:47-22:35 PST). **warmy passes, certx fails.**

| arm | nodes / fathomed | cut phases (fathomed) | DNN iterations (first + phase) | DNN CPU | wall (3 slots) |
|---|---|---|---|---|---|
| zone (Nov 4 afternoon) | 366 / 358 | 8 (0) | 51,150 (40,650 + 10,500) | 7,985 s | 3,192 s |
| warmy | 232 / 227 | 6 (1) | 30,260 (24,910 + 5,350) | 4,811 s | 2,036 s |
| certx | 322 / 315 | 10 (3) | 49,150 (34,900 + 14,250) | 6,871 s | 2,833 s |

warmy (Y kept between cut rounds, 12 rounds of 200, quiet at most 1,000,
first step damped 0.3 on inherited pools):
- Kill line 1.05x on the matched nodes: **1.159x** (186 matched nodes,
  27,510 -> 23,730 iterations); end to end **1.690x**, tree 366 -> 232.
- The root's phase reaches 48,725.6 in its 1,200 iterations against
  48,685.1 for zone's cold rounds: +40 more from the same budget, and the
  children inherit the stronger pool. First solves 40,650 -> 24,910.
- Two depth-5 nodes flip to fathoming (one in its phase at 48,814.4, one
  in the first solve at 48,814.2); the 134 nodes under them are the whole
  tree difference. No node is only in warmy. The 6 shared phases gain +40.1
  over the first solve against +17.1.
- Stacked on zone over base: 83,640 -> 30,260 DNN iterations, 2.76x end to
  end, 454 -> 232 nodes.

certx (the certificate polish within 10 units of the cutoff):
- Kill line 1.03x on the matched nodes: **0.983x** (321 matched, 48,190 ->
  49,030), FAIL. End to end 1.041x and 1.16x DNN CPU, all from one depth-6
  node whose phase reaches 48,814.1 at 1,780 iterations where zone's
  stopped at 48,808.5 after 2,400; the 44 nodes under it are the tree
  difference. Per node the polish costs slightly more than it saves.
- The one flip is a near-miss closed by the polish, which is the effect it
  was built for, but it does not pay on the matched nodes. Retired from the
  production candidates; the near-miss mechanism, if wanted, belongs in the
  cut phase's quiet round rather than a separate pass.

Production candidate after this pair: zone + warmy
(`QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12
QAP_DNN_CUT_ITERS=200 QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3`).
Replicates (`zone2 base2`) and the `ship` arm are the remaining pairs.

### Replicates and ship (laptop, Nov 5 early morning; `runs/zone2`, `runs/base2`, `runs/ship`)

All at 08ff251. `zone2 base2` ran together (03:44-04:57 PST, 3 slots
each); `ship` ran alone afterwards at the same 3 slots (04:57-05:22).

| arm | nodes / fathomed | cut phases (fathomed) | warm / pooled nodes | DNN iterations (first + phase) | DNN CPU | wall (3 slots) |
|---|---|---|---|---|---|---|
| base | 454 / 444 | 15 (6) | 222 / 406 | 83,640 | 11,318 s | 4,431 s |
| base2 | 454 / 444 | 15 (6) | 229 / 406 | 83,030 | 11,102 s | 4,373 s |
| zone | 366 / 358 | 8 (0) | 248 / 365 | 51,150 | 7,985 s | 3,192 s |
| zone2 | 366 / 358 | 8 (0) | 248 / 365 | 51,150 | 7,895 s | 3,111 s |
| ship | 322 / 315 | 7 (0) | 321 / 321 | 30,810 (22,710 + 8,100) | 3,769 s | 1,473 s |

Replicates: the search is reproducible through the coordinator. zone2 is
zone to the iteration (366 nodes, 51,150 iterations, every per-depth count
equal); base2 is base to within 0.7% in iterations (the 7 nodes whose warm
status depends on unit boundaries: 229 warm against 222) and identical on
the 447 matched nodes. CPU repeats within 1-2%. The 1.635x of zone over
base is therefore not noise, and differences of a few percent between arms
on this subtree are real.

(Cloud box, Nov 5:) One caution on "real". The replicates show that one
configuration repeats. Any change to the configuration also moves the
handful of marginal nodes of the 49:45 subtree, which fathom or branch
into 44-134 nodes. certx is the example: a change that does nothing on
the 321 matched nodes (0.983x) flipped one node and came out 1.041x end
to end. So an arm's end-to-end factor is its effect plus how its
perturbation lands on those few nodes. The per-depth figures carry the
mechanism: depth 4 (node 10's children) is 1.65x for zone and 1.92x more
for ship, with no marginal flips involved. A second subtree would put an
error bar on the end-to-end numbers (`../roadmap/README.md`, the sampling
plan).

ship (zone + the parent's DNN state carried across coordinator units,
`QAP_DNN_SHIP=32`), against zone:
- DNN iterations 51,150 -> 30,810: **1.660x**; DNN CPU 7,985 -> 3,769 s:
  2.12x; wall 3,192 -> 1,473 s: 2.17x. Tree 366 -> 322.
- Every node is warm (321 of 321 against 248 of 365): the 118 nodes that
  started cold at unit boundaries in zone are the whole first-solve
  difference (40,650 -> 22,710). Cut phases 8 -> 7, phase iterations
  10,500 -> 8,100.
- The CPU ratio (2.12x) exceeds the iteration ratio (1.66x) because ship's
  iterations were cheaper at every depth: 139 against 173 ms at depth 4,
  124 against 157 at depth 5, 102 against 136 at depth 6 (first solves; zone
  and zone2 agree within 2 ms). That is the machine, not the mix of
  depths: ship ran alone (3 solver processes, not 6). Warm iterates might
  also make each iteration cheaper, but a run alone cannot tell. The
  `full ship2` pair, with ship2 beside full, settles it. An earlier version
  of this bullet put the CPU gap down to cold first solves running longer;
  that is in the iteration ratio, not on top of it.
- The tree shrinks by 44 nodes, the same depth-7 subtree certx removed:
  warm first solves certify a node zone's cold first solve only branched.
  No node is only in ship.
- On the 205 matched nodes (same warm status in both) ship is 0.978x:
  carrying the state costs nothing and gains nothing there; the whole
  effect is the cold starts. The 7 shared phases are identical (8,100
  iterations in both). Kill line 1.10x end to end: PASS (1.660x).
- The CPU figure is the one caveat: ship ran alone after the pair, not
  beside zone, so its slots had the machine's full memory bandwidth;
  the iteration ratio is machine-independent and is the one to quote.

Production candidate after all seven arms: zone + warmy + ship. The two
are independent levers (warmy strengthens the pool, ship removes cold
starts at unit boundaries) and both were measured against zone; their
combination has not been run. Against base on this subtree, zone + warmy
is 2.76x in iterations and zone + ship 2.71x; if they multiply, the pair
is near 4.5x, which is the next and last pair to run.

Next (cloud box, Nov 5): `run_pair.sh full ship2`.
- `full` is zone + warmy's phase settings + ship. It is judged against ship
  and against warmy (kill line 1.1x end to end against each, so the second
  lever must add something on top of the first), and reported against zone
  and base.
- `ship2`, beside it, repeats ship under the same load. That gives the CPU
  per iteration with the parent state carried, which the ship run alone
  could not.
- Both arms were smoke-tested on nug20 (every node warm in both). Each arm
  needs about 1-2 GB of disk for its store.

### full and ship2 (laptop, Nov 5; `runs/full`, `runs/ship2`)

Both at f763734, 3 slots each, run together (06:35-07:06 PST). **full
passes both gates: the levers stack.** ship2 repeats ship to the iteration
and gives ship a fair CPU figure.

| arm | nodes / fathomed | cut phases (fathomed) | warm nodes | DNN iterations (first + phase) | DNN CPU | wall (3 slots) |
|---|---|---|---|---|---|---|
| zone | 366 / 358 | 8 (0) | 248 / 365 | 51,150 (40,650 + 10,500) | 7,985 s | 3,192 s |
| warmy | 232 / 227 | 6 (1) | 163 / 231 | 30,260 (24,910 + 5,350) | 4,811 s | 2,036 s |
| ship | 322 / 315 | 7 (0) | 321 / 321 | 30,810 (22,710 + 8,100) | 3,769 s | 1,473 s |
| ship2 (beside full) | 322 / 315 | 7 (0) | 321 / 321 | 30,820 (22,720 + 8,100) | 4,562 s | 1,847 s |
| **full** = zone + warmy + ship | 232 / 227 | 6 (1) | 231 / 231 | **19,640** (14,840 + 4,800) | 3,239 s | 1,405 s |

full:
- Against warmy 1.541x in iterations (gate 1.10x), against ship 1.569x
  (gate 1.10x): PASS on both. Against zone 2.604x, against base **4.259x**
  (83,640 -> 19,640; DNN CPU 11,318 -> 3,239 s, 3.49x; wall 3.15x).
- The tree is warmy's (232 nodes, the same two depth-5 flips, every node
  warm). Against warmy the gain is per node: 1.130x on the 164 matched
  nodes and the 68 cold nodes made warm; depth 4 goes 12,870 -> 6,560
  iterations, depth 5 15,890 -> 11,580. Against ship the gain is the
  pool: the two depth-5 flips remove 90 nodes, and the 6 shared phases
  gain +36.9 over the first solve against +17.1.
- Levers multiply to within noise: warmy 1.690x times ship 1.660x
  predicts 2.81x over zone, measured 2.60x; the overlap is the depth-4
  first solves, which both levers shorten.

ship2: 30,820 iterations against ship's 30,810, same 322 nodes, same 7
phases. Run beside full rather than alone its CPU is 4,562 s, so ship's
fair figure against zone is 1.75x in CPU (7,985 -> 4,562) for 1.66x in
iterations, not the 2.1x of the solo run; the 0.826x between ship and
ship2 is the machine's load, as the cloud's check said.

Production configuration after the eight arms: `QAP_DNN_CUT_ZONE=0.005
QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200
QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=32` on top of
base. On this subtree it is 4.26x fewer DNN iterations than the Nov 4
production configuration; the caveat that it is one subtree and that the
49:45 flips carry 44-134 nodes each stands, so the per-depth figures
(depth 4 3.5x, depth 5 2.7x against base) are the ones to carry to the
wil50 estimate.

### fast vs full (cloud box, Nov 5 evening; `runs_cloud/full`, `runs_cloud/fast`, `runs_cloud/summary.txt`)

`fast` is `full` plus the two Nov 5 speed switches that passed their tests
(`../speed/`): `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7` and
`DNN_R32_NONCERT=1`. Both arms ran at the same time on the cloud box,
with 2 slots each, `DNN_AMX=1`, and binary 1ebffcb. The pair was killed at
19:57 UTC (a detached launch did not survive the end of a session turn) and
resumed at 20:45 from its ledger. The iterations and CPU are unaffected;
the wall and busy columns include the gap.

| | full | fast | full / fast |
|---|---|---|---|
| DNN nodes (fathomed) | 277 (271) | 232 (227) | 1.19x |
| DNN iterations (first solves + phases) | 23,530 (17,530 + 6,000) | 18,820 (13,390 + 5,430) | **1.250x** |
| DNN CPU s | 6,340 | 5,088 | 1.246x |
| the 232 shared nodes: iterations | 21,490 | 18,820 | **1.142x** |
| the 232 shared nodes: CPU s | 5,881 | 5,088 | 1.156x |
| ms per first-solve iteration (shared) | 268.5 | 268.7 | 1.00x |
| ms per phase iteration (shared) | 287.1 | 274.4 | 1.05x |

- **The β schedule is the gain: 1.14x fewer iterations on the matched
  nodes.** 231 of 232 nodes start warm, so this is almost all the `w=0.7`
  key (β 0.21·m/3 on warm calls). The 8-child test in
  `../speed/beta_sched/` gave 1.18x (450 -> 380 iterations); this is the
  same effect across a whole subtree.
- **The rest of the end-to-end 1.25x is one marginal node.**
  `[0:25 25:49 30:48 39:27 49:45]` at depth 5: full's phase ended at 48,813.5
  and branched into 45 depth-6 nodes (2,040 iterations); fast's ended at
  48,814.4, 0.35 above the cutoff, and fathomed. Such flips go either way
  (`base vs zone` above). Count on the matched 1.14x, not the 1.25x.
- **The fp32 rank term buys about nothing with AMX.** First-solve
  iterations cost the same in both arms (268.5 against 268.7 ms); phase
  iterations 5% less. The 4.2% of `../speed/r32_noncert/` was measured
  without AMX, where the R-step's dsyrk is a larger share of the iteration.
  It stays useful on CPUs without AMX (c7a, Graviton), which is what the
  fleet may run.
- No outcome changed on the shared nodes except that one flip; both arms
  found 48,816 and proved nothing better.

**Verdict:** `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7` joins the production
candidate (about 1.14x in iterations and CPU on this subtree). `DNN_R32_NONCERT=1`
is neutral with AMX and 4% without; it is safe to include, since it changes no
outcome and no certificate.
