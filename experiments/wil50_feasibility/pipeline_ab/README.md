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
  marginal nodes fathom and which nodes start cold. The end-to-end factor
  of a single pair is therefore uncertain by about 1.4x. Zone passes its
  kill line even against the luckier base (1.15x). The 1.05x and 1.03x
  lines of warmy and certx are below that noise, so `summarize.py` judges
  them on the matched nodes (shared, same start) and, for warmy, on the
  phases both arms ran. The replicates `zone2` and `base2` measure the noise
  on one machine.
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
