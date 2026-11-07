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

So a week on 1,000 cores (6.05e8 core-seconds) needs about **41x** at the
center and **17-100x** across the range, after the symmetry file and the cut
phase (corrected Nov 3: this line said 12x and 3-45x, which do not follow
from the range above).

## 2. What the three angles offer (factors on the total)

| angle | item | factor | effort | first test |
|---|---|---|---|---|
| tree | symmetry file for wil50 (transposed instance) | 3.98x against today's runs, already in the estimate | **done** (`sym/README.md`: `latest/SYM/wil50t.sym`, `qaplib/wil50t.dat`) | passed: 15 and 635 children at depths 1 and 2 (651 nodes against 2,501 plain); node 10's bound agrees to 0.02 after a cut round with `-side 1` |
| tree | branch on the DNN marginal line (rule 9) at every DNN depth, DNN from depth 1 | **measured on the node-10 subtree: 1.41x less CPU, 1.13x fewer nodes** (rule 8 gives the same tree at 1.23x; `stronger_bounds/README.md` step 12); the 2-5x estimate was too high because rule 2 already finds decent lines | done | the node-3 subtree would confirm; DNN from depth 1 untested |
| bound | children inherit the parent's cut pool and multipliers (and a phase at clear misses to equip their children) | 1.5-4x on the dominant class (bound memo); 1.1-1.3x (tree memo). **Gate passed Nov 3** (Astra, `research/pool-inheritance`): one frozen solve on the inherited pool gains a median +54 over the child's own uncut value, d1 and s5 fathom from it; continuing with Polyak rounds adds only +10 to +21 before the quiet round, at 3.3x the iterations, so the frozen solve is the policy. **Measured in the search Nov 4: 1.37x end to end on top of the stacked arm** (`pool_search/README.md`): the same 322-node tree in 23,905 CPU-s against 32,740; all 228 children that inherited a pool fathomed on their first solve (101 iterations on average), cut phases fell from about 14,400 to 7,000 CPU-s | done (`QAP_DNN_POOL=1`, the pool travels in state files) | production: `QAP_DNN_POOL=1` |
| bound | children start from the parent's final Y as well as its Z (`warm_y/README.md`) | **measured Nov 3 on 16 depth-4 children: 1.17x per child** (1.29x on fathoms, 24% fewer iterations; +67 at iteration 100 on misses; no outcome changes) | done; in the search as `QAP_DNN_WARM_Y=2` | in the stacked node-10 A/B (`run_ab6.sh`) |
| speed | AMX bf16 for the filter and Rayleigh-Ritz products, plus CPU trims | **measured Nov 3 on the 32-node sample: 1.71x on fathoms, 1.75x on misses, same 18 fathoms, same iterations** (bundle D: AMX kernel for the five W*Y products, degree 3 with tight lambda, one-pass CholQR, refresh every 40; `amx/README.md`). Closed as not worth it: guard K/40, fp32 R (certificate margin), the small products on AMX, the row fusion; neutral: the symmetrisation skip. Left: the kernel's last third (ceiling 2 TF against 1 TF reached) | done | production: `DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40` |
| speed | cut-phase trims: shorter inner solves with Y/Z restoration (Astra's Y12s: `-warmY`, 12 rounds x 200 iterations, quiet cap 1,000) | measured Nov 3 on the phase: 1.56x. Nov 4: the search's rounds restarted Y from e0e0', the cause of most of the round-end lag (`convergence2/REPORT.md`, F1). **Measured Nov 5 in the search** (`pipeline_ab/`, warmy against zone on the laptop, with the first step of an inherited phase damped to 0.3): **1.69x fewer DNN iterations end to end** (51,150 -> 30,260, tree 366 -> 232), 1.16x on the matched nodes, its phases gain +40 over the first solve against +17 | built (`QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3`); Nathan's Nov 3 call reversed Nov 4 | with ship: the `full` arm |
| speed | coordinator instead of the checkpoint barrier and node-end steals | 1.15-1.35x of the fleet [e]; **built Nov 4** (`coord/`, plan.md phase 1, local and AWS backends, exact totals under worker kills, coordinator restarts and drains); on the 4-slot node-10 subtree the same tree with **1.09x fewer DNN iterations** (warm starts kept: 23 units, no barrier) and the same 81% busy (the subtree's tail) (`coord_search/README.md`) | phase 2 of plan.md (Terraform, Spot handler) to run it on a fleet | the fleet factor needs plan.md phase 2 (100 Spot cores) |
| tree | confirm give-ups uncut before the phase (QAP_DNN_GU_ZONE was off in run_ab3), primal gate in the phase, zone 0.0035 | **replayed, about 1.1x; closed** (`trims/README.md`): the deferral has nothing to save (the extrapolation rule already lets slow uncut fathoms finish), the primal gate 4:20 is in the code (default off, no false stop on 88 logged phases) and worth about 6%, the zone a modelled 4%; the search A/B was abandoned as not worth four hours | done | switch `QAP_DNN_CUT_PGATE=4:20` on in the production configuration; no further test |
| bound | skew-symmetric certificate polishing at non-converged iterates (the Hahn-Grant RLT1 dual ascent on the Lagrangian matrix L + Z') | measured Nov 4 per node (`certificate/REPORT.md`: 19/19 fathoming nodes cross 1-4 checks earlier at equal iterates), then **in the search Nov 5** (`pipeline_ab/`, certx against zone on the laptop): **0.983x on the 321 matched nodes**, 1.041x end to end from one flipped near miss. The rest of the certificate is closed (margins 4e-6 units, free multipliers neutral or blocked, child bounds flat at convergence) | built, opt-in (`DNN_CERTX`); **closed** | - |
| bound | equip phase: a short cut phase at out-of-zone misses that will branch, so their children inherit a pool (`convergence2/REPORT.md` #1) | **measured Nov 4-5 in the search** (`pipeline_ab/README.md`, laptop, base and zone at the same time): `QAP_DNN_CUT_ZONE=0.005` gives node 10 a 1,200-iteration phase, and its 47 children need 1.65x fewer iterations (43 first-solve fathoms against 38, 4 phases against 9). End to end **1.64x** in DNN iterations against the same-hour base (1.45x on the shared tree), **1.15x** against the luckier cloud base run of the same configuration. Single runs vary by about 1.4x (marginal nodes, cold starts at unit boundaries); its phases from inherited pools still gain little (+7.7 on average) | done; **adopted** (`QAP_DNN_CUT_ZONE=0.005`) | warm-Y with first-step damping on top (pipeline_ab `warmy`) |
| speed | keep the parent's DNN state (Y, Z) when its children cross a unit boundary (`convergence2/REPORT.md` #2) | **measured Nov 5 in the search** (`pipeline_ab/`, ship against zone on the laptop): every node warm (321 of 321 against 248 of 365), **1.66x fewer DNN iterations end to end** (51,150 -> 30,810), node 10's children 1.92x; 0.978x on the matched nodes, so the whole effect is the cold starts. bf16 = fp32 = exact on 8 wil50 children; fp32 is the tested setting. nug30 on the laptop: 15% off the Nov 3 record | built (`QAP_DNN_SHIP=32`, state files v6, `make shipcheck`) | with warm-Y: the `full` arm |
| bound | regional M2 block as a Lagrangian decomposition (Astra's 10b) | +20 units (10-15% chance of +65); **downgraded Nov 4 to about 0** (`new_cuts/SYNTHESIS.md`: after the triangle phase every single threshold is already saturated, so per-region blocks have nothing left) | 2-3 weeks | not worth running |
| bound | new cut families, three research agents (literature, grid, facility; `new_cuts/SYNTHESIS.md`) | **none meets +40; no path to +130 through linear inequalities on Y.** The residual gap is the global coupling of the middle thresholds (117-143 at depth 3), which only branching resolves. The one side-bet, large-Q box inequalities in the pool (+10 to +25 estimated), **failed its gate Nov 4: -32.3 on node 10, -22.6 on node 13** (`new_cuts/box_pool/`): no Lagrangian gain inside the Polyak rounds, a wider certificate lag | done | closed |
| speed | GPU port, 8-16 nodes per device | wall clock: 100 H100s = about 7x the 1,000 cores | 3-6 weeks | one-day CuPy prototype; dropped above 10 ms per node-iteration on an A100 |

Closed this week, do not revisit without new evidence: 5-point counts, triple
consistency, the dual method (bundle), the local-7 family as a pool cut type,
Anderson/Halpern, momentum, atomic triangles, eliminations from the
certificate. Another full level from new inequalities at depths 3-5 is rated
under 10% within a month (`bound_strength.md`). Closed Nov 4 by the three-agent study (`new_cuts/`): every facility-side family, the location-side geometric families (parity, dispersion, free regions, rectangles, l1/hypermetric), the row/column decomposition, objective splitting with a coarse optimum, the three local RLT2 projections, RLT2/RLT3/Lasserre-2, disjunctive and lift-and-project cuts, and the per-region moment block; and, after its own test, the large-Q box inequalities (`new_cuts/box_pool/`).

## 3. Multiplying out (revised Nov 3, evening)

The Nov 2 centers were rule-9 branching 3x, inheritance 2x, per-node speed
2x, phase trims and fleet efficiency 1.3x: about 15x. What has been measured
since (`[m]`) and what is still an estimate (`[e]`):

| lever | factor | status |
|---|---|---|
| rule 9 branching | 1.41x [m] | end to end on the node-10 subtree (`stronger_bounds/README.md` step 12) |
| per-node speed, bundle D | 1.71x [m] | per node on the 32-node sample, same decisions (`amx/README.md`) |
| parent-Y child warm start | 1.17x [m] | per child on 16 depth-4 children (`warm_y/README.md`); 25% on the laptop's full nug30 search with one-pass CholQR and two give-up settings |
| **all three in the search** | **2.34x [m]** | **end to end on the node-10 subtree** (`stronger_bounds/README.md` step 13): the same tree as rule 9 at 32,740 CPU-s against 76,755 for the baseline configuration; the per-node levers give 1.67x in the search against 2.0x measured separately (lost warm starts at checkpoints and restarts, and cut phases not measured per node) |
| pool inheritance in the search | 1.37x [m] | end to end on the node-10 subtree on top of the stacked arm (`pool_search/README.md`): same tree, 23,905 against 32,740 CPU-s. **Same-host check through the coordinator: 1.40x on the 322 shared nodes, 1.56x end to end** (without the pool two marginal nodes branch into 89 more) |
| **all four in the search** | **3.21x [m]** | **against the baseline configuration's 76,755 CPU-s on the same subtree (5.8x against the plain search's 139,149)** |
| coordinator (fleet busy 70-85% -> 93-97%) | 1.15-1.35x [e] of wall clock; 1.09x [m] in iterations on the 4-slot subtree | built (`coord/`); the fleet factor is measured at plan.md phase 2 |
| shorter cut phase (Y12s), phase trims | 1.1-1.2x each [m on the phase] | not pursued; Y12s is now the switch `QAP_DNN_CUT_WARMY=1` with 200-iteration rounds |
| Nov 4-5 levers (`convergence2/`, `pipeline_ab/`): equip phase, warm-Y rounds with damping, parent state across units | in DNN iterations on the node-10 subtree, laptop: zone 1.64x over the Nov 4 production configuration; on top of zone, warm-Y 1.69x and ship 1.66x. **All three together (`full`): 4.26x in iterations and 3.49x in DNN CPU over the same-hour base** (83,640 -> 19,640 iterations; 11,318 -> 3,239 s); 3.0x in iterations over the cloud's run of the Nov 4 configuration. The levers multiply to within 8% (1.69 x 1.66 = 2.81x predicted over zone, 2.60x measured). Per depth against base: depth 4 3.5x, depth 5 2.7x | production: `QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=32` on top of the Nov 4 setting; certificate polish closed |
| Nov 5 speed agents (`speed/`): two-phase β, fp32 rank term off certificate iterations, restricted PRSM pieces, complementary-subspace Rayleigh-Ritz | **β schedule: 1.14x fewer DNN iterations on 232 matched nodes of the node-10 subtree on top of `full`** (1.25x end to end, the rest one marginal node; `pipeline_ab/` arm `fast`, cloud box); the gain is the warm-call key `w=0.7`. fp32 rank term: 4.2% per iteration without AMX, about 1% with it. The restricted PRSM pieces (-4% CPU) and complementary-subspace Rayleigh-Ritz (7 ms of 356 per iteration) were killed | production candidate: `full` + `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7` (+ `DNN_R32_NONCERT=1` on CPUs without AMX) |
| Nov 6 levers (`speed/nov6_admm/`, `speed/nov6_core/`, `pipeline_ab/`): warm-call β 0.25 for children with cut rounds kept at 0.7, the certificate at the predicted cutoff crossing, three core switches | on the node-10 subtree, laptop: **warm-call β (`wlowc`, after the Nov 6 `wc` fix): 1.24x fewer DNN iterations and 1.12x less DNN CPU on 232 matched nodes** (1.155x and 1.057x end to end, one marginal phase lost its fathom; depth-5 first solves 1.60x; each iteration 10-15% dearer at the lower β); crossing certificate 1.035x in iterations and 1.04x in DNN CPU; core switches 1.08-1.10x in DNN CPU on the same tree. **Stacked on the fixed binary (`bestcore` against `fast`): 1.23x less DNN CPU on the 232 matched nodes, 1.16x end to end**; on node 5's tail the warm β gives 1.32x per matched level-5 node with an identical tree. The two-phase warm β (`wtwo`) ties `wlowc` per node and kept node 10's marginal fathom; it is the fallback | production candidate (`bestcore`): `fast` with `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1` |
| Nov 7: CP-ADMM, the cut pool inside the PRSM's Y-step by Dykstra sweeps (`admm/cpadmm/`; memo item 2 of `admm/alternatives.md`) | single nodes, cloud box: **node 13 certifies 48,753 at 1,070 iterations after the give-up point with cuts separated inside the solve every 100 iterations, against 4,860 for the Polyak protocol (4.5x in iterations and in CPU: 589 s against 2,659)**; on a frozen 13,000-cut pool 1,040; with only the 1,000 cuts of the give-up point it still reaches the Polyak protocol's 12-round value by 770 against 3,600. Node 14 fathoms at 790 against 4,390. Production settings: 48,753 at 1,210 against 4,290. +4-14% per iteration. d2 does not fathom: 48,808.0 at 4,500, 6.1 short (the regional triangles' cap there). Estimate if the phases gain alike in the search: about 1.3x on the node-10 subtree [e] | opt-in (`DNNParams::ycuts`, `dnncuts -cpadmm/-cpsep`), not in bnbdrv; validity by brute force (19,138 certified values) and the switch-off identity; next: a CP cut phase in the search (pooled children from mu_D, per-check stopping rules, a lazy certificate), then the node-10 A/B against `bestcore` |
| Nov 7: the CP cut phase in the search (`QAP_DNN_CUT_CP=k:K`; `admm/cpadmm/` "In the search", `pipeline_ab/runs_cp/`) | node-10 subtree, cloud box, both arms at once with 2 slots each: **`cpcut` (bestcore + `QAP_DNN_CUT_CP=100:1000`) against `bestcore`: 1.307x in DNN CPU on the 94 matched nodes (kill line 1.10x); end to end 1.86x in DNN CPU (4,821 s to 2,586), 2.03x in DNN iterations (15,941 to 7,868) and 1.77x in wall time, with the tree going from 277 nodes to 94**. Three depth-4 phases fathom. Under bestcore, three phases from an inherited pool gain nothing in 600 iterations, and the CP phases at those nodes end 18-23 higher. One depth-5 node fathoms in its first solve on the CP pool (bestcore: 0.8 short after a 1,600-iteration phase). +13% per phase iteration. `cpcut2` (+ `QAP_DNN_POOL_CP=1`, CP first solves for pooled children): a further 1.05x in DNN CPU (1.089x in iterations) in a same-session pair. `cpcutl` (+ `DNN_CP_LAZYCERT=1`): the same tree, about 2% less CPU | opt-in, all three off by default; validity by brute force through bnbdrv (17 instances, 97 phases, four switch sets, no bound above its node's optimum) and the switch-off identity (nug20; a node-10 slice at bestcore settings); gates pass on the same solver tree (`admm/cpadmm/logs/gates_search.txt`). One subtree, and the end-to-end factor rests on 4 changed outcomes. Next: the laptop pair (3 slots), then cpcut against bestcore on a second subtree (the node-5 tail probe) |
| Nov 7: the CP projector on side 1, the production configuration (`wil50t.dat`, `QAP_DNN_CUT_SIDE=1`; `admm/cpadmm/` "Side 1", `pipeline_ab/runs_cpt/`, `tail_probe/runs_node5_cpt/`) | Until this change, `cpcut` ran the rounds in the production configuration, because the projector was side 0 only. Side 1 costs +0.6% per iteration on node 13, and node 13 certifies 48,753 at 1,070 iterations on both sides. **Node-10 subtree on `wil50t.dat`, both arms at once with 2 slots each: `cpcutt` against `bestcoret` 1.760x in DNN CPU on the 94 matched nodes, 2.20x end to end (4,446 s to 2,021), 2.37x in DNN iterations and 1.96x in wall time; the tree goes from 231 nodes to 94.** `bestcoret`'s tree is 231 nodes against `bestcore`'s 277 on `wil50.dat`, because node 10 branches on a different line in each orientation (rule 9 is not invariant under transposition). Iterations agree on the shared nodes. **Second subtree, node 5 of the tail probe (matched probe, estimated): 1.15x less CPU and 1.21x fewer iterations, inside one standard error.** At level 4 the CP phases cost more (0.71x in CPU on 47 matched nodes) and fathom 2 more children; level 5 is 1.37x | opt-in as before. Validity: the self-test, with side 1 bit-identical to side 0 on the transposed iterate; the transpose check (certified values within 0.19 from t = 300, less than the plain PRSM's own difference); brute force through bnbdrv with side-1 phases forced; the switch-off identity; gates pass on the same solver tree (`admm/cpadmm/logs/gates_side1.txt`). Atomic cuts still run the rounds (bnbdrv never sets them). **Decided Nov 7: `cpcutt` is the full run's setting** (`fullrun/settings.sh`, the same environment as the measured arm); laptop replicate optional |

Against the baseline of section 1 (120-700 days on 1,000 Xeon cores, center
290; a week needs about 41x at the center):

- **Measured end to end (3.21x, Nov 4, with inheritance):** center about
  **90 days**, range about 40-220 (it was 125 days at 2.34x before
  inheritance was measured).
- **Plus the coordinator at 1.25x (about 4.0x):** center about **72 days**,
  range about 30-175. A week is reached only at the low end of the tree
  range. (The Nov 3 projection assumed inheritance at 2x; it measured 1.37x.)
- **Plus the equip phase (Nov 5, 1.15-1.64x in DNN iterations,
  `pipeline_ab/`):** center about **55-80 days** on the measured levers,
  about 45-63 with the coordinator.
- **Plus warm-Y rounds or the shipped parent state (Nov 5, each about 1.9x
  over the cloud base and 2.7x over the laptop base together with zone):**
  center about **33-48 days** on the measured levers, about 26-38 with the
  coordinator.
- **Plus all three together (`full`, Nov 5: 3.0-4.3x in iterations, 3.5x in
  CPU over the Nov 4 configuration):** center about **21-30 days** on the
  measured levers, about **17-24 with the coordinator**.
- **Plus the β schedule (Nov 5, 1.14x on matched nodes in the search):**
  center about **18-26 days** on 1,000 Xeon cores, about **15-21 with the
  coordinator**. Measured on node 10 only; most of the search's nodes start
  warm, as there.
- **The tail class, measured Nov 5 (`../tail_probe/`, matched probe of node
  1, 405 short, on the laptop): the levers hold there.** The production
  candidate (`full`) against the Nov 4 setting: 4.9x in estimated DNN
  iterations and 3.1x in CPU on node 1's subtree (node 10: 4.26x and
  3.49x), 2.68x per node on 87 matched level-5 nodes (83 fathomed against
  79). Node 1's subtree is about 27,000 DNN nodes under the Nov 4 setting
  and 14,000 under `full` (SE about 70%; level 7, 84% of the estimate,
  rests on 6 bounded nodes). In Xeon CPU (2.0x the laptop per core) that is
  within 1.3-1.8x of what this section's centers assume for a tail node,
  well inside the probe's own error: **the 21-30-day center stands**, now
  with the tail class measured on one node instead of assumed. The 329-short
  node of the class is not probed.
- **The caveat that remains.** Every lever was measured on two subtrees,
  node 10 and node 1. The total still rests on the tree model's count of
  tail nodes per orbit, which only the sampling plan measures. The zone
  factor is measured in iterations; its CPU ratio on the laptop was 1.42x against 1.64x, and part of that gap
  may be the laptop's mixed cores.
- **Per core, an M2 Pro is 2.0x the cloud's AMX Xeon** on identical DNN
  nodes (`pipeline_ab/README.md`). This is not in the estimates, which are
  in Xeon-core days.
- **100 H100-class GPUs in place of the cores** (about 7x the wall-clock
  throughput at similar dollars, `node_speed.md`): center about 13 days on
  the measured levers alone, about 10 days with the coordinator.

**A week on 1,000 Xeon-class cores is not reachable at the center of the tree
range with the levers we have.** It comes into reach if the tree is at the
small end of its range (the 16 sampled depth-3 nodes branch on random lines,
the production search on rule-9 lines that fathom far more children), or on
GPUs. Caveat on the arithmetic: the 3.21x is measured end to end in the
search, but on one subtree (node 10); how much of the full tree's work sits
in the children of cut-phase misses, where inheritance acts, is for the
sampling plan to show.

The two measurements that decide between these outcomes are cheap:

1. **The sampling plan** (`tree_size.md` section 6): build the production
   frontier at depths 0-2 (1.6e5 CPU-s), take 400 frontier nodes to the
   give-up point, run stratified subtrees end to end, dive into the ten
   worst tail nodes, and run a cutoff ladder. About 1e7-2e7 CPU-s, 3-6 hours
   on 1,000 cores, 2-3% of the week, and every subtree finished is a piece of
   the proof. It pins the total within 1.5x instead of 15x.
2. **The subtree A/Bs** on the cloud box: rule 9 is done (1.41x); the
   stacked arm is done (2.34x end to end); inheritance in the search is
   done (1.37x on top, 3.21x in all, Nov 4).

## 4. Order of work

1. ~~Symmetry file for wil50~~ done Nov 3 (`sym/README.md`): every wil50 run
   from here uses `qaplib/wil50t.dat` with `QAP_DNN_CUT_SIDE=1` and `-F` specs
   as j:i; the node-10 subtree A/Bs are unaffected (no symmetry below node 10).
2. ~~Inheritance test and the rule-9 subtree A/B~~ done Nov 3: rule 9 1.41x;
   frozen-pool inheritance passed its gate; the stacked node-10 A/B gives
   2.34x end to end (Nov 4, step 13); inheritance in the search adds 1.37x
   on top (3.21x, Nov 4, `pool_search/`).
3. ~~Confirm give-ups uncut, the primal gate, zone 0.0035~~ closed Nov 3 by
   replay (`trims/README.md`, about 1.1x; the gate ships default off). The
   Y12s phase budget was not pursued.
4. ~~AMX and the CPU trims~~ done Nov 3 (1.71x per node); ~~the coordinator
   from docs/distributed/plan.md~~ built Nov 4 (phase 1, `coord/`; on the
   node-10 subtree the same tree with 1.09x fewer iterations,
   `coord_search/`). Next for it: plan.md phases 2-3 (Terraform, the Spot
   handler, 100 cores with induced interruptions) once an AWS account is
   connected.
5. The Nov 4 per-node candidates (`convergence2/REPORT.md`,
   `certificate/REPORT.md`), cheapest first, each a node-10 subtree run
   through the coordinator against the 18,496 CPU-s pool arm: the equip
   phase (`QAP_DNN_CUT_ZONE=0.005`, no code), warm-Y rounds with it
   (`QAP_DNN_CUT_WARMY=1`, 200-iteration rounds, first-step damping),
   the certificate polish; then the parent state in state files (1-2
   days) if the first two hold.
6. The sampling plan on the first 1,000 cores, before the full run.
7. In parallel, gated: the M2 projection test and the GPU prototype.
8. Laptop items (Nov 5), cheapest settlement first:
   - **The 329-short tail node**: the matched probe of `../tail_probe/`
     (base against full) on node 5 of `../tuning/giveup/all32.txt`
     (`23:49,37:15,6:5`, 316 short converged, 329 at 300 iterations), the
     one tail-class node not yet probed; checks that node 1 stands for the
     class. **Done Nov 5 evening** (`../tail_probe/README.md`, "Matched
     results, node 5"): the levers hold (1.82x CPU, 4.04x per matched
     level-5 node), but node 5's subtree is 21x smaller than node 1's under
     the same setting (1,290 against 27,000 DNN nodes), in line with the
     per-level model for an 89-unit difference in converged shortfall. Node
     1 is the worst of the class, not its representative; the tail weight
     has to be spread over converged shortfall, which the sampling plan
     measures.
   - **RLT3 violation sample**: before building any RLT3 cut type, sample a
     few hundred random RLT3 inequalities (products of three assignment
     constraints) at the converged iterates of nodes 10, 13 and d2 and
     report the violation distribution, as step 8 of
     `../stronger_bounds/README.md` did for the RLT2 triples (median about
     0.01, a few units of gain). A median in that range closes selective
     RLT3 on wil50 without a cut type; only a median several times the
     regional triangles' violations (0.1-0.2) would justify one. Prompted by
     Jeff Linderoth's RLT-n plus PDLP suggestion (Nov 5 thread).
     **Tabled Nov 6** (`../stronger_bounds/rlt3/README.md`): the quadruple
     separation LP has 4.3 M placement rows and no valid shortcut (merging
     locations or dropping zero-mass potentials loses the violation), so one
     quadruple costs 10+ minutes on the laptop; a cloud job if picked up.
     By-product: node 10's first triple sample (30 triples, median 0.0175,
     max 0.049), the same order as nodes 13 and d2.
   - **sko49 or wil50, settled Nov 6** (`../../sko49_feasibility/vs_wil50.md`):
     matched probes on two sko49 tail nodes; sko49 nodes cost 2-4x a wil50
     node at equal shortfall, cancelling its narrower tree: about 15 against
     18 days on 1,000 cores, sko49's range four times tighter. Stay on
     wil50; switch only if the sampling plan puts wil50 above about 1.5e9
     laptop CPU s.
