## Follow-up (Nov 16): sizing tai50b's tree with `cut.regions`

Stage samples in the search (`size50/`: `par/`, `runs/TAG/{cmd,nodes,dnn,log}.txt`, `root.sh`,
`kid.sh`, `pass1.sh`, `cutset.py`, `collect.py`, the queues `q*.sh`), warm from the root with
its pool, one solver process at a time on the shared box. The Knuth dives and the root's
plateau run elsewhere (coordinator's scope change); the estimate here is the cross-check from
the samples. [m] measured (run named), [e] estimated.

### What the samples say

1. **The top of the tree is one path, as on tai30b.** The root branches on location 31, then
   33, then 32: the three locations of the isolated 3-point cluster. At depths 2 and 3,
   below the node on the path to the best known solution, **20 of 20 random children are
   fathomable at once**: 11 fathom in the pooled first solve (30-229 iterations), the other 9
   end it 0.15M-2.2M short, and the furthest tested fathom after 144-227 phase iterations.
   Only the on-path child survives each level.
2. **From depth 4 it is not one path.** The fourth branching location (34) lies in a
   16-point cluster. Of 10 random children there, 2 fathom in the first solve and 8 end it
   1.5M-3.2M short (the on-path child 2.9M: no different from its siblings). A child whose
   facility belongs to *another* cluster in the best known solution fathoms after 200 phase
   iterations (`34:36`, 2.36M short); one whose facility belongs to the *same* cluster does
   not (`34:28`: 0.55M short after 1,000 phase iterations and flattening). In the best
   known solution 15 other facilities share that cluster, so about 15 of the 46 children are
   of the surviving kind [e].
3. **A survivor's children die quickly**: under `34:28` (branched 0.55M short), 9 of 10 random
   depth-5 children fathom in the first solve (40-91 iterations, 9-21 CPU-s); the tenth,
   again a same-cluster facility (`42:35`), ends it 0.47M short and **does not fathom either**:
   106k short after 1,000 phase iterations, gaining 15k per 100. So the same-cluster fringe
   compounds for at least one more level, with gaps falling about 5x per level (3.1M at
   depth 3, 0.55M at depth 4, 0.11M at depth 5).
4. **The on-path nodes close slowly**: 4.82M short at depth 1, 3.88M at depth 2, 3.10M at
   depth 3 when the stall rule stops their phases after 600-760 iterations (0.68-1.05% of
   the best known value; on tai30b the same nodes were within 0.025%). **Run long they close much further**: the on-path depth-4 node with 3,000 phase
   iterations, `wc=2.88` and no stall rule reaches 458,412,325, 0.41M (0.089%) under the
   best known value (`runs/d4_34_6path_long`); its stalled parent was 3.10M short. And its children
   are then far easier: of 10 random depth-5 children of that node, **8 fathom in the first
   solve** (45-141 iterations), 2 end 0.74M short, against 2 of 10 and 1.5M-3.2M below its
   stalled parent one level up (`runs/d5p1_under_longpath`, `runs/d4p1_under_path`).
5. **The phase penalty matters here too** (the coordinator's lead, `../cpfast`): on the median
   depth-1 child `31:23` (5.54M short at its first solve) the phase needs
   **589 iterations with `dnn.beta_sched = w=0.36,wc=2.88`, 743 with `wc=1.0`, 640 with
   `wc=5.76`, and more than 1,000 without `wc`** (0.36M short at 1,000, about 1,300 by its
   trend [e]): 2.2x. Every phase number in the first report (sections 5-10) and the runs
   here without "wc" in their name predate it.

### Estimate (cross-check from the samples; the dives are the number to use)

- **The part I measured, to depth 3: about 55,000 iterations, 7 core-hours** [e from m]:
  root 3,000; depth 1, 50 children at 155-1,110 iterations each (about 700 mean with
  `wc=2.88`), the far ones branching into children that fathom in about 35 iterations each;
  depths 2-3, 48 children at about 150 iterations plus the on-path node's 700.
- **Below depth 3 the tree is a path with a fringe of same-cluster swaps.** Per on-path
  level, with k surviving siblings (k about 15 at depth 4 [e], from 1 phase and 4 of 10
  first solves) each costing its own phase (1,100 iterations) and 46 children (9 in 10 at
  about 55 iterations, 1 in 10 a phase or one more level): about 4,000-9,000 iterations per
  survivor, so 60,000-140,000 iterations per on-path level at depth 4, less below as the
  clusters fill and m falls (0.30 s per iteration at m = 46, 0.15 at m = 40 [e], 0.04 at
  m = 30 [m, tai30b]).
- **If every survivor's children died at once: about 1e6-3e6 iterations, 10,000-60,000 nodes,
  50-200 core-hours on this box** [e]. The samples say this is too low: one depth-5 child in
  ten survives again (1 of 2 same-cluster ones sampled).
- **With the fringe compounding** at rho surviving children per survivor for L more levels,
  the cost multiplies by about rho^L. The samples give rho between 0.5 and 20 (1 of 10;
  point value 4.5, i.e. half of the same-cluster children) and gaps falling 5x per level,
  which should end the compounding within 2 to 4 levels [e]. rho = 3, L = 3: x27, **about
  1,500-5,000 core-hours here (2,000-7,000 c7a core-hours), of the order of 1e6 nodes**: my
  central reading. rho = 2, L = 2: 200-800 core-hours. rho = 5, L = 4: 3e4-1e5 core-hours.
- So: **a run of thousands of core-hours, with a range from several hundred to several
  tens of thousands**; hundreds of thousands only if the same-cluster fringe keeps
  compounding below depth 6, where I have no sample. For comparison the current figures are
  4.1e9 CPU-s (1.1e6 core-hours; tuned, 10 dives) to 1e11. The dives are the number to use:
  everything here hangs on rho and L, measured on 12 nodes.
- Two things should pull it down, neither measured in a tree: on-path phases run to
  convergence instead of to the stall rule (0.41M short instead of 3.1M at depth 4), so that
  their children start from a tighter state (8 of 10 children then fathom in the first solve, not 2 of 10); and a branching order that fills the
  small clusters first and leaves near-equivalent placements to the end.
- What it rests on: one root (2,970 iterations, not converged), 14 depth-1 first solves and
  4 depth-1 phases, 21 + 11 + 11 + 21 first solves at depths 2-5, 8 phases at depths 2-5,
  all below one branching sequence, most with phases capped at 1,000 and the stall rule on. The depth-4 survivor count is the weakest number.

### Stage tables [m]

**Root** (`runs/rootL`, `par/rootL.par`): first solve "hopeless" at 470 iterations
(417,845,879), phase 2,500 iterations at the cold penalty (no `w`): **450,196,995 certified**
(8,624,521 short), not converged; 2,042 CPU-s. (The first report's `dnncuts` root reached
449,469,933 at 1,900 iterations; neither is a plateau. The other container reports that.)

**Depth 1** (children `31:x`, warm, pooled). First solves of 13 random children + the on-path
one (`runs/d1p1`; stopped by hand after 14 of 21, none was fathoming): all give up at 110
iterations (49-59 CPU-s), shortfalls 1.44M, 4.91M, 5.30M, 5.54M, 5.62M, 5.76M, 5.83M, 6.02M,
6.22M, 6.49M, 6.70M, 7.29M, 7.49M; on-path `31:7` 8.36M. Phases:

| child | first solve: short | settings | phase iterations | outcome | total iterations | CPU-s |
|---|---|---|---|---|---|---|
| `31:1` (nearest) | 1,439,949 | no wc | 45 | fathomed | 155 | 85 |
| `31:23` (median) | 5,535,206 | no wc | 1,000 (cap) | 359,412 short, gaining 164k per 100 | 1,110 | 640 |
| `31:23` | | wc=1.0 | 743 | fathomed | 853 | 449 |
| `31:23` | | **wc=2.88** | **589** | fathomed | 699 | 338 |
| `31:23` | | wc=5.76 | 640 | fathomed | 750 | 349 |
| `31:16` (furthest random) | 7,489,213 | wc=2.88 | 1,000 (cap) | 1,393,444 short, gaining 156k per 100 and slowing: branches [e] | 1,110 | 562 |
| `31:7` (on-path, incumbent + 1) | 8,362,751 | no wc | 760 | stalled 4,824,476 short | 870 | 487 |

**Depth 2.** First solves (`runs/d2p1_under_31_7`, `runs/d2p1_under_31_23`; children `33:x`):

| parent | sample | fathom in the first solve | the rest: short at give-up (110 iterations) |
|---|---|---|---|
| `31:7` (on-path, 4.82M short) | 10 random | 5 (30-85 iterations, 15-40 CPU-s) | 0.20M, 0.53M, 1.24M, 1.34M, 2.20M |
| `31:7` | on-path `33:32` | | 4.97M |
| `31:23` (median, branched 0.36M short) | 10 random | **10** (20-61 iterations, 10-28 CPU-s) | - |

Phases under `31:7` (wc=2.88): `33:28` (2.20M) fathoms after 227 iterations (337 in all, 138
CPU-s); `33:9` (1.34M) after 144 (254, 102 CPU-s); on-path `33:32` (no wc, incumbent + 1)
stalls 3,882,511 short after 600 (710, 338 CPU-s).

**Depth 3** (`runs/d3p1_under_path`; children `32:x` of `31:7,33:32`): of 10 random, 6 fathom in
the first solve (32-229 iterations), 4 end 0.15M, 0.79M, 2.02M, 2.16M short; on-path `32:0`
4.30M, phase stalls 3,095,422 short after 620 (`runs/d3_32_0path`, 730 iterations, 307 CPU-s).

**Depth 4** (`runs/d4p1_under_path`; children `34:x` of `31:7,33:32,32:0`): of 10 random, 2
fathom in the first solve (120, 133 iterations); 8 end 1.55M-3.16M short; on-path `34:6`
2.94M. Phases (wc=2.88): `34:36` (2.36M; facility of another cluster) fathoms after 200
iterations (310, 100 CPU-s); `34:28` (3.16M; same cluster) 546,053 short at the 1,000 cap,
gaining 43k per 100: branches (1,110 iterations, 337 CPU-s).

**Depth 5** (`runs/d5p1_under_34_28`; children `42:x` of the survivor `34:28`): 9 of 10 random
fathom in the first solve (40-91 iterations, 9-21 CPU-s); `42:35` (same cluster) 474,550
short. Its phase (`runs/d5_42_35`, wc=2.88): 106,016 short at the 1,000 cap (303 CPU-s).

**Depth 4 run long, and depth 5 below it** (`runs/d4_34_6path_long`, `par/kidLong.par`: wc=2.88,
`cut.cp_iters = 3000`, `cut.stall = 0`, incumbent + 1): on-path `34:6`, 2.94M short at its
first solve, 409,192 short after 3,000 phase iterations (3,110 in all, 1,012 CPU-s). Its
children `37:x` (`runs/d5p1_under_longpath`): of 10 random, 8 fathom in the first solve (45-141
iterations, 11-33 CPU-s), 2 (same-cluster facilities) end 0.74M short; on-path `37:28` 1.29M.

### Validity [m]

With the incumbent at best known + 1 (cutoff 458,821,517.05) the nodes on the path to
`QAPSOLN/tai50b.sln` never fathom and never exceed the best known value: depth 1 `31:7`
453,997,041; depth 2 `33:32` 454,939,006; depth 3 `32:0` 455,726,095; depth 4 `34:6` run to 3,000 phase
iterations 458,412,325 (0.41M under); its on-path child `37:28` 1.29M short at its first solve. Their off-path
siblings fathom around them, as above.

### Calibration

- This box, loaded (load 5-9 on 4 cores, 00:50-03:00): 0.69 s per iteration at m = 50
  (root, with a second process of mine for the first 12 minutes), 0.50-0.58 at m = 49.
- This box nearly idle (load 2, from 03:00): **0.40-0.42 s at m = 47-48, 0.30-0.32 at
  m = 46**, 0.21-0.28 at m = 45 (first solves). My best figure for an idle core at m = 48:
  **0.40 s** [e]. `tai50b.par` has no `dnn.amx`, so these are without the AMX kernels.
- c7a: `../tho40_search` measured this box at 0.70 of c7a's CPU per matched node (loaded,
  m <= 38), so about 0.55-0.6 s per iteration at m = 48 on c7a [e]. Iterations and nodes are
  the numbers to carry over.

### For the dive runner

**Build.** Apply `regions.patch` to `cognac` HEAD 6e993cc7 (`git apply`), then
`make -C latest -j4 && make -C latest dnntest dnncuts`. It adds `cut.regions` to bnbdrv and
`-regions` to dnncuts; without it `cut.regions` is an unknown key.

**Regions file.** `regions/tai50b.comb.txt` (61 regions, 0-based location indices, the first
matrix = `cut.side = 0`). Rebuild with
`python3 -I gen_regions.py qaplib/tai50b.dat line:1,top:5,linegap:4:0.03 -o FILE`.
`cut.regions` takes a path; my par files carry an absolute one
(`/home/user/qap/experiments/bound_lab/regions/regions/tai50b.comb.txt`): change it.

**Settings block** (`size50/par/kidLong.par`; =
`experiments/estimates/tai50b/tai50b.par` with these changes):

```
dnn.depths = 0:48        # was 1:48. Pitfall 1
dnn.giveup = 100         # was 600: hand over to the phase early (it is a floor, pitfall 2)
dnn.giveup_warm = 100
dnn.gu_zone = 0          # was 0.003: no 100-iteration wait before a near-miss give-up
dnn.beta_sched = w=0.36,wc=2.88   # was w=0.36: the phase's penalty (2.2x on 31:23)
pool.on = 1
cut.regions = <path>/tai50b.comb.txt
cut.side = 0
cut.zone = 1             # a phase at every DNN miss (the plain bound is 6% off at the root)
cut.warmy = 1
cut.rounds = 12
cut.iters = 200          # with cut.stall = 3: the stall window is 600 iterations
cut.freeze = 0
cut.cp = 100:1000
cut.cp_iters = 3000      # most of my samples: 1000 with cut.stall = 3. Pitfall 3
cut.stall = 0
```

The root's par is the same with `cut.cp_iters` large and `cut.stall = 0` (`par/rootL.par`
used 2,500 and no `wc`: at the cold root the phase runs at penalty factor 1, not `w`).

**Pitfalls.**
1. **The root gets no DNN bound under `dnn.depths = 1:48`** (the estimates' setting), hence no
   phase, no pool and no warm state: all 50 depth-1 children start cold and uncut. Use
   `0:48`. (In the tai30b solve of section 10 the root was not bounded, which is why depth 1
   took 57% of that run.)
2. **Start the phase early, and `dnn.giveup` is only a floor.** The uncut first solve is
   nearly worthless here. tai30b root, same total iterations, phase after 100 uncut
   iterations against 300 (`logs/r30_line1_I100.log`) [m]: 601.2M against 598.4M at 400,
   616.9M against 612.9M at 600, 633.7M against 633.1M at 1,600. But the cold tai50b root's
   first solve ran 470 iterations before the extrapolation rule called it hopeless; warm
   pooled children stop at 110.
3. **Phase length: long, and the stall rule is too eager for this instance.** A phase that
   fathoms is far cheaper than branching (a branch costs at least 46 x 40 iterations even
   when every child fathoms at once): `31:23` needs 589 iterations, cold depth-3 node 10
   needed 2,480 without `wc`. And a node that cannot fathom (on the path, or a same-cluster
   swap) should still run long, because its children inherit its state: below the on-path
   depth-3 node stopped by `cut.stall = 3` after 620 iterations (3.1M short) 2 of 10 children
   fathom in the first solve; below the on-path depth-4 node run 3,000 iterations (0.41M
   short) 8 of 10 do. 3,000 iterations cost less than the 6 x 46 x 400 they save. I would
   run the dives with `cut.cp_iters = 3000` and `cut.stall = 0` (or a much wider window),
   not with the 1,000 / stall 3 that most of my samples used; the price is 3,000 iterations at
   one or more nodes per level of every dive. Say which limit you used.
4. **A warm root for samples or dives**: `bnbdrv -p ROOT.par -u 458821517 -d 0 -o root.state
   qaplib/tai50b.dat` bounds the root (first solve + phase) and writes its 50 children, set
   aside by `-d 0`, with the root's shipped DNN state (`ship.bits = 64`) and pool (134 MB).
   One child: `python3 ../tho40_search/cutstate.py root.state kid.state 1 F:L`, then
   `bnbdrv -p KID.par --allow-settings-change -r kid.state -d 1 -o out.state ...` (that node
   only, its children written for the next level), or without `-d` for its subtree.
   **`-r` refuses a state written under other settings unless `--allow-settings-change`** is
   given. `size50/root.sh`, `kid.sh` (one node), `pass1.sh` + `cutset.py` (first solves of a
   random subset of a state's open nodes, `par/kid0.par` = no phase) do exactly this. The
   incumbent travels in the state file's `incumbent` line (`kid.sh` rewrites it for the
   best known + 1 runs).
5. **`bnbdrv -e -r root.state` runs but I do not trust what it did**: on tai25b
   (`size50/logs/est25/`) two dives from the saved root each reported 22-27 level-1 nodes,
   all fathomed, 22 of them "by QPB", and no survivor, although the on-path child cannot
   fathom. Check how the estimator treats a resumed stack before using it; otherwise the
   estimator recomputes the root (2,970+ iterations) and bounds 50 depth-1 children with a
   phase each (about 35,000 iterations) in every dive.
6. `QAP_NODE_TRACE`, `QAP_DNN_TRACE`, `QAP_DNN_TRACE_HIST=1` (environment, diagnostic keys):
   one line per node / per DNN node with first-solve and phase iterations and status, and the
   phase's running best every 10 iterations (`collect.py` tabulates them). The history is
   empty when `cut.stall = 0`.
7. **Where the dives will be decided**: depth 4 onward, the same-cluster children. Watch the
   number of unfathomed children per node there (I expect about 15 of 46 at depth 4 below
   the on-path node, and 0-5 below a survivor) and how many levels that lasts.
8. The search's `dnn.beta` / `dnn.cheb` are tuned for plain nug-type solves; `dnncuts`'
   defaults did at least as well on the cold root at equal iterations. One A/B on the idle
   box may pay. `pool.cp = 1` (a pooled child's first solve in CP form) was prepared
   (`par/kidQ.par`) and not run.
9. Validity guard for any run: with the incumbent at best known + 1 (`-u 458821518`) the
   path to `QAPSOLN/tai50b.sln` must never fathom (table above, and section 7).

