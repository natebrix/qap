# tai60b and tai80b with the cluster-region cuts (bound_lab agent "tai60b", Nov 16 2026)

Question: do the cluster-region cuts that transformed tai50b (`../regions/README.md`) do the same for
the two larger open instances of the family? Is either within reach of a few thousand cloud
core-hours, and does one root run give a new best lower bound?

**[m]** measured here (the run is named; `runs/TAG/`), **[e]** estimated. Every bound quoted is a
certified value of the solver (a dual certificate at each check). Iterations and nodes are the
metrics; CPU seconds are from a 4-core box that ran 3 to 4 solver processes at once.

## Follow-up (Nov 16, afternoon): Knuth dives for tai60b

Asked by the coordinating session: dives from the saved states, to replace the model's "400 to more than
10,000 core-hours" below with a measured estimate. [m]/[e] as in the rest of this file.

### Update (Nov 16, evening): more dives, the mirror's path, tai80b

**In one line: tai60b 279 core-hours (202 to 363) [e], up 8%; the mirror behaves as a copy of the path for the
five levels measured; tai80b's root levels off at 801.4M, 46% of its gap, not tai60b's 81%.**

Asked by the coordinating session after the section below was merged. **The numbers in this block supersede
the ones in "Answer" below where they differ.**

**tai60b: E2 = 29,357 nodes (21,502 to 38,340), 1.52e6 ADMM iterations (1.12e6 to 1.97e6), 1,003,935 CPU-s =
279 core-hours on this box (202 to 363)** [e]; E1 = 15,207 nodes, 147 core-hours (106 to 192). **It moved from
258 to 279 (+8%) with 2.75 times the dives at the three thin levels (20 to 55), and the interval narrowed** (it was 180 to 344).
My range stays 150 to 600 for the reason given below (paths under survivors).

| set | dives before | dives now | estimate before: CPU-s | now |
|---|---|---|---|---|
| `near3` (19 start nodes) | 4 | 11 | 142,436 ± 45,477 | 197,799 ± 66,001 |
| `near9` (7) | 8 | 22 | 126,370 ± 77,318 | 102,446 ± 52,251 |
| `near10` (8) | 8 | 22 | 54,677 ± 39,728 | 22,289 ± 14,850 |
| `nearx1` (14) | 6 | 8 | 18,232 ± 13,246 | 77,522 ± 62,476 |
| `near4` (11) | 5 | 7 | 24,391 ± 14,330 | 43,612 ± 18,114 |
| `near6` (10) | 4 | 6 | 29,487 ± 10,280 | 19,808 ± 8,930 |
| `near8` (10) | 4 | 5 | 15,803 ± 4,116 | 16,031 ± 3,197 |

- `near3` went up by 39%, as the cut-off jobs suggested it would; it is now 38% of E1's cost, depth 9 20%,
  depth 1 13%, depth 4 9%. `nearx1` quadrupled on two more dives (one of them went two levels down): depth 1 is the thinnest level left (8 dives over 14 start nodes).
- 141 dive jobs, 284 finished dives in all (`logs/dive_est.txt`). 0 pool or cold retries.
- **The survive rule, now on 286 distinct off-path dive nodes: right on 270. All 39 nodes with excess below
  s(d) branched**; the 16 errors are again all nodes that branched above it (8 of 20 between 1 and 1.6 s(d),
  6 of 39 between 1.6 and 3, 2 of 23 between 3 and 5, 0 of 165 above 5):

| excess / s(depth) | nodes | branched | fathomed by the phase | fathomed by the first solve | mean iterations of the fathomed |
|---|---|---|---|---|---|
| 0 to 0.5 | 22 | 22 | 0 | 0 | nan |
| 0.5 to 1 | 17 | 17 | 0 | 0 | nan |
| 1 to 1.6 | 20 | 8 | 9 | 3 | 482 |
| 1.6 to 3 | 39 | 6 | 14 | 19 | 256 |
| 3 to 5 | 23 | 2 | 7 | 14 | 209 |
| 5 to inf | 165 | 0 | 4 | 161 | 42 |

**The mirror's path, measured** (`mirror.txt`: the best known permutation with the facilities of locations 43
and 44 exchanged, cost 608,217,840; `chain.sh` with `ONPATH=mirror.txt`, runs `mir*`): the on-mirror node at
each depth below `44:58`, warm from its parent, incumbent at best known + 1.

| depth | on-mirror node | first solve: short | phase end: short | iterations | on-path node at that depth: phase end, iterations |
|---|---|---|---|---|---|
| 1 | `44:58` | 3,846,034 | 2,293,578 | 930 | 2,722,281, 890 |
| 2 | `43:45` | 1,966,505 | 1,554,923 | 710 | 1,576,818, 880 |
| 3 | `12:0` | 1,883,425 | 1,043,301 | 1,070 | 1,065,357, 710 |
| 4 | `1:33` | 1,041,300 | 845,817 | 710 | 733,869, 940 |
| 5 | `6:24` | 1,302,931 | 815,158 | 1,020 | 609,108, 930 |

So E2's "copy" holds for the five levels measured: the mirror's nodes cost what the path's do (4,440 against
4,350 iterations) and end as far or a little further from the cutoff (0.82M against 0.61M at depth 5; the
mirror branched on location 12 at depth 3 where the path took the column of facility 26). If anything its
subtree is a little larger than a copy. Its off-mirror children were not sampled.

**Second memory kill.** The restarted tai80b root was killed at 16:12 at phase iteration 640 (7.5 GB, beside
three dive processes one of which had grown to 2.6 GB; the container allows 14.3 GB and counts the page
cache). I saw it at 16:41. `runs/root80_killed640` has its trace, identical to the first run's to that point.

**tai80b root, run on: 801,397,334 certified at about 470 [e] + 2,320 iterations: 46% of the published gap
closed, gap 2.08%. It levels off; it does not do what tai60b's root did (81%).** [m `runs/root80`, the third
run: `root.sh root80 tai80b 818415043` with `THREADS=4`, alone on the box, 3.8 s of wall per phase iteration
against 8.1 on one thread; stopped by hand at 19:48, values from `live.txt`. `threads.patch`: three lines in
`latest/driver.C`, `QAP_BLAS_THREADS=N`, default 1 as before. The multi-threaded run's values differ from the
single-threaded one's in the seventh digit: 798,331,438 against 798,336,758 at 1,000.]

| phase iterations | certified running best | short of best known | share of the published gap closed | gain over the previous 250 |
|---|---|---|---|---|
| 250 | 787,033,428 | 31,381,615 | 0% |  |
| 500 | 793,523,956 | 24,891,087 | 21% | 6,490,528 |
| 750 | 796,394,398 | 22,020,645 | 30% | 2,870,442 |
| 1,000 | 798,331,437 | 20,083,606 | 36% | 1,937,039 |
| 1,250 | 799,481,341 | 18,933,702 | 40% | 1,149,904 |
| 1,500 | 800,214,468 | 18,200,575 | 42% | 733,127 |
| 1,750 | 800,711,844 | 17,703,199 | 44% | 497,376 |
| 2,000 | 801,055,717 | 17,359,326 | 45% | 343,873 |
| 2,250 | 801,316,552 | 17,098,491 | 46% | 260,835 |
| 2,320 | 801,397,334 | 17,017,709 | 46% |  |

- Gains per 250 iterations: 1.94M, 1.15M, 0.73M, 0.50M, 0.34M, 0.26M; per 100 over the last 800: 238k, 224k,
  162k, 123k, 95k, 104k, 123k, 115k. The geometric decay (about 0.68 per 250) stops at about iteration 1,900
  and what is left is a creep of about 110k per 100 iterations, 0.35% of the published gap. By the decay alone
  the plateau is near 801.9M (48%); with the creep, 3,000 more iterations would add about 3M (56%) if it held,
  which tai50b's and tai60b's roots suggest it would not [e].
- **So tai80b's root closes about half of its gap with this region set, against 81% for tai60b**, and its
  shortfall is 17.0M (2.08%) where tai60b's is 2.9M (0.48%). With the survive rule (a node branches when it
  holds a solution within the slack) that is a tree of another order: tai60b's on-path slack at depth 1 was
  already 0.45%. **On this evidence tai80b is not worth sizing with `line:1,top:7,linegap:4:0.03`.**
- What is not settled: whether a region set made for two dimensions closes more. tai80b's clusters are spread
  in the plane (first-axis share 0.78; tai60b 0.92), this set has the prefixes of one axis only, and
  `../tai_solves` solved tai35b and tai40b with two families each. One root with `line:2` or `line:4` is the
  next thing to try: 2.5 hours on four threads for 2,300 iterations, 7.3 GB, alone on a box of this size.
- The first restart (single-threaded, killed at 640) and the first run agree with each other to the last digit.


### Answer

**tai60b with these cuts is a run of about 250 core-hours on this box: 31,551 nodes (18,585 to 46,675),
1.5e6 ADMM iterations (1.0e6 to 2.1e6), 927,005 CPU-s = 258 core-hours (180 to 344)** [e, estimate E2].
About 350 c7a core-hours if this box's 0.70 of a c7a carries over. **My range is 150 to 600 core-hours**: the
bootstrap interval is too narrow (it does not know what four dives per level cannot see), the model below
says 463, and the three levels that carry 63% of the cost rest on 4 to 8 dives each. The first report's
central 500-2,500 was too high and its ">10,000" is not supported by anything the dives met.

| estimate | nodes | ADMM iterations | CPU-s | core-hours here | what it is |
|---|---|---|---|---|---|
| E1 | 16,399 (9,869 to 24,066) | 808,518 | 519,990 | **144** (101 to 192) | the sampling design alone: root and the 42 on-path nodes measured, the off-path children of every path level by Knuth dives |
| E2 | 31,551 (18,585 to 46,675) | 1,531,296 | 927,005 | **258** (180 to 344) | the same, with the depth-1 mirror `44:58` counted as a copy of the on-path child and its subtree |
| model (`logs/model_m42_x1.0.txt`) | 47,143 | 3.8e6 | 1.67e6 | 463 | the survive rule with the slack now measured at every depth 1-42 (no extrapolation), priced with the first report's costs per node class |
| model, slack x 0.5 / x 1.5 | 13,889 / 139,670 | 1.1e6 / 1.2e7 | | 95 / 1,764 | |

(Intervals: 5th to 95th bootstrap percentile, dives resampled within sets.)

- **Why E2 and not E1.** The mirror `44:58` (the two facilities of the 2-point cluster exchanged; its best
  completion costs 2,786 more than the best known solution) must be a second copy of the whole path: its
  phase ends 2.29M short where the on-path node's ends 2.72M short, and nothing can separate the two until
  the bottom. A dive under it does not see that: the one that drew it branched there and then drew an ordinary
  depth-2 child, which fathomed in 64 iterations (`runs/near1/s201`). The dives are unbiased in expectation, but
  a path one node wide under a survivor is met with probability about 1/50 per level, so in practice they
  miss it. E2 puts it back for the one node where it is certain. **I did not bound a node on the mirror's own
  path below depth 1**, so "copy" is an assumption resting on that one phase [e].
- **Where the cost is** (E1 by path level): the subtrees hanging off the path at depth 3 28%, depth 9 25%,
  depth 1 13%, depth 10 11%, depth 6 6%, depth 4 5%; everything from depth 12 down 6%. The first report's
  model put a second bulge at depths 17-30; **the dives find nothing there**: 48 near dives and 69 far dives at
  depths 16-43 estimate 4% of the cost between them, and the slack measured on the path keeps falling
  (61k at depth 18, 37k at 30, 4k at 42).
- **The shape.** Off the path a near node usually branches (it holds a solution within the bound's slack) and
  the child the dive then draws fathoms: 42 of the 113 finished near dives met a branching node, 5 of them
  two levels running, none three; no far dive did (114 node bounds: 111 fathom in the first solve, 3 in a
  phase). The deepest any dive went below its start is two levels. A branching dive node costs 1,023
  iterations on average (32 of the 47 ran to the 1,000-iteration phase limit).
- **The survive rule against the dives, node by node** (235 distinct off-path nodes; s(d) = the on-path
  shortfall at that depth): right on 223. **All 24 nodes with excess below s(d) branched.** The 12 errors are
  all nodes that branched with excess above s(d): 6 of 16 between 1 and 1.6 s(d), 5 of 37 between 1.6 and 3,
  1 of 18 between 3 and 5, none of 140 above 5. Eleven of the 12 (depths 4-16) are near misses that reached the
  1,000-iteration phase limit 26k to 198k short; the twelfth is at depth 40. So the rule's threshold is right for "cannot fathom" and about
  1.5 times low for "does not fathom under `cut.cp_iters = 1000`"; a longer limit may move those eleven (not tried).
- **Path check**: 42 on-path nodes bounded (depths 1-42, incumbent at best known + 1), none fathomed, the
  largest certified value 608,210,961 at depth 42, against 608,215,054. No dive bounded an on-path node (they
  start off the path), so the chain is the whole check.
- **Solver faults**: 0 pool retries and 0 cold retries in the 274 traced dive node bounds and the 24 new chain nodes
  (`pool_retry.patch` in the build; counted from the `retry:` lines of every `err.txt`).
- **tai80b: not continued.** No core was free: four cores ran the dives and the chain throughout, and under the
  container's 14.3 GB a tai80b root (7.3 GB) would have displaced two of them. Its bound stands at
  798,336,758 at about 470 + 1,000 iterations.

### What was run

- **Build**: the first report's worktree plus `../tai50b_size/est_resume.patch` and
  `../tai_solves/pool_retry.patch` (both applied cleanly).
- **The chain** (`chain.sh`, as before): the on-path node of every depth from 19 to 42 (m = 18), each warm
  from its parent with the pool, `par/kidW.par` (`cut.cp_iters = 1000`, `cut.stall = 3`). Shortfall at the end
  of the phase: 52k at depth 19, 70k, 70k, 90k, 67k, 69k, 45k, 37k, 41k, 55k, 43k, 37k (depth 30), 31k, 30k,
  27k, 18k, 20k, 12k, 22k, 11k, 18k, 18k, 11k, 4k at depth 42 (`logs/pathtable.md`). Not monotone, and not exact
  by depth 42: tai30b's path was within 1 from depth 15.
- **Strata** (`strata.py`): at every path level D = 1..43 the off-path children of the on-path node are
  split by the excess of the best completion a 2-swap descent finds: **near** = below twice the on-path
  first-solve shortfall at that depth, **far** = the rest (`runs/strata/offD.txt` lists every child with its
  excess; not committed, rebuilt by `strata.py`). 0 to 19 near nodes and 9 to 58 far nodes per level.
- **Dives** (`run_dive.sh`: `bnbdrv -e -r STATE`, start level = the stratum's depth, exponent 1.5, `est.xbar =
  1`, `est.xbar_mix = 0.3`, one seed per process; `worker.sh` + `jobs.txt`: four workers on a queue, each job a
  set, a seed and a number of dives, so that no dive is cut by a clock). 99 jobs, 242 finished dives, 274 node
  bounds in the DNN traces. Far sets: 3 dives each. Near sets: 2 to 8 dives each. 15 far dives at depths 33-42
  have an estimate and no trace line (their node fathomed before the DNN bound); in those eight processes
  `dive_est.py` pairs estimates with traces in order and may be off by one, and `far36`, `far39`, `far42`
  (no trace at all) are counted at one node per start node. Together they are under 0.1% of the cost.
- **The estimate** (`dive_est.py` -> `logs/dive_est.txt`): root (`rootM`: 1,470 iterations, 2,170 CPU-s) +
  for every level D the on-path node (measured) + the mean over the dives of the far and of the near set
  (each finished dive is an unbiased estimate of the nodes and CPU-s of all subtrees of its set). Iterations
  are the dive's CPU estimate times the iterations per CPU-second of its own nodes [e]. Below depth 43
  nothing is counted (m <= 17, milliseconds per node).

The near sets (all sets: `logs/dive_est.txt`):

| set | start nodes | dives | node bounds: first-solve fathom / phase fathom / branched | deepest dive | estimate: nodes | iterations | CPU-s |
|---|---|---|---|---|---|---|---|
| `near1` | 15 | 2 | 1 / 1 / 1 | 2 | 836 ± 822 | 47,576 | 62,533 ± 57,759 |
| `nearx1` | 14 | 6 | 3 / 3 / 1 | 2 | 468 ± 454 | 13,871 | 18,232 ± 13,246 |
| `near3` | 19 | 4 | 3 / 1 / 5 | 5 | 3,444 ± 2,035 | 146,398 | 142,436 ± 45,477 |
| `near4` | 11 | 5 | 5 / 0 / 2 | 5 | 551 ± 381 | 27,084 | 24,391 ± 14,330 |
| `near5` | 12 | 4 | 4 / 0 / 0 | 5 | 12 ± 0 | 1,403 | 978 ± 134 |
| `near6` | 10 | 4 | 2 / 2 / 3 | 7 | 945 ± 534 | 38,689 | 29,487 ± 10,280 |
| `near7` | 5 | 4 | 2 / 2 / 1 | 8 | 66 ± 62 | 4,758 | 3,349 ± 2,064 |
| `near8` | 10 | 4 | 3 / 1 / 3 | 9 | 378 ± 123 | 24,398 | 15,803 ± 4,116 |
| `near9` | 7 | 8 | 8 / 0 / 10 | 11 | 5,012 ± 3,067 | 230,644 | 126,370 ± 77,318 |
| `near10` | 8 | 8 | 7 / 1 / 5 | 12 | 2,947 ± 2,759 | 106,681 | 54,677 ± 39,728 |
| `near11` | 8 | 4 | 3 / 1 / 1 | 12 | 99 ± 92 | 6,821 | 3,105 ± 2,153 |
| `near12` | 4 | 4 | 3 / 1 / 2 | 13 | 102 ± 58 | 6,257 | 2,775 ± 1,515 |
| `near13` | 4 | 4 | 0 / 4 / 0 | 13 | 4 ± 0 | 1,676 | 689 ± 4 |
| `near14` | 6 | 4 | 3 / 1 / 1 | 15 | 68 ± 62 | 5,779 | 2,286 ± 1,519 |
| `near16` | 2 | 2 | 1 / 1 / 2 | 17 | 81 ± 0 | 19,612 | 5,916 ± 4,194 |
| `near17` | 1 | 2 | 0 / 2 / 0 | 17 | 1 ± 0 | 471 | 145 ± 2 |
| `near18` | 1 | 2 | 0 / 2 / 0 | 18 | 1 ± 0 | 561 | 156 ± 3 |
| `near20` | 6 | 3 | 1 / 2 / 0 | 20 | 6 ± 0 | 2,563 | 705 ± 310 |
| `near21` | 2 | 3 | 1 / 2 / 1 | 22 | 25 ± 23 | 1,712 | 651 ± 458 |
| `near22` | 3 | 3 | 0 / 3 / 0 | 22 | 2 ± 0 | 1,362 | 340 ± 47 |
| `near23` | 1 | 3 | 3 / 0 / 0 | 23 | 1 ± 0 | 111 | 16 ± 0 |
| `near24` | 2 | 3 | 2 / 1 / 2 | 25 | 61 ± 33 | 2,393 | 790 ± 347 |
| `near26` | 2 | 3 | 1 / 2 / 1 | 27 | 23 ± 21 | 1,848 | 424 ± 95 |
| `near27` | 1 | 3 | 0 / 3 / 0 | 27 | 1 ± 0 | 421 | 76 ± 0 |
| `near29` | 7 | 3 | 2 / 1 / 2 | 30 | 156 ± 76 | 6,341 | 981 ± 364 |
| `near30` | 2 | 3 | 3 / 0 / 0 | 30 | 2 ± 0 | 211 | 18 ± 2 |
| `near34` | 1 | 3 | 0 / 3 / 0 | 34 | 1 ± 0 | 391 | 21 ± 1 |
| `near39` | 5 | 3 | 3 / 0 / 0 | 39 | 5 ± 0 | 540 | 7 ± 2 |
| `near40` | 4 | 3 | 3 / 0 / 4 | 42 | 145 ± 98 | 25,912 | 353 ± 274 |
| `near41` | 1 | 3 | 3 / 0 / 0 | 41 | 1 ± 0 | 112 | 1 ± 0 |
| `near42` | 3 | 3 | 3 / 0 / 0 | 42 | 2 ± 0 | 258 | 2 ± 0 |

Far sets: 43 (40 with traced dives), 3 dives each; of their 114 traced node bounds 111 fathom in the first
solve and 3 in a phase, none branched; each set's estimate is about its number of nodes (10 to 58), 2 to
4,000 CPU-s per set.

By level, E1, the levels above 4,000 CPU-s (all levels: `logs/dive_est.txt`):

| D | near: start nodes, dives | nodes | iterations | CPU-s | share of E1's CPU |
|---|---|---|---|---|---|
| 1 | 15, 2 | 878 | 51,455 | 67,690 | 13.0% |
| 3 | 19, 4 | 3,483 | 148,990 | 144,869 | 27.9% |
| 4 | 11, 5 | 601 | 29,966 | 26,810 | 5.2% |
| 6 | 10, 4 | 1,000 | 40,766 | 31,056 | 6.0% |
| 7 | 5, 4 | 111 | 6,982 | 4,825 | 0.9% |
| 8 | 10, 4 | 417 | 27,103 | 17,394 | 3.3% |
| 9 | 7, 8 | 5,058 | 233,074 | 127,715 | 24.6% |
| 10 | 8, 8 | 2,986 | 110,258 | 56,473 | 10.9% |
| 11 | 8, 4 | 138 | 9,763 | 4,394 | 0.8% |
| 16 | 2, 2 | 127 | 21,567 | 6,639 | 1.3% |

The model check:

| excess / s(depth) | nodes | branched | fathomed by the phase | fathomed by the first solve | mean iterations of the fathomed |
|---|---|---|---|---|---|
| 0 to 0.5 | 14 | 14 | 0 | 0 | nan |
| 0.5 to 1 | 10 | 10 | 0 | 0 | nan |
| 1 to 1.6 | 16 | 6 | 7 | 3 | 459 |
| 1.6 to 3 | 37 | 5 | 14 | 18 | 261 |
| 3 to 5 | 18 | 1 | 7 | 10 | 235 |
| 5 to inf | 140 | 0 | 3 | 137 | 42 |

### What these numbers rest on, and what went wrong

- **Few dives where the cost is.** `near3` (19 start nodes, 28% of E1) has 4 dives, `near9` 8, `near10` 8,
  `near1`/`nearx1` 2 and 6, `near4` 5, `near6` 4. Single dives of one set range from 7 to 22,254 nodes (`near10`)
  and from 6,009 to 544,984 CPU-s (`near9`). The standard errors in the table are 30% to 95% of the means.
- **Paths under survivors are under-sampled** (the mirror's case, above). E2 corrects it at depth 1 only. The
  14 nodes below 0.5 s(d) that the dives met at depths 3-10 are the same kind in smaller: each has its own
  narrow path for the few levels until the slack falls below its excess. A dive samples that path one child
  in about 50, so the estimate has the right expectation and a long upper tail that 4 to 8 dives do not show.
  This is the main reason my range's upper end is 600 and not 344.
- **I lost 47 minutes of dives to my own mistake.** At 12:26 a cleanup command killed a scheduled-stop timer
  in a way that fired the stop: all four dive processes were terminated then. Four jobs were cut (`near3/s303`
  after its first dive, which is counted; `near16/s316`, `near17/s317`, `near18/s318` before any). I noticed at
  13:13 and ran the three expensive sets again until 13:49; `near3/s703` was cut then, in its first dive. A
  cut dive is a long dive, so those cuts bias `near3` low.
- **The dives' cutoff is 608,215,053.05, not best known + 1**: `-u 608215055` does not raise the incumbent
  saved in a state file. Off the path that changes nothing.
- **Settings as in the first report** (`par/kidW.par`), in particular the 1,000-iteration phase limit that ten
  of the twelve rule errors ran into; the tai50b dives used 2,400. Samples start from the 1,000-iteration
  root; a run would start from a longer one.
- **Node cost against the model's**: the model's 463 core-hours is 1.8 times E2 at 1.5 times the nodes; its
  prices (910 iterations for a branching node, 510 for a phase fathom, 45 for a first-solve fathom) are about
  right for the first two (branching dive nodes: 1,023; phase fathoms: 235-459) and its count of branching
  nodes is higher than what the dives saw below depth 12.

### Not done

More dives at depths 1, 3, 9 and 10; any node on the mirror's path below depth 1; a longer phase limit for
near misses; the tai80b root; a run of any subtree end to end as a check on the dives.

## The first report (morning): roots, tree samples, model

(Its estimate, item 3 below, is superseded by the dives above.)

### Answer (first report)


1. **Both roots give new best lower bounds** (against Kim, Kojima, Toh's plain DNN bounds, the
   published ones given in the brief):

   | instance | best known (checked from `QAPSOLN/*.sln`) | published bound | gap | this root, certified | gap | share of the published gap closed |
   |---|---|---|---|---|---|---|
   | tai60b | 608,215,054 | 592,829,061 | 2.53% | **605,315,696** at 470 + 3,000 iterations [m `rootL`] | **0.48%** | **81%** |
   | tai80b | 818,415,043 | 786,900,529 | 3.85% | **798,336,758** at about 470 [e] + 1,000 iterations, still rising [m `root80`] | **2.45%** | **36%** |

   tai60b's root passes the published bound 250 iterations into the phase (592,871,409 at
   470 + 250 [m `ab_c6`]); tai80b's at phase iteration 250 (787,039,986 [m `root80`]). Neither has
   plateaued: tai60b gains 124k over its last 250 iterations (207k, 167k, 143k, 124k per 250 from
   2,000 on; by that trend about 0.8M more, a plateau near 606.1M [e]); tai80b was stopped by the
   clock while gaining 486k per 100 (802k, 668k, 486k over its last three hundreds).

2. **tai60b's tree is closer to tai30b's than to tai50b's at the top, and then it is a path with a
   fringe**, like tai50b's:
   - **Depths 1 and 2: every off-path child is fathomable at once, except the mirror.** The root
     branches on location 44, then 43: the isolated 2-point cluster. Of 12 random depth-1 children
     10 fathom in the first (pooled, warm) solve in 10-91 iterations, one fathoms after 482 phase
     iterations, and one does not fathom: `44:58`, the child that swaps the two facilities of that
     cluster (its best completion costs 2,786 more than the best known solution; no bound can
     separate that before the bottom of the tree). Under the on-path node 12 of 12 random depth-2
     children fathom in the first solve (10-60 iterations) [m `d1p1`, `d1_44_39`, `d1_44_58`,
     `d2p1_under_path`].
   - **The fringe appears with the first location in a large cluster** (depth 4, location 1 of the
     14-point cluster; depth 3 is a column branch on a facility). There 5 of 5 random children
     fathom in the first solve, and of 9 children whose facility belongs to the same cluster in the
     best known solution 5 fathom in the first solve and 4 do not; 3 of those 4 were given a phase
     and all three end it unfathomed (0.23M, 0.25M and 0.40M short at the 1,000-iteration limit)
     [m `d4p1_under_path`, `d4_1_39`, `d4_1_59`, `d4_1_29`].
   - **A survivor's children mostly die at once, but not all, and the ones that do not are again
     same-cluster ones**: under the depth-4 survivor `1:39`, 8 of 10 random children fathom in the
     first solve (27-90 iterations); the other two (0.24M and 0.62M short) carry a facility of their
     location's cluster; their phases end 0.17M and 0.02M short at the 1,000-iteration limit
     [m `d5p1_under_1_39`, `d5_1_39_22_48`, `d5_1_39_22_9`]. tai50b's survivor had 1 of 10. So the fringe does
     compound for at least two levels.
   - **Which children survive is predictable from the near-optimal solutions**: a child survives
     when the best completion a 2-swap descent finds for it costs less than best known + the
     slack of the bound at that depth (the on-path node's shortfall). That rule gets 26 of the 27
     depth-1 and depth-4 first solves above right (`nearopt.py`, section 4). The best known
     solution has 8 pair swaps within 10,000 of it, 40 within 100,000 and 213 within 1M, all but
     two inside one cluster [m `nearopt.py swaps`]: that is the fringe, and it compounds only as
     far as the bound's slack exceeds those excesses.
   - **The on-path nodes close faster than tai50b's**: 2.72M short at depth 1, 1.58M, 1.07M, 0.73M,
     0.61M, 0.46M, 0.41M, 0.37M, 0.34M, 0.22M, 0.14M, 0.13M at depth 12, then 0.10M, 0.10M, 0.09M,
     0.08M, 0.08M, 0.06M at depths 13-18 (0.45% of the best known value at depth 1, 0.02% at depth
     12, 0.01% at depth 18; tai50b: 1.05%, 0.85%, 0.67% at depths 1-3) [m chain `d*_path`]. That is
     far closer than tai50b and far less close than tai30b (0.01% at depth 3, exact from depth 15).

3. **Estimate for tai60b: 500-2,500 core-hours on this box (5e6-4e7 ADMM iterations, 6e4-4e5
   nodes; about 700-3,500 c7a core-hours), with a range of 400 to more than 10,000** [e]. It is a
   model, not a dive: the tree of all prefixes whose best completion lies within the measured
   slack of the bound, costed with the measured iterations per node class (section 5, with what it
   rests on). About 440 of those core-hours are depths 0-16, where the slack is measured; the rest
   is the extrapolation. So: **tai60b looks within reach of a few thousand core-hours, but not
   safely: the number to get next is the bound's slack at depths 17-30 and off the path**, which
   decides between 400 and 10,000. Knuth dives from the saved root would give it. tai80b's tree was
   not sampled (its root alone costs 8 s per iteration here).

## 1. Build, regions, cost per iteration

**Build.** A clean worktree of this branch's base (26894b2) with `../regions/regions.patch`
(`git apply`, no conflict, no fix-up to the patch) and `make -C latest -j4 && make -C latest dnntest
dnncuts`. Two things about the container, not the patch: without the `scipy-openblas32` wheel the
makefile links the reference BLAS (7x slower, RUNNING.md), so `pip install scipy-openblas32` and a
rebuild came first; `gen_regions.py -stats` needs `scipy`. One diagnostic addition of mine,
`live.patch` (10 lines in `latest/dnncut.C`): `QAP_CP_LIVE=FILE` appends the phase's certified
running best at every check, so that a long root leaves a record if it is killed (`runs/*/live.txt`).
It does not touch the numbers.

**Best known values** (`structure.py`): `sum A[i,k] B[p(i),p(k)]`, A the first matrix, gives
608,215,054 for `tai60b.sln` and 818,415,043 for `tai80b.sln`, the listed values [m].

**Structure** (`gen_regions.py -stats`, `structure.py`) [m]:

| instance | 2-D MDS stress | share of the first axis | natural clusters (sizes) | merge height jump | mean distance within / between | share of sum(D) between |
|---|---|---|---|---|---|---|
| tai50b (for comparison) | 0.0003 | 0.971 | 5 (3, 4, 16, 15, 12) | 100 -> 457 | 71 / 1070 | 0.979 |
| tai60b | 0.0003 | 0.918 | 4 (14, 2, 8, 36); at the next levels the 36 split into 20, 12, 4 | 341 -> 791 | 195 / 1188 | 0.892 |
| tai80b | 0.0005 | 0.781 | 7 (10, 10, 21, 10, 6, 20, 3) | 137 -> 248 | 51 / 640 | 0.983 |

tai60b is the collinear kind, less cleanly than tai50b: one 36-point blob (diameter 546 against 1,188
between clusters) holds 11% of the distance sum inside the "natural" clusters, and the 2-point cluster
{43, 44} lies off the axis. tai80b is two-dimensional (as tai35b and tai40b, where the same recipe
closed 66-71%), with 7 tight clusters carrying 98% of the distance.

**Regions chosen** (location side, `cut.side = 0`; files in `regions/`), on the tai50b pattern
`line:1,top:K,linegap:4:0.03`:
- tai60b: **K = 6** (`tai60b.comb6.txt`, 86 regions: 59 prefixes along the principal axis, the 31
  unions of the six clusters 14, 2, 8, 20, 12, 4, and the between-cluster half-plane cuts in four
  directions). K = 6 and not the "natural" 4 because the 36-point blob is not tight: its three parts
  are as far apart as tai50b's clusters are wide. Checked against K = 4 (`tai60b.comb4.txt`, 64
  regions) at the root [m `rootM`, `c4M`, same settings]: equal to 250 phase iterations (592.87M
  against 592.94M), then K = 6 ahead: 601.51M against 601.37M at 700, **602,989,748 against
  602,786,441 at 1,000**, at 0.4% more CPU per iteration.
- tai80b: **K = 7**, the natural clusters (`tai80b.comb7.txt`, 138 regions). Not compared with
  anything: no time at 8 s per iteration. With clusters spread in two dimensions `line:2` or `line:4`
  may do better (not tried here or on tai35b/40b).

**Phase penalty**: `dnn.beta_sched = w=0.36,wc=2.88` from the start. At the tai60b root [m `ab_c6`,
`ab_c6_nowc`, 470 + 250 iterations]: 592,871,409 with `wc`, 587,781,768 without; the run without
reaches at 250 what the run with reaches at about 115.

**Cost** [m]:

| | m = 50 (tai50b, `../regions/size50`) | m = 60 | m = 80 |
|---|---|---|---|
| s per iteration, first solve (cold root) | 0.53 | 1.36 (`rootL`: 470 iterations, 647 CPU-s) | about 6.4 (`root80`: about 3,010 CPU-s for about 470 [e]) |
| s per iteration, phase | 0.72 loaded, 0.4-0.5 idle | 1.55 (`rootL`: 3,000 in 4,649 CPU-s), 1.1-1.3 at depth 1 | 7.9-8.3 (`root80`) |
| peak memory | | 2.5 GB | 7.3 GB in the phase |
| root state file (`ship.bits = 64`) | 134 MB | 277 MB | not written (run stopped) |

So 3x and 16x tai50b, as expected from (m-1)^6. Cores: every solver process had one core
(`OPENBLAS_NUM_THREADS=1`; a tai80b run started with 2 BLAS threads still used one). The container
has 4 cores and a 14.3 GB memory limit: **the first tai80b root was killed for memory** after 51
minutes when three tai60b processes ran beside it (`runs/root80_oomkilled`); after that at most two
tai60b processes ran beside tai80b, and one core idled. `dnn.amx` was not tried (`par/ab_c6_amx.par`
prepared; the CPU has AMX-BF16).

## 2. tai60b root [m]

`root.sh NAME`: `bnbdrv -p par/NAME.par -u 608215054 -d 0 -o runs/NAME/root.state` (the settings of
`../regions/size50/par/rootL.par` with `wc=2.88`, `cut.stall = 0`, `cut.cp_iters` as below). The
first solve stops as "hopeless" at 470 iterations at **581,315,106** (26.9M short): that is the plain
bound at the give-up, not converged; Kim, Kojima and Toh's 592,829,061 is the converged plain DNN
bound, which I did not rerun.

| phase iterations | certified (`rootL`, comb6, wc=2.88) | short of best known | share of the published gap closed |
|---|---|---|---|
| 0 (470 plain) | 581,315,106 | 26,899,948 | (below the published bound) |
| 250 | 592,871,409 | 15,343,645 | 0.3% |
| 500 | 599,190,046 | 9,025,008 | 41% |
| 1,000 | 602,989,748 (= `rootM`) | 5,225,306 | 66% |
| 1,500 | 604,106,475 | 4,108,579 | 73% |
| 2,000 | 604,674,348 | 3,540,706 | 77% |
| 2,500 | 605,048,662 | 3,166,392 | 79% |
| 3,000 | **605,315,696** | **2,899,358** | **81%** |

5,296 CPU-s in all. Not a plateau (gains per 250: 207k, 167k, 143k, 124k); the on-path depth-1 child
reaches 605.49M and the depth-2 one 606.64M, so a longer root has at most about 1M to give.
`rootM` (1,000 phase iterations, 5.23M short) is the root every tree sample below starts from: it was
ready 50 minutes earlier. A dive runner should start from the longer one.

## 3. tai60b tree samples [m]

Protocol of `../regions/size50`: `kid.sh` (one node warm from its parent's state with the pool, first
solve then phase, `par/kidW.par`: `cut.cp_iters = 1000`, `cut.stall = 3`), `pass1.sh` + `cutset.py`
(first solves only of a random subset, `par/kid0.par`), `chain.sh` (the on-path node of each level,
incumbent at best known + 1, cutoff 608,215,054.05). `F:L` as in the node trace: entry F (a location,
the first matrix) of the permutation holds facility L; the on-path child of a branch on row r is
`r:p[r]` with p from `tai60b.sln` (`onpath.txt`). `collect.py` prints the table of all runs.

**The path** (`runs/d*path`; first solves are 110 iterations each):

| depth | child (m) | first solve: certified, short | phase end: certified | short | % of best known | phase iterations (end) | CPU-s | s/it |
|---|---|---|---|---|---|---|---|---|
| 1 | `44:45` (59) | 603,840,108, 4,374,946 | 605,492,773 | 2,722,281 | 0.448% | 780 (stalled) | 1,185 | 1.33 |
| 2 | `43:58` (58) | 606,006,543, 2,208,511 | 606,638,236 | 1,576,818 | 0.259% | 770 (stalled) | 968 | 1.10 |
| 3 | `45:26` (57) | 606,838,380, 1,376,674 | 607,149,697 | 1,065,357 | 0.175% | 600 (stalled) | 678 | 0.95 |
| 4 | `1:33` (56) | 607,144,783, 1,070,271 | 607,481,185 | 733,869 | 0.121% | 830 (stalled) | 803 | 0.85 |
| 5 | `6:24` (55) | 607,330,940, 884,114 | 607,605,946 | 609,108 | 0.100% | 820 (stalled) | 707 | 0.76 |
| 6 | `7:19` (54) | 607,184,172, 1,030,882 | 607,752,706 | 462,348 | 0.076% | 1000 (limit) | 805 | 0.72 |
| 7 | `0:32` (53) | 607,757,734, 457,320 | 607,805,120 | 409,934 | 0.067% | 600 (stalled) | 503 | 0.71 |
| 8 | `14:54` (52) | 607,767,805, 447,249 | 607,841,134 | 373,920 | 0.061% | 600 (stalled) | 463 | 0.65 |
| 9 | `24:2` (51) | 607,818,783, 396,271 | 607,878,851 | 336,203 | 0.055% | 600 (stalled) | 401 | 0.56 |
| 10 | `28:50` (50) | 607,774,646, 440,408 | 607,995,616 | 219,438 | 0.036% | 1000 (limit) | 603 | 0.54 |
| 11 | `18:43` (49) | 607,914,004, 301,050 | 608,071,159 | 143,895 | 0.024% | 1000 (limit) | 543 | 0.49 |
| 12 | `23:18` (48) | 607,998,639, 216,415 | 608,087,137 | 127,917 | 0.021% | 1000 (limit) | 515 | 0.46 |
| 13 | `38:41` (47) | 608,033,437, 181,617 | 608,111,140 | 103,914 | 0.017% | 1000 (limit) | 468 | 0.42 |
| 14 | `15:12` (46) | 607,879,688, 335,366 | 608,113,739 | 101,315 | 0.017% | 1000 (limit) | 478 | 0.42 |
| 15 | `31:10` (45) | 608,097,002, 118,052 | 608,121,522 | 93,532 | 0.015% | 600 (stalled) | 314 | 0.44 |
| 16 | `20:21` (44) | 608,084,974, 130,080 | 608,130,965 | 84,089 | 0.014% | 1000 (limit) | 426 | 0.38 |
| 17 | `16:17` (43) | 608,086,740, 128,314 | 608,136,872 | 78,182 | 0.013% | 900 (stalled) | 357 | 0.35 |
| 18 | `19:13` (42) | 608,113,127, 101,927 | 608,154,280 | 60,774 | 0.010% | 920 (stalled) | 383 | 0.37 |

Shortfalls are against the best known value 608,215,054 (the cutoff is that + 0.05). "limit": the phase
reached `cut.cp_iters = 1000` still gaining; "stalled": stopped by `cut.stall = 3`.

The branching order is the search's own (rule 9): the 2-point cluster's locations 44 and 43, a column
branch on facility 26 (its children are the 58 free locations; on-path 45), then locations of the
14-point cluster (1, 6, 7, 0), then others.

**Depth 1** (`d1p1`: 12 random children + the on-path one, first solves; then phases):

| child | first solve | phase | outcome |
|---|---|---|---|
| `44:2`, `44:5`, `44:9`, `44:22`, `44:23`, `44:25`, `44:40`, `44:42`, `44:54`, `44:59` | fathomed in 10-91 iterations (13-100 CPU-s) | | 10 of 12 fathom at once |
| `44:39` | 1,958,421 short at 110 | fathomed after 482 iterations | fathomed (592 iterations, 741 CPU-s) |
| `44:58` (the swap of the cluster's two facilities; excess 2,786) | 3,846,034 short | stalled 2,293,578 short after 820 | **branches**: the mirror of the path |
| `44:45` (on-path) | 4,374,946 short | stalled 2,722,281 short after 780 | branches |

**Depth 2** (`d2p1_under_path`: 12 random children `43:x` of `44:45`): **12 of 12 fathom in the first
solve**, 10-60 iterations, 13-60 CPU-s. On-path `43:58`: 2.21M short, phase stalls 1.58M short.

**Depth 3**: a column branch (facility 26 over 58 locations); only the on-path child was run.

**Depth 4** (`d4p1_under_path`: children `1:x` of `44:45,43:58,45:26`; location 1 is in the 14-point
cluster, whose facilities in the best known solution are 0, 19, 23, 24, 29, 32, 33, 35, 36, 39, 49,
55, 57, 59):

| child | kind | excess of its best completion (`nearopt.py`) | first solve | phase |
|---|---|---|---|---|
| `1:8`, `1:13`, `1:17`, `1:50`, `1:53` | random (5 of 5), other clusters | 13.6M-42.9M (`1:53`, `1:8`, `1:13`, `1:17`, `1:50`: 13.6M, 15.9M, 30.3M, 33.0M, 42.9M) | fathomed, 23-50 iterations | |
| `1:23`, `1:0` | same cluster | 2.86M, 2.47M | fathomed, 56 and 70 iterations | |
| `1:32`, `1:35`, `1:36` | same cluster | 1.57M, 1.33M, 1.23M | fathomed, 87-117 iterations (by 2k-14k) | |
| `1:29` | same cluster | 0.96M | 705,728 short | 247,676 short at the 1,000 limit, gaining 25k per 100: **branches** |
| `1:24` | same cluster | 0.51M | 1,124,476 short | not run |
| `1:59` | same cluster | 0.51M | 798,843 short | 231,441 short at the 1,000 limit: **branches** |
| `1:39` | same cluster | 0.25M | 922,686 short | 396,970 short at the 1,000 limit: **branches** |
| `1:33` (on-path) | | 0 | 1,070,271 short | stalled 733,869 short after 830 |

So at depth 4 the 43 other-cluster children fathom at once [e from 5 of 5], and of the 13 same-cluster
ones 8 fathom in the first solve and **about 5 survive besides the on-path one** [e: 4 of the 9 run
miss the first solve, `1:19` (excess 0.96M) is predicted to, and 3 of 3 phases end unfathomed at the limit]. tai50b had about 15 of 46 surviving at its
depth 4. The phases of these survivors end 0.74M to 1.2M below the upper bounds on their subtrees'
optima (shortfall + excess), against 0.73M for the on-path node: off the path the slack is up to 1.6
times larger.

**A survivor's children** (`d5p1_under_1_39`: first solves of 10 random children of `1:39`, which
branches on location 22, in the 12-point cluster {14..25}, whose facilities in the best known solution
are 2, 9, 12, 13, 17, 18, 21, 25, 31, 43, 48, 54): 8 fathom in the first solve (27-90 iterations,
20-62 CPU-s), among them one same-cluster facility (`22:17`); `22:9` ends 237,006 short and `22:48`
624,649 short, both same-cluster. Phases: `22:48` 174,189 short at the 1,000 limit: branches
(`d5_1_39_22_48`); `22:9` 22,694 short at the limit, gaining 7k per 100: counted as branching here; about 300 more
iterations would fathom it [e] (`d5_1_39_22_9`). Of the 4 survivors' phases at depths 4-5, all were still
gaining at 1,000 iterations: a longer limit trades iterations at the node for a smaller subtree.

## 4. Why the fringe is what it is: near-optimal solutions (`nearopt.py`) [m for the counts]

- Pair swaps of the best known permutation: 3 cost less than 1,000 more, 8 less than 10,000, 18 less
  than 50,000, 40 less than 100,000, 77 less than 250,000, 123 less than 500,000, 213 less than 1M;
  211 of those 213 exchange two locations of the same cluster. The 8-point cluster {52..59} (diameter
  14) alone has 22 of the 40 cheapest: almost any arrangement of its eight facilities is within 0.1M.
- A node whose subtree holds a solution with excess e cannot fathom before the bound's slack at its
  depth falls below e. Measured at depth 1 and 4: children with excess up to 0.96M did not fathom in
  the first solve (first-solve slack of the on-path node there: 4.37M and 1.07M), children from 1.23M
  up did at depth 4; at depth 1 the first solve fathomed children with excess 6.28M and more (the one
  at 6.28M, `44:40`, by 22k after 91 iterations) and missed one at 6.75M (`44:39`, 1.96M short, then
  fathomed by its phase). Of the 27 sampled first solves at depths 1 and 4, the rule "fathoms in the
  first solve iff excess > the on-path node's first-solve shortfall" is wrong on that one only: the
  slack off the path is somewhat larger than on it. The same shows in the phases: the three depth-4
  survivors end 0.74M-1.2M below their subtrees' upper bounds, the on-path node 0.73M below its optimum.
- This is the tai50b agent's "same-cluster children survive" seen from the solutions' side, and it
  says when the fringe stops: when the slack is below the swap excesses, i.e. gradually from about
  0.5M (123 swaps) down to 10k (8 swaps).

## 5. Tree-size estimate for tai60b (`model.py`) [e]

**Model.** Fix a branching order of the rows. A depth-d prefix *branches* when the best completion
found for it (2-swap best-improvement descent from its parent's completion) has excess below s(d),
the phase's slack at depth d; it *fathoms in a phase* when its excess is between s(d) and 1.6 s(d)
(the first-solve slack: on-path first-solve shortfalls are 1.3-2.2 times the phase-end ones); otherwise
it *fathoms in the first solve*. s(d) = the on-path node's shortfall at the end of its phase, measured
for d = 1..16 (or 18), continued below by a ratio r per level. Costs from section 3: a branching node
110 + 800 iterations, a phase fathom 110 + 400, a first-solve fathom 45; 1.3 s per iteration at
m = 59 scaled by (m/59)^6 (measured: 0.93-0.95 s at m = 56-57, 0.72 at m = 54, 0.46 at m = 48).
Order: the search's measured rows down the path, then the remaining locations cluster by cluster
(`path`) or at random (`rand`).

s(d) measured to depth 16 (2.72M ... 0.084M, the table in section 3); `logs/model16_*.txt` have the
per-depth profiles:

| order | s(d) for d = 1..16 | ratio r below depth 16 | branching nodes | nodes bounded | ADMM iterations | core-hours on this box |
|---|---|---|---|---|---|---|
| path | measured | 0.8 | 987 | 43,057 | 3.4e6 | 475 |
| path | measured | 0.9 | 4,090 | 137,417 | 1.1e7 | 661 |
| random (seed 1) | measured | 0.9 | 1,805 | 60,347 | 4.9e6 | 492 |
| path | measured | 0.95 | 22,888 | 635,972 | 5.3e7 | 1,378 |
| path | measured x 1.5 | 0.9 | 11,574 | 407,171 | 3.4e7 | 2,362 |
| path | measured x 1.5 | 0.95 | 76,303 | 2.2e6 | 1.8e8 | 4,925 |
| path | measured | 1.0 (the slack never falls below 0.084M) | 1.7e6 | 1.8e7 | 2.3e9 | 7,918 |
| path, with depths 17-18 measured too (0.078M, 0.061M; rows 16, 19) | measured | 0.9 below depth 18 | 2,482 | 83,063 | 6.7e6 | 524 |
| same | measured | 0.95 | 10,458 | 271,358 | 2.3e7 | 723 |

**Reading.**
- **Depths 0-16, where s(d) is measured, cost about 440 core-hours in the model** (560 branching
  nodes and 2.2M iterations; most of it at depths 4-10, where a node costs 0.5-1 s per iteration and
  the same-cluster children of the large clusters survive). That part does not depend on r, only on
  the survive rule; with the slack x 1.5 it is about 1,700 (2,060 branching nodes, 9.0M iterations).
- **Below depth 16 everything depends on how fast the slack keeps falling.** Measured ratios per
  level: 0.58-0.83 at depths 2-6, 0.89-0.91 at depths 7-9, 0.65 at depths 10-11, then 0.89, 0.81,
  0.97, 0.92, 0.90 at depths 12-16: about 0.9 lately. At r = 0.9 the lower tree adds 220 core-hours;
  at 0.95, 940; if the slack stopped at 0.084M (40 pair swaps are cheaper than that) the fringe would
  run to the bottom: 7,500, and with the larger slack as well, several times that.
- The last two rows came in at the end (depths 17 and 18: ratios 0.93 and 0.78, in line with 0.9) and
  are lower than the rows above them at the same r, because two more rows of the order are the
  search's own: the model moves by a factor of 2 with the order alone.
- **Central reading: 500-2,500 core-hours on this box (5e6-4e7 iterations, 6e4-4e5 nodes), about
  700-3,500 c7a core-hours; range 400 to more than 10,000** [e]. The x 1.5 rows are not a
  pessimistic decoration: the three depth-4 survivors I ran sit 1.0-1.6 times the on-path slack below
  their subtrees' upper bounds.
- For comparison, the tai50b agent's cross-check for tai50b was 50-200 core-hours if its fringe dies
  in one or two levels and 1e4-1e5 if it compounds; here the fringe is seen compounding for two
  levels (depth 4 and 5) and the model says it goes on at a slowly falling rate.

**What it rests on.** One root; the on-path chain to depth 18 (one node per level); 13 + 12 + 14 + 10
first solves at depths 1, 2, 4, 5 and 7 off-path phases; the survive rule checked on those nodes
only at depths 1 and 4; a descent that finds upper bounds on subtree optima (so it can only
under-count survivors where it is poor); one fixed branching order where the search picks its own in
every subtree (the survivor `1:39` branched on location 22, the path on location 6); phases stopped by
`cut.stall = 3` or the 1,000 limit, so s(d) is what these settings leave, not the relaxation's limit
(the root: 5.2M short at 1,000 iterations, 2.9M at 3,000). The cost side assumes off-path survivors
cost what on-path nodes do.

**What would move it.** (a) The slack below depth 18 and off the path: r = 0.9 against 1.0 is a factor
of 12, the slack x 1.5 a factor of 3.5. On tai30b the on-path bound was within 0.01% from depth 3 and
within 1 from depth 15; here it is at 0.01% at depth 18. (b) Longer phases at survivors and a longer root lower
s(d) at a cost per node the model does not trade off. (c) Symmetric-looking arrangements of the
8-point cluster: if the bound cannot tell them apart until they are fixed, that cluster alone is up
to 8! = 40,320 leaves; the model has this through the swap excesses, but only as far as the descent
finds them. (d) c7a against this box: `../tho40_search` measured this box at 0.70 of c7a per matched
node, so core-hours here are about 1.4 c7a core-hours when the box is loaded [e].

## 6. Validity [m]

With the incumbent at best known + 1 (cutoff 608,215,054.05) no node on the path to
`QAPSOLN/tai60b.sln` fathoms and no certified value there exceeds 608,215,054: 18 nodes, depths 1-18,
the largest 608,154,280 at depth 18 (table in section 3). Off the path, the mirror `44:58` (subtree optimum at most 608,217,840)
ends at 605,921,475 and the depth-4 same-cluster survivors at 607,818,084 and 607,983,613, below
their subtrees' upper bounds 608,463k and 608,727k. The root values are below the best known value.
No brute-force check was rerun here; `../regions/README.md` section 3 has those for the patch.

## 7. tai80b root [m `root80`]

(Morning run, now `runs/root80_first1000`; the evening run to 2,320 iterations is in the update block at the top.)
`root.sh root80 tai80b 818415043` (`par/root80.par`: the same settings, `tai80b.comb7.txt`, `dnn.depths
= 0:78`, `cut.cp_iters = 5000`). One core; two tai60b processes beside it. The run was stopped by
hand at phase iteration 1,000, at the end of the session, so its values are the certified running best at each check
(`runs/root80/live.txt`), and the first solve's iteration count was never printed (about 470 [e]
from about 3,010 CPU-s; the cold tai50b and tai60b roots both stopped at 470).

| phase iterations | certified running best | short of best known | share of the published gap closed | CPU-s |
|---|---|---|---|---|
| 0 (first solve, about 470 [e]) | 774,516,079 | 43,898,964 | (below the published bound) | about 3,010 |
| 100 | 781,173,418 | 37,241,625 | (below) | 3,803 |
| 200 | 785,588,035 | 32,827,008 | (below) | 4,584 |
| 250 | 787,039,986 | 31,375,057 | 0% | 4,996 |
| 300 | 789,369,685 | 29,045,358 | 8% | 5,386 |
| 400 | 791,804,155 | 26,610,888 | 16% | 6,195 |
| 500 | 793,499,663 | 24,915,380 | 21% | 7,023 |
| 600 | 795,282,477 | 23,132,566 | 27% | 7,856 |
| 700 | 796,379,920 | 22,035,123 | 30% | 8,664 |
| 800 | 797,182,337 | 21,232,706 | 33% | 9,500 |
| 900 | 797,850,656 | 20,564,387 | 35% | 10,317 |
| 1000 | 798,336,758 | 20,078,285 | 36% | 11,125 |

The first solve gave up at 774,516,079 (43.9M short; the published converged plain bound is
786,900,529). The phase passes the published bound at iteration 250 and
ends **11.4M above it: 798,336,758, a gap of 2.45% against the published 3.85%**. Gains per 100
iterations: 2.4M-4.4M up to 400, then 1.70M, 1.78M, 1.10M, 0.80M, 0.67M, 0.49M: slowing by about 0.8 per
hundred, so perhaps 2M more by that trend, an end near 800M [e]; tai50b and tai60b both kept creeping
past such a reading. The share closed (36% at 1,000, against 66% for tai60b at 1,000) is what the
two-dimensional siblings showed in `../regions` (tai35b 71%, tai40b 66% at 2,300 with one axis of
prefixes): more half-plane directions are the thing to try. 11,125 CPU-s (3.1 core-hours) in all, 8.1 s
per phase iteration. **A root run to 3,000 phase iterations is about 7.5 core-hours on this box and
needs 7.3 GB** [e from m]; a depth-1 sample would cost about 15 minutes per first solve (110 iterations)
and 2 hours per phase, which is why the tree was not touched.

## 8. Not measured

- tai60b: a search or a Knuth dive; anything off the path below depth 5; a phase for `1:24`; the survivors' phases beyond 1,000 iterations (all four were still gaining); any
  depth-3 (column-branch) sibling; the root beyond 3,000 phase iterations; the converged plain bound;
  samples warm from the 3,000-iteration root (all samples start from the 1,000-iteration one);
  `dnn.amx`; three Dykstra sweeps / extrapolation (`../cpfast`).
- tai80b: anything but the root; other region sets; the root beyond 1,000 phase iterations; the first
  solve's iteration count (the run was stopped before it printed its trace line).
- Whether 608,215,054 and 818,415,043 are optimal.

## Files

- `root.sh`, `kid.sh`, `pass1.sh`, `cutset.py`, `chain.sh`, `collect.py`, `mkpar.sh`, `rss.sh`, the
  queues `q*.sh`; `structure.py`, `nearopt.py`, `model.py`, `pathrows.txt`, `onpath.txt`; `live.patch`.
- `par/`: `rootL` (3,000), `rootM` (1,000), `c4M` (comb4), `ab_*` (250-iteration A/B), `kidW`, `kid0`,
  `root80`. They carry absolute paths to `regions/`.
- `regions/`: `tai60b.comb6.txt` (used), `tai60b.comb4.txt`, `tai80b.comb7.txt`.
- `runs/TAG/`: `cmd.txt`, `log.txt`, `nodes.txt`, `dnn.txt` (with the phase history), `live.txt`.
  State files are not committed (277 MB each).
- `logs/`: timing runs (`time60.log`, `time80.log`), queue outputs, `model*.txt`.
