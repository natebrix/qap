# tho40: the whole tree under the cut phase, by Knuth dives (bound_lab agent `tho40_dives`, Nov 16 2026)

Question: what is tho40's whole tree under the cut-phase settings, in nodes and CPU, by Knuth-style
dives, with an interval, and how does it compare with the same estimator under production settings?
**[m]** measured here (file named), **[e]** estimated. "c7a CPU-days" = CPU s here x 1.29 / 86,400
(calibration below).

## Answer

1. **Zone 0.06 (`tho40.par` + the recommended block): about 226 c7a CPU-days, standard error 39
   (bootstrap 80% band 179-273), and 1.10 million nodes (standard error 0.27 million)** [m, 51 dives
   over all 12 depth-1 strata, `results/final.txt`]. In CPU here: 15.2 million s (se 2.6 million).
2. **Against base: 5.8x less CPU and 31x fewer nodes** than the laptop study's 1,308 c7a CPU-days and
   34.4 million nodes (ratio 0.173, se 0.042; nodes 0.032, se 0.011) [m]. Against base dives run here
   with the same estimator on the same strata (48 dives: 1,036 c7a CPU-days, se 463): ratio 0.22,
   se 0.105 [m]. The two base figures agree within their errors.
3. **The two hand-over arms take about half of zone 0.06's CPU on the same dive paths, but their own
   whole-tree totals are too noisy to show it** [m, `results/paired_*.txt`, `results/final.txt`]:

   | arm | dives | whole tree, c7a CPU-days (se; bootstrap 80%) | nodes (se) | CPU against zone 0.06, same-path pairs (80% band) | against zone 0.06, all pairs |
   |---|---|---|---|---|---|
   | zone 0.06 | 51 | 226 (39; 179-273) | 1.10M (0.27M) | 1 | 1 |
   | hand-over `20:2:3` | 31 | 160 (55; 106-214) | 1.34M (0.67M) | 0.55 (0.49-0.62), 28 pairs | 0.67 (0.48-0.89), 31 pairs |
   | combined policy | 39 | 215 (82; 132-303) | 1.87M (0.74M) | 0.49 (0.43-0.55), 34 pairs | 0.95 (0.65-1.27), 39 pairs |

   Combined against hand-over, same-path pairs: 0.86 (0.76-0.93), 25 pairs; all pairs 1.33 (0.93-1.78).
4. **What the dives say about "about 165 c7a CPU-days (140-210)" for the combined policy.** They do
   not contradict it and they do not confirm it. The combined arm's own total is 215 with a standard
   error of 82, so anything from about 110 to 320 is within 1.3 standard errors. The better-determined
   route is zone 0.06's total times the same-path ratio: 226 x 0.49 = **about 110 c7a CPU-days [e]**
   (roughly 80-150 from the two errors combined), which is below the study's range. The gap between
   110 and 215 is three dives: with the same seeds the hand-over arms drew a different depth-3 or
   depth-4 node than zone 0.06 in 5 of 39 pairs, and three of those landed on depth-4 nodes whose
   phase stalls 1,000-1,600 short and whose 36 children all need work (single dives worth 5.7, 6.0 and
   13.3 million CPU-s here; `results/paired_comb_cut06.txt`). Such nodes are in every arm's tree; zone
   0.06's 51 dives met fewer of them. So: the per-node saving (about 2x) carries over the whole tree
   on 34 paired dives; the absolute figure is uncertain by a factor of about two either way, and the
   study's interval (140-210) is narrower than these dives can support.
5. **Where the cost is under the cut arm (zone 0.06)** [m, `results/analyze.txt`]: depth 4 carries
   52% of the CPU, depth 5 30%, depth 3 9%, depth 6 8%; depths 1-2 under 1%. Base here: depth 6 43%,
   depth 5 25%, deeper the rest. Phase iterations are 57% of all DNN iterations (hand-over 50%,
   combined 43%).
6. **Phases at their 3,400 limit**: 3 of 150 phases seen under zone 0.06 (one at depth 3, two at depth
   5), which the estimator weights to about 2,000 of 89,000 phases in the tree (se 1,300; 52% of it
   from one dive). The hand-over arms hit the limit only at depths 1-2 (all 5 depth-1 phases, 11 of 39
   depth-2 phases seen under the combined policy), under 1% of the CPU. No phase ended by give-up.
7. **Sanity: PASS** [m, `results/sanity.txt`]. No dive run (57 runs) found any solution under the
   incumbent 240,516. All 50 nodes on the path to the best known solution met in any run have bounds
   at most 240,516. A forced run down that path (zone 0.06, incumbent 240,518, `runs/path_s2`)
   branched at every one of its 36 nodes (depths 1-36), no bound above 240,516.0000, and reached 240,516.

### What it rests on, and what could make it wrong

- **Few dives, heavy tail.** 2-6 dives per stratum and arm. The largest single dive is 12% of zone
   0.06's CPU total, 26% of hand-over's, 31% of combined's. Standard errors from 2-6 values are
   themselves poorly known, and the bootstrap bands are too narrow for the same reason: read each
   interval as a lower bound on the uncertainty. A stratum whose few dives all missed its heavy
   depth-4 nodes is estimated low; the laptop study called its own base figure "likely low" for that
   reason, and the same applies here to every arm.
- **Same-path ratios condition on the path**; they measure the saving at the nodes both arms visit,
   not differences in which children survive. The all-pairs ratios include those and are noisy.
- **CPU calibration** (1.29) comes from one 40-node comparison under 4-process load; per-iteration
   times under load varied by a few percent between arms (`results/analyze.txt`, last table per arm).
- **Incumbent**: 240,516 here; the laptop study used 240,518. The base dives here use 240,516.
- **Cold start**: each stratum starts with `-F 26:c`, so the depth-1 node is bounded without the
   root's warm start; depths 1-2 are under 1% of the cost.
- **Not run**: the `cut.zone = 0.03` arm (time went to the two arms the coordinator added).
- The hand-over and combined arms ran 2 dives in the six light strata, 3-6 in the heavy ones.

## Method

- **Estimator.** The solver's own (`bnbdrv -e`, `latest/estimate.C`): a dive draws one child per
   level with probability proportional to (child's QPB bound score)^10 mixed with 10% uniform
   (`QAP_EST_CHILD_QPB=1`), and weights what it sees by 1/probability. All probabilities are positive,
   so the estimate of any additive quantity (nodes, CPU, iterations, phases) is unbiased (Knuth 1975).
- **Exact frontier instead of recomputing the top.** The root (1 node, 0.47 CPU-s, QPB only, branches
   on facility 26 into 12 children [m, `logs/root/`]) is not rerun. Each process takes one depth-1 node
   (`-F 26:c`, start level 2), bounds it exactly once, and dives from its depth-2 children. Tree =
   root + sum over the 12 strata of (depth-1 node, exact + mean of that stratum's dives). The strata
   partition the tree and each mean is unbiased for its subtree, so the sum is unbiased; the variance
   is the sum of the strata's s^2/n.
- **Full enumeration below depth 4** (`DIVE_FULL_DEPTH=4`, cut arms only): draws are made at depths
   2, 3 and 4; under the sampled depth-4 node every child is taken, so that subtree is solved exactly,
   pools and all. Inclusion probabilities stay positive and known, so the estimate stays unbiased, and
   the pooled children run as they would in the search. Base dives sample at every level, as the
   laptop study did.
- **Fixed batches.** Every process ran a fixed number of dives (3 or 2; base 4); only finished batches
   are counted, and the order of the queue never depended on results, so there is no bias from
   stopping. No batch was cut off.
- **Pairing.** Hand-over batch `a`/`b` uses the seed of zone 0.06's batch `a`/`b` in the same stratum,
   and the combined arm the same again. Draws depend on the seed and the QPB scores, so the arms
   follow the same path unless their bounds leave different children.
- **Patch** `dive_trace.patch` (against `latest/`): per-node trace lines, `DIVE_FULL_DEPTH`,
   `DIVE_PATH` (follow a permutation), and resetting the importance weight per node. The hand-over arm
   adds `handover.patch`, the combined arm `combined.patch` (the coordinator's, from
   `tho40_search/stop/`; both applied and built cleanly on top of the dive patch).
- **Arms** (`par/`): `base` = `experiments/estimates/tho40/tho40.par`; `cut06` = base + the block
   with `cut.zone = 0.06`; `hand` = cut06 with `dnn.gu_zone = 0`, `dnn.giveup_warm = 100`,
   `cut.handover = 20:2:3`; `comb` = hand with `dnn.beta_sched = w=0.36,wc=2.88`,
   `cut.handover = 20:0:0.22:1.5:100:3:1.25`, `cut.cp_wc_gap = 300`; `repro` = the laptop study's
   `infra/phase1/recipes/tho40.par`. All `-u 240516` except `repro` and the path run (240,518).

## Step 1: the existing estimate's method reproduced [m, `results/repro.txt`]

`runs/repro_s0`: stratum `26:0` exactly as the laptop study ran it (start level 3, 20 dives).

| the top (40 nodes, deterministic) | laptop | here |
|---|---|---|
| depth-3 start nodes | 1,482 | 1,482 |
| sum of their scores | 4.01095e-07 | 4.01095e-07 |
| FW iterations | 150,150 | 150,150 |
| DNN iterations, depth 1 / depth 2 | 250 / 8,980 | 250 / 8,970 |
| CPU s | 856.9 | 1,485.8 |

Dives (369 on the laptop, 20 here): median node estimate 21,755 / 19,256; mean 204,668 / 408,521;
last depth 3: 27% / 30%, 4: 53% / 50%, 5: 18% / 20%. Same estimator, same top, same distribution
within what 20 dives can show.

## CPU calibration

- Machine [m, `logs/calib/machine.txt`]: `nproc` 4; Intel Xeon 2.10 GHz (family 6 model 207), 4
   cores, 1 thread per core, KVM, AVX-512 and AMX.
- Idle core, DNN at m = 36: 300 iterations in 31.87 s = **0.106 s per iteration** [m,
   `logs/calib/idle_m36.txt`]. Under the 4-process load of the runs: 0.119 at m = 36, 0.139 at m = 37,
   0.166 at m = 38, 0.103 at m = 35 [m, `results/analyze.txt`].
- Here against the laptop: the reproduction's top, identical work, 1,485.8 / 856.9 = 1.73x the
   laptop's CPU [m]. The laptop study's factor is c7a = 2.23 x laptop (its README), so
   **c7a s = 2.23 / 1.73 = 1.29 x this box's loaded s** [e]. Base's reference: 1,300 c7a CPU-days,
   about 3.9 c7a CPU-s per node; zone 0.06 here: about 17.8 c7a CPU-s per node [e].

## Tables

### Per stratum, zone 0.06 against base [m, `results/final.txt`]

| stratum | cut dives | cut nodes (se) | cut c7a CPU-days (se) | base here, c7a CPU-days (se), 4 dives | laptop study, c7a CPU-days (se) | CPU cut / base here | CPU cut / laptop |
|---|---|---|---|---|---|---|---|
| 26:0 | 3 | 22,819 (11,323) | 1.9 (0.4) | 0.2 (0.0) | 8.9 (2.9) | 9.61 | 0.22 |
| 26:1 | 3 | 14,218 (6,635) | 8.3 (3.5) | 15.6 (6.8) | 55.4 (16.8) | 0.53 | 0.15 |
| 26:2 | 6 | 73,700 (54,984) | 17.7 (8.9) | 443.7 (428.6) | 53.0 (11.9) | 0.04 | 0.33 |
| 26:3 | 3 | 16,715 (8,647) | 4.4 (1.3) | 89.3 (88.6) | 65.5 (10.7) | 0.05 | 0.07 |
| 26:8 | 3 | 7,636 (7,075) | 1.3 (0.7) | 11.7 (5.2) | 32.8 (8.4) | 0.11 | 0.04 |
| 26:9 | 6 | 207,850 (101,230) | 25.8 (13.6) | 11.2 (7.2) | 92.4 (29.1) | 2.30 | 0.28 |
| 26:10 | 6 | 257,184 (145,951) | 61.1 (18.9) | 71.2 (58.8) | 179.8 (44.7) | 0.86 | 0.34 |
| 26:11 | 3 | 33,759 (17,769) | 5.3 (0.8) | 32.0 (13.3) | 277.8 (192.7) | 0.17 | 0.02 |
| 26:16 | 3 | 12,969 (7,766) | 3.7 (2.2) | 13.7 (5.0) | 44.5 (10.8) | 0.27 | 0.08 |
| 26:17 | 3 | 33,172 (19,114) | 2.8 (0.4) | 123.3 (62.7) | 252.9 (81.8) | 0.02 | 0.01 |
| 26:18 | 6 | 266,717 (170,126) | 56.0 (27.0) | 57.9 (37.7) | 110.1 (24.0) | 0.97 | 0.51 |
| 26:19 | 6 | 148,298 (82,608) | 37.6 (13.3) | 165.9 (116.1) | 134.5 (39.9) | 0.23 | 0.28 |

The per-stratum ratios against 4 base dives are noise-dominated; the laptop column (241-369 dives per
stratum) is the steadier denominator. Strata 26:10, 26:18 and 26:19 hold two thirds of the cut arm's
cost. Per stratum for the hand-over and combined arms: `results/final.txt`.

### Per depth, zone 0.06 (tree estimates) [m, `results/analyze.txt`]

| depth | nodes | CPU s here | share | first-solve it | phase it | phases | fathomed | stalled | at the limit | phases seen: fathomed / stalled / limit |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 12 | 1,249 | 0.0% | 2,990 | 3,340 | 5 | 0 | 5 | 0 | 0 / 5 / 0 |
| 2 | 401 | 54,851 | 0.4% | 87,234 | 230,905 | 401 | 50 | 351 | 0 | 1 / 50 / 0 |
| 3 | 10,802 | 1,363,636 | 9.0% | 2,206,351 | 7,393,571 | 9,050 | 2,833 | 5,915 | 302 | 6 / 39 / 1 |
| 4 | 240,799 | 7,920,274 | 52.2% | 17,521,697 | 47,661,781 | 43,428 | 26,826 | 16,602 | 0 | 8 / 9 / 0 |
| 5 | 597,672 | 4,590,519 | 30.2% | 27,144,214 | 15,929,518 | 29,622 | 22,612 | 5,271 | 1,739 | 15 / 6 / 2 |
| 6 | 245,350 | 1,249,649 | 8.2% | 9,323,534 | 4,238,014 | 6,628 | 6,628 | 0 | 0 | 8 / 0 / 0 |

Depth 1 is exact (5 of the 12 depth-1 nodes are within the zone and ran a phase; none fathomed).
Nothing was seen below depth 6 in any cut arm.

### Whole-tree iterations [m, `results/analyze.txt`]

| arm | first-solve it (se) | phase it (se) | phase share | phases (se) |
|---|---|---|---|---|
| zone 0.06 | 56.3M (12.8M) | 75.5M (16.1M) | 57% | 89,135 (18,980) |
| hand-over | 47.3M (22.3M) | 47.6M (14.9M) | 50% | 102,700 (39,300) |
| combined | 76.1M (30.4M) | 56.8M (23.6M) | 43% | 151,644 (64,715) |
| base here | 728M (359M) | 0 | 0 | 0 |

### Per dive

`results/dives.csv`: every dive of every arm (169): last depth, nodes visited, process CPU, and its
estimates of nodes, first-solve iterations, phase iterations and CPU. Mean process CPU per dive [m]:
zone 0.06 about 420 s, hand-over 620 s, combined 650 s, base 80 s. `results/paired_*.txt` list each
pair.

## Files

- `README.md`; `results/` (the scripts' output: `final.txt`, `analyze.txt`, `paired_*.txt`,
   `sanity.txt`, `repro.txt`, `dives.csv`).
- `scripts/`: `run.sh` (one batch), `worker.sh` (queue), `analyze.py`, `final.py`, `paired.py`,
   `repro.py`, `sanity.py`, `dives_table.py`. `queue.txt`, `logs/queue.log`: what ran, when. The queue
   line `path_s2` taken at 01:52 was a first path run that left the path at depth 2 (the first draw
   was not forced); it was deleted and rerun by hand at 02:06 with the fix now in the patch.
- `runs/<arm>_s<c>_<batch>/`: `cmd.txt` (command, exit code), `estimate.par`, `est.out`, `err.txt`,
   `dives.txt` (the solver's history line per dive), `trace.txt` (per-node and DNN trace).
- `par/`, `dive_trace.patch`, `handover.patch`, `combined.patch`, `logs/calib/`, `logs/root/`.
