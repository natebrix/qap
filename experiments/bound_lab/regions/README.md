# Regional triangle cuts on non-grid regions: tai*b (bound_lab agent "regions", Nov 15 2026)

Question: tai50b has no grid, so it runs with no cut phase. Do the regional triangle
inequalities of `latest/dnncut.h`, which are valid for any region, bite when the regions
are built from the structure of Taillard's type-b distances (clusters, half-plane prefixes)?

**Answer: yes, very strongly, with regions on the location side that separate the
clusters.** Numbers below are certified values; [m] = measured (log named), [e] = estimated.

## Follow-up (Nov 16): sizing tai50b's tree with `cut.regions` -- IN PROGRESS

Files in `size50/` (`par/`, `runs/`, scripts). Stages: root with a saved state, settings,
depth-1/2/3 samples warm from the root with its pool, validity. The Knuth dives and the
root's plateau are run elsewhere (coordinator's scope change); this section gives the
settings and the stage samples, and a rough cross-check estimate at the end.

### For the dive runner

(first version, written before the tai50b samples; revised below as they come in)

**Build.** Apply `regions.patch` to `cognac` HEAD 6e993cc7 (`git apply`), then
`make -C latest -j4 && make -C latest dnntest dnncuts`. It adds `cut.regions` to bnbdrv and
`-regions` to dnncuts; without it `cut.regions` is an unknown key.

**Regions file.** `regions/tai50b.comb.txt` (61 regions, 0-based location indices, the first
matrix = `cut.side = 0`). Rebuild with
`python3 -I gen_regions.py qaplib/tai50b.dat line:1,top:5,linegap:4:0.03 -o FILE`.
`cut.regions` takes a path; my par files carry an absolute one
(`/home/user/qap/experiments/bound_lab/regions/regions/tai50b.comb.txt`): change it.

**Settings block** (`size50/par/kidA.par` = `experiments/estimates/tai50b/tai50b.par` with
these changes; the root's own par differs only in `cut.cp_iters` and `cut.stall = 0`):

```
dnn.depths = 0:48        # was 1:48. Pitfall 1 below
dnn.giveup = 100         # was 600: hand over to the phase early
dnn.giveup_warm = 100
dnn.gu_zone = 0          # was 0.003: no 100-iteration wait before a near-miss give-up
pool.on = 1
cut.regions = <path>/tai50b.comb.txt
cut.side = 0
cut.zone = 1             # a phase at every DNN miss (the plain bound is 6% off at the root)
cut.warmy = 1
cut.rounds = 12
cut.iters = 200          # with cut.stall = 3: the stall window is 600 iterations
cut.freeze = 0
cut.cp = 100:1000
cut.cp_iters = 3000      # node 10 (depth 3, cold) needed 2,480 phase iterations to fathom
cut.stall = 3
```

**Pitfalls.**
1. **The root gets no DNN bound under `dnn.depths = 1:48`** (the estimates' setting), hence no
   phase, no pool and no warm state: all 50 depth-1 children start cold and uncut. Use
   `0:48`. In my tai30b solve (section 10) the root was not bounded, which is why depth 1
   took 57% of that run.
2. **Start the phase early.** The first (uncut) solve is nearly worthless here: the phase
   passes the plain bound's plateau within about 100 of its own iterations. tai30b root,
   same total iterations, phase started after 100 uncut iterations against 300
   (`logs/r30_line1_I100.log`, `logs/r30_line1.log`) [m]: 601.2M against 598.4M at 400 total,
   616.9M against 612.9M at 600, 627.8M against 625.7M at 1,000, 633.7M against 633.1M at
   1,600. With `dnn.giveup = 100` the first solve stops as "hopeless" at 110 iterations.
3. **Phases run long and are worth it**, and the default limit (rounds x iters + freeze) cuts
   them: give `cut.cp_iters` explicitly. The stall rule (after 600 iterations: stop when the
   gap is more than 3 x the last 600 iterations' gain) does not stop a node on the optimal
   path, whose gap and gains shrink together: on tai25b the on-path depth-1 child ran to the
   3,000 limit (8,256 short, `size50/logs/t25kid.*`). Such nodes are one per level, so this
   costs 48 x 3,000 iterations in a whole run, but in a *dive* every level's on-path-like
   child pays it: consider `cut.giveup = T:G` or a smaller `cut.cp_iters` for the dives and
   state which you used.
4. **A warm root for samples or dives**: `bnbdrv -p ROOT.par -u 458821517 -d 0 -o root.state
   qaplib/tai50b.dat` bounds the root (first solve + phase) and writes its 50 children,
   set aside by `-d 0`, with the root's shipped DNN state (`ship.bits = 64`) and pool.
   One child: `python3 ../tho40_search/cutstate.py root.state kid.state 1 F:L` (keeps the
   depth-1 node line whose permutation has entry F equal to L; the node trace prints the same
   `F:L`), then `bnbdrv -p KID.par --allow-settings-change -r kid.state -d 1 -o out.state ...`
   (bound that node only and write its children) or without `-d` for its whole subtree.
   **`-r` refuses a state written under other settings unless `--allow-settings-change`** is
   given (the root's par differs from the children's in `cut.cp_iters`/`cut.stall`).
   `size50/root.sh` and `size50/kid.sh` do exactly this.
5. `QAP_NODE_TRACE`, `QAP_DNN_TRACE`, `QAP_DNN_TRACE_HIST=1` (environment, diagnostic keys) give
   one line per node / per DNN node with first-solve and phase iterations, status and the
   phase's running best every 10 iterations.
6. A dive (`bnbdrv -e`) bounds every child of each node on its way down. At depth 1 that is
   50 warm pooled children; expect most of a dive's cost there and at depth 2 (numbers in
   the tables below when measured).
8. **Measured on tai50b so far (Nov 16, 02:00; `size50/runs/rootL`, `size50/runs/d1p1`)** [m]:
   - Root in the search (`par/rootL.par`: the block above with `cut.cp_iters = 2500`,
     `cut.stall = 0`): the first solve did **not** stop at 110 iterations but at 470
     ("hopeless", 417.85M): `dnn.giveup = 100` is only the floor of the extrapolation rule.
     Then 2,500 phase iterations: **450,196,995** (8.62M short), 2,970 iterations, 2,042 CPU-s
     (0.69 s per iteration at m = 50 on the loaded box). That is slower per iteration count
     than `dnncuts`' defaults on the same root (449.47M at 300 + 1,600; section 5): the
     search's `dnn.beta = 0.21` / `dnn.beta_sched = w=0.36` / `dnn.cheb = 3` are tuned for
     plain nug-type solves, not for this phase. Worth one A/B on the idle box before the dives
     (dnncuts' defaults: `-g 1 -G 1 -d 4`, no beta schedule).
   - It branches on row 31 (children `31:x`). The root state file is 134 MB (64-bit).
   - Depth-1 children, warm with the root's pool, first solve only: the first 13 of a
     random 20 all give up at 110 iterations (55-59 CPU-s each), **none fathoms**, 1.4M to
     8.4M short (median 6.0M; the on-path child `31:7` is the furthest, 8.36M). So at depth 1
     every child needs a phase, and a dive pays 50 x (110 + a phase) there unless the root is
     run much longer first (a converged root is cheap relative to that: run it to its
     plateau and start the dives from its state).
9. **`-e` and a saved root.** I have not checked whether `bnbdrv -e` accepts `-r`; if it does
   not, `estimate.par`'s start level (`start_dive_level`) makes the estimator bound the top
   levels itself, which at 50 x phases is hours per dive on tai50b. Check this first.
7. Validity guard for any run: with the incumbent at best known + 1 (`-u 458821518`) the
   path to `QAPSOLN/tai50b.sln` must never fathom; my on-path checks are in section 7 and
   below.

## Headline numbers

- **tai30b root: 99.7% of the DNN gap closed.** Plain DNN 598.3M at 3,300 iterations
  (optimum 637,117,113, gap 38.9M); with the 29 prefixes of the locations along their
  principal axis, 635.9M at 2,300 iterations and **636,993,730 at 6,300** (0.019% under the
  optimum) [m: `logs/r30_line1_long.log`]. Random regions and facility-side regions: 0%.
- **Roots, DNN gap closed at equal iterations (2,300): tai25b 96%, tai30b 97%, tai35b 71%,
  tai40b 66%**; tai20b's plain bound is already exact [m: section 6].
- **tai50b root: 449,469,933 certified at 1,900 iterations and rising**, against the plain
  430,147,161 at 1,500 and QAPLIB's 431,090,700: the gap to the best known value falls from
  6.25% to at most 2.04% [m: `logs/t50_root_comb.log`].
- **tai50b depth 3: 3 of 3 sampled misses now fathom**: the two nearest (2.7M and 3.8M short
  at 1,500 plain iterations) at 628 and 680 total iterations, and the furthest of the
  12-node sample (14.4M short) at 2,780 [m: `logs/t50_d3_node*_comb.log`].
- **Validity**: 0 violations in 240 brute-force node checks and 32 small-root checks;
  nodes on the optimal path of tai30b (depths 3, 6, 10) and tai50b (depth 3) stay below the
  optimum; the search of tai30b from an incumbent one above the optimum finds it, with every
  on-path bound at or below it [m: sections 3, 7, 10].
- **In the search** (`cut.regions`, section 10): **tai30b solved to proven optimality in 363
  nodes and 1,520 CPU-s** from an incumbent one above the optimum. 29 of the 30 depth-1
  children fathom (8 by the plain bound, 21 by the cut phase); the one that does not is the
  optimal child, and below it the tree is a single path [m].

## 0. What was changed (`regions.patch`, worktree `/home/user/wt-regions`)

- `DNNCutParams::regions` (`latest/dnncut.h`): arbitrary regions as lists of original
  indices of `side`. `cut_regions()` (`latest/dnncut.C`) maps them to the node's free
  elements exactly as `grid_regions()` does (drops empty, full, duplicate and complementary
  sets); with a grid given too they are added to its prefixes. Everything downstream
  (separation, `assemble`, the CP projector's region-pair classes, pools) already worked on
  general `Region`s, so both forms (`dnn_cut_rounds`, `dnn_cut_cp`) and `dnn_cp_bound` use
  them with no other change.
- `dnncuts -regions FILE [-side 0|1]`: one region per line, original 0-based indices, `#`
  comments.
- `cut.regions = FILE` for the search (`latest/bnb.C`, `settings.C`, `RUNNING.md` row;
  `tests/settings_doc.py` passes). Used for the tai30b solve of section 10.
- In `tai*b.dat` the first matrix is the symmetric distance matrix (side 0 = locations) and
  the second the asymmetric flow matrix (side 1 = facilities); checked for n = 12..50: the
  `.sln` permutations evaluate to the listed optima with `sum A[i,k] B[p(i),p(k)]`, A the
  first matrix. Optima verified from `QAPSOLN/*.sln`: tai12b 39,464,925, tai15b 51,765,268,
  tai20b 122,455,319, tai25b 344,355,646, tai30b 637,117,113, tai35b 283,315,445, tai40b
  637,250,948; tai50b best known 458,821,517.

## 1. Structure of the tai*b distance matrices (`gen_regions.py -stats`, `logs/structure.txt`) [m]

| instance | 2-D MDS stress | share of the first axis | natural clusters (sizes) | merge height jump | mean distance within / between | share of sum(D) between clusters |
|---|---|---|---|---|---|---|
| tai20b | 0.0025 | 0.994 | 2 (16, 4) | 72 -> 771 | 36 / 771 | 0.916 |
| tai25b | 0.0009 | 0.971 | 3 (15, 4, 6) | 120 -> 369 | 55 / 739 | 0.949 |
| tai30b | 0.0005 | 0.970 | 2 (16, 14); 3 at the next level (16, 4, 10) | 313 -> 1106 | 106 / 1106 | 0.917 |
| tai35b | 0.0007 | 0.836 | 3 (15, 2, 18) | 168 -> 618 | 72 / 660 | 0.922 |
| tai40b | 0.0004 | 0.636 | 3 (12, 9, 19); 4 at the next level (12, 9, 12, 7) | 257 -> 897 | 115 / 1020 | 0.943 |
| tai50b | 0.0003 | 0.971 | **5 (3, 4, 16, 15, 12)** | 100 -> 457 | 71 / 1070 | **0.979** |

- The distances are Euclidean distances of points in the plane, rounded (stress 0.03-0.25%).
- The points sit in a few tight clusters: the dendrogram's merge heights jump 3-10x at the
  cluster level, and 92-98% of the total distance is between clusters.
- On tai20b/25b/30b/50b the clusters are nearly collinear (97-99% of the MDS variance on one
  axis): the distance matrix is close to a *line metric*, |x_a - x_b|, which decomposes
  exactly into the cut semimetrics of the prefixes of the sorted order, just as a Manhattan
  grid distance decomposes into row and column prefixes. tai35b and tai40b are genuinely 2-D.
- Flows: 43-57% nonzero, heavy-tailed (the top 10% of nonzero pairs carry 53-66%).

So the objective is close to "sum over thresholds t of (gap_t) x (flow cut by the set of
facilities placed left of t)": a stack of weighted graph-partitioning problems with
prescribed part sizes, on a heavy-tailed flow graph. The membership vector of a prefix or
cluster region S (z_i = facility i sits in S) is exactly the partition indicator, and the
regional triangle inequalities are the triangle inequalities on it.

## 2. Region families (`gen_regions.py`)

Location side (side 0): `hier` every cluster of the average-linkage dendrogram; `top:K` the
K top clusters and all their unions; `line:K` half-plane prefixes, K directions over the
2-D MDS coordinates, direction 0 = first principal axis (n-1 prefixes per direction);
`linegap:K:G` only the prefixes whose gap to the next point is at least G x diameter (the
cuts that fall between clusters); `ball:K` balls around every location; `random:R`
(control). Facility side (side 1): `fhier` (dendrogram of the flow graph), `fhub:K` (top
hubs by total flow), `fspec:K` (prefixes along Laplacian eigenvectors). `comb` =
`line:1,top:K,linegap:4:0.03` (K the natural number of clusters) is the set used on tai50b.
Region files are in `regions/`.

## 3. Validity [m]

Every number is a certified value of `dnncuts` (worktree build). None exceeded a known optimum.

| check | runs | result | log |
|---|---|---|---|
| tai12b, tai15b roots, 12 random regions x 4 seeds x both sides x both forms (CP form, Polyak rounds with `-delta` = 1% of the optimum), cutoff 2 x optimum | 32 | all <= optimum; largest 39,464,924.9992 (tai12b, optimum 39,464,925) and 51,570,574 (tai15b, optimum 51,765,268) | `logs/val_*.log` |
| brute force (`validity.py`): tai12b, tai15b and three synthetic clustered instances (n = 8-10), 3 random nodes each with 8 free elements, node optimum by enumeration; families random (both sides), line:2, hier, ball:3, fhier, fspec:1, a mixed set; both forms; cutoff 2 x node optimum | 240 | **0 violations**; 87 runs end within 0.05 of the node optimum, largest (certified - optimum)/optimum = -1.8e-11 | `logs/validity_seed1.txt` |
| tai30b root, `line:1`, 6,300 iterations, cutoff 2 x optimum | 1 | 636,993,730 <= 637,117,113: the bound climbs to 0.019% under the optimum and flattens there (636,931,896 at CP iteration 4,000, 636,991,811 at 6,000) | `logs/r30_line1_long.log` |
| nodes on the path to the optimum (tai30b depths 3, 6, 10; tai50b depth 3), cutoff = optimum + 1000 | see section 7 | | `logs/onpath_*.log` |

A warning that applies to everything below (the one `dnn_cuts/README.md` gives): on tai*b the
plain bound converges slowly, and on tai12b it is exact, so "gain over the 300-iteration
uncut bound" means nothing by itself. tai12b: every run "gains" 29,231, and so does the plain
solve (`logs/plain_tai12b_cpplain.log`: 39,464,924.9993 at 2,000 more iterations). tai15b: the
plain bound plateaus at 51,584,127 (0.35% under the optimum) and random regions do not beat
it. All comparisons below are against the plain solve continued from the same 300-iteration
point for the same number of iterations (`dnncuts -cpplain`), and against its plateau.

## 4. tai30b root: which regions bite (`logs/r30_*.log`, `logs/plain_tai30b.log`) [m]

300 uncut iterations (583,806,468), then N more: plain (`-cpplain`) or the CP form
(`-I 300 -k 1000 -cpsep 100:1000`). Optimum 637,117,113. "CP it" counts iterations after
the first 300. Last column: (cut bound - plain bound) / (optimum - plain bound), both at the
run's last iteration.

| regions file | side | regions | violated at round 0 (max) | CP it 100 | CP it 300 | CP it 500 | CP it 1000 | CP it 2000 | final | gap closed vs plain at equal iterations |
|---|---|---|---|---|---|---|---|---|---|---|
| plain (no cuts) | | | | 588,971,576 | 592,608,458 | 594,182,914 | 596,084,588 | 597,633,906 | 598,265,483 (it 3000) | |
| line1 | 0 | 29 | 9016 (0.2185) | 598,393,201 | 612,856,374 | 620,678,742 | 630,337,061 | 635,876,664 | 635,876,664 (it 2000) | 96.9% |
| line2 | 0 | 57 | 10777 (0.2185) | 598,740,205 | 613,942,046 | 622,487,312 | 632,184,740 |  | 632,184,740 (it 1000) | 88.0% |
| line4 | 0 | 101 | 22088 (0.2185) | 598,897,011 | 614,272,302 | 623,762,138 | 633,256,214 |  | 633,256,214 (it 1000) | 90.6% |
| linegap40.03 | 0 | 5 | 4138 (0.2185) | 599,873,018 | 614,913,020 | 622,879,321 | 632,254,269 |  | 632,254,269 (it 1000) | 88.1% |
| hier | 0 | 27 | 5486 (0.2106) | 598,573,709 | 614,097,218 | 622,039,390 | 631,397,350 |  | 631,397,350 (it 1000) | 86.1% |
| top | 0 | 3 | 2790 (0.2106) | 596,872,234 | 605,953,755 | 614,279,630 | 626,628,320 |  | 626,628,320 (it 1000) | 74.4% |
| ball4 | 0 | 39 | 5758 (0.1360) | 593,960,175 | 602,405,927 | 605,752,758 | 609,026,885 |  | 609,026,885 (it 1000) | 31.5% |
| random30 | 0 | 30 | 705 (0.0642) | 589,023,942 | 592,680,062 | 594,225,168 | 596,119,457 |  | 596,119,457 (it 1000) | 0.1% |
| fhier | 1 | 27 | 9 (0.0093) | 588,968,360 | 592,607,191 | 594,182,416 | 596,086,876 |  | 596,086,876 (it 1000) | 0.0% |
| fhub10 | 1 | 10 | 0 (0.0000) | 588,971,576 | 592,608,458 | 594,182,914 | 596,084,588 |  | 596,084,588 (it 1000) | -0.0% |
| fspec2 | 1 | 56 | 0 (0.0000) | 588,971,576 | 592,610,188 | 594,184,366 | 596,085,764 |  | 596,085,764 (it 1000) | 0.0% |

- **Location-side regions that separate the clusters close 74-97% of the DNN gap.** The
  three top clusters alone (3 regions) give 74% at 1,000 iterations; five between-cluster
  half-plane cuts (`linegap`) 88%; the 29 prefixes along the principal axis 88% at 1,000
  iterations, 97% at 2,000 and 99.7% at 6,000 (636,993,730; plain plateau about 598.3M to
  598.8M, so the gap goes from about 38.5M to 0.12M).
- **Controls give nothing**: 30 random location sets 0.1%; facility-side regions (flow
  clusters, hubs, spectral prefixes) 0.0%, with 0 to 9 violated triangles at the plain
  iterate. The effect is the geometry of the distance matrix, not "more iterations" and not
  "any regions".
- Violations at the plain 300-iteration iterate: up to 0.22 (2,800-22,000 violated
  triangles) for the structured location regions, against 0.004-0.02 for atomic triangles on
  nug30 (`dnn_cuts/README.md`) and 0.06 for random sets here.
- Cost: the CP iterations cost the same as plain ones with up to 60 regions (49 s per 1,000
  at m = 30 either way); 101 regions (`line:4`) cost 45% more per iteration.
- The cut bound passes the plain bound's 3,300-iteration value after 100 CP iterations.

## 5. tai50b [m]

Regions: `regions/tai50b.comb.txt` = `line:1,top:5,linegap:4:0.03`, 61 regions on the
locations (the 49 prefixes along the principal axis, the unions of the 5 clusters, the
between-cluster half-plane cuts in 4 directions). True cutoff 458,821,516.05. Machine shared
(5-6 busy processes on 4 cores): 0.5-0.8 s per iteration at m = 47-50.

### Root (`logs/t50_root_comb.log`; plain: `../../tai50b_feasibility/root/root.log`)

`dnncuts -c 458821516.05 -regions regions/tai50b.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 3000 -T 1`,
stopped by a 26-minute cap at 1,900 iterations. At the plain 300-iteration iterate 89,148
triangles are violated, the largest by 0.347.

| total iterations | plain DNN | short | with the regions | short | share of the plain shortfall closed |
|---|---|---|---|---|---|
| 400 | 421,577,777 | 37,243,739 | 419,179,071 | 39,642,445 | -6% |
| 600 | 426,923,126 | 31,898,391 | 430,765,110 | 28,056,406 | 12% |
| 800 | 428,620,522 | 30,200,994 | 438,094,599 | 20,726,917 | 31% |
| 1,000 | 429,428,065 | 29,393,451 | 442,541,276 | 16,280,240 | 45% |
| 1,200 | 429,637,492 | 29,184,024 | 445,188,225 | 13,633,291 | 53% |
| 1,500 | 430,147,161 | 28,674,355 | 447,585,096 | 11,236,420 | **61%** |
| 1,900 | (431,090,700 is QAPLIB's best bound, BBCPOP) | (27,730,816) | **449,469,933** | **9,351,583** | 66% of BBCPOP's gap |

- The root bound is 18.4M above QAPLIB's best published bound for tai50b and still rising
  by 380k per 100 iterations when the cap hit (increments 496k, 434k, 382k per 100; their
  ratio is drifting up from 0.80 to 0.88, so no plateau can be read off: at least 2-3M more
  [e]). The root gap is at most 2.0% of the best known value, against 6.0-6.25%.
- In branching terms: the plain shortfall falls 1.8x per level (NOTES.md), so the root with
  cuts at 1,900 iterations sits where plain depth-2 nodes do (28.7M / 1.8^2 = 8.9M).

### Depth-3 nodes of `tai50b_feasibility/NOTES.md` section 2 (`logs/t50_d3_node*_comb.log`)

Same command with `-x SPEC`, 25-minute cap. Plain columns from `depth3/nodeN.log`.

| node | spec | plain: short at 400 / 1,000 / 1,500 iterations | plain status | with the regions | total iterations (300 uncut + CP) | CPU-s | violated at the plain iterate (max) |
|---|---|---|---|---|---|---|---|
| 7 | `0:44,49:49,25:31` | 7,642,007 / 3,272,870 / 2,679,320 | miss | **fathomed** (458,840,171) | **628** (300 + 328) | 301 | 52,684 (0.212) |
| 1 | `20:19,28:24,14:35` | 9,281,971 / 4,506,776 / 3,847,935 | miss | **fathomed** (458,876,867) | **680** (300 + 380) | 332 | 60,094 (0.254) |
| 10 | `5:40,40:22,44:3` (the worst of the 12) | 18,492,881 / 15,127,655 / 14,437,225 | miss | **fathomed** (458,866,243) | **2,780** (300 + 2,480) | 1,389 | 59,366 (0.175) |

- **3 of 3 fathom**, including the node that was furthest from the cutoff in the 12-node
  sample (14.4M short at 1,500 plain iterations). Without cuts 2 of the 12 fathom (nodes 4
  and 6, at 600-630 iterations) and the other 10 end 2.7M-14.4M short.
- The two near misses fathom in about the iterations (and CPU) the two plain fathomers
  took: 628-680 against 600-630.
- Node 10's shortfall with cuts: 10.4M at 600 total iterations, 6.3M at 800, 3.2M at 1,100,
  1.8M at 1,500, 0.66M at 2,000. A plain depth-3 node this far away needs about four more
  levels of branching (14.4M / 1.8^4 = 1.4M); here it costs 2,780 iterations at one node.

## 6. Roots of the solved siblings (`logs/plain_tai*b.log`, `logs/r*b_comb.log`) [m]

300 uncut iterations, then 2,000 more: plain (`-cpplain`) or CP form with the `comb`
regions (`line:1,top:K,linegap:4:0.03`). Cutoff 2 x optimum so the runs go on.

| instance | optimum | plain DNN at 2,300 | plain gap | regions | violated at the plain 300-iteration iterate (max) | with cuts at 2,300 | gap left | **DNN gap closed** (equal iterations) | cut bound passes the plain 2,300-iteration value after |
|---|---|---|---|---|---|---|---|---|---|
| tai20b | 122,455,319 | 122,455,316 | 3 (exact) | 20 | 1,028 (0.062) | 122,455,319 (.59 under) | 0.4 | nothing to close | - |
| tai25b | 344,355,646 | 337,285,170 | 7,070,476 (2.1%) | 27 | 4,572 (0.098) | 344,084,304 | 271,342 (0.08%) | **96%** | 200 CP iterations |
| tai30b (`line:1`) | 637,117,113 | 597,633,906 | 39,483,207 (6.2%) | 29 | 9,016 (0.219) | 635,876,664 | 1,240,449 (0.19%) | **97%** (99.7% at 6,300: 636,993,730) | 100 |
| tai35b | 283,315,445 | 269,414,422 | 13,901,023 (4.9%) | 42 | 14,511 (0.123) | 279,325,031 | 3,990,414 (1.4%) | **71%** | 100 |
| tai40b | 637,250,948 | 608,514,446 | 28,736,502 (4.5%) | 48 | 24,230 (0.281) | 627,372,144 | 9,878,804 (1.6%) | **66%** (still rising 340k per 100) | 300 |
| tai50b (best known) | 458,821,517 | 430,147,161 at 1,500 | 28,674,356 (6.25%) | 61 | 89,148 (0.347) | 449,469,933 at 1,900 | 9,351,584 (2.0%) | **61%** at 1,500, 66% at 1,900 (still rising 380k per 100) | 300 |

- At equal iteration counts and against the plain plateau (the plain bound gains 0.05-0.1M
  per 100 iterations at the end on tai30b-40b, the cut bound is 10M-38M above it), so this is
  not the "more iterations" artefact.
- tai25b and tai30b (nearly collinear clusters) become almost exact. tai35b and tai40b, whose
  clusters are spread in two dimensions (first-axis share 0.84 and 0.64), close 66-71% with a
  region set built around one axis; tai35b has nearly flattened (279.33M), tai40b has not.
  `line:4` or the full dendrogram were not tried on them (no time): on tai30b `line:4` was
  the best set at 1,000 iterations. tai50b is of the collinear kind (first-axis share 0.97).

## 7. Nodes on the path to the best known solution (`onpath.txt`, `logs/onpath_*.log`) [m]

A node whose fixes agree with an optimal (best known) permutation has that value as its
optimum, so no certified value there may exceed it. Cutoff = optimum + 1,000 so the runs go
on; `comb` regions; 300 uncut + up to 2,500 CP iterations (tai50b: 15-minute cap).

| instance | depth, spec | certified at the end | optimum / best known | left |
|---|---|---|---|---|
| tai30b | 3, `18:27,20:11,26:24` | 637,055,611 (2,800 it) | 637,117,113 | 61,502 (0.010%) |
| tai30b | 6, `1:7,6:20,8:13,19:21,21:8,29:2` | 637,105,977 | 637,117,113 | 11,136 (0.002%) |
| tai30b | 10, `2:10,7:4,9:29,11:12,17:23,20:11,22:6,27:15,28:0,29:2` | 637,077,476 | 637,117,113 | 39,637 (0.006%) |
| tai50b | 3, `7:40,40:39,42:43` | 454,411,867 at 1,400 it (cap; +350k per 100) | 458,821,517 | 4,409,650 (0.96%) |

None exceeds the optimum; the tai30b ones flatten within 0.01% under it. The tai50b on-path
node is at 452.2M after the 1,000 iterations at which the off-path nodes 7 and 1 had long
fathomed (they pass 458.82M at 628-680), so the cuts are separating nodes, not lifting
everything.

## 8. Why it works, and why only on the location side

- A type-b distance matrix is, to 0.1%, a planar Euclidean metric concentrated on a few
  clusters, and on tai25b/30b/50b nearly a line metric. A line metric is a nonnegative
  combination of the cut semimetrics of its prefixes; a planar Euclidean one of half-plane
  cuts (Crofton). With z^S_i = [facility i is placed in S], the objective is then
  sum_S w_S sum_{i<j} f_ij (z^S_i + z^S_j - 2 z^S_i z^S_j): one weighted graph-partition
  objective per cut S, in exactly the moments (p, W) the regional triangles constrain.
  This is the same reason the grid prefixes work on nug/sko/wil/tho; type b just has far
  fewer cuts carrying the weight (5 clusters carry 98% of tai50b's distance sum) and
  heavy-tailed flows, so the DNN relaxation's "fractional partition" is far from Boolean:
  violations 0.2-0.35 at the plain iterate, against 0.004-0.02 for atomic triangles on nug30.
- Random location sets are not cuts of the metric (0.1%), balls are poor ones (32%), and
  the flow matrix has no such decomposition (facility-side families: 0 to 9 violated
  triangles, 0.0%).

## 9. What this implies for tai50b, and what it does not show

- The 12-node depth-3 sample of NOTES.md was 2 fathomed / 10 missed (2.7M-14.4M short at
  1,500 iterations) and underlies the "floor 1e11 CPU-s / tuned 4.1e9 +- 3.0e9" estimates,
  whose cost sits at depths 5-8. With the regions, the nearest two misses and the furthest
  one fathom at depth 3. If the other seven behave like these (they lie between them in
  shortfall; not run) [e], a uniformly drawn depth-3 node fathoms at depth 3, in 600-2,800
  iterations (300-1,400 CPU-s on this shared box).
- A frame, not an estimate [e]: bounding every one of the 50 x 49 x 48 = 117,600 depth-3
  nodes at about 700 CPU-s each would be 8e7 CPU-s (about one day on 1,000 cores), 50-1,000x
  under the current figures. The real tree is smaller above (depth-1/2 nodes that fathom:
  the root itself is at 2.0% and rising) and larger below (subtrees near optimal solutions
  cannot fathom at depth 3: the on-path depth-3 node is still 4.4M short at 1,400
  iterations). Only Knuth dives with the cut phase can size it; that is the next step.
- Not shown: any search run on tai50b; the root's plateau (the 26-minute cap stopped it at
  1,900 iterations while rising); how far below the best known value the strengthened
  relaxation sits near optimal solutions on tai50b (on tai30b: 0.02% at the root); whether
  458,821,517 is optimal.
- By-product: the tai50b root value 449,469,933 is a certified lower bound 18.4M above the
  431,090,700 that NOTES.md quotes as QAPLIB's best (gap 6.04% -> at most 2.04%). The same
  construction applies unchanged to the other open type-b instances (tai60b, tai80b, tai100b,
  tai150b).

## 10. In the search: tai30b solved with `cut.regions` (`bnb/tai30b_cut.par`, `logs/bnb_tai30b_cut.*`) [m]

`bnbdrv -p bnb/tai30b_cut.par -u 637117114` (incumbent = optimum + 1, so the run has to find
the optimum and prove it; cutoff 637,117,113.05): tai50b.par's settings with give-up 300,
DNN at every depth, and the cut phase at every DNN miss (`cut.zone = 1`, `cut.regions =
regions/tai30b.line1.txt`, `cut.cp = 100:1000`, `pool.on = 1`). The root is not bounded by
the DNN in the search (depths 1-28). Run in three pieces because of my time caps (`-o/-t`,
then `-r` twice; logs `.log`, `.log2`, `.log3`, node traces `.trace*`).

**Result: tai30b solved, 637,117,113 found and proved optimal in 363 nodes and 1,520 CPU-s**
(one core of the shared box; DNN bound 299 nodes, 250 fathomed, 709 s; cut phase 49 nodes,
22 fathomed, 801 s).

Depth 1 (the root branches on column 19, 30 children):

| fate | children | detail |
|---|---|---|
| fathomed by the plain DNN bound | 8 | 5-45 s each |
| plain DNN gave up at 300 iterations (623.4M-636.1M, i.e. 1.0M-13.7M short), **cut phase fathomed** | 21 | 22-46 s per node in all |
| not fathomed | 1 | `5:19`, the optimal child (p(5) = 19 in `QAPSOLN/tai30b.sln`): plain 623,179,768, cut phase 636,960,944, branched |

Below it **the tree is one path**: at each depth exactly one node branches, always the pair
of the optimal permutation (`11:12`, `22:6`, `25:17`, `8:13`, `2:10`, `0:3`, `9:29`, ...),
and the siblings fathom (levels 2-8: 23-29 nodes each, 91-97% fathomed, by the DNN
bound with the inherited pool in about 60 iterations, 1-2 s; QPB takes over from depth 9).
The optimum is found at depth 27. The on-path nodes' certified cut bounds are 636,960,944,
637,023,811, 637,047,851, 637,055,843, 637,094,293, 637,088,196, 637,110,195, 637,109,018
at depths 1-8 and 637,117,112.55 to 637,117,112.9999 from depth 15 on: never above the
optimum 637,117,113, and within 0.025% of it from depth 1.

- 57% of the CPU is depth 1 (865 s): 300 uncut iterations and then the phase at 22 nodes.
  Starting the phase earlier is the obvious saving. The on-path nodes' phases (35-80 s
  each at depths 1-9) cannot fathom and run to the stall rule; a give-up for nodes that end
  within 0.03% of the cutoff, or an incumbent heuristic, would remove most of that.
- No plain-search comparison was run (no time). For scale: the 22 depth-1 children that
  needed the cut phase were 1.0M-13.7M short at the plain give-up, where the plain bound
  gains about 0.05M per 100 iterations; without cuts each of them is a subtree.

## 11. Next (in order)

1. **Review the patch and rerun the validity set** (`validity.py`, the on-path nodes, the
   tai30b search from `-u 637117114`; then tai35b and tai40b solved by the
   search as end-to-end checks: they are solved instances, the search must return 283,315,445
   and 637,250,948).
2. **tai50b Knuth dives with the cut phase** (`cut.regions`, `cut.zone = 1`, CP form, pools):
   that is the estimate the current 4e9-1e11 CPU-s figures must be replaced by. Settings to
   tune first, all specific to this regime (bound far from converged at the give-up, gaps of
   1e6-1e7): start the cut phase early (the first 300 uncut iterations are mostly wasted:
   at the tai50b root the cut bound passes the plain 1,500-iteration value 300 iterations
   into the phase), the phase's iteration limit and stall/give-up rules (node 10 needed
   2,480 CP iterations and was still worth it), and a give-up for on-path nodes.
3. **The root of tai50b to its plateau** (a few hours of one core): how close the relaxation
   gets to 458,821,517 decides whether the tree is "one path" as on tai30b.
4. **Region sets**: fewer regions (`linegap` + `top`, 16 on tai50b) against `comb` (61) for
   speed; `line:4` / `hier` for the 2-D instances tai35b, tai40b (66-71% with one axis).
5. **Same idea elsewhere**: any instance whose distance matrix is (near) l1-embeddable has
   such cuts: the other open type-b instances (tai60b-tai150b), and Hamming or tree
   metrics. On the grid instances the row/column prefixes are already this construction.

## Files

- `regions.patch`: `git diff` of the worktree `/home/user/wt-regions` (dnncut.h/.C, dnncuts.C,
  bnb.C, settings.C, RUNNING.md).
- `gen_regions.py` (families, `-stats`), `regions/*.txt` (region files), `validity.py`
  (brute force), `run.sh` (one `dnncuts` run, filtered log), `tab.py` (tables),
  `sweep30.sh`, `nodes50.sh`, `tail.sh`, `rest.sh` (the queues, one solver process at a
  time), `onpath.txt`, `bnb/*.par`.
- `logs/`: every run; `structure.txt`, `validity_seed1.txt`, `table_*.md`.

Not done: the rounds (Polyak) form beyond the validity runs (it needs `-delta` scaled to
the objective, about 1% of it; the CP form was used for every measurement); families on
tai35b/40b beyond `comb`; tai50b depth-3 nodes other than 7, 1, 10; a plain-search
comparison on tai30b; no search or dive on tai50b.
