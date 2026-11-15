# Regional triangle cuts on non-grid regions: tai*b (bound_lab agent "regions", Nov 15 2026)

Question: tai50b has no grid, so it runs with no cut phase. Do the regional triangle
inequalities of `latest/dnncut.h`, which are valid for any region, bite when the regions
are built from the structure of Taillard's type-b distances (clusters, half-plane prefixes)?

**Answer: yes, very strongly, with regions on the location side that separate the
clusters.** Numbers below are certified values; [m] = measured (log named), [e] = estimated.
Status: DRAFT, being filled as runs finish (see "Status" at the end).

## Headline numbers

(filled at the end)

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
  `tests/settings_doc.py` passes). Compiled, **not yet run in a search** (see "Next").
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
