# Is sko49 reachable with the certified DNN bound? A numbers-based feasibility note (Nov 2 2026)

Companion to `../wil50_feasibility/NOTES.md`, written on a 4-vCPU cloud Xeon (OpenBLAS
`scipy-openblas32`, 1 thread per process, 4 processes at a time) in about 2.5 hours of
compute. Every bound below is a certified value printed by `dnntest` or `dnncuts` unless
it is marked *estimate*. Scripts and logs: `root/`, `depth3/`, `depth4/`, `cuts/`, `sym/`;
`summarize.py` builds the tables from the logs.

## 0. The instance

- `qaplib/sko49.dat` (already in the repo): n = 49, **the first matrix is the rectilinear
  distance on a 7 x 7 grid** (checked entry by entry: `d((r,c),(r',c')) = |r-r'| + |c-c'|`
  with row-major cells), the second is a symmetric random flow matrix (entries 0-10, 68%
  nonzero, zero diagonal). Same layout as wil50 (grid first), so the regional cuts take
  `-grid 7:7 -side 0`.
- **Best known 23,386**, confirmed from QAPLIB's `sko49.sln` (fetched from the Lehigh
  mirror; the permutation evaluates to 23,386 with the first matrix as the distance).
  The file is now in `QAPSOLN/sko49.sln`. Never proved optimal.
- `bnbdrv` prints `objective values are even integers: fathoming tolerance 1.95`, so the
  **cutoff is 23,384.05**. Root QPB bound 21,003.7 (gap 2,382, 10.2%).

## 5. Tree structure and symmetry (section numbers follow the task's measurement order)

Branching fixes one facility per level, so a depth-d node has 49 - d children (fewer where
the QPB reduced costs eliminate a location; `dnntest` runs on the plain node). The 7 x 7
grid distance is invariant under the dihedral group of order 8 (rotations and reflections),
against order 4 for wil50's 5 x 10 grid. Counting ordered d-tuples of distinct cells up to
symmetry (Burnside; `sym/make_sym.py` prints this):

| depth | plain nodes | inequivalent nodes | ratio | wil50 (order 4, estimate) |
|---|---|---|---|---|
| 1 | 49 | 10 | 4.9 | 15 |
| 2 | 2,352 | 315 | 7.5 | 635 |
| 3 | 110,544 | 13,923 | 7.9 | 29,580 |
| 4 | 5,085,024 | 636,048 | 8.0 | 1.4e6 |
| 5 | 2.29e8 | 2.86e7 | 8.0 | 6.4e7 |

So sko49's tree is **2.1 x narrower than wil50's at every level from depth 3**: the
symmetry (8 against 4) is worth 2 x, and n = 49 against 50 another 6%.

No `SYM/sko49.sym` existed. The format (`latest/branch.C`, `read_symmetry`) is J1 = the
locations to branch on at the root, then J2/J3 pairs: when every assigned location lies in
J2 (the fixed cells of a reflection), the children need only cover J3 (a fundamental domain
of that reflection); `SYM/nug30.sym` does exactly this for the 5 x 6 grid (J1 = the 3 x 3
corner block, J2 = the middle row, J3 = the top three rows). `sym/make_sym.py` writes the
7 x 7 version in a few lines (`sym/sko49.sym`: J1 = the 10 cells with row <= col <= 3; four
J2/J3 pairs, one per reflection, |J2| = 7, |J3| = 28). One caveat: **the solver's symmetry
acts on the columns (locations, the second matrix) and sko49's grid is its first matrix**,
so the file applies to the transposed instance `sym/sko49t.dat` (B then A; the optimal
value is unchanged, q = p^-1), on which the regional cuts would take `-side 1`. The
mechanism uses only the first matching J2 set at a node and gives nothing once a location
off every axis is assigned, so it realises the order-8 reduction at depth 1 (10 children)
and most of it at depth 2; below that the tree is the plain one under each inequivalent
depth-2 node, which is where the factor 8 of the table comes from.

## 1. Root DNN bound (`root/root.log`)

`dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u 23386 -v` (PRSM step, degree-4
filter, subspace tracker, one thread, four processes sharing the machine):

| | sko49 | wil50 (`../wil50_feasibility/NOTES.md`) |
|---|---|---|
| certified root DNN bound | **22,652.2** (primal 22,640.0 at 1,500 iterations; 22,626.7 at 400, 22,649.6 at 800: converged to within about 3) | about 48,130 (root DNN "1.4% below the optimum", gap about 680) |
| best known | 23,386 | 48,816 |
| gap, units | **733.8** | about 680 |
| gap, percent | **3.14%** | about 1.4% |
| root QPB bound | 21,003.7 (10.2%) | – |
| ms per iteration (m = 49, 1 thread, 4 processes running) | 823 | wil50 m = 50: 530-600 on the M2 Pro with `-f` |
| CPU for 1,500 iterations | 1,234 s | – |

The root gap in units is the same as wil50's, but sko49's objective is half the size, so
**relative to the cutoff the DNN relaxation is 2.2 x weaker on sko49 at the root**. (Fujii
et al. 2021 say that DNN-type lower bounds for sko42-sko81 were computed by Oliveira et al.
and by the Newton-bracketing method but do not print the sko49 values; not checked here.)

## 2. Depth-3 sample (`depth3/`, 12 random nodes, seed 20261102)

Drawn as the wil50 sample was (three random facilities paired with three random locations,
`depth3/draw_nodes.py`; 16 drawn, the queue stopped after 12 to keep the budget), each run
with `dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u 23386 -v` (PRSM, degree 4, one
thread, no give-up, three or four processes at a time), cutoff 23,384.05. "At 400" is the
certified bound at the search's give-up point; every value is certified. Table by
`summarize.py depth3`:

| node | spec | lb at 400 | short at 400 | lb at 1,500 | short | status | iters | CPU s |
|---|---|---|---|---|---|---|---|---|
| 1 | `36:40,4:41,24:10` | 23,234.8 | 149.3 | 23,249.6 | 134.4 | iter_limit | 1500 | 838 |
| 2 | `7:48,26:41,2:12` | 23,068.8 | 315.2 | 23,079.7 | 304.3 | iter_limit | 1500 | 855 |
| 3 | `25:29,33:12,39:38` | 23,061.0 | 323.1 | 23,084.2 | 299.9 | iter_limit | 1500 | 868 |
| 4 | `11:26,45:36,35:28` | 23,106.5 | 277.6 | 23,120.9 | 263.2 | iter_limit | 1500 | 846 |
| 5 | `37:34,20:42,12:32` | 23,036.8 | 347.2 | 23,047.2 | 336.9 | iter_limit | 1500 | 854 |
| 6 | `3:6,40:40,5:4` | 23,191.3 | 192.8 | 23,208.7 | 175.4 | iter_limit | 1500 | 843 |
| 7 | `33:0,43:1,34:29` | 23,101.1 | 282.9 | 23,124.7 | 259.3 | iter_limit | 1500 | 782 |
| 8 | `22:4,36:27,34:40` | 23,106.4 | 277.6 | 23,125.7 | 258.3 | iter_limit | 1500 | 792 |
| 9 | `46:7,48:42,38:14` | 23,163.7 | 220.3 | 23,173.7 | 210.4 | iter_limit | 1500 | 791 |
| 10 | `43:6,32:18,38:36` | 23,091.1 | 292.9 | 23,104.3 | 279.7 | iter_limit | 1500 | 777 |
| 11 | `33:7,34:18,27:41` | 23,086.9 | 297.2 | 23,105.8 | 278.2 | iter_limit | 1500 | 777 |
| 12 | `40:46,36:17,10:20` | 23,051.4 | 332.7 | 23,069.6 | 314.4 | iter_limit | 1500 | 776 |

12 done: fathomed 0, short values 134, 175, 210, 258, 259, 263, 278, 280, 300, 304, 314, 337, median 271

Primal and certified values agree within 3-20 at 1,500 iterations, so these are the
relaxation values to within a few units (the bound gains only 10-24 between iteration 400
and 1,500; wil50's depth-3 nodes gained 26-36 cold). Per iteration 520-580 ms at m = 46 with
four processes running; 780-870 CPU-s per converged node.

| | sko49 (12 nodes) | wil50 (16 nodes, `../wil50_feasibility/NOTES.md` sections 2, 7a) |
|---|---|---|
| fathomed at depth 3 | **0 of 12** | 2 of 16 |
| converged near misses (within about 120 of the cutoff) | **0** | 5 (59-116 short) |
| clear misses | **12 of 12: 134-337 short, median 271** | 9: 185-405 short, median about 240 |
| shortfall as a fraction of the cutoff | **0.57-1.44%, median 1.16%** | 0.12-0.83%, clear misses 0.4-0.8% |
| give-up point (400 iterations) | 149-347 short | 300-515 short at 300 iterations (`-f`, ADMM) |

In units the sko49 depth-3 nodes look like wil50's clear misses (median 271 against 240),
but there is no fathomed or near-miss class at all, and relative to the cutoff the
shortfall is twice wil50's: every one of these nodes would be given up by the search and
branched.

## 6b. Fujii, Ito, Kim, Kojima, Shinano, Toh (arXiv 2101.09629, ZIB report 21-02)

Fetched as data (`curl` through the proxy, read with `pdftotext`; not in the repo). Their
solver is a parallel branch and bound (UG framework) over a DNN-type relaxation solved by
their Newton-bracketing / BBCPOP methods, with three branching rules M (fix a facility to
each location), P and D.

- **Table 4 (solved instances; branching rule M):** sko42, optimum 15,812, **6,019,419
  nodes, 5.12e5 s of execution in parallel on 5,184 cores** (about 2.65e9 core-seconds, 84
  core-years, 5.9 days of wall time); nug30 took 26,181 nodes and 3.14e3 s on 1,728 cores,
  tai30a 34,000,579 nodes and 5.81e5 s on 1,728 cores.
- **Table 5 (estimated total nodes, by sampling the tree the way our model does):**

| instance | rule M | rule P | rule D | their comment |
|---|---|---|---|---|
| sko42 | 3,575,067 (6,019,419 generated) | 3,071,294 | 2,764,606 | solved |
| **sko49** | **1.6e8** | **1.3e9** | **1.2e9** | "we may challenge sko49 (M)" |
| wil50 | 32,963,810 | 20,523,849 | 39,330,160 | "could be solved" (P) |
| tho40 | 7.1e11 | 37,382,192 | 97,834,698 | "could be solved" (P) |
| tai35a | 2.1e11 | 3.9e10 | 1.3e10 | too difficult |

They note that the nodes actually generated were at most 3 x the estimate on the solved
instances. So by their own count **sko49 is 5 x (rule M) to 60 x (rules P, D) more nodes
than wil50**, and 27 x sko42, which already took 5,184 cores for six days. Their relaxation
is the Lagrangian DNN relaxation (the same relaxation our certified bound solves, so the
same strength up to solver accuracy) and the paper does not mention symmetry, so their
tree is the plain one; their node estimates are an independent check on section 6 once
divided by the symmetry factor 8 (sko49) or 4 (wil50).

## 3. Per-level closure: depth-4 children of two misses (`depth4/`)

Eight random children each of node 1 (134.4 short at convergence) and node 2 (304.3 short):
one random free facility per parent, eight random free locations
(`depth4/draw_children.py`, seed 20261102), same `dnntest` settings as depth 3. A child
that fathoms stops when its certificate clears the cutoff, so its gain is a lower bound.
Table by `summarize.py depth4`:

| child | spec | lb at 400 | short at 400 | lb at end | short | status | iters | CPU s |
|---|---|---|---|---|---|---|---|---|
| 1 | `36:40,4:41,24:10,39:4` | 23,338.7 | 45.4 | 23,352.9 | 31.1 | iter_limit | 1500 | 742 |
| 2 | `36:40,4:41,24:10,39:25` | 23,357.2 | 26.9 | 23,373.1 | 10.9 | iter_limit | 1500 | 706 |
| 3 | `36:40,4:41,24:10,39:43` | 23,372.3 | 11.8 | 23,384.1 | -0.0 | fathomed | 1260 | 591 |
| 4 | `36:40,4:41,24:10,39:44` | 23,335.3 | 48.8 | 23,348.3 | 35.7 | iter_limit | 1500 | 698 |
| 5 | `36:40,4:41,24:10,39:11` | 23,347.1 | 36.9 | 23,359.3 | 24.8 | iter_limit | 1500 | 693 |
| 6 | `36:40,4:41,24:10,39:7` | 23,382.3 | 1.8 | 23,384.2 | -0.1 | fathomed | 420 | 193 |
| 7 | `36:40,4:41,24:10,39:27` | 23,385.2 | -1.1 | 23,385.2 | -1.1 | fathomed | 240 | 105 |
| 8 | `36:40,4:41,24:10,39:2` | 23,353.4 | 30.7 | 23,367.2 | 16.9 | iter_limit | 1500 | 697 |
| 9 | `7:48,26:41,2:12,44:13` | 23,218.9 | 165.1 | 23,232.4 | 151.7 | iter_limit | 1500 | 700 |
| 10 | `7:48,26:41,2:12,44:26` | 23,210.5 | 173.6 | 23,222.8 | 161.3 | iter_limit | 1500 | 691 |
| 11 | `7:48,26:41,2:12,44:34` | 23,214.5 | 169.6 | 23,223.3 | 160.8 | iter_limit | 1500 | 696 |
| 12 | `7:48,26:41,2:12,44:40` | 23,332.0 | 52.1 | 23,347.8 | 36.3 | running | 1040 | nan |

11 done: fathomed 3, short values -1, -0, -0, 11, 17, 25, 31, 36, 152, 161, 161, median 25
DEPTH4_SUMMARY

## 4. Regional cuts on the near misses (`cuts/`)

The wil50 step-6/7 protocol from the give-up point, `dnncuts -x <node> -u 23386 -grid 7:7
-side 0 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500`: 400 uncut
PRSM iterations, 12 rounds of 300 with 1,000 new regional triangles each (the row and
column prefixes of the 7 x 7 grid, 12 regions), uncapped Polyak step toward best + 50, then
a quiet round of up to 1,500 iterations from the best round's multipliers. The `five` arm
adds `-five 500` (500 five-point count inequalities a round). Every value certified.

| node (short at convergence) | arm | uncut lb at 400 | best Polyak round | after the quiet round | gain over the give-up point | fathomed | CPU s |
|---|---|---|---|---|---|---|---|
| 1 `36:40,4:41,24:10` (134.4) | tri | 23,234.8 (149.3 short) | 23,379.3 (+144.5, round 12) | **23,384.23** | **+149.5** | **yes** (380 quiet iterations) | 2,310 |
| 1 | five | 23,234.8 | 23,378.2 (+143.5, round 8) | **23,384.21** | **+149.4** | **yes** (420 quiet iterations) | 2,320 |
CUTS_6_9

Node 1 is a depth-3 node that converges 134 short: its uncut subtree is 46 depth-4
children, five of eight of which miss (section 3) and would each be branched into 45
depth-5 children. The cut phase fathoms it at depth 3 in 2,300 CPU-s, about 3 x the uncut
converged solve and 11 x the give-up budget. Compared with wil50: there the same protocol
gained +117 to +134 on depth-4 children 62-99 short and +122 on depth-3 node 14 (84
short), and fathomed nothing more than about 85-100 short at convergence; here it gains
**+149.5 on a node 134 short**, i.e. the cuts are worth 10-25% more on the 7 x 7 grid,
and a full level or more (sections 3 and 5). The five-point family is again a wash (its
rounds are noisier, with a -94 overshoot in round 9, and the quiet round lands at the same
value). Rounds cost 150-170 s (300 iterations at m = 46), the same as the uncut solve
per iteration.

## 6. Putting it together: the tree model (`model.py`, estimates)

The `../wil50_feasibility/NOTES.md` section-5 model, applied to the depth-3 sample: a
depth-3 node that converges s short fathoms k = round(s / g) levels lower (at least one when
s > 0), g the gain per level; its subtree is 46 + 46 x 45 + ... DNN nodes; 13,818
inequivalent depth-3 nodes (symmetry 8). Extra bound strength X per node (a cut phase at
every DNN node) is s -> s - X. Costs on this 2.1 GHz Xeon, one thread: 250 CPU-s per plain
DNN node (the give-up at 400 iterations costs 210; fathoms 105-590), 2,300 with the measured
cut phase. All estimates; the inputs are the 12 converged shortfalls of section 2.

`python3 model.py 49 8 <g> 250 2300 X=0,130,149,260 -- 134 175 210 258 259 263 278 280 300 304 314 337`

| g per level | extra bound per node | DNN nodes | CPU-s | CPU-years | 1,000 cores | this 4-core box |
|---|---|---|---|---|---|---|
| 130 (wil50's) | none | 1.3e8 | 3.2e10 | 1,025 | 1.0 years | 250 years |
| 130 | +130 (a full level, cut cost) | 3.0e6 | 7.0e9 | 221 | 81 days | 55 years |
| 130 | **+149 (the measured cut gain)** | **6.0e5** | **1.4e9** | **43** | **16 days** | 11 years |
| 130 | +260 | 3.8e5 | 8.8e8 | 28 | 10 days | 7 years |
| 145 (measured here, section 3) | none | 2.2e7 | 5.5e9 | 175 | 64 days | 44 years |
| 145 | +149 | 6.0e5 | 1.4e9 | 43 | 16 days | 11 years |
| 120 (pessimistic) | none | 3.4e8 | 8.5e10 | 2,686 | 2.7 years | 670 years |
| 120 | +149 | 3.0e6 | 6.9e9 | 217 | 79 days | 54 years |

The same code on wil50's 16-node sample (`model.py 50 4 130 250 2300 -- -1 -2 59 74 79 83
116 405 202 242 244 329 198 233 182 184`, the six unconverged clear misses entered at their
300-iteration value minus 120 as NOTES.md section 7a does) gives 3.9e8 DNN nodes plain
(NOTES.md: 4e7 for the median subtree, "nearer 1e8" with the depth-6 tail; the rounding
puts two more sampled nodes in the depth-6 class), 8.8e6 at +130 and 2.0e5 at +260.

Reading the table against Fujii et al.: their 1.6e8 (rule M, no symmetry) is 2e7 after the
symmetry factor 8, inside the range of the plain rows here (2e7-3e8). Their wil50 estimate,
3.3e7 plain or 8e6 after symmetry 4, is far below this model's plain wil50 row (3.9e8) and
below NOTES.md's 4e7-1e8 as well, so the two models agree on sko49 and disagree on wil50;
the reason is not established (different branching, a tree sampled rather than modelled from
a depth-3 shortfall). Either way every row without a cut phase is a 1,000-core job of months
to years, and the only rows under a month on 1,000 cores need the regional cuts at every
DNN node.

## 7. Verdict

VERDICT
