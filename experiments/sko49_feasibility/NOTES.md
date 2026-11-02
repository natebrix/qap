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
