# Is tai50b reachable with the certified DNN bound? (Nov 8 2026)

Companion to `../sko49_feasibility/NOTES.md` and `../wil50_feasibility/NOTES.md`,
run on the 2-core cloud Xeon @ 2.1 GHz (AMX; OpenBLAS `scipy-openblas32`, one thread
per process, two processes at a time). Every bound is a certified value printed by
`dnntest`. Scripts and logs: `run.sh`, `root/`, `depth3/`; `../sko49_feasibility/
summarize.py depth3 458821516.05` builds the depth-3 table.

## Summary

- Root DNN bound **430,147,161 certified, 28.7M (6.25%) below the best known
  458,821,517**. QAPLIB's best bound is BBCPOP's 431,090,700 (6.04%); ours is 0.2% under
  it at 1,500 iterations, still rising slowly.
- **2 of 12 random depth-3 nodes fathom** (at 600 and 630 iterations, 332-349 CPU-s); the
  other 10 stop 2.7M-14.4M short (0.58-3.1% of the cutoff, median 6.8M).
- From the root to depth 3 the shortfall falls 5.6x (28.7M to a median of 5.1M over all
  12), about **1.8x per level**. wil50's falls about 1.25x per level; tho40's halves. On
  sko49 none of 12 depth-3 nodes fathomed.
- So the relative root gap overstates tai50b, as it did tai40b (4.5% at the root, solved
  in 278,465 nodes by Fujii et al. 2021). Fujii et al. estimate 5.5e7 nodes for tai50b
  with their rule M. The tree size is what `estimate/` measures next.

## 0. The instance

- `qaplib/tai50b.dat`: n = 50, Taillard's "real-life-like" type b. The first matrix is
  symmetric (distances, 0-1,843, 98% nonzero); the second is asymmetric (flows, 0-8,071,
  43% nonzero, heavy-tailed: median nonzero flow 58, and the top 10% of nonzero flows carry
  61% of the total). bnbdrv and dnntest symmetrize the second matrix (objective
  values unchanged) and use fathoming tolerance 0.95 (integer, not even), **cutoff
  458,821,516.05**.
- **Best known 458,821,517** (`QAPSOLN/tai50b.sln`, Ro-TS); the permutation evaluates to
  458,821,517 with the first matrix as the distance (757,932,044 the other way round).
  Never proved optimal.
- Objectives are doubles in the solver (`solution.h`), so the 4.6e8 magnitude is no
  overflow risk; the largest possible sum, n² · 1,843 · 8,071 = 3.7e10, also fits.
- No grid and no symmetry file, so neither the symmetry file nor the regional cuts of
  wil50 and sko49 apply. Rule 9 branching, warm starts, pool inheritance and the
  per-node speed work do.

## 1. Root (`root/root.log`)

`dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u 458821517 -v`, the sko49 note's
settings: 421,577,777 at 400 iterations, 429,428,065 at 1,000, **430,147,161 at
1,500** (gaining about 1.4M per 1,000 iterations at the end); 833 ms per iteration
with both cores busy, 1,250 CPU-s.

## 2. Depth-3 sample (`depth3/`, 12 random nodes, seed 20261108)

Drawn as the sko49 and wil50 samples were (three random facilities paired with three
random locations, `depth3/draw_nodes.py`), same `dnntest` settings, cutoff
458,821,516.05. "At 400" is the certified bound at the search's give-up point.

| node | spec | lb at 400 | short at 400 | final lb | short | status | iters | CPU s |
|---|---|---|---|---|---|---|---|---|
| 1 | `20:19,28:24,14:35` | 449,539,545 | 9,281,971 | 454,973,581 | 3,847,935 | iter_limit | 1500 | 847 |
| 2 | `0:42,35:27,48:5` | 442,863,825 | 15,957,691 | 449,238,914 | 9,582,602 | iter_limit | 1500 | 902 |
| 3 | `27:0,44:23,6:19` | 449,930,773 | 8,890,743 | 455,582,587 | 3,238,929 | iter_limit | 1500 | 857 |
| 4 | `42:5,25:16,29:23` | 456,140,025 | 2,681,492 | 458,849,521 | fathomed | fathomed | 630 | 349 |
| 5 | `44:27,27:36,0:21` | 441,358,403 | 17,463,113 | 448,190,032 | 10,631,484 | iter_limit | 1500 | 818 |
| 6 | `12:31,27:41,25:11` | 455,943,022 | 2,878,494 | 458,823,193 | fathomed | fathomed | 600 | 332 |
| 7 | `0:44,49:49,25:31` | 451,179,509 | 7,642,007 | 456,142,196 | 2,679,320 | iter_limit | 1500 | 854 |
| 8 | `43:29,6:45,28:31` | 447,419,407 | 11,402,109 | 452,625,028 | 6,196,488 | iter_limit | 1500 | 845 |
| 9 | `4:41,35:32,24:5` | 448,550,160 | 10,271,356 | 454,731,218 | 4,090,298 | iter_limit | 1500 | 821 |
| 10 | `5:40,40:22,44:3` | 440,328,635 | 18,492,881 | 444,384,291 | 14,437,225 | iter_limit | 1500 | 807 |
| 11 | `17:0,31:16,9:33` | 443,884,812 | 14,936,704 | 451,515,658 | 7,305,858 | iter_limit | 1500 | 811 |
| 12 | `25:25,6:40,16:13` | 443,972,299 | 14,849,217 | 451,078,758 | 7,742,758 | iter_limit | 1500 | 872 |

- 538-601 ms per iteration at m = 47 with both cores busy; 810-900 CPU-s per unfathomed
  node at 1,500 iterations.
- The bound gains 4M-8M between iteration 400 and 1,500 at the unfathomed nodes, a larger
  share of the remaining shortfall than on sko49 (10-24 units of 134-337): these nodes are
  further from convergence at the give-up point, so the give-up iteration matters more here.
- Both fathomed nodes were 2.7M-2.9M short at 400 iterations, about the same as the
  closest unfathomed node at the end (node 7). The two classes overlap at the give-up point.

| | tai50b (12 nodes) | sko49 (12 nodes) | wil50 (16 nodes) |
|---|---|---|---|
| root gap, % of cutoff | 6.25% | 3.14% | about 1.4% |
| fathomed at depth 3 | **2 of 12** | 0 of 12 | 2 of 16 |
| unfathomed shortfall, % of cutoff | **0.58-3.1%, median 1.5%** | 0.57-1.44%, median 1.16% | 0.12-0.83% |
| shortfall, root to depth-3 median | **5.6x (about 1.8x per level)** | 2.7x (1.4x) | about 2x (1.25x) |

## 3. What this says, and what is next

- tai50b's bound closes much faster with depth than sko49's or wil50's. That fits its type:
  a minority of large flows carries most of the objective, and fixing the facilities that carry them settles most of
  the objective. Which three facilities a node fixes matters a lot (shortfalls 2.7M to
  14.4M at depth 3), so branching order should matter more here than on the grid instances.
- Per node, a tai50b DNN solve costs what a wil50 one does (same m, 0.55-0.6 s per
  iteration here).
- The tree size is not settled by 12 nodes. Next: the Knuth estimate in `estimate/` with
  the nug30 record recipe (`../nug30_solve/README.md`: rule 9 through depth 6, DNN at
  depths 1-7, give-up 200) and the DNN-marginal sampler with a 30% uniform share,
  per-level table after every dive.
