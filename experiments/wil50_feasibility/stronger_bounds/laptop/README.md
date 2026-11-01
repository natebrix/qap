# Laptop runs, Nov 1 2026: the quiet round fathoms the depth-4 near misses

M2 Pro (Accelerate, `VECLIB_MAXIMUM_THREADS=1`), 5-8 processes at a time, so
CPU seconds here are 1.3-1.6x what one process alone would take (the machine
was swapping for part of the afternoon); comparisons within a table are fair,
comparisons to the cloud's Xeon numbers only roughly (laptop ≈ 0.5x Xeon CPU).
Scripts: `dual/run_laptop_dual*.sh` (dnncuts), `search/run_laptop_ab.sh`
(bnbdrv -F, the cut phase in the search). Node ids as in `../dual/run_dual.sh`:
all32 line numbers for depth 3, `d<k>` = row k of `../depth4/nodes.txt`.

## 1. The cut phase in the search on the four depth-3 near misses (`search/`)

Same settings as `../search/run_ab.sh` (PRSM + degree 4, `nug_deep.par`,
`-F node -x 0:4:1:300:1`, zone 0.002, grid 5:10). Base: the node gives up and
its 47 children fathom. `nocap`: `QAP_DNN_CUT_NOCAP=1 QAP_DNN_CUT_FREEZE=600`.

| node (short at convergence) | base, 48 nodes | cut phase, capped | cut phase, uncapped + quiet 600 |
|---|---|---|---|
| 6 (56) | 1,790 s | 713 s (2.5x) | **394 s (4.5x)** |
| 9 (44) | 1,902 s | 572 s (3.3x) | **268 s (7.1x)** |
| 15 (60) | 2,279 s | 683 s (3.3x) | **321 s (7.1x)** |
| 16 (0.003 short in step 2) | 2,136 s | 978 s (2.2x) | 1,200 s (1.8x): the phase did **not** fathom (48 nodes; phase ~510 s, children 692 s) |

On 6, 9 and 15 every phase fathoms its node. Uncapping the Polyak target
removes the crawl at 48,814.0x: the capped phases spend 4-6 rounds there, the
uncapped ones fathom in 2-3 rounds, **4.5-7x** against the capped 2.3-3.3x.
Node 16 is the exception: capped it fathoms (978 s), uncapped it does not, and
the 600-iteration quiet round (from the last round's multipliers, pre-fix
binary) did not recover it, so its 47 children ran: 1,200 s, still 1.8x over
the base. The larger uncapped steps overshoot more; the quiet round from the
best round's multipliers (section 2) is the intended remedy and node 16 is
the node to re-test it on. A depth-4 child costs 35-46 CPU s here (base DNN time
over 47 children), a depth-3 give-up 70-130 s.

## 2. The quiet round at depth 4: four of four "unreachable" children fathom (`dual/`)

`research_bound_quality.md` section 1.2 predicted that the depth-4 children
62-109 short at convergence (seven of the sixteen in `../depth4/`, none
fathomed by the cut phase) are held back by the certificate's lag behind the
Lagrangian value at 300-iteration rounds, not by the inequalities. Protocol
from the give-up point as in step 3 (`-I 400 -i 300 -delta 50 -k 1000 -nocap`,
no stall rule), then a quiet round (multipliers frozen, no new cuts).

| child (short at convergence / at the give-up point) | cut rounds | quiet round | certified at the end | iterations, CPU | result |
|---|---|---|---|---|---|
| d1 (62 / 117) | 12 | 860 of 2,000, from round 12 | 48,814.34 | 4,860, 1,215 s | **fathomed**; before the quiet round 48,764 (50 short), primal 48,847 |
| d1 | 4, quiet from the *last* round (an overshoot, −72) | 1,500 | 48,724 | 3,100, 729 s | not fathomed: the overshoot's multipliers |
| d1 | 6, quiet from the last round | 1,500 | 48,797.6 | 3,700, 850 s | 16 short; the round-6 multipliers' value is below the cutoff |
| d1 | 4, quiet from the **best** round (3) | 1,500 | 48,814.04 | 3,100, 742 s | 0.009 short (primal 48,821) |
| d1 | 3, quiet from the best round (3) | 1,510 of 2,500 | 48,814.12 | 2,810, **689 s** | **fathomed** |
| d4 (66 / 123) | 6, quiet from the best round | 1,570 of 2,500 | 48,814.06 | 3,770, 856 s | **fathomed** |
| d6 (99 / 134) | 6, quiet from the best round | 960 of 2,500 | 48,814.25 | 3,160, 761 s | **fathomed**, gain +134 = one level |
| d11 (66 / 124) | 12 (two overshoots), quiet from the best round (9) | 1,020 of 2,000 | 48,814.62 | 5,020, 1,048 s | **fathomed** |

* The inner lag at depth 4 with 300-iteration rounds is **50-60 units**: the
  certificate sits that far below the Lagrangian value of the multipliers it
  was computed with, and 900-1,600 quiet iterations recover it. Together with
  the cuts' own +50-95 that is +117 to +134 on these children, a full level.
* **The quiet round must start from the best round's multipliers.** From the
  last round's, an overshoot (which happens on most nodes, −35 to −90) freezes
  a bad multiplier and the quiet round converges to its value (d1 rows 2-3,
  and the two depth-3 step-A runs below). `DNNCutParams::freeze_best`
  (default 1; `dnncuts -freeze-last` for the old behaviour) keeps the
  multipliers and iterate of the best round and restores them before the quiet
  round. Only the search's `QAP_DNN_CUT_FREEZE` runs after this change use it
  (bnbdrv rebuilt after the search runs above finished; those used the last
  round, which happened not to be an overshoot on nodes 6, 9, 15).
* Cost: 690-860 CPU s for a fathom with 3-6 rounds against about 1,500 s for
  the child's 46 uncut depth-5 children (35 s each): **about 2x on the
  dominant class**, against the 10-25% of `../README.md` step 4. Fewer rounds
  is better as long as the best round's multipliers are kept: 3 rounds
  fathomed d1 at 689 s, 12 rounds at 1,215 s.
* Zone: these children were 117-134 short at the give-up point, outside the
  search's `QAP_DNN_CUT_ZONE=0.002` (98 units). Reaching them needs about
  0.003, and then the stall rule has to carry the failures in that band.

## 3. What a failed phase costs (`dual/clear_*`)

Clear misses from the give-up point with the stall rule (3 rounds, factor 3)
and a 600-iteration quiet round: depth-3 nodes 10 (158 short at convergence)
and 13 (182), depth-4 child d2 (109).

| node | stalled after round | best gain | end | iterations, CPU |
|---|---|---|---|---|
| 10 | 4 | +70 | 131 short | 2,200, 617 s (quiet from the last round, pre-fix) |
| 13 | 5 | +106 | 112 short | 2,500, 683 s |
| d2 | 4 | +77 | 80 short | 2,200, 576 s |

About 600-700 CPU s each, a third of a clear miss's 47-child subtree. Nodes
10 and 13 were 201 and 217 short at the give-up point, d2 157: none would
enter a 0.003 zone (146 units), so in the search these costs are not paid.
The band that pays is 100-146 short at the give-up point, where this sample
has only fathoms (d1, d4, d6, d11); a larger depth-4 sample in that band is
the thing to run next.

## 4. Step A at depth 3: spoiled by the round-8 overshoot

`stepA_2`, `stepA_14` (step-2 protocol, `-freeze 8:2000 -nocap`, quiet from
the last round): on both nodes round 8 overshot (−125, −19) and the quiet
round converged to the overshoot's value (48,577 and 48,802 against bests of
48,700.5 and 48,808.1). The reruns with `freeze_best` were stopped when the
machine ran out of memory for other reasons; the two-hour protocol is parked.
At 700-iteration rounds the inner lag is small anyway (primal minus certified
11 and 17 at the best round); the quiet round is a 300-iteration-round tool.
