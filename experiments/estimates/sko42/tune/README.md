# sko42: settings for an estimate (cloud box, Nov 11-12)

sko42 (6 x 7 grid, QAPLIB optimum 15,812, solved by Fujii et al. with 6.0e6
nodes) has one search on record, Nov 8 (`../../../sko42_solve/`). It used
the nug30 record recipe with rule 8 and was stopped after 5.5 h with its
frontier still widening. That recipe has no cut phase, although sko42 is a
grid instance like wil50 and sko49, whose regional cuts are the strongest
lever measured on either.

## The arms

| arm | files | what it is |
|---|---|---|
| A | `A.par`, `A.sh` | `../sko42.par` until Nov 12, when it became C: the nug30 recipe with rule 9, DNN at depths 1-7, give-up 200, no cut phase |
| B | `B.par`, `B.sh` | wil50's full run (`../../../wil50_feasibility/fullrun/fullrun.par`) with `cut.grid = 6:7`: CP cut phase, `cut.zone = 0.0075` (119 short on sko42), DNN at depths 1-4 then 20 levels below a miss, give-up 300, the two-phase beta schedule |
| C | `C.par`, `C.sh` | B with `cut.zone = 0.010` (158 short): **the best of the three** (below) |

The rule blocks above `[settings]` are the same in both (and in
`fullrun.par`). B differs from `fullrun.par` only in the grid line. The grid:
`sko42t.dat`'s second matrix is the rectilinear distance on 6 rows of 7
cells, row-major (checked entry by entry), the same convention as wil50t's
5:10 and sko49t's 7:7, so `cut.side = 1`. Both arms pick the same root line,
facility 36 across the 12 cells of the symmetry's fundamental domain
(cells 0-3, 7-10, 14-17).

## Depth 1: the first solve (`d1/`)

Each depth-1 node bounded alone (`bnbdrv -F 36:c -d 1`, one process per
node, four at a time). The cutoff is 15,810.05; "short" is the cutoff minus
the first solve's bound.

| cell | A short | A iterations | A CPU s | B short | B iterations | B CPU s | B tighter by |
|---|---|---|---|---|---|---|---|
| 0 | 318.6 | 260 | 51.6 | 231.4 | 300 | 54.1 | 87.2 |
| 1 | 440.5 | 230 | 44.2 | 307.5 | 300 | 53.5 | 133.1 |
| 2 | 496.3 | 220 | 41.7 | 334.5 | 300 | 54.4 | 161.7 |
| 3 | 504.8 | 220 | 42.3 | 341.7 | 300 | 54.1 | 163.1 |
| 7 | 457.3 | 220 | 41.8 | 303.2 | 300 | 53.2 | 154.1 |
| 8 | 519.9 | 220 | 41.9 | 357.3 | 300 | 52.9 | 162.5 |
| 9 | 551.2 | 220 | 41.7 | 379.9 | 300 | 52.8 | 171.3 |
| 10 | 560.5 | 220 | 41.6 | 385.6 | 300 | 52.2 | 174.9 |
| 14 | 483.5 | 220 | 41.3 | 326.1 | 300 | 52.9 | 157.4 |
| 15 | 543.2 | 220 | 42.1 | 375.0 | 300 | 52.8 | 168.2 |
| 16 | 546.0 | 230 | 44.0 | 391.0 | 300 | 52.4 | 155.0 |
| 17 | 552.7 | 230 | 45.0 | 393.0 | 300 | 52.3 | 159.7 |

B's first solve is 87-175 tighter (median shortfall 349 against 512) for 23%
more CPU per node (53.1 s against 43.3 s on the cloud box). No depth-1 node
is within the cut zone, so this is the first solve alone: B's beta schedule
and its 300 iterations against A's 200. Whether B's cut phase pays for
itself deeper in the tree is what complete subtrees measure.

## Stages A and B under B (`sampling/`)

The wil50 sampling plan (`../../../wil50_feasibility/sampling/README.md`)
with `B.sh`, on the cloud box, 4 slots:

```
export H=experiments/estimates/sko42 R=<run dir>
SAMPLING_SETTINGS=$H/tune/B.sh SAMPLING_PARAMS=$H/tune/B.par \
  python3 experiments/wil50_feasibility/sampling/stage_a.py --run $R/stageA --slots 4 --cpus 0,1,2,3
SAMPLING_SETTINGS=$H/tune/B.sh SAMPLING_PARAMS=$H/tune/B.par \
  python3 experiments/wil50_feasibility/sampling/stage_b.py --stage-a $R/stageA --run $R/stageB --n 250 --converge 40
```

**Stage A** (`sampling/stageA_summary.md`; 2 h 53 min, 40,762 CPU s): every
node at depths 0-2.

| depth | nodes | fathomed | of which by the cut phase | CPU s per node | children kept |
|---|---|---|---|---|---|
| 0 | 1 | 0 | | 0.6 | 12 |
| 1 | 12 | 0 | | 53.6 | 438 |
| 2 | 438 | 164 (37%) | 123 | 91.6 | 10,726 |

- **The cut phase is what fathoms at depth 2.** The first solve fathomed 41 of
  the 438 nodes. The cut phase ran on the 138 within the zone (119 short) and
  fathomed 123 of them (89%), from a median first-solve shortfall of 65. It
  cost 13,300 CPU s for those 123, 108 s a node, against branching each into
  40 children at about 28 s each before any grandchild.
- **The rest:** 6 hit the phase's iteration limit at 0.6-4.3 short, and 9
  stalled at 11-35 short.
- **Survivors:** 274 depth-2 nodes, median 181 short, none farther than 292,
  leave a frontier of 10,726 depth-3 nodes.
- Three lines of the summary are wil50's constants, not sko42's: "the cut zone
  is 244" (wil50's 0.005 zone), "expected 15 orbits" and the "10 x 49 + 5 x 29"
  check.

**Stage B** (`sampling/stageB_summary.md`): 250 frontier nodes, seed
20261107, each to its first solve's give-up (300 iterations, warm from its
parent, the cut phase held off at depth 3):

| fathomed in the first solve | 0-50 short | 50-100 | 100-150 | 150-210 | >= 210 |
|---|---|---|---|---|---|
| 123 (49%) | 55 | 44 | 23 | 5 | 0 |

- **The misses are close.** Median 62.9 short, 90th percentile 124.0, the
  farthest 188.8. wil50's farthest was 338.4.
- **The cut phase covers almost all of them.** 114 of the 127 misses are within
  its zone, where production runs it at depth 3.
- **Converging the farthest four** (157.7-188.8 short, 4,000 iterations, 527 s
  each) gained only 3.1-5.3 units, under 3% of the shortfall, so the stage
  was stopped there: the estimate indexes
  subtrees by the give-up shortfall, which is what production sees at a node.
- **The run.** The first launch bounded none of the 250: the solver's settings
  check on `-r`, added after wil50's stage B, refused the intended
  `cut.maxdepth 2`. Fixed in `stage_b.py`, which passes
  `--allow-settings-change`; the solver log lists `cut.maxdepth=2` as the only
  difference.

## Anchors: complete subtrees (`subtrees/`)

Ten frontier nodes from stage B's sample, spread over the give-up shortfalls of
its 127 misses, each run to the end with `../../subtree.sh` (cold from `-F`,
through the coordinator, on the cloud box, one or two slots each):

```
EST_RUNS=<run dir> experiments/estimates/subtree.sh experiments/estimates/sko42/tune/B.sh Bs189 7:35,15:0,36:15 2
```

DNN CPU counts first solves and cut phases. "Cold root" is the subtree root's
own first solve, without its parent's warm start. Stage B's shortfall is warm.

| node | spec | stage B short | cold root (B) | A: nodes, CPU s | B: nodes, cut phases, CPU s | C: nodes, cut phases, CPU s |
|---|---|---|---|---|---|---|
| s20 | `9:41,12:0,36:15` | 20.2 | 20.2 | 40, 217 | 1, 1, 83 | |
| s63 | `21:0,22:23,36:9` | 62.9 | 69.2 | 40, 200 | 1, 1, 93 | |
| s94 | `0:0,36:16,40:41` | 93.5 | 99.7 | 116, 588 | 1, 1, 511 | |
| s108 | `22:9,26:35,36:10` | 108.4 | 129.4 | | 40, 4, 494 | |
| s124 | `7:41,12:0,36:15` | 124.0 | 136.4 | 417, 1,722 | 40, 7, 609 | 40, 1, 373 |
| s136 | `22:0,36:16,41:41` | 136.4 | 145.7 | | 154, 13, 2,476 | 116, 6, 1,871 |
| s145 | `7:41,22:0,36:8` | 145.4 | 153.1 | | 78, 14, 1,744 | 78, 7, 1,109 |
| s155 | `15:41,22:0,36:17` | 155.2 | 167.5 | | 40, 12, 992 | |
| s168 | `7:0,22:13,36:9` | 168.2 | 188.3 | | 229, 25, 4,746 | 229, 21, 3,628 |
| s189 | `7:35,15:0,36:15` | 188.8 | 205.9 | stopped open at 2,100, 8,619 | 191, 34, 4,474 | 191, 25, 3,549 |

Under A, the cold roots are looser everywhere except s20 (A gives up at 200
iterations and uses its own beta): 89.7, 163.5, 184.0 and 264.3 short on s63,
s94, s124 and s189.

**A against B: B wins on every node.** B's subtree takes 1.15x (s94) to 2.8x
(s124) less CPU. On the farthest node, A was stopped with its tree still open
at depth 9, already at 1.9x B's whole subtree. A's nodes are cheap (4-5 s,
about 30 iterations each) but numerous. B's cut phase prunes a node that A
would branch. In s189, 31 of B's 34 phases prune their node, at a median of
34 s each. Stage A adds the top of the tree, where B's advantage is wider
still: at depth 1 B is 87-175 tighter, and at depth 2 its cut phase fathoms
123 of 438 nodes.

**C against B: zone 0.010 is 21-39% cheaper on all five nodes tried.**
- **s124, s136 and s145.** Their cold roots (136-153 short) are outside B's
  zone, so B branches them without a phase. C runs the phase, which stalls
  (at 34, 65 and 69 short). But the children start from the tighter bound
  and the inherited pool, and need far fewer phases of their own: 1, 6 and 7
  against B's 7, 13 and 14.
- **s168 and s189.** Their roots (188 and 206 short) are outside both zones,
  and C gives the same tree as B, node for node by depth. It runs fewer
  phases (21 and 25 against 25 and 34) for 24% and 21% less CPU. A depth-4
  node 119-158 short that runs a phase under C passes a pool to its
  children, which then fathom in their first solve where under B they
  needed a phase.
- **Below 119 short** the arms differ only in descendants 119-158 short;
  s20, s63 and s94 have none, so their trees under C are B's. s108 and s155
  were not run under C: the estimate under C uses B's, an upper bound if C
  does there what it did on the five.
- **Depth 2** is in stage A under C (next section).

## The estimate under B (`estimate.py`)

```
python3 experiments/estimates/sko42/tune/estimate.py subtrees/Bs20=20.2 ... subtrees/Bs189=188.8
```

Each frontier node in stage B's sample takes the anchors' size at its give-up
shortfall: 1 node and its own first-solve CPU if fathomed there; otherwise
log-linear interpolation between the anchors, or a log-linear fit over them
with the mean correction for their scatter. Frontier 10,726; above it, stage
A's 451 nodes and 40,762 CPU s; 95% bootstrap over the 250 sampled nodes.

| quantity | model | mean per frontier node | total (95%) |
|---|---|---|---|
| DNN nodes | interpolation | 7.9 | 84,845 (52,347-122,068) |
| DNN nodes | fit (x49.6 per 100 short, scatter x3.0) | 11.7 | 126,224 (73,468-195,249) |
| CPU s | interpolation | 198.6 | 2.17e6 (1.54e6-2.91e6) = 603 core-hours |
| CPU s | fit (x14.1 per 100 short, scatter x1.7) | 220.3 | 2.40e6 (1.75e6-3.15e6) = 668 core-hours |

**sko42 under B: about 1e5 DNN nodes and 600-670 core-hours on the cloud
box's Xeon (430-870 at the ends of the intervals).**
- **On the laptop.** The cloud box takes about twice the laptop's time per
  DNN iteration (0.33 s against 0.16 s on wil50's a220, `../../README.md`,
  "Cold against warm").
  That makes it about 300-335 core-hours there: some 2-2.5 days on the 6
  performance cores, if the coordinator keeps them busy as it did here (99%).
  On this box's 4 cores it is about a week.
- **Against the yardstick.** Fujii et al. solved sko42 with 6.0e6 nodes. The
  Nov 8 run (the nug30 recipe, rule 8) bounded 9,651 nodes in 5.5 h with its
  frontier still widening.
- **How the work splits.** By the interpolation, the 5 sampled nodes 150-189
  short (2% of the frontier) carry 26% of the frontier's CPU, the 23 at
  100-150 41%, the 44 at 50-100 21%, the 55 under 50 9%, and the 123
  fathomed in their first solve 3%.

What this estimate does not cover:
- **One subtree per shortfall.** Two subtrees at nearly the same shortfall
  can differ 2-6x: s136 has 154 nodes and s145 has 78; s155 has 40 and s168
  has 229. The fit's scatter (x1.7 in CPU) is in the "fit" row, not in the
  bootstrap.
- **Cold roots.** Every anchor starts cold. On wil50's a220, warm and cold
  gave the same tree. Here, though, a cold root can fall outside the zone
  where the warm one is inside it (s108: 108.4 warm, 129.4 cold), and then it
  branches instead of running a phase. Which way that moves the total is not
  measured: s94 (phase, 511 s) and s108 (branch, 494 s) cost about the same.
- **The incumbent is the optimum,** as in every estimate here.

## Stages A and B under C, and the estimate under C (`sampling_C/`)

The same runs as under B, with `C.sh` and `C.par`: stage A in 4 h 9 min, stage
B in 53 min, both on 4 slots (`sampling_C/`, the same files as `sampling/`).

| | B | C |
|---|---|---|
| depth-2 nodes bounded, fathomed | 438, 164 | 438, 165 |
| cut phases at depth 2 | 138 | 218 |
| CPU s at depths 0-2 | 40,762 | 58,946 |
| frontier (depth 3) | 10,726 | 10,704 |
| stage B: fathomed in the first solve, of 250 | 123 | 140 |
| stage B misses: median, 95th percentile, largest (short) | 62.9, 136.4, 188.8 | 67.7, 140.1, 188.8 |

- **At depth 2, C's extra phases do not pay for themselves.** 80 more phases
  cost 18,184 CPU s, fathom one more node and remove 22 frontier nodes (about
  4,000 s of subtree). That is about 5 core-hours, under 1% of the total.
- **At depth 3 they do.** The pools those phases pass down let 17 of the
  250 sampled nodes fathom in their first solve where under B they missed
  (249 of the 250 nodes are the same: the sample is fixed by the seed). The
  other misses keep B's shortfall (median change 0; the largest, 67.7, a
  node that moved closer). The anchors' shortfalls are the same in both.

```
python3 experiments/estimates/sko42/tune/estimate.py --table sampling_C/stageB_summary.md --top 451 \
  --top-cpu 58946 --front 10704 subtrees/Bs20=20.2 subtrees/Bs63=62.9 subtrees/Bs94=93.5 \
  subtrees/Bs108=108.4 subtrees/Cs124=124.0 subtrees/Cs136=136.4 subtrees/Cs145=145.4 \
  subtrees/Bs155=155.2 subtrees/Cs168=168.2 subtrees/Cs189=188.8
```

C's five subtrees, B's s20, s63 and s94 (the same trees under C), and B's s108
and s155 in place of C's (upper bounds):

| quantity | model | mean per frontier node | total (95%) |
|---|---|---|---|
| DNN nodes | interpolation | 7.4 | 79,854 (50,141-115,535) |
| DNN nodes | fit (x48.7 per 100 short, scatter x2.9) | 10.9 | 116,853 (67,647-178,818) |
| CPU s | interpolation | 164.6 | 1.82e6 (1.33e6-2.40e6) = 506 core-hours |
| CPU s | fit (x11.5 per 100 short, scatter x1.6) | 179.2 | 1.98e6 (1.46e6-2.55e6) = 549 core-hours |

**sko42 under C: about 8e4-1.2e5 DNN nodes and 506-549 core-hours on the
cloud box (370-710 at the ends of the intervals), 16-18% under B's 603-668.**
On the laptop that is about 255-275 core-hours, some 1.8-1.9 days on the 6
performance cores. C is the setting to solve sko42 with. What the estimate
under B does not cover (above) holds here too.
