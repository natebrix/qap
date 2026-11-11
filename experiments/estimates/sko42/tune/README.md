# sko42: settings for an estimate (cloud box, Nov 11)

sko42 (6 x 7 grid, QAPLIB optimum 15,812, solved by Fujii et al. with 6.0e6
nodes) has one search on record, Nov 8 (`../../../sko42_solve/`). It used
the nug30 record recipe with rule 8 and was stopped after 5.5 h with its
frontier still widening. That recipe has no cut phase, although sko42 is a
grid instance like wil50 and sko49, whose regional cuts are the strongest
lever measured on either.

## The arms

| arm | files | what it is |
|---|---|---|
| A | `A.par`, `A.sh` | `../sko42.par` as it stands: the nug30 recipe with rule 9, DNN at depths 1-7, give-up 200, no cut phase |
| B | `B.par`, `B.sh` | wil50's full run (`../../../wil50_feasibility/fullrun/fullrun.par`) with `cut.grid = 6:7`: CP cut phase, `cut.zone = 0.0075` (119 short on sko42), DNN at depths 1-4 then 20 levels below a miss, give-up 300, the two-phase beta schedule |

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

## Anchors (running)

Complete subtrees (`../../subtree.sh`, cold from `-F`, the coordinator on the
cloud box) at points of stage B's give-up shortfalls:

| name | node | give-up short | where in the 127 misses |
|---|---|---|---|
| s189 | `7:35,15:0,36:15` | 188.8 | the farthest |
| s124 | `7:41,12:0,36:15` | 124.0 | 90th percentile, just outside the zone |
| s94 | `0:0,36:16,40:41` | 93.5 | 75th percentile |
| s63 | `21:0,22:23,36:9` | 62.9 | the median |
| s20 | `9:41,12:0,36:15` | 20.2 | 25th percentile |

Each runs under B; s94 also runs under A, as the settings comparison.
