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

## Next

Stage A of the wil50 sampling plan under B (every node at depths 0-2, the
depth-3 frontier saved): the depth-2 shortfalls, and the frontier from which
stage B samples, and the A/B subtrees are drawn.
