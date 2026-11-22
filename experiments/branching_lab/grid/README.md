# Does the top-of-tree line matter on a grid? tho30 with the DNN lookahead at depth 1

Nov 21 2026. [m] = measured here, [e] = estimate or inference.

**Yes. The lookahead at tho30's ten depth-1 nodes (150 iterations a child, four lines a side) cuts the whole tree
from 14,815 to 11,100 nodes (25%). With 30 iterations a child most of that is lost (14,221, 4%).**

## The runs (`run.sh`)

Whole trees through the coordinator's local mode (`infra/laptop/solve.sh` of the experiment build: cognac
0e6c2ecf + `experiments/line_signals/experiment_build.patch` + `../lsla/lsla.patch`), the nug30 record recipe
(`infra/phase1/recipes/nug30.par`: rule 9, no cut phase, no DNN at the root) with tho30's symmetry file, the
incumbent the best known + 2, 3 slots, the root not from a cache. `QAP_LA=K:CAP:DMIN:DMAX`: at a node of depth
DMIN..DMAX the K rows and K columns with the largest marginal scores are each scored by bounding their children
with CAP warm ADMM iterations (no certificate), and the node branches on the line with the smallest sum of
shortfalls. The recipe has no DNN at the root, so depth 1 is the first depth where it can fire.

| arm | solver's nodes | trace lines | CPU s | lookahead CPU s | CPU net of the lookahead | branched at depths 2 / 3 / 4 |
|---|---|---|---|---|---|---|
| `prod` (rule 9) | 14,815 | 14,900 | 9,375 | 0 | 9,375 | 108 / 279 / 165 |
| `la150_d1` (`QAP_LA=4:150:1:1`) | **11,100** | 11,165 | 13,115 | 5,535 | 7,580 | 89 / 180 / 125 |
| `la30_d1` (`QAP_LA=4:30:1:1`) | 14,221 | 14,285 | 10,836 | 1,514 | 9,322 | 92 / 238 / 183 |

All three closed at 149,936 [m]. Production's count here (14,815) is the reference run's (14,822, cognac
`experiments/tho30_solve`) to 0.05%.

## What it says

- **The line at ten nodes moves a quarter of the tree** [m]. The 150-iteration lookahead chose a line other
  than production's at all ten depth-1 nodes (`la.log` against the depth-1 lines of `prod`'s trace), and often
  not one of rule 9's two candidates: c9 at four nodes, r17 and r18 at three (the fourth row by marginal score).
  So on a grid, too, the high-stakes tier exists and the shortlist has to be wider than rule 9's.
- **30 iterations are not enough** [m]: the cheap lookahead agrees with the 150-iteration one at 3 of 10 nodes
  and keeps 4 points of the 25. The same thing was measured at tai60b's chain nodes (`../lsla/README.md`: CAP
  30 mis-orders the lines, CAP 150 orders them).
- **On a tree this small the lookahead loses on CPU** [m]: 5,535 s for ten firings against a 9,375 s tree. Its
  cost per node is fixed by the node's size (about 550 s at m = 29); what it saves scales with the subtree. On
  an instance whose depth-1 subtrees cost thousands of CPU hours the same trade would be far in its favour, if
  the 25% carries over [e: not measured on a larger grid].
- CPU: this container ran 3 slots on 4 cores, partly beside another run, so the node counts are the clean
  comparison [e].
- The `la150_d1` arm was stopped once by the two-hour limit on background commands, 65 minutes in, and resumed
  from its store (`cmd.first.txt`, `cmd.txt`). Its ten firings are ten distinct nodes [m `la.log`].

## Later (Nov 21, `../root/README.md`)

With the lookahead restricted to a row's real children at nodes with symmetry left (`../root/root.patch`), the
same arm gives 10,083 nodes (32% fewer than production), and nug30 with the same lookahead 4,037 against 4,285
to 4,365.

## Not run

- The lookahead at depths 1 and 2 (about a hundred more firings: roughly 50,000 CPU s here [e]).
- The rule "tail" (`../tail`) on a grid: it may choose a column, so it needs the rows-only case for instances
  with symmetry left.
- A second grid (nug30, sko42).
