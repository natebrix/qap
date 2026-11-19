# sko42: cut.zone 0.025 against 0.010 on two whole-tree runs (Nov 19)

The two arms of `infra/phase1/specs/sko42-zone.spec`, `sko42-zone-261119-z10` (`sko42.par`, then `cut.zone =
0.010`) and `sko42-zone-261119-z25` (`sko42_z25.par`: that line at 0.025, nothing else changed). Same commit,
incumbent 15,814, unit sizing and slots (a c7a.2xlarge coordinator with 7, one 64-slot worker from 08:48 UTC).
Both traced, both stopped part-way at 20:48 UTC (`qapctl stop`) after 18.3 h, to give the Spot quota to the
tho40 runs. Their records and node traces are in `infra/phase1/runs/`.

| | z10 | z25 |
|---|---|---|
| CPU | 872.4 h | 842.0 h |
| nodes | 99,377 | 224,765 |
| depth-1 subtrees closed | 2 of 12 | 5 of 12 |
| depth-2 branched nodes closed | 224 of 270 | 258 of 272 |
| depth-3 | 970 of 1,083 | 1,041 of 1,062 |
| depth-4 | 638 of 784 | 931 of 963 |
| open nodes | 9,585 | 16,652 (mostly at depths 8-11, QPB-fathomed) |

**From Nov 19, `sko42.par` has `cut.zone = 0.025`.** `../sko42_z10.par` is the file z10 ran, byte for byte (the
solver's settings hash covers the file's bytes: f0f67bfef2d59e5b, which it checks on a resume); `../sko42_z25.par`
is z25's. The spec's z10 arm points at `sko42_z10.par` now.

## Paired subtrees: the same trees for 0.65 of the CPU (`closed_pairs.py`)

The subtrees both arms closed, paired by their root's fixes: 3,885 maximal pairs (closed in both, their
parent not), 71% of z10's CPU and 48% of z25's.

| root depth | subtrees | z10 nodes | z25 nodes | z10 CPU h | z25 CPU h | z25/z10 CPU (5-95%) |
|---|---|---|---|---|---|---|
| 1 | 2 | 3,245 | 3,436 | 27.0 | 20.4 | 0.76 (0.62-0.78) |
| 2 | 299 | 32,783 | 32,968 | 355.0 | 216.6 | 0.61 (0.58-0.64) |
| 3 | 1,138 | 18,540 | 18,826 | 176.9 | 120.6 | 0.68 (0.65-0.71) |
| 4 | 1,461 | 6,739 | 7,365 | 53.1 | 37.5 | 0.71 (0.66-0.76) |
| all (depths 1-8) | 3,885 | 62,913 | 64,128 | 618.7 | 400.9 | 0.65 (0.63-0.67) |

The trees are the same size (1.02x the nodes); on the 452 branched roots whose subtrees have exactly the same
node count in both, the ratio is 0.65 too.

**Why (`outcomes.py`).** The wider zone runs the cut phase at more parents (those whose uncut bound ends
within 395 of the cutoff, against 158), and their children inherit the cuts with the parent's state, so the
DNN bound alone fathoms them. Of the 39,032 nodes both arms bounded, 2,628 were fathomed by their own cut
phase in z10 and by the DNN bound alone in z25 (median 122 against 18 CPU s); z25 ran the cut phase at
fewer nodes overall (2.7% against 9.1%). On nodes decided alike: DNN-fathomed in both 0.83 of the CPU, cut-phase
fathomed in both 0.55, branched in both 0.86.

**The caveat: the machines.** z10's worker was an r7a.16xlarge throughout; z25's were two c7a.16xlarge, then
two m7a.16xlarge after Spot interruptions. All three are the same Zen 4 cores, but tai50b showed 13-20% from
run to run on the same nodes (`experiments/estimates_nov19`). That can be some of the 0.83 on DNN-fathomed
nodes; not the outcome shift or the 0.55 on cut-phase nodes. The depth-1 units ran first, on the two arms'
identical c7a.2xlarge coordinators: at 36:0 (282 nodes in both) z25 took 2.19 CPU h against 3.56.

## The tree's size (`rest.py`, `backtest.py`)

`rest.py` imputes each open node (a child made, not traced) the mean of the closed subtrees rooted at its
depth, or of its closed siblings'. Closed subtrees are the small ones, so both run low. `backtest.py` cuts
z25's trace at 30-85% of its lines and compares with the depth-1 subtrees finished by the end: the siblings'
imputation gave 0.55-0.79 of the true rest (z10: 0.38-1.03), the depth means 0.20-0.35. That backtest is
itself favourable (the subtrees that finished are the ones without a long tail).

| | z10 (0.010) | z25 (0.025) |
|---|---|---|
| done | 872 CPU h | 843 CPU h |
| rest, siblings' imputation | 128 CPU h | 467 CPU h, 313,000 nodes |
| rest, corrected by the backtest (/0.79 to /0.55) | - | 590-850 CPU h |
| tree | - | about 1,400-1,700 CPU h, 0.6-0.8 million nodes; at least 1,030 |

z25's rest is nearly all in two depth-1 subtrees whose deep parts it has reached and z10 has not: 36:3 (134
CPU h done, 9,268 open nodes, to depth 14) and 36:16 (200 h, 6,720 open, to depth 16). z10's imputation gives
36:3 one hour more after 32: its 186 open nodes hide the tail z25 is in. So z10's own figure (938-1,000 CPU h
for the tree) is not its tree; at the paired 0.65, 0.010's tree would be about 2,200-2,600 CPU h, four to five
times the 506-549 core-hours estimated for it on Nov 12 (`../../README.md`). Any depth-1 subtree here could
hide such a tail; only a finished run, or Knuth's dives from a paused frontier (`coord/dives.py`), would say.

On 71 slots (the coordinator and one worker) 1,400-1,700 CPU h is about 20-24 h; on 135, 10.5-12.5 h. Fujii
et al. needed 6.0 million nodes.

## Open

- **Wider still.** tho40 runs at 0.06 (`infra/phase1/recipes/tho40c.par`); on sko42 0.025 is the widest tried.
- **sko49** uses sko42's settings with `cut.zone = 0.010` (`../../sko49/`); the same question applies.
