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

## Since (Nov 19)

- **Wider still: probably not worth it.** The bound lab found sko49's zone, 0.018 (421 short), already wide
  enough that every node below depth 2 which its first solve does not fathom runs a phase
  (`../../../bound_lab/wil50_policy/README.md`, "What did not carry over"). sko42's 0.025 is 395 short: about
  the same. tho40's 0.06 is a different instance's need (its far phases still close 20 units per 100
  iterations).
- **The sko49 policy block, in `sko42.par` and `sko49.par`.** The bound lab's arm `hoe2w`: early entry
  (`dnn.giveup_warm = 100`), the hand-over rule fitted on sko49 (`cut.handover = 27:0:5`) and the phase's own
  penalty `wc=1.4` (0.7 before), 1.21x and 1.34x fewer iterations on two sko49 subtrees. The next sko42 run
  (`infra/phase1/specs/sko42.spec`) runs it; paired with z25, which differs from it in the block only, it
  measures the block on a whole tree, in DNN iterations (z25's DNN trace is in its export, not in git).
- **A check on sko42 before the run [m].** Two of the subtrees both arms closed, solved cold from `-F` on the
  cloud box under `sko42_z25.par` and under `sko42.par` with the block, side by side on one 4-core machine, incumbent
  15,814. The same trees node for node, both to the end, no new objective:

  | subtree | nodes (B / D / C) | iterations z25 / block | CPU s z25 / block |
  |---|---|---|---|
  | `36:9,15:0,22:6` | 116 (3 / 110 / 3) | 15,135 / 9,880 (**1.53x**) | 1,455 / 975 (1.49x) |
  | `36:16,21:0,22:2` | 40 (1 / 34 / 5) | 6,263 / 5,446 (**1.15x**) | 769 / 595 (1.29x) |

  The bound lab's 1.21x and 1.34x on sko49 bracket these.

## Failed pooled solves (Nov 20)

A child's first solve starts from the cut pool its parent handed down. When the pool no longer fits the child the
solve can fail outright (`latest/bnb.C`, "tai_solves (pool_retry.patch)"): its certified bound ends far below the
parent's, the node looks far outside `cut.zone`, gets no phase and branches, and its children inherit the failed
state and fail too, down to depths where the QPB alone fathoms. `dnn.retry = 1` solves such a node again without the
pool, then cold; tho40c has it since Nov 16, `sko42.par` and `sko49.par` did not until Nov 20.

Unfathomed DNN solves that certified more than 1% of the incumbent (158) below their parent's bound, the outermost
of each such chain:

| run | failed roots | their subtrees | CPU h |
|---|---|---|---|
| z10 | 4 (in 36:15 and 36:17) | 12,113 nodes | 14 |
| z25 | 3: `36:3,22:21,9:6` (-8,073 against 15,758), `36:16,7:0,33:41,1:6` (-2.4 million), `36:16,40:0,3:41,12:11` (3,950) | 117,254 nodes, the first still open | 189 |
| sko42-261120 (the block, 08:46 UTC) | `36:7,28:0` at depth 2 (-36,269 against 15,640; 41 nodes and 0.19 CPU h in z25) | 235,764 nodes, open, to depth 16 | 387 of 955 |

So z25's two deep tails, which carried nearly all of its imputed rest, were failures, not sko42's tree: the
"hidden tail" read above for z10's 36:3 is wrong (z10 did not fail there). Without the failed subtrees (`rest.py`'s
imputation on the rest of each tree, corrected as above): **z25 about 680-690 CPU h, z10 about 960-1,000** (0.70,
close to the paired 0.65). With the policy block's 1.2-1.35x, sko42 under `sko42.par` with `dnn.retry = 1` should
take about 500-700 CPU h. The table under "The tree's size" above stands as written for the traces as they are.
