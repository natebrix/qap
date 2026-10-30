# nug30, proven optimal again (October 2026)

nug30's optimum, 6124, was proven again with this solver: QPB branch and
bound plus the certified DNN bound (latest/dnn.C) near the top of the tree.
Anstreicher, Brixius, Goux and Linderoth proved it first in 2000.

## Result: clean run (one setting throughout)

| | |
|---|---|
| nodes | 6,532,419 (5,455,913 fathomed) |
| CPU | 18,614 s (5.2 h), on 4 cores |
| DNN bound | 5,764 nodes, 5,551 fathomed, 17,192 s of the CPU |
| wall | about 1 h 20 min of solving (1 h 27 min elapsed, including a container restart) |
| root QPB | 5332.33 |

Solver at f8107dd, run with:

```sh
python3 bench/parallel.py qaplib/nug30.dat --jobs 4 --ub 6124 \
    --params PARAM/nug_deep.par --extra "-x 2:6:1:300" \
    --checkpoint nug30.ckpt --checkpoint-every 600
```

1. `clean_run1.log`: 00:32 UTC on a 2.1 GHz Xeon. A container restart killed
   it at 01:48; the last checkpoint (01:42) was at 6,531,809 nodes with 93
   open nodes. The 6 minutes after that checkpoint are lost and not counted.
2. `clean_run2.log`: rebuilt on the new machine (2.8 GHz Xeon) and resumed;
   the 93 open nodes took 9.5 minutes (2,215 CPU s).

Compared with the first solve below, the DNN bound at depth 6 from the start
cuts the tree 7x and the CPU 2.2x: almost all of the first solve's extra work
was QPB search below depth-6 nodes the depth-5 setting left to QPB. 92% of the
CPU is now in the DNN bound, and 96% of the nodes it bounds are fathomed.

## Improved settings (October 29, solver at ff1f044)

The DNN bound to depth 7, with each node's ADMM warm-started from its
parent's multiplier when the parent's DNN bound missed:

```sh
python3 bench/parallel.py qaplib/nug30.dat --jobs 4 --ub 6124 \
    --params PARAM/nug_deep.par --extra "-x 2:7:1:300:1" \
    --checkpoint nug30.ckpt --checkpoint-every 600
```

This container (2.8 GHz Xeon) is about 1.45x slower per DNN iteration
than the one that ran most of the clean run (the same code and node take
73 ms per iteration here), so the clean settings were rerun on it as an
A/B comparison with the same code:

| on the same machine | `-x 2:6:1:300` (clean settings) | `-x 2:7:1:300:1` |
|---|---|---|
| nodes | 6,532,419 | 5,817 |
| DNN nodes / fathomed | 5,764 / 5,551 | 5,807 / 5,595 |
| DNN CPU | 26,726 s | 23,537 s |
| CPU | 28,428 s | 23,618 s |
| wall | about 2.1 h | 1 h 41 min |
| logs | `ab_run1.log`, `ab_run2.log` (resumed once, 2-hour task limit) | `v3_run1.log` |

The A/B run reproduces the clean run's tree exactly (same nodes, DNN nodes
and fathoms), so the fixes since f8107dd change nothing on nug30. The new
settings take 17% less CPU. Depth 7 fathoms all 96 depth-7 nodes with the
DNN bound (178 s) and removes the QPB search below them. Warm starts cut
the DNN time per node at depths 3 and 4 by 16% and 13% (5.3 against 6.3 s,
3.4 against 3.8 s); they cut iterations per fathomed node by about 40%,
but a warm-started ADMM's iterations cost more.

A second run with the same settings (`v4_run1.log`, final per-level
statistics in `v4_final.state`) took 5,819 nodes, 24,768 CPU s and
1 h 46 min: which nodes are warm-started depends on scheduling (a node handed
to another worker loses its parent's multiplier). Its nodes per level, next to
the 2000 run's (`nug30_2000_levels.txt`, from Table 5.12 of Nathan Brixius's
thesis), are Table 6 of docs/dnn_bound/dnn_bound.pdf.

## On a laptop (October 29, solver at 2d7579c)

The improved settings on an Apple M2 Pro (6 performance + 4 efficiency
cores), built against Accelerate, 8 workers:

```sh
VECLIB_MAXIMUM_THREADS=1 python3 bench/parallel.py qaplib/nug30.dat --jobs 8 --ub 6124 \
    --params PARAM/nug_deep.par --extra "-x 2:7:1:300:1" \
    --checkpoint nug30.ckpt --checkpoint-every 300
```

| | laptop, 8 workers | 2.8 GHz Xeon, 4 workers (above) |
|---|---|---|
| nodes | 585,550 | 5,817 |
| DNN nodes / fathomed | 5,834 / 5,619 | 5,807 / 5,595 |
| DNN CPU | 8,240 s | 23,537 s |
| CPU | 8,355 s | 23,618 s |
| wall | 19 min 42 s | 1 h 41 min |
| log | `laptop_run1.log` | `v3_run1.log` |

Per DNN node the laptop is about 3x faster (depth 3: 1.9 s against 5.8;
depth 4: 1.2 against 3.5), and 8 workers gave 7.1x CPU over wall. The DNN
part of the tree is the Xeon's within 0.5%. The large node count is two of
the 96 depth-7 nodes that the Xeon runs fathomed and this run did not
(Accelerate rounds differently, and which nodes are warm-started depends on
scheduling, as above): with `-x 2:7` their subtrees went to QPB alone,
580,000 nodes for 84 CPU s.

## On an Apple M1 desktop (October 30, solver at a77c88e)

The same command on an Apple M1 (4 performance + 4 efficiency cores,
Accelerate, 8 workers), run from a fresh clone of the branch after
`make test` and `make dnncheck` passed on the machine:

| | M1, 8 workers | M2 Pro, 8 workers (above) |
|---|---|---|
| nodes | 536,075 | 585,550 |
| DNN nodes / fathomed | 5,858 / 5,643 | 5,834 / 5,619 |
| DNN CPU | 18,231 s | 8,240 s |
| CPU | 18,383 s | 8,355 s |
| wall | 42 min 54 s | 19 min 42 s |
| log | `m1_run1.log` | `laptop_run1.log` |

The DNN part of the tree again matches the Xeon runs within 0.5%; the
QPB-only nodes below a few depth-7 misses are the rest. Parallel efficiency
was 7.14x CPU over wall. The M1 is about 2.2x slower than the M2 Pro per
DNN node, in line with four performance cores against six. Its root QPB
bound is 5336.58 against 5332.33 on the Xeon: the QPB path rounds
differently on Apple silicon (see RUNNING.md, "Correctness checks").

## First solve (mixed settings)


| | |
|---|---|
| nodes | 45,486,332 (38,587,819 fathomed) |
| CPU | 40,555 s (11.3 h), on 4 cores |
| DNN bound | 5,615 nodes, 5,403 fathomed, 26,770 s of the CPU |
| wall | about 3.0 h of solving (3 h 11 min elapsed, including a container restart) |
| root QPB | 5332.33 |

## How it was run

`bench/parallel.py` with 4 workers, solver at d431c39. The known optimum
6124 was the incumbent (`--ub 6124`), as in the 2000 solve:

```sh
python3 bench/parallel.py qaplib/nug30.dat --jobs 4 --ub 6124 \
    --params PARAM/nug_deep.par --extra "-x 2:5:1:300" \
    --checkpoint nug30.ckpt --checkpoint-every 600
```

It ran in four sessions; the logs are here:

1. `run1.log`: 21:04-23:04 UTC, `-x 2:5:1:300`, 2.8 GHz Xeon. Stopped by the
   2-hour limit on background tasks; checkpoint at 14.6M nodes.
2. `run2.log`: resumed, same settings. A container restart killed it; the
   last checkpoint (23:44) was at 38.2M nodes and 35,579 CPU s. The 8
   minutes after that checkpoint are lost and not counted.
3. `run3.log`: resumed for a minute on the new machine (2.1 GHz Xeon,
   rebuilt), then stopped to change the setting.
4. `run4.log`: resumed with the DNN bound down to depth 6
   (`-x 2:6:1:300`). It finished the remaining 506 open nodes, including
   a whole depth-1 subtree, in 21 minutes (4,976 CPU s).

The totals therefore mix two machines and two settings (DNN to depth 5,
then 6). Every DNN bound is certified, so the mix does not affect the
proof. In sessions 1-2, QPB search below the depth-5 DNN misses cost about
12,000 CPU s; running the bound to depth 6 from the start should be
noticeably faster.

## Settings

- `PARAM/nug_deep.par`: default.par (lookahead at depths <= 2, rule 7 at
  depths 3-4) plus rule 6 (eccentric location) down to depth 6, and 50
  "can't fathom" FW iterations below the top levels.
- `-x D1:D2:1:300`: the DNN bound on every node QPB does not fathom at
  depths D1-D2, giving up on a node after 300 iterations when its
  certified bound extrapolates to below the cutoff. A give-up only means
  the node is branched on, and its children are bounded again.
- The DNN bound is used only to fathom; branching (and elimination with
  bound + U) uses the QPB bound.
