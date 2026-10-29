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
