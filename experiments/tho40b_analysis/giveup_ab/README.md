# An earlier DNN give-up: A/B on a complete tho40 subtree (Nov 16)

**Result (laptop, night of Nov 15):** giving up at 100 iterations without the zone saves
10.1% of the CPU on the tho40ab subtree, for 10.4% more nodes.

| arm | nodes | fathomed | wall s | CPU s (M2 Pro) | nodes / g200 | CPU / g200 |
|---|---|---|---|---|---|---|
| g200 (tho40-b) | 22,895 | 22,244 | 4,173 | 31,594 | 1.000 | 1.000 |
| g100 | 24,336 | 23,643 | 3,686 | 29,358 | 1.063 | 0.929 |
| g100z0 | 25,273 | 24,553 | 3,529 | 28,407 | 1.104 | 0.899 |

- **The baseline checks the setup:** g200 is within 1.1% of the A/B's 22,650 nodes
  (Accelerate rounds differently from OpenBLAS), at 1/2.49 of the c7a's CPU.
- **Less than the ceiling:** the earlier give-up saves about half of the 19% its
  iterations past the 100th would suggest. The rest goes to the nodes it branches
  that a later give-up would have fathomed.
- **Next:** whether 75 or 50 does better (the trend had not flattened at 100 without
  the zone), and whether the saving holds on a second subtree: `26:3,38:0`, a complete
  tho40-b depth-2 subtree of 7,100 nodes and 28,123 c7a s, whose root starts 5,827 short
  of the cutoff against the tho40ab subtree's 8,634.


## The cut-phase arms (Nov 16)

From `experiments/bound_lab` (README.md, section 2): its Knuth dives put tho40's tree at 226
c7a CPU-days with the grid cut phase at `cut.zone = 0.06`, against 1,308 as tho40-b runs, and
about half that again with a hand-over rule and the phase's own penalty. Three arms put that on
this subtree, on top of tho40-b's recipe (`ship.bits = 32` as in production, where the dives ran
single-process at 64):

| arm | settings | build |
|---|---|---|
| cut06 | the dives' zone 0.06 arm exactly (`tho40_dives/par/cut06.par`): `pool.on = 1`, `cut.grid = 5:8`, `cut.side = 1`, `cut.cp = 100:1000`, `cut.warmy = 1`, `cut.rounds = 12`, `cut.iters = 200`, `cut.freeze = 1000`, `cut.zone = 0.06` | cognac |
| np10 | the README's ten-line block: cut06 with `dnn.giveup_warm = 100`, and `dnn.gu_zone = 0`, `dnn.beta_sched = w=0.36,wc=2.88` | cognac |
| full12 | np10 with `cut.handover = 20:0:0.22:1.5:100:3:1.25` and `cut.cp_wc_gap = 300`: the dives' combined arm exactly (`comb.par`) | integration b86352e, built in its own worktree (`../qap-b86352e`, or `INTEGRATION_WT`) |

**Through the coordinator** (checked on this box, the depth-3 subtree `26:10,8:0,30:32`: 571
nodes and 2,122 c7a s under tho40-b's settings; 2 slots, 60-s units so that pools cross unit
boundaries, 38 units):
- **The same tree as one process:** full12 gave the same 38 nodes and the same outcome at every
  node (1 branched, 37 fathomed) as a single `bnbdrv` run on the same build. The branched node's
  bounds matched to the last digit. The children's DNN bounds differ by 0.06-0.6 units, from the
  fp32 parent states the coordinator ships.
- **CPU here:** np10 344 s, full12 267 s (tho40-b's settings: about 1,540 s on this box).

`BUDGET` sets the units' CPU budget (default 300 s).

In tho40-b's part-way export, the nodes whose DNN bound gave up as hopeless are 3% of
the nodes but 30% of the DNN CPU, and 19% of all of it is their iterations past the
100th (`../giveup_probe/README.md`). Past depth 4, at most 1.6% of the nodes that fathom need
more than 100 iterations. The probe showed that giving up later does not pay. This tests
giving up sooner, on a whole subtree, so the nodes it branches instead of fathoming are
counted too.

`run.sh` searches the tho40ab subtree (`-F 26:10,8:0`, `-u 240516`: 22,650 nodes and 78,804
c7a CPU s in the shipping A/B, `infra/phase1/runs/tho40ab-32`) three times through the
coordinator's local mode. Each arm uses tho40-b's settings (`infra/phase1/recipes/tho40.par`;
for this subtree, whose deepest node is at depth 8, they bound as the A/B's did) except
the give-up:

| arm | `dnn.giveup` | `dnn.gu_zone` | a hopeless node gives up at |
|---|---|---|---|
| g200 | 200 | 0.003 | 200; 300 within 0.3% of the cutoff (tho40-b) |
| g100 | 100 | 0.003 | 100; 200 within the zone |
| g100z0 | 100 | 0 | 100 |

The incumbent is fixed and parent states are shipped, so the trees are deterministic on
one machine: the arms' node counts differ only by the give-up. g200 checks the setup
against 22,650 (close, not exact: Accelerate rounds differently from OpenBLAS).

```
cd ~/dev/qap && git pull origin cognac
caffeinate -i experiments/tho40b_analysis/giveup_ab/run.sh
```

Any arm `gN` (`dnn.giveup = N` with the zone) or `gNz0` (without it) runs by name,
e.g. `ARMS="g200 g100z0 g75z0 g50z0"`; arms already in `OUT` are skipped.

About an hour an arm on the M2 Pro's 12 cores (78,804 c7a CPU s at the laptop's 2.23x
under load). Keep the laptop plugged in with the lid open. The table is in
`runs/summary.txt`. Each arm's results, log and traces (`runs/ARM.export/`: node and
DNN traces, about 5 MB) stay in `runs/` for a commit. An arm already done is skipped,
so after an interruption the same command goes on. Tested Nov 16 on the cloud box on a
37-node subtree (`SUB=26:10,31:0,34:32,2:7 SLOTS=4`, 6.5 min for the three arms).
