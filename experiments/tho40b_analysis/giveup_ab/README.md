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
