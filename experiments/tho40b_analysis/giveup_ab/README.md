# An earlier DNN give-up: A/B on a complete tho40 subtree (Nov 16)

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

About an hour an arm on the M2 Pro's 12 cores (78,804 c7a CPU s at the laptop's 2.23x
under load). Keep the laptop plugged in with the lid open. The table is in
`runs/summary.txt`. Each arm's results, log and traces (`runs/ARM.export/`: node and
DNN traces, about 5 MB) stay in `runs/` for a commit. An arm already done is skipped,
so after an interruption the same command goes on. Tested Nov 16 on the cloud box on a
37-node subtree (`SUB=26:10,31:0,34:32,2:7 SLOTS=4`, 6.5 min for the three arms).
