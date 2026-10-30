# tho40: DNN bound a few levels down (Oct 30 2026, cloud Xeon 2.8 GHz)

36 random nodes (seed 20261030; `nodes.txt`: depth and assignments), certified
DNN bound with `dnntest -f -T 1 -H 0 -i 1000 -u 240516` (best known 240516,
even objective, cutoff 240514.05); `out_q.txt` hold nodes q, q+4, ...

| depth | nodes | fathomed | gap to the cutoff min / median / max | median iterations | median CPU s |
|---|---|---|---|---|---|
| 0 (screen) | 1 | 0 | 14,003 | | |
| 1-2 (screen) | 16 | 0 | about 4,200-12,000 | | |
| 3 | 10 | 0 | 606 / 4,515 / 7,252 | 1000 (limit) | 298 |
| 4 | 10 | 1 | 0 / 2,272 / 5,939 | 1000 (limit) | 249 |
| 5 | 10 | 6 | -110 / -10 / 1,397 | 375 | 81 |
| 6 | 6 | 6 | all fathomed | 150 | 29 |

The gap halves per level (wil50: about 20% per level), but fathoming starts
only at depth 5-6. Assuming symmetry of order 4, a level has about
40!/(40-d)!/4 nodes: 1.5e4 at depth 3, 5.5e5 at 4, 2e7 at 5, 7e8 at 6. With
almost nothing fathomed by depth 4 and about 40% of depth-5 nodes surviving,
that is about 1e8-3e8 DNN nodes at 30-300 s: about 1e10 CPU-s here, about 100
CPU-years on an M2 Pro. Random nodes are easier than the ones a tree keeps
(those are children of unfathomed parents), so this is if anything
optimistic. Same verdict as wil50.
