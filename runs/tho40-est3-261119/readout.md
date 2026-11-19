# tho40-est3-261119

A second estimate of tho40-est2-261119: its paused frontier dived again (its own readout has its first).

**PART-WAY**: the estimate's dives run tho40-est3-261119 is running, 252 of 548 dives units merged

- PART-WAY: 296 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (74% of the rest's CPU)
- SLOW PENDING: 36 pending units have run longer than 90% of the merged ones took (2234 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- GAP: 22 open nodes of pending units are at depths no merged unit has open nodes at yet (depth 7: 22): not counted, so the estimate is low by their share until one is merged
- WIDE: the jackknife over units gives a standard error of 47% of the rest's CPU

## The run

| | |
|---|---|
| problem | tho40t.dat -s SYM/tho40t.sym |
| params | tho40c.par |
| incumbent | 240518.0 |
| slots | 135; units of 1200.0 CPU s, traced |
| created | 2026-11-19 02:26 UTC |
| status | paused |
| estimate | dives run tho40-est3-261119: k1.5, sample 1.0 (unif), QAP_EST_MIN_DIVES=5, QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 |

## The work

| | |
|---|---|
| wall | 1 h 03 min |
| nodes | 18,764 (18,763 DNN, 16,788 DNN-fathomed) |
| CPU | 116.5 h |
| units | merged 596, queued 548 |
| status written | 2026-11-19 03:30 UTC |

## The estimate

```
PART-WAY: 296 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (74% of the rest's CPU)
SLOW PENDING: 36 pending units have run longer than 90% of the merged ones took (2234 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
GAP: 22 open nodes of pending units are at depths no merged unit has open nodes at yet (depth 7: 22): not counted, so the estimate is low by their share until one is merged
WIDE: the jackknife over units gives a standard error of 47% of the rest's CPU
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 560,221 nodes (5-95%: 518,585 to 611,979), 2,778.3 CPU h (2,458.0 CPU h to 3,214.3 CPU h)
estimated tree: 578,988 nodes (5-95%: 537,352 to 630,746), 2,895.5 CPU h (2,575.1 CPU h to 3,331.5 CPU h)
2055 open nodes (501 with children that branch), 14380 dives in 501 of them, deepest dive node at depth 7; estimator 71.1 CPU h
2,055 of the 8,158 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 86,017 nodes, 1,292.2 CPU h
the largest dive makes 9% of the rest's CPU, the 10 largest 14%; 0 of 501 dived open nodes got one dive; estimator CPU: setup 31.4 CPU h, dives 39.8 CPU h
dives going deeper, by depth (reached, went on): 3: 99/23, 4: 7124/215, 5: 7383/68, 6: 68/0, 7: 12/0


(nodes at depth 8 are in the solve's own part; the dives reached depth 7)
```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 19 | 45.12 | 1.64 |
| 3 | 2,924 | 2,715 | 28.66 | 50.23 |
| 4 | 8,475 | 32,239 | 30.50 | 174.53 |
| 5 | 5,851 | 79,039 | 10.21 | 453.78 |
| 6 | 1,004 | 51,636 | 1.27 | 50.91 |
| 7 | 148 | 1,172 | 0.16 | 1.25 |
| 8 | 33 | 0 | 0.03 | 0.00 |
