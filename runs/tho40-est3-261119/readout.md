# tho40-est3-261119

A second estimate of tho40-est2-261119: its paused frontier dived again (its own readout has its first).

**PART-WAY**: the estimate's dives run tho40-est3-261119 is running, 174 of 548 dives units merged

- PART-WAY: 374 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (96% of the rest's CPU)
- SLOW PENDING: 49 pending units have run longer than 90% of the merged ones took (1999 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- GAP: 22 open nodes of pending units are at depths no merged unit has open nodes at yet (depth 7: 22): not counted, so the estimate is low by their share until one is merged
- WIDE: the jackknife over units gives a standard error of 121% of the rest's CPU

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
PART-WAY: 374 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (96% of the rest's CPU)
SLOW PENDING: 49 pending units have run longer than 90% of the merged ones took (1999 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
GAP: 22 open nodes of pending units are at depths no merged unit has open nodes at yet (depth 7: 22): not counted, so the estimate is low by their share until one is merged
WIDE: the jackknife over units gives a standard error of 121% of the rest's CPU
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 780,885 nodes (5-95%: 751,054 to 817,820), 11,685.8 CPU h (11,407.5 CPU h to 12,113.6 CPU h)
estimated tree: 799,652 nodes (5-95%: 769,821 to 836,587), 11,802.9 CPU h (11,524.7 CPU h to 12,230.7 CPU h)
1289 open nodes (246 with children that branch), 9472 dives in 246 of them, deepest dive node at depth 6; estimator 39.6 CPU h
1,289 of the 5,611 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 495,371 nodes, 14,167.5 CPU h
the largest dive makes 2% of the rest's CPU, the 10 largest 3%; 0 of 246 dived open nodes got one dive; estimator CPU: setup 16.2 CPU h, dives 23.4 CPU h
dives going deeper, by depth (reached, went on): 3: 13/7, 4: 4413/109, 5: 5162/43, 6: 43/0


(nodes at depth 8 are in the solve's own part; the dives reached depth 6)
```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 2 | 45.12 | 0.35 |
| 3 | 2,924 | 1,231 | 28.66 | 21.46 |
| 4 | 8,475 | 14,872 | 30.50 | 84.49 |
| 5 | 5,851 | 40,215 | 10.21 | 317.26 |
| 6 | 1,004 | 30,109 | 1.27 | 31.78 |
| 7 | 148 | 0 | 0.16 | 0.00 |
| 8 | 33 | 0 | 0.03 | 0.00 |
