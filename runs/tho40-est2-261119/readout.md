# tho40-est2-261119

**PART-WAY**: the estimate's dives run tho40-est2-261119-dives is running, 152 of 548 dives units merged

- STALE: the coordinator last wrote its status 34 minutes ago (an interrupted coordinator? collect.sh shows the instances and why they stopped)
- PART-WAY: 396 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (83% of the rest's CPU)
- GAP: 107 open nodes of pending units are at depths no merged unit has: not counted

## The run

| | |
|---|---|
| problem | tho40t.dat -s SYM/tho40t.sym |
| params | tho40c.par |
| incumbent | 240518.0 |
| slots | 135; units of 1200.0 CPU s, traced |
| created | 2026-11-19 02:26 UTC |
| status | paused |
| estimate | dives run tho40-est2-261119-dives: k1.5, sample 1.0 (unif), QAP_EST_MIN_DIVES=5, QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 |

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
PART-WAY: 396 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (83% of the rest's CPU)
GAP: 107 open nodes of pending units are at depths no merged unit has: not counted
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 394,438 nodes (5-95%: 361,709 to 442,722), 1,191.3 CPU h (1,141.1 CPU h to 1,251.3 CPU h)
estimated tree: 413,205 nodes (5-95%: 380,476 to 461,489), 1,308.5 CPU h (1,258.3 CPU h to 1,368.5 CPU h)
1103 open nodes (160 with children that branch), 4931 dives in 160 of them, deepest dive node at depth 7; estimator 29.5 CPU h
1,103 of the 4,992 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 161,332 nodes, 262.7 CPU h
the largest dive makes 2% of the rest's CPU, the 10 largest 4%; 0 of 160 dived open nodes got one dive; estimator CPU: setup 11.1 CPU h, dives 18.5 CPU h
dives going deeper, by depth (reached, went on): 4: 1772/56, 5: 3138/47, 6: 124/2, 7: 2/0


(nodes at depth 8 are in the solve's own part; the dives reached depth 7)
```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 0 | 45.12 | 0.00 |
| 3 | 2,924 | 867 | 28.66 | 11.92 |
| 4 | 8,475 | 9,199 | 30.50 | 78.30 |
| 5 | 5,851 | 20,829 | 10.21 | 68.62 |
| 6 | 1,004 | 31,570 | 1.27 | 40.24 |
| 7 | 148 | 11,853 | 0.16 | 6.28 |
| 8 | 33 | 0 | 0.03 | 0.00 |
