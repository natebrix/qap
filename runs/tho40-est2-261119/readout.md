# tho40-est2-261119

**PART-WAY**: the estimate's dives run tho40-est2-261119-dives is running, 437 of 548 dives units merged

- PART-WAY: 111 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (39% of the rest's CPU)

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
PART-WAY: 111 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (39% of the rest's CPU)
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 1,173,461 nodes (5-95%: 952,957 to 1,480,864), 3,892.6 CPU h (3,477.2 CPU h to 4,351.4 CPU h)
estimated tree: 1,192,228 nodes (5-95%: 971,724 to 1,499,631), 4,009.8 CPU h (3,594.4 CPU h to 4,468.6 CPU h)
4052 open nodes (1164 with children that branch), 16307 dives in 1164 of them, deepest dive node at depth 7; estimator 179.6 CPU h
4,052 of the 14,148 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 306,190 nodes, 577.0 CPU h
the largest dive makes 5% of the rest's CPU, the 10 largest 16%; 0 of 1164 dived open nodes got one dive; estimator CPU: setup 74.6 CPU h, dives 105.1 CPU h
dives going deeper, by depth (reached, went on): 3: 240/96, 4: 8518/519, 5: 8070/116, 6: 210/2, 7: 2/0


(nodes at depth 8 are in the solve's own part; the dives reached depth 7)
```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 53 | 45.12 | 4.68 |
| 3 | 2,924 | 6,048 | 28.66 | 150.44 |
| 4 | 8,475 | 85,927 | 30.50 | 885.02 |
| 5 | 5,851 | 276,758 | 10.21 | 941.64 |
| 6 | 1,004 | 344,624 | 1.27 | 405.39 |
| 7 | 148 | 11,887 | 0.16 | 6.53 |
| 8 | 33 | 0 | 0.03 | 0.00 |
