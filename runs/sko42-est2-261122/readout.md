# sko42-est2-261122

A second estimate of sko42-est-261121b: its paused frontier dived again (its own readout has its first).

**PART-WAY**: the estimate's dives run sko42-est2-261122 is running, 685 of 689 dives units merged

- PART-WAY: 4 of 689 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
- SLOW PENDING: 4 pending units have run longer than 90% of the merged ones took (1182 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- HEAVY TAIL: the 10 largest dives make 76% of the rest's nodes and 94% of its CPU (the largest 55%, 87%): more dives where they come from
- WIDE: the jackknife over units gives a standard error of 87% of the rest's CPU

## The run

| | |
|---|---|
| problem | sko42t.dat -s SYM/sko42t.sym |
| params | sko42.par |
| incumbent | 15814.0 |
| slots | 263; units of 1200.0 CPU s, traced |
| created | 2026-11-21 04:50 UTC |
| status | paused |
| estimate | dives run sko42-est2-261122: k1.5, sample 1.0 (unif), QAP_EST_DNN_TRACE=1, QAP_EST_MIN_DIVES=5, QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 |

## The work

| | |
|---|---|
| wall | 0 h 27 min |
| nodes | 7,792 (7,791 DNN, 7,030 DNN-fathomed) |
| CPU | 46.5 h |
| units | merged 2340, queued 689 |
| status written | 2026-11-21 05:17 UTC |

## The estimate

```
PART-WAY: 4 of 689 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
SLOW PENDING: 4 pending units have run longer than 90% of the merged ones took (1182 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
HEAVY TAIL: the 10 largest dives make 76% of the rest's nodes and 94% of its CPU (the largest 55%, 87%): more dives where they come from
WIDE: the jackknife over units gives a standard error of 87% of the rest's CPU
done so far: 7,793 nodes, 46.7 CPU h
estimated rest: 244,447 nodes (5-95%: 90,081 to 509,676), 3,237.2 CPU h (262.0 CPU h to 8,846.5 CPU h)
estimated tree: 252,240 nodes (5-95%: 97,874 to 517,469), 3,283.9 CPU h (308.7 CPU h to 8,893.2 CPU h)
2050 open nodes (326 with children that branch), 12898 dives in 326 of them, deepest dive node at depth 6; estimator 65.3 CPU h
2,050 of the 6,428 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 139,024 nodes, 2,808.2 CPU h
the largest dive makes 87% of the rest's CPU, the 10 largest 94%; 0 of 326 dived open nodes got one dive; estimator CPU: setup 26.3 CPU h, dives 39.0 CPU h
dives going deeper, by depth (reached, went on): 3: 4304/168, 4: 7691/75, 5: 1146/9, 6: 9/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 2.05 | 0.00 |
| 2 | 269 | 169 | 9.54 | 7.67 |
| 3 | 3,859 | 6,550 | 25.62 | 67.36 |
| 4 | 3,112 | 36,046 | 8.33 | 139.86 |
| 5 | 493 | 67,018 | 1.03 | 293.05 |
| 6 | 47 | 132,366 | 0.08 | 2,722.37 |
