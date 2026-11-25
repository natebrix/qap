# wil50-est-261124

**PART-WAY**: the estimate's dives run wil50-est-261124-dives is running, 997 of 998 dives units merged

- PART-WAY: 1 of 998 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
- SLOW PENDING: 1 pending units have run longer than 90% of the merged ones took (3710 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- SETUP-BOUND: bounding the sampled open nodes took 196.0 CPU h, the dives 141.6 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
- WIDE: the jackknife over units gives a standard error of 27% of the rest's CPU

## The run

| | |
|---|---|
| problem | wil50t.dat -s SYM/wil50t.sym |
| params | wil50_d2r10_g150xl2.par |
| incumbent | 48818.0 |
| slots | 263; units of 1200.0 CPU s, traced |
| created | 2026-11-24 16:47 UTC |
| status | paused |
| estimate | dives run wil50-est-261124-dives: k1.5, sample 1.0 (unif), QAP_EST_MIN_DIVES=5, QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 |

## The work

| | |
|---|---|
| wall | 2 h 00 min |
| nodes | 17,309 (17,308 DNN, 15,630 DNN-fathomed) |
| CPU | 298.0 h |
| units | merged 3313, queued 998 |
| status written | 2026-11-24 18:47 UTC |

## The estimate

```
PART-WAY: 1 of 998 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
SLOW PENDING: 1 pending units have run longer than 90% of the merged ones took (3710 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
SETUP-BOUND: bounding the sampled open nodes took 196.0 CPU h, the dives 141.6 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
WIDE: the jackknife over units gives a standard error of 27% of the rest's CPU
done so far: 17,309 nodes, 298.0 CPU h
estimated rest: 984,894 nodes (5-95%: 756,860 to 1,260,763), 8,620.0 CPU h (6,054.2 CPU h to 11,561.0 CPU h)
estimated tree: 1,002,203 nodes (5-95%: 774,169 to 1,278,072), 8,918.0 CPU h (6,352.2 CPU h to 11,859.0 CPU h)
4515 open nodes (950 with children that branch), 10224 dives in 950 of them, deepest dive node at depth 7; estimator 337.6 CPU h
4,515 of the 18,576 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 142,251 nodes, 2,342.6 CPU h
the largest dive makes 16% of the rest's CPU, the 10 largest 41%; 0 of 950 dived open nodes got one dive; estimator CPU: setup 196.0 CPU h, dives 141.6 CPU h
dives going deeper, by depth (reached, went on): 3: 776/99, 4: 6526/426, 5: 3356/61, 6: 152/2, 7: 2/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 15 | 0 | 5.93 | 0.00 |
| 2 | 476 | 159 | 51.46 | 19.65 |
| 3 | 6,486 | 7,669 | 160.21 | 395.47 |
| 4 | 9,590 | 106,958 | 76.60 | 1,707.10 |
| 5 | 741 | 439,732 | 3.79 | 4,901.75 |
| 6 | 0 | 363,390 | 0.00 | 1,359.96 |
| 7 | 0 | 63,164 | 0.00 | 209.22 |
