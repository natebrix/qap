# tai50b-est-261119

**PART-WAY**: the estimate's dives run tai50b-est-261119-dives is running, 258 of 259 dives units merged

- PART-WAY: 1 of 259 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
- SLOW PENDING: 1 pending units have run longer than 90% of the merged ones took (743 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- FEW DIVES: 6 of 6 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
- SETUP-BOUND: bounding the sampled open nodes took 25.8 CPU h, the dives 6.2 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

## The run

| | |
|---|---|
| problem | tai50b.dat |
| params | tai50b.par |
| incumbent | 458821518.0 |
| slots | 135; units of 600.0 CPU s, traced |
| created | 2026-11-19 03:37 UTC |
| status | paused |
| estimate | dives run tai50b-est-261119-dives: 0.4, sample 1.0 (unif) |

## The work

| | |
|---|---|
| wall | 1 h 09 min |
| nodes | 1,116 (1,116 DNN, 1,004 DNN-fathomed) |
| CPU | 32.5 h |
| units | merged 1116, queued 259 |
| status written | 2026-11-19 04:47 UTC |

## The estimate

```
PART-WAY: 1 of 259 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
SLOW PENDING: 1 pending units have run longer than 90% of the merged ones took (743 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
FEW DIVES: 6 of 6 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
SETUP-BOUND: bounding the sampled open nodes took 25.8 CPU h, the dives 6.2 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 1,117 nodes, 32.8 CPU h
estimated rest: 1,688 nodes (5-95%: 1,688 to 1,688), 31.3 CPU h (31.3 CPU h to 31.3 CPU h)
estimated tree: 2,805 nodes (5-95%: 2,805 to 2,805), 64.1 CPU h (64.1 CPU h to 64.1 CPU h)
281 open nodes (6 with children that branch), 6 dives in 6 of them, deepest dive node at depth 8; estimator 32.1 CPU h
jackknife over units: standard error 395 nodes, 5.1 CPU h
the largest dive makes 7% of the rest's CPU, the 10 largest 24%; 6 of 6 dived open nodes got one dive; estimator CPU: setup 25.8 CPU h, dives 6.2 CPU h

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 1.20 | 0.00 |
| 1 | 50 | 0 | 7.25 | 0.00 |
| 2 | 490 | 0 | 16.60 | 0.00 |
| 3 | 576 | 0 | 7.72 | 0.00 |
| 4 | 0 | 281 | 0.00 | 16.36 |
| 5 | 0 | 737 | 0.00 | 9.65 |
| 6 | 0 | 446 | 0.00 | 3.55 |
| 7 | 0 | 132 | 0.00 | 1.24 |
| 8 | 0 | 86 | 0.00 | 0.43 |
