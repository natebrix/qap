# tai50b-est-261119

**PART-WAY**: the estimate's dives run tai50b-est-261119-dives is running, 257 of 259 dives units merged

- PART-WAY: 2 of 259 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (1% of the rest's CPU)
- SLOW PENDING: 2 pending units have run longer than 90% of the merged ones took (742 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- FEW DIVES: 5 of 5 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
- SETUP-BOUND: bounding the sampled open nodes took 24.6 CPU h, the dives 4.1 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

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
PART-WAY: 2 of 259 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (1% of the rest's CPU)
SLOW PENDING: 2 pending units have run longer than 90% of the merged ones took (742 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
FEW DIVES: 5 of 5 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
SETUP-BOUND: bounding the sampled open nodes took 24.6 CPU h, the dives 4.1 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 1,117 nodes, 32.8 CPU h
estimated rest: 1,514 nodes (5-95%: 1,514 to 1,514), 28.7 CPU h (28.7 CPU h to 28.7 CPU h)
estimated tree: 2,631 nodes (5-95%: 2,631 to 2,631), 61.5 CPU h (61.5 CPU h to 61.5 CPU h)
280 open nodes (5 with children that branch), 5 dives in 5 of them, deepest dive node at depth 8; estimator 28.8 CPU h
jackknife over units: standard error 444 nodes, 5.5 CPU h
the largest dive makes 8% of the rest's CPU, the 10 largest 20%; 5 of 5 dived open nodes got one dive; estimator CPU: setup 24.6 CPU h, dives 4.1 CPU h

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 1.20 | 0.00 |
| 1 | 50 | 0 | 7.25 | 0.00 |
| 2 | 490 | 0 | 16.60 | 0.00 |
| 3 | 576 | 0 | 7.72 | 0.00 |
| 4 | 0 | 280 | 0.00 | 16.09 |
| 5 | 0 | 691 | 0.00 | 8.74 |
| 6 | 0 | 401 | 0.00 | 2.82 |
| 7 | 0 | 88 | 0.00 | 0.68 |
| 8 | 0 | 43 | 0.00 | 0.19 |
