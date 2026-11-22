# sko42-est2-261122

A second estimate of sko42-est-261121b: its paused frontier dived again (its own readout has its first).

**PART-WAY**: the estimate's dives run sko42-est2-261122 is running, 541 of 689 dives units merged

- PART-WAY: 148 of 689 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (56% of the rest's CPU)
- SLOW PENDING: 148 pending units have run longer than 90% of the merged ones took (106 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- SETUP-BOUND: bounding the sampled open nodes took 5.8 CPU h, the dives 1.3 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
- DEEP TAIL: 6% of the dives at depth 5 went deeper, where few dives reached: the deep levels are thinly sampled

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
PART-WAY: 148 of 689 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (56% of the rest's CPU)
SLOW PENDING: 148 pending units have run longer than 90% of the merged ones took (106 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
SETUP-BOUND: bounding the sampled open nodes took 5.8 CPU h, the dives 1.3 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
DEEP TAIL: 6% of the dives at depth 5 went deeper, where few dives reached: the deep levels are thinly sampled
done so far: 7,793 nodes, 46.7 CPU h
estimated rest: 12,070 nodes (5-95%: 10,820 to 14,422), 46.0 CPU h (42.4 CPU h to 52.4 CPU h)
estimated tree: 19,863 nodes (5-95%: 18,613 to 22,215), 92.7 CPU h (89.1 CPU h to 99.0 CPU h)
999 open nodes (6 with children that branch), 429 dives in 6 of them, deepest dive node at depth 6; estimator 7.1 CPU h
999 of the 4,093 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 3,323 nodes, 7.9 CPU h
the largest dive makes 3% of the rest's CPU, the 10 largest 7%; 0 of 6 dived open nodes got one dive; estimator CPU: setup 5.8 CPU h, dives 1.3 CPU h
dives going deeper, by depth (reached, went on): 3: 277/10, 4: 86/7, 5: 83/5, 6: 5/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 2.05 | 0.00 |
| 2 | 269 | 35 | 9.54 | 1.07 |
| 3 | 3,859 | 794 | 25.62 | 6.07 |
| 4 | 3,112 | 3,374 | 8.33 | 9.30 |
| 5 | 493 | 961 | 1.03 | 2.09 |
| 6 | 47 | 802 | 0.08 | 1.77 |
