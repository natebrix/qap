# tho40-est3-261119

A second estimate of tho40-est2-261119: its paused frontier dived again (its own readout has its first).

**PART-WAY**: the estimate's dives run tho40-est3-261119 is running, 533 of 548 dives units merged

- PART-WAY: 15 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (5% of the rest's CPU)
- SLOW PENDING: 15 pending units have run longer than 90% of the merged ones took (4205 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
- SETUP-BOUND: bounding the sampled open nodes took 135.7 CPU h, the dives 131.1 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

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
PART-WAY: 15 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (5% of the rest's CPU)
SLOW PENDING: 15 pending units have run longer than 90% of the merged ones took (4205 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
SETUP-BOUND: bounding the sampled open nodes took 135.7 CPU h, the dives 131.1 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 2,828,534 nodes (5-95%: 1,666,014 to 4,287,436), 5,939.8 CPU h (4,571.9 CPU h to 7,716.6 CPU h)
estimated tree: 2,847,301 nodes (5-95%: 1,684,781 to 4,306,203), 6,056.9 CPU h (4,689.1 CPU h to 7,833.8 CPU h)
5872 open nodes (2299 with children that branch), 32280 dives in 2299 of them, deepest dive node at depth 8; estimator 266.8 CPU h
5,872 of the 17,916 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 724,847 nodes, 1,070.7 CPU h
the largest dive makes 7% of the rest's CPU, the 10 largest 33%; 0 of 2299 dived open nodes got one dive; estimator CPU: setup 135.7 CPU h, dives 131.1 CPU h
dives going deeper, by depth (reached, went on): 3: 418/130, 4: 20811/1352, 5: 12349/208, 6: 380/5, 7: 17/1, 8: 1/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 82 | 45.12 | 6.81 |
| 3 | 2,924 | 9,786 | 28.66 | 243.94 |
| 4 | 8,475 | 158,778 | 30.50 | 1,021.53 |
| 5 | 5,851 | 643,556 | 10.21 | 2,344.33 |
| 6 | 1,004 | 993,844 | 1.27 | 1,291.47 |
| 7 | 148 | 408,611 | 0.16 | 462.49 |
| 8 | 33 | 470,460 | 0.03 | 257.74 |
