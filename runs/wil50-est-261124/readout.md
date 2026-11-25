# wil50-est-261124

**COMPLETE**: the estimate (wil50-est-261124-dives finished)

- SETUP-BOUND: bounding the sampled open nodes took 197.0 CPU h, the dives 142.7 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

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
the estimate's cost: 339.7 CPU h of estimator CPU, 3% of the estimated rest of the tree (9,847.1 CPU h)
SETUP-BOUND: bounding the sampled open nodes took 197.0 CPU h, the dives 142.7 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 17,309 nodes, 298.0 CPU h
estimated rest: 1,005,804 nodes (5-95%: 757,015 to 1,291,987), 9,847.1 CPU h (6,714.6 CPU h to 14,144.7 CPU h)
estimated tree: 1,023,113 nodes (5-95%: 774,324 to 1,309,296), 10,145.1 CPU h (7,012.6 CPU h to 14,442.8 CPU h)
4525 open nodes (958 with children that branch), 10264 dives in 958 of them, deepest dive node at depth 7; estimator 339.7 CPU h
4,525 of the 18,597 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 167,302 nodes, 2,100.8 CPU h
the largest dive makes 14% of the rest's CPU, the 10 largest 47%; 0 of 958 dived open nodes got one dive; estimator CPU: setup 197.0 CPU h, dives 142.7 CPU h
dives going deeper, by depth (reached, went on): 3: 776/99, 4: 6566/433, 5: 3363/62, 6: 153/2, 7: 2/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 15 | 0 | 5.93 | 0.00 |
| 2 | 476 | 159 | 51.46 | 19.65 |
| 3 | 6,486 | 7,686 | 160.21 | 397.11 |
| 4 | 9,590 | 107,568 | 76.60 | 1,720.56 |
| 5 | 741 | 444,299 | 3.79 | 4,955.67 |
| 6 | 0 | 382,929 | 0.00 | 2,544.89 |
| 7 | 0 | 63,164 | 0.00 | 209.22 |
