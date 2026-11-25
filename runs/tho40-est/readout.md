# tho40-est

**PART-WAY**: the estimate's dives run tho40-est-dives is running, 524 of 552 dives units merged

- STALE: the coordinator last wrote its status 8771 minutes ago (an interrupted coordinator? collect.sh shows the instances and why they stopped)
- PART-WAY: 28 of 552 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
- HEAVY TAIL: the 10 largest dives make 58% of the rest's nodes and 37% of its CPU (the largest 10%, 6%): more dives where they come from
- FEW DIVES: 828 of 1108 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
- SETUP-BOUND: bounding the sampled open nodes took 68.3 CPU h, the dives 26.3 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

## The run

| | |
|---|---|
| problem | tho40t.dat -s SYM/tho40t.sym |
| params | tho40c.par |
| incumbent | 240518.0 |
| slots | 128; units of 1200.0 CPU s, traced |
| created | 2026-11-18 20:08 UTC |
| status | paused |
| estimate | dives run tho40-est-dives: k1.5, sample 0.2 (qpb) |

## The work

| | |
|---|---|
| wall | 1 h 04 min |
| nodes | 17,579 (17,578 DNN, 15,704 DNN-fathomed) |
| CPU | 114.0 h |
| units | merged 560, queued 552 |
| status written | 2026-11-18 21:13 UTC |

## The estimate

```
PART-WAY: 28 of 552 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (0% of the rest's CPU)
HEAVY TAIL: the 10 largest dives make 58% of the rest's nodes and 37% of its CPU (the largest 10%, 6%): more dives where they come from
FEW DIVES: 828 of 1108 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
SETUP-BOUND: bounding the sampled open nodes took 68.3 CPU h, the dives 26.3 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 17,619 nodes, 114.3 CPU h
estimated rest: 1,850,462 nodes (5-95%: 1,230,914 to 2,588,069), 4,391.9 CPU h (3,389.3 CPU h to 5,440.2 CPU h)
estimated tree: 1,868,081 nodes (5-95%: 1,248,533 to 2,605,688), 4,506.2 CPU h (3,503.6 CPU h to 5,554.5 CPU h)
3479 open nodes (1108 with children that branch), 7165 dives in 1108 of them, deepest dive node at depth 7; estimator 94.6 CPU h
3,479 of the 16,889 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 344,370 nodes, 559.3 CPU h
the largest dive makes 6% of the rest's CPU, the 10 largest 37%; 828 of 1108 dived open nodes got one dive; estimator CPU: setup 68.3 CPU h, dives 26.3 CPU h
dives going deeper, by depth (reached, went on): 3: 158/12, 4: 3868/169, 5: 3180/30, 6: 169/0, 7: 1/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.03 | 0.00 |
| 2 | 297 | 101 | 44.13 | 9.76 |
| 3 | 2,625 | 10,731 | 27.29 | 287.63 |
| 4 | 8,099 | 175,520 | 30.38 | 1,417.00 |
| 5 | 5,556 | 652,489 | 10.19 | 1,762.26 |
| 6 | 927 | 1,010,511 | 1.21 | 913.37 |
| 7 | 102 | 182 | 0.09 | 0.19 |
