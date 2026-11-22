# sko42-est2-261122

A second estimate of sko42-est-261121b: its paused frontier dived again (its own readout has its first).

**COMPLETE**: the estimate (sko42-est2-261122 finished)

- HEAVY TAIL: the 10 largest dives make 76% of the rest's nodes and 93% of its CPU (the largest 49%, 85%): more dives where they come from
- WIDE: the jackknife over units gives a standard error of 85% of the rest's CPU

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
the estimate's cost: 69.2 CPU h of estimator CPU, 2% of the estimated rest of the tree (3,310.0 CPU h)
HEAVY TAIL: the 10 largest dives make 76% of the rest's nodes and 93% of its CPU (the largest 49%, 85%): more dives where they come from
WIDE: the jackknife over units gives a standard error of 85% of the rest's CPU
done so far: 7,793 nodes, 46.7 CPU h
estimated rest: 278,243 nodes (5-95%: 105,506 to 515,527), 3,310.0 CPU h (320.3 CPU h to 8,743.1 CPU h)
estimated tree: 286,036 nodes (5-95%: 113,299 to 523,320), 3,356.7 CPU h (366.9 CPU h to 8,789.7 CPU h)
2169 open nodes (377 with children that branch), 13183 dives in 377 of them, deepest dive node at depth 6; estimator 69.2 CPU h
2,169 of the 6,701 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 143,198 nodes, 2,823.3 CPU h
the largest dive makes 85% of the rest's CPU, the 10 largest 93%; 0 of 377 dived open nodes got one dive; estimator CPU: setup 28.8 CPU h, dives 40.4 CPU h
dives going deeper, by depth (reached, went on): 3: 4304/168, 4: 7971/90, 5: 1166/11, 6: 11/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 2.05 | 0.00 |
| 2 | 269 | 169 | 9.54 | 7.67 |
| 3 | 3,859 | 6,710 | 25.62 | 71.35 |
| 4 | 3,112 | 39,384 | 8.33 | 154.91 |
| 5 | 493 | 73,103 | 1.03 | 335.15 |
| 6 | 47 | 158,878 | 0.08 | 2,740.95 |
