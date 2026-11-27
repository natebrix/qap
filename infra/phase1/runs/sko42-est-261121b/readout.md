# sko42-est-261121b

**COMPLETE**: the estimate (sko42-est-261121b-dives finished)

## The run

| | |
|---|---|
| problem | sko42t.dat -s SYM/sko42t.sym |
| params | sko42.par |
| incumbent | 15814.0 |
| slots | 263; units of 1200.0 CPU s, traced |
| created | 2026-11-21 04:50 UTC |
| status | paused |
| estimate | dives run sko42-est-261121b-dives: k1.5, sample 1.0 (unif), QAP_EST_MIN_DIVES=5, QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 |

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
the estimate's cost: 70.0 CPU h of estimator CPU, 25% of the estimated rest of the tree (277.5 CPU h)
done so far: 7,793 nodes, 46.7 CPU h
estimated rest: 94,047 nodes (5-95%: 74,827 to 125,846), 277.5 CPU h (226.4 CPU h to 357.1 CPU h)
estimated tree: 101,840 nodes (5-95%: 82,620 to 133,639), 324.1 CPU h (273.1 CPU h to 403.7 CPU h)
2246 open nodes (382 with children that branch), 13949 dives in 382 of them, deepest dive node at depth 6; estimator 70.0 CPU h
2,246 of the 6,701 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 13,764 nodes, 39.5 CPU h
the largest dive makes 6% of the rest's CPU, the 10 largest 22%; 0 of 382 dived open nodes got one dive; estimator CPU: setup 28.7 CPU h, dives 41.2 CPU h
dives going deeper, by depth (reached, went on): 3: 4469/179, 4: 7883/85, 5: 1809/18, 6: 70/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 2.05 | 0.00 |
| 2 | 269 | 169 | 9.54 | 7.67 |
| 3 | 3,859 | 6,771 | 25.62 | 72.06 |
| 4 | 3,112 | 40,976 | 8.33 | 117.19 |
| 5 | 493 | 30,788 | 1.03 | 66.59 |
| 6 | 47 | 15,343 | 0.08 | 13.99 |
