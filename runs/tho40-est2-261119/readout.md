# tho40-est2-261119

**COMPLETE**: the estimate (tho40-est2-261119-dives finished)

## The run

| | |
|---|---|
| problem | tho40t.dat -s SYM/tho40t.sym |
| params | tho40c.par |
| incumbent | 240518.0 |
| slots | 135; units of 1200.0 CPU s, traced |
| created | 2026-11-19 02:26 UTC |
| status | paused |
| estimate | dives run tho40-est2-261119-dives: k1.5, sample 1.0 (unif), QAP_EST_MIN_DIVES=5, QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 |

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
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 2,053,682 nodes (5-95%: 1,735,371 to 2,430,901), 5,610.5 CPU h (5,064.9 CPU h to 6,138.1 CPU h)
estimated tree: 2,072,449 nodes (5-95%: 1,754,138 to 2,449,668), 5,727.7 CPU h (5,182.1 CPU h to 6,255.2 CPU h)
6307 open nodes (2575 with children that branch), 24106 dives in 2575 of them, deepest dive node at depth 7; estimator 373.2 CPU h
6,307 of the 18,578 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 199,950 nodes, 416.5 CPU h
the largest dive makes 3% of the rest's CPU, the 10 largest 13%; 0 of 2575 dived open nodes got one dive; estimator CPU: setup 156.4 CPU h, dives 216.8 CPU h
dives going deeper, by depth (reached, went on): 3: 390/161, 4: 15437/1609, 5: 9955/181, 6: 275/3, 7: 3/0


(nodes at depth 8 are in the solve's own part; the dives reached depth 7)
```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 85 | 45.12 | 7.29 |
| 3 | 2,924 | 10,535 | 28.66 | 316.08 |
| 4 | 8,475 | 184,741 | 30.50 | 2,143.98 |
| 5 | 5,851 | 773,167 | 10.21 | 2,108.26 |
| 6 | 1,004 | 1,043,469 | 1.27 | 1,014.49 |
| 7 | 148 | 41,686 | 0.16 | 20.40 |
| 8 | 33 | 0 | 0.03 | 0.00 |
