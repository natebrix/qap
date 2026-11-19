# tai50b-est-261119

**COMPLETE**: the estimate (tai50b-est-261119-dives finished)

- HEAVY TAIL: the 10 largest dives make 65% of the rest's nodes and 44% of its CPU (the largest 42%, 27%): more dives where they come from
- FEW DIVES: 7 of 7 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
- SETUP-BOUND: bounding the sampled open nodes took 28.0 CPU h, the dives 9.5 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
- WIDE: the jackknife over units gives a standard error of 31% of the rest's CPU

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
HEAVY TAIL: the 10 largest dives make 65% of the rest's nodes and 44% of its CPU (the largest 42%, 27%): more dives where they come from
FEW DIVES: 7 of 7 dived open nodes got one dive (QAP_EST_MIN_DIVES, a longer dive budget)
SETUP-BOUND: bounding the sampled open nodes took 28.0 CPU h, the dives 9.5 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
WIDE: the jackknife over units gives a standard error of 31% of the rest's CPU
done so far: 1,117 nodes, 32.8 CPU h
estimated rest: 2,959 nodes (5-95%: 2,959 to 2,959), 44.4 CPU h (44.4 CPU h to 44.4 CPU h)
estimated tree: 4,076 nodes (5-95%: 4,076 to 4,076), 77.2 CPU h (77.2 CPU h to 77.2 CPU h)
282 open nodes (7 with children that branch), 7 dives in 7 of them, deepest dive node at depth 11; estimator 37.5 CPU h
jackknife over units: standard error 1,292 nodes, 13.7 CPU h
the largest dive makes 27% of the rest's CPU, the 10 largest 44%; 7 of 7 dived open nodes got one dive; estimator CPU: setup 28.0 CPU h, dives 9.5 CPU h

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 1.20 | 0.00 |
| 1 | 50 | 0 | 7.25 | 0.00 |
| 2 | 490 | 0 | 16.60 | 0.00 |
| 3 | 576 | 0 | 7.72 | 0.00 |
| 4 | 0 | 282 | 0.00 | 16.62 |
| 5 | 0 | 783 | 0.00 | 11.53 |
| 6 | 0 | 663 | 0.00 | 7.13 |
| 7 | 0 | 344 | 0.00 | 3.75 |
| 8 | 0 | 293 | 0.00 | 2.23 |
| 9 | 0 | 203 | 0.00 | 1.20 |
| 10 | 0 | 198 | 0.00 | 1.26 |
| 11 | 0 | 193 | 0.00 | 0.72 |
