# tai50b-est-261119

**PART-WAY**: the estimate's dives run tai50b-est-261119-dives is running, 252 of 259 dives units merged

- PART-WAY: 7 of 259 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (2% of the rest's CPU)
- SETUP-BOUND: bounding the sampled open nodes took 18.7 CPU h, the dives 0.0 CPU s (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

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
PART-WAY: 7 of 259 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (2% of the rest's CPU)
SETUP-BOUND: bounding the sampled open nodes took 18.7 CPU h, the dives 0.0 CPU s (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 1,117 nodes, 32.8 CPU h
estimated rest: 754 nodes (5-95%: 754 to 754), 19.2 CPU h (19.2 CPU h to 19.2 CPU h)
estimated tree: 1,871 nodes (5-95%: 1,871 to 1,871), 51.9 CPU h (51.9 CPU h to 51.9 CPU h)
275 open nodes (0 with children that branch), 0 dives in 0 of them, deepest dive node at depth 0; estimator 18.7 CPU h
jackknife over units: standard error 130 nodes, 2.0 CPU h
the largest dive makes 0% of the rest's CPU, the 10 largest 0%; 0 of 0 dived open nodes got one dive; estimator CPU: setup 18.7 CPU h, dives 0.0 CPU s

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 1.20 | 0.00 |
| 1 | 50 | 0 | 7.25 | 0.00 |
| 2 | 490 | 0 | 16.60 | 0.00 |
| 3 | 576 | 0 | 7.72 | 0.00 |
| 4 | 0 | 275 | 0.00 | 14.75 |
| 5 | 0 | 460 | 0.00 | 3.93 |
