# tho40-est3-261119

A second estimate of tho40-est2-261119: its paused frontier dived again (its own readout has its first).

**PART-WAY**: the estimate's dives run tho40-est3-261119 is running, 458 of 548 dives units merged

- PART-WAY: 90 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (31% of the rest's CPU)
- SLOW PENDING: 22 pending units have run longer than 90% of the merged ones took (3092 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge

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
PART-WAY: 90 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (31% of the rest's CPU)
SLOW PENDING: 22 pending units have run longer than 90% of the merged ones took (3092 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 2,104,526 nodes (5-95%: 1,176,845 to 3,922,885), 4,506.6 CPU h (3,434.3 CPU h to 6,445.1 CPU h)
estimated tree: 2,123,293 nodes (5-95%: 1,195,612 to 3,941,652), 4,623.8 CPU h (3,551.5 CPU h to 6,562.2 CPU h)
4457 open nodes (1429 with children that branch), 26910 dives in 1429 of them, deepest dive node at depth 8; estimator 181.3 CPU h
4,457 of the 15,155 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 778,470 nodes, 964.3 CPU h
the largest dive makes 10% of the rest's CPU, the 10 largest 31%; 0 of 1429 dived open nodes got one dive; estimator CPU: setup 87.6 CPU h, dives 93.7 CPU h
dives going deeper, by depth (reached, went on): 3: 283/84, 4: 15703/691, 5: 11515/156, 6: 328/4, 7: 16/1, 8: 1/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 52 | 45.12 | 4.56 |
| 3 | 2,924 | 6,899 | 28.66 | 154.09 |
| 4 | 8,475 | 98,538 | 30.50 | 607.53 |
| 5 | 5,851 | 320,555 | 10.21 | 1,146.06 |
| 6 | 1,004 | 250,800 | 1.27 | 530.30 |
| 7 | 148 | 398,165 | 0.16 | 424.36 |
| 8 | 33 | 470,460 | 0.03 | 257.74 |
