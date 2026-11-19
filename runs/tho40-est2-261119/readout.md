# tho40-est2-261119

**PART-WAY**: the estimate's dives run tho40-est2-261119-dives is running, 527 of 548 dives units merged

- PART-WAY: 21 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (8% of the rest's CPU)
- SLOW PENDING: 19 pending units have run longer than 90% of the merged ones took (5890 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge

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
PART-WAY: 21 of 548 dives units not merged; their open nodes are imputed, depth by depth, from the merged units' (8% of the rest's CPU)
SLOW PENDING: 19 pending units have run longer than 90% of the merged ones took (5890 s): a unit runs long when its open nodes branch, so these likely hold the large subtrees, and the imputation stands in for them with small ones: read the estimate as a floor until they merge
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 1,611,542 nodes (5-95%: 1,308,186 to 1,943,911), 4,936.6 CPU h (4,392.6 CPU h to 5,512.9 CPU h)
estimated tree: 1,630,309 nodes (5-95%: 1,326,953 to 1,962,678), 5,053.8 CPU h (4,509.7 CPU h to 5,630.0 CPU h)
5771 open nodes (2175 with children that branch), 22106 dives in 2175 of them, deepest dive node at depth 7; estimator 316.1 CPU h
5,771 of the 17,587 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 246,229 nodes, 371.3 CPU h
the largest dive makes 4% of the rest's CPU, the 10 largest 14%; 0 of 2175 dived open nodes got one dive; estimator CPU: setup 132.9 CPU h, dives 183.2 CPU h
dives going deeper, by depth (reached, went on): 3: 370/147, 4: 13538/1136, 5: 9387/149, 6: 243/2, 7: 2/0


(nodes at depth 8 are in the solve's own part; the dives reached depth 7)
```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 81 | 45.12 | 6.99 |
| 3 | 2,924 | 9,608 | 28.66 | 274.03 |
| 4 | 8,475 | 158,564 | 30.50 | 1,754.66 |
| 5 | 5,851 | 561,570 | 10.21 | 1,667.20 |
| 6 | 1,004 | 737,868 | 1.27 | 821.14 |
| 7 | 148 | 11,887 | 0.16 | 6.53 |
| 8 | 33 | 0 | 0.03 | 0.00 |
