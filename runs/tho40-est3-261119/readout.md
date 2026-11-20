# tho40-est3-261119

A second estimate of tho40-est2-261119: its paused frontier dived again (its own readout has its first).

**COMPLETE**: the estimate (tho40-est3-261119 finished)

- HEAVY TAIL: the 10 largest dives make 52% of the rest's nodes and 35% of its CPU (the largest 13%, 7%): more dives where they come from
- SETUP-BOUND: bounding the sampled open nodes took 149.9 CPU h, the dives 143.1 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)

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
the estimate's cost: 293.0 CPU h of estimator CPU, 4% of the estimated rest of the tree (6,726.9 CPU h)
HEAVY TAIL: the 10 largest dives make 52% of the rest's nodes and 35% of its CPU (the largest 13%, 7%): more dives where they come from
SETUP-BOUND: bounding the sampled open nodes took 149.9 CPU h, the dives 143.1 CPU h (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)
done so far: 18,767 nodes, 117.2 CPU h
estimated rest: 3,685,908 nodes (5-95%: 2,199,210 to 5,662,799), 6,726.9 CPU h (5,291.7 CPU h to 8,645.7 CPU h)
estimated tree: 3,704,675 nodes (5-95%: 2,217,977 to 5,681,566), 6,844.1 CPU h (5,408.9 CPU h to 8,762.9 CPU h)
6240 open nodes (2574 with children that branch), 33655 dives in 2574 of them, deepest dive node at depth 8; estimator 293.0 CPU h
6,240 of the 18,578 open nodes sampled (each weighted 1/pi)
jackknife over units: standard error 779,429 nodes, 859.9 CPU h
the largest dive makes 7% of the rest's CPU, the 10 largest 35%; 0 of 2574 dived open nodes got one dive; estimator CPU: setup 149.9 CPU h, dives 143.1 CPU h
dives going deeper, by depth (reached, went on): 3: 433/144, 4: 22130/1640, 5: 12692/227, 6: 399/7, 7: 19/1, 8: 1/0

```

| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.00 | 0.00 |
| 1 | 12 | 0 | 1.23 | 0.00 |
| 2 | 319 | 85 | 45.12 | 7.10 |
| 3 | 2,924 | 10,453 | 28.66 | 269.45 |
| 4 | 8,475 | 178,941 | 30.50 | 1,150.31 |
| 5 | 5,851 | 762,000 | 10.21 | 2,628.83 |
| 6 | 1,004 | 1,192,671 | 1.27 | 1,657.77 |
| 7 | 148 | 1,071,298 | 0.16 | 755.74 |
| 8 | 33 | 470,460 | 0.03 | 257.74 |
