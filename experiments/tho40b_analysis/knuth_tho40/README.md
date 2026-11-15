# tho40-b's whole tree by Knuth's estimator (Nov 15)

The method of Brixius's thesis, chapter 4 (`bnbdrv -e`, `latest/estimate.C`), on the M2
Pro laptop for about three hours on Nov 15, with tho40-b's recipe
(`infra/phase1/recipes/tho40.par`, `SYM/tho40t.sym`, `-u 240518`).
- **Strata:** one process per depth-1 node, `-F 26:c`, 12 in all.
- **Each process:** bounds its node and the node's 39 depth-2 children exactly, then
  dives from depth 3 (`estimate.par`: start level 3, q = 10).
- **Data:** `s<c>/dives.txt` holds one line per dive (last depth, node estimate, laptop
  CPU estimate); `s<c>/est.out` the run's output.

`python3 strata.py` (root shard of tho40-b's live tree read at 20:45 UTC):

| depth-1 node | dives | nodes: mean (80% band) | one dive's share | merged by tho40-b | CPU, c7a days |
|---|---|---|---|---|---|
| 26:0 | 369 | 204,668 (116,291-311,234) | 0.37 | 41,555 | 8.9 |
| 26:1 | 310 | 1,250,608 (700,439-1,893,499) | 0.35 | 50,787 | 55.5 |
| 26:2 | 283 | 1,325,691 (905,301-1,806,210) | 0.20 | 57,780 | 53.1 |
| 26:3 | 262 | 1,680,991 (1,195,381-2,179,727) | 0.15 | 422,789 | 65.6 |
| 26:8 | 310 | 755,288 (514,706-1,009,602) | 0.19 | 57,092 | 32.9 |
| 26:9 | 265 | 2,189,907 (1,250,470-3,301,589) | 0.33 | 276,903 | 92.5 |
| 26:10 | 242 | 4,320,435 (2,843,236-5,899,654) | 0.19 | 80,409 | 180.2 |
| 26:11 | 268 | 9,428,431 (1,940,750-17,189,320) | 0.76 | 139,375 | 278.4 |
| 26:16 | 275 | 1,104,175 (742,185-1,458,388) | 0.18 | 244,560 | 44.6 |
| 26:17 | 247 | 6,311,404 (3,255,688-9,603,642) | 0.25 | 35,472 | 253.4 |
| 26:18 | 241 | 2,617,824 (1,669,310-3,698,461) | 0.26 | 337,210 | 110.3 |
| 26:19 | 250 | 3,243,881 (2,012,375-4,565,526) | 0.20 | 23,152 | 134.8 |
| **tree** | 3,322 | **3.4e7 (2.5e7-4.5e7)** | | 1.77e6 | **1,310 (1,049-1,614)** |

Laptop CPU is put in c7a CPU by the same nodes in tho40-b: the 12 depth-1 nodes and
the 142 depth-2 nodes it had merged, a factor of 2.23. Every stratum's estimate is
above what tho40-b had merged under it. 26:11's estimate is carried by one dive (76%
of its sum), so it is the least certain.

**Read it as likely low.** The bands come from the dives, and a heavy tail makes them
too narrow. On the A/B subtree, whose size is known, the same estimator read 0.53 of
the nodes and about 0.47 of the CPU after 173 dives (`../README.md`; the simulation's median at that
count is 0.79 and 0.88). The probable cause is how a dive weights children: by a
predicted bound (the parent's bound plus the parent's QPB increment `U_ij`,
`latest/branch.C`), not by each child's own bound, as the simulation does.

**What it means for the run.** tho40-b did 73.3 CPU days in 28.0 h of wall (63 CPU days per
day on its c7a.16xlarge). 1,310 CPU days is about 21 days of wall on one instance,
20 more from Nov 15 (80% band 16 to 25), and more if the estimate reads low. The
72-hour limit ends the first instance on Nov 17; the 7-day resume planned then covers a
third of what remains.

**Compared with the earlier estimates:**

| estimate | nodes |
|---|---|
| the Nov 8 estimate | 6e5 to 3e6 |
| tho40-b's rates projection, 27.4 h in | 9.5e6, still rising |
| this one | 3.4e7 |
