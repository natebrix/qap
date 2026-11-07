# Stage B: 16 of 16 sampled nodes bounded (frontier 552 nodes, complete; seed 20261107; sample of 16)
# binary 4d8693d0d9, give-up -x 1:4:1:300:1, converge -x 1:4:1:0:1, phase off at depth 3

| node | parent | inherited | first lb | short | status | iters | CPU s | warm | pool | children | converged lb | converged short | conv iters | conv CPU s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0:8,4:24,6:22 | r.5.23 | 6082.1 | 6163.15 | -41.1 | fathomed | 20 | 1 | 1 | 0 | 0 |  | | | |
| 0:7,4:21,14:24 | r.4.13 | 6066.5 | 6136.71 | -14.7 | fathomed | 10 | 0 | 1 | 0 | 0 |  | | | |
| 0:1,9:24,16:0 | r.1.15 | 6081.9 | 6153.69 | -31.6 | fathomed | 10 | 0 | 1 | 0 | 0 |  | | | |
| 0:8,4:0,6:14 | r.5.0 | 6081.5 | 6100.47 | 21.6 | hopeless | 300 | 6 | 1 | 0 | 27 | 6101.61 | 20.4 | 4000 | 73 |
| 0:12,8:24,23:6 | r.6.6 | 6055.3 | 6131.33 | -9.3 | fathomed | 20 | 1 | 1 | 0 | 0 |  | | | |
| 0:13,5:24,6:8 | r.7.8 | 6074.6 | 6123.31 | -1.3 | fathomed | 20 | 1 | 1 | 0 | 0 |  | | | |
| 0:2,4:5,21:14 | r.2.4 | 6085.1 | 6112.46 | 9.6 | hopeless | 300 | 6 | 1 | 0 | 27 | 6113.54 | 8.5 | 4000 | 70 |
| 0:2,4:29,27:0 | r.2.28 | 6084.0 | 6122.10 | -0.1 | fathomed | 30 | 1 | 1 | 0 | 0 |  | | | |
| 0:1,4:16,23:0 | r.1.22 | 6074.0 | 6123.31 | -1.3 | fathomed | 10 | 0 | 1 | 0 | 0 |  | | | |
| 0:13,6:8,22:24 | r.7.8 | 6074.6 | 6122.67 | -0.6 | fathomed | 40 | 1 | 1 | 0 | 0 |  | | | |
| 0:14,6:8,8:24 | r.8.8 | 6088.7 | 6136.75 | -14.7 | fathomed | 10 | 0 | 1 | 0 | 0 |  | | | |
| 0:7,4:24,24:5 | r.4.3 | 6074.5 | 6122.89 | -0.8 | fathomed | 35 | 1 | 1 | 0 | 0 |  | | | |
| 0:2,4:29,24:0 | r.2.28 | 6084.0 | 6122.24 | -0.2 | fathomed | 26 | 1 | 1 | 0 | 0 |  | | | |
| 0:1,4:18,23:0 | r.1.22 | 6074.0 | 6122.08 | -0.0 | fathomed | 48 | 1 | 1 | 0 | 0 |  | | | |
| 0:8,4:24,6:14 | r.5.23 | 6082.1 | 6103.49 | 18.6 | hopeless | 300 | 6 | 1 | 0 | 27 | 6104.62 | 17.4 | 4000 | 72 |
| 0:8,4:0,6:13 | r.5.0 | 6081.5 | 6122.58 | -0.5 | fathomed | 35 | 1 | 1 | 0 | 0 |  | | | |

Cutoff 6122.05. Classes (converged shortfall for the 3 converged nodes, the give-up shortfall, an upper bound on it, for the others):

| class | nodes | fraction | 1-sigma interval |
|---|---|---|---|
| fathomed | 13 | 0.812 | 0.698-0.891 |
| near (< 150) | 3 | 0.188 | 0.109-0.302 |
| 150-210 | 0 | 0.000 | 0.000-0.059 |
| 210-300 | 0 | 0.000 | 0.000-0.059 |
| tail (>= 300) | 0 | 0.000 | 0.000-0.059 |

Give-up shortfall of the 3 misses: min 9.6, 10% 9.6, median 18.6, 90% 21.6, max 21.6

CPU: give-up 26 s on 16 nodes (2 s each), converge 215 s on 3 nodes
Restart check 0:8,4:24,6:22: no warm start 6125.6529 (66 it), -F 6125.6529 (66 it): same
Calibration 0:7,13:24,22:5 Fgiveup: lb 6132.07, -10.0 short, 55 iterations, fathomed, CPU 1 s
Calibration 0:7,13:24,22:5 Fconverge: lb 6132.07, -10.0 short, 55 iterations, fathomed, CPU 1 s
