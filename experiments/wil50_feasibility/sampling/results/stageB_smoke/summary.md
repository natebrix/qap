# Stage B: 4 of 4 sampled nodes bounded (frontier 14192 nodes, complete; seed 20261107; sample of 400)
# binary a667727601, give-up -x 1:4:1:300:1, converge -x 1:4:1:0:1, phase off at depth 3

| node | parent | inherited | first lb | short | status | iters | CPU s | warm | pool | children | converged lb | converged short | conv iters | conv CPU s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 21:2,27:0,45:45 | r.2.26 | 48678.2 | 48944.63 | -130.6 | fathomed | 20 | 8 | 1 | 1 | 0 |  | | | |
| 21:13,24:0,34:40 | r.8.23 | 48788.5 | 48815.63 | -1.6 | fathomed | 30 | 11 | 1 | 1 | 0 |  | | | |
| 21:13,28:0,45:14 | r.8.27 | 48734.4 | 48816.07 | -2.0 | fathomed | 20 | 8 | 1 | 1 | 0 |  | | | |
| 21:2,25:0,39:49 | r.2.24 | 48524.1 | 48740.74 | 73.3 | hopeless | 300 | 103 | 1 | 0 | 47 | 48745.54 | 68.5 | 4000 | 1267 |

Cutoff 48814.05. Classes (converged shortfall for the 1 converged nodes, the give-up shortfall, an upper bound on it, for the others):

| class | nodes | fraction | 1-sigma interval |
|---|---|---|---|
| fathomed | 3 | 0.750 | 0.500-0.900 |
| near (< 150) | 1 | 0.250 | 0.100-0.500 |
| 150-210 | 0 | 0.000 | 0.000-0.200 |
| 210-300 | 0 | 0.000 | 0.000-0.200 |
| tail (>= 300) | 0 | 0.000 | 0.000-0.200 |

Give-up shortfall of the 1 misses: min 73.3, 10% 73.3, median 73.3, 90% 73.3, max 73.3

CPU: give-up 129 s on 4 nodes (32 s each), converge 1267 s on 1 nodes
Restart check 21:2,27:0,45:45: no warm start 48828.0264 (98 it), -F 48828.0264 (98 it): same
Calibration 49:23,15:37,5:6 Fgiveup: lb 48450.48, 363.6 short, 300 iterations, hopeless, CPU 88 s
Calibration 49:23,15:37,5:6 Fconverge: lb 48498.65, 315.4 short, 4000 iterations, iter_limit, CPU 1174 s

Projection for the full stage B [e]: 400 give-ups x 32 s (mean of 4) = 12,932 s; 80 converges x 1267 s (measured on 1) = 101,361 s; 9 parents to rebuild, 2,408 s at their stage-A CPU; total 116,701 CPU s, 8.2 h of wall on 4 slots; if every give-up cost what the hopeless ones did (103 s): 144,915 CPU s, 10.2 h
