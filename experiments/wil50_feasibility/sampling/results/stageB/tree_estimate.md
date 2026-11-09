# Tree estimate from stages A and B (400 of 400 sampled frontier nodes; 80 converged)

Frontier 14,192 depth-3 nodes (stage A), 651 nodes bounded at depths 0-2; 29,580 depth-3 orbits (Burnside), so 0.480 frontier nodes per orbit.

| size model | mean subtree per frontier node (95%) | total DNN nodes (95%) | per depth-3 orbit | tail (>= 300) share (95%) | CPU s at 22 s per DNN node [e] |
|---|---|---|---|---|---|
| anchors (central) | 34.6 (19.5-52.6) | 491,127 (278,033-747,795) | 16.6 | 0.44 (0.12-0.66) | 1.10e+07 |
| per-level model 47^(s/130) | 209.2 (74.3-386.4) | 2,969,982 (1,055,541-5,484,852) | 100.4 | 0.71 (0.29-0.89) | 6.68e+07 |

By class (anchors):

| class | sampled nodes | shortfall range | mean subtree | share of the total |
|---|---|---|---|---|
| fathomed | 297 | -1.0--1.0 | 1.0 | 0.021 |
| near (< 150) | 69 | 0.4-149.4 | 11.9 | 0.059 |
| 150-210 | 19 | 151.4-207.6 | 128.4 | 0.176 |
| 210-300 | 11 | 212.2-296.6 | 378.4 | 0.301 |
| tail (>= 300) | 4 | 312.2-338.4 | 1,526.7 | 0.442 |

Anchors (shortfall, DNN nodes): (0.0, 1), (158.2, 94), (315.7, 1,198), (404.8, 13,245); CPU per DNN node 22.5 s (node 10 and node 5 under cpcutt).
Against roadmap/tree_size.md: 1.0e4-6e4 DNN nodes per depth-3 orbit, 3e8-1.8e9 in all, 85-91% from the tail (the Oct 29 settings, the random-line 16-node sample).
