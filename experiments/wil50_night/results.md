## a220: base 232 nodes, 13,063 DNN iterations, 3,860 CPU s

| kind of node (base) | nodes | first-solve iterations, CPU s | phase iterations, CPU s | share of DNN CPU |
|---|---|---|---|---|
| phase misses | 5 | 740, 217 | 4,140, 1,342 | 41% |
| first solve fathoms | 224 | 6,583, 1,815 | 0, 0 | 47% |
| phase fathoms | 3 | 330, 94 | 1,270, 374 | 12% |

First solves 55% of the DNN CPU, phases 45%; DNN 100% of the run's CPU.

| arm | nodes | DNN iterations (first + phase) | CPU s | iterations / base | CPU / base | nodes both bounded: iterations / base |
|---|---|---|---|---|---|---|
| base | 232 | 13,063 (7,653 + 5,410) | 3,860 | 1.000 | 1.000 | 232 nodes: 1.000 |
| extrap | 232 | 12,865 (7,755 + 5,110) | 3,874 | 0.985 | 1.004 | 232 nodes: 0.985 |
| g150 | 232 | 12,732 (7,740 + 4,992) | 4,133 | 0.975 | 1.071 | 232 nodes: 0.975 |
| g150x | 232 | 12,407 (7,864 + 4,543) | 3,880 | 0.950 | 1.005 | 232 nodes: 0.950 |
| gw150 | 232 | 13,122 (7,982 + 5,140) | 3,966 | 1.005 | 1.027 | 232 nodes: 1.005 |
| refresh80 | 232 | 12,987 (7,464 + 5,523) | 3,877 | 0.994 | 1.005 | 232 nodes: 0.994 |
| wc025 | 232 | 14,563 (7,732 + 6,831) | 4,753 | 1.115 | 1.231 | 232 nodes: 1.115 |
| wc2g150 | 232 | 12,200 (7,877 + 4,323) | 4,004 | 0.934 | 1.037 | 186 nodes: 0.936 |
| wcg150 | 232 | 12,122 (7,662 + 4,460) | 3,660 | 0.928 | 0.948 | 186 nodes: 0.932 |
| wcg250 | 232 | 12,122 (7,662 + 4,460) | 4,028 | 0.928 | 1.043 | 186 nodes: 0.932 |
| wcg80 | 232 | 12,122 (7,662 + 4,460) | 4,063 | 0.928 | 1.052 | 186 nodes: 0.932 |
| xl2 | 232 | 12,959 (7,549 + 5,410) | 3,967 | 0.992 | 1.028 | 186 nodes: 1.000 |

## b83: base 181 nodes, 12,609 DNN iterations, 3,018 CPU s

| kind of node (base) | nodes | first-solve iterations, CPU s | phase iterations, CPU s | share of DNN CPU |
|---|---|---|---|---|
| phase misses | 4 | 640, 151 | 3,160, 861 | 34% |
| first solve fathoms | 169 | 5,578, 1,222 | 0, 0 | 41% |
| phase fathoms | 8 | 880, 183 | 2,351, 587 | 26% |

First solves 52% of the DNN CPU, phases 48%; DNN 100% of the run's CPU.

| arm | nodes | DNN iterations (first + phase) | CPU s | iterations / base | CPU / base | nodes both bounded: iterations / base |
|---|---|---|---|---|---|---|
| base | 181 | 12,609 (7,098 + 5,511) | 3,018 | 1.000 | 1.000 | 181 nodes: 1.000 |
| extrap | 181 | 12,519 (7,099 + 5,420) | 3,047 | 0.993 | 1.010 | 181 nodes: 0.993 |
| g150 | 181 | 11,252 (7,042 + 4,210) | 2,847 | 0.892 | 0.943 | 181 nodes: 0.892 |
| g150x | 181 | 11,032 (7,038 + 3,994) | 2,724 | 0.875 | 0.902 | 181 nodes: 0.875 |
| gw150 | 181 | 12,257 (9,177 + 3,080) | 2,926 | 0.972 | 0.969 | 181 nodes: 0.972 |
| refresh80 | 181 | 13,242 (6,935 + 6,307) | 3,180 | 1.050 | 1.054 | 181 nodes: 1.050 |
| wc025 | 181 | 11,833 (7,118 + 4,715) | 3,074 | 0.938 | 1.018 | 181 nodes: 0.938 |
| wc2g150 | 181 | 11,328 (7,153 + 4,175) | 2,964 | 0.898 | 0.982 | 181 nodes: 0.898 |
| wcg150 | 181 | 11,277 (7,081 + 4,196) | 2,699 | 0.894 | 0.894 | 181 nodes: 0.894 |
| wcg250 | 181 | 11,833 (7,118 + 4,715) | 3,138 | 0.938 | 1.040 | 181 nodes: 0.938 |
| wcg80 | 181 | 11,277 (7,081 + 4,196) | 3,006 | 0.894 | 0.996 | 181 nodes: 0.894 |
| xl2 | 137 | 9,162 (5,377 + 3,785) | 2,304 | 0.727 | 0.763 | 47 nodes: 1.000 |

## c135: base 28 nodes, 2,264 DNN iterations, 770 CPU s

| kind of node (base) | nodes | first-solve iterations, CPU s | phase iterations, CPU s | share of DNN CPU |
|---|---|---|---|---|
| phase misses | 1 | 300, 85 | 970, 362 | 58% |
| first solve fathoms | 27 | 994, 320 | 0, 0 | 42% |

First solves 53% of the DNN CPU, phases 47%; DNN 100% of the run's CPU.

| arm | nodes | DNN iterations (first + phase) | CPU s | iterations / base | CPU / base | nodes both bounded: iterations / base |
|---|---|---|---|---|---|---|
| base | 28 | 2,264 (1,294 + 970) | 770 | 1.000 | 1.000 | 28 nodes: 1.000 |
| xl2 | 28 | 2,264 (1,294 + 970) | 811 | 1.000 | 1.053 | 28 nodes: 1.000 |

## d178: base 140 nodes, 10,817 DNN iterations, 3,688 CPU s

| kind of node (base) | nodes | first-solve iterations, CPU s | phase iterations, CPU s | share of DNN CPU |
|---|---|---|---|---|
| phase misses | 3 | 520, 161 | 1,830, 688 | 23% |
| first solve fathoms | 131 | 4,296, 1,384 | 0, 0 | 38% |
| phase fathoms | 6 | 660, 211 | 3,511, 1,229 | 39% |

First solves 48% of the DNN CPU, phases 52%; DNN 100% of the run's CPU.

| arm | nodes | DNN iterations (first + phase) | CPU s | iterations / base | CPU / base | nodes both bounded: iterations / base |
|---|---|---|---|---|---|---|
| base | 140 | 10,817 (5,476 + 5,341) | 3,688 | 1.000 | 1.000 | 140 nodes: 1.000 |
| xl2 | 48 | 4,179 (2,665 + 1,514) | 1,450 | 0.386 | 0.393 | 1 nodes: 1.000 |

## e212: base 186 nodes, 11,220 DNN iterations, 3,681 CPU s

| kind of node (base) | nodes | first-solve iterations, CPU s | phase iterations, CPU s | share of DNN CPU |
|---|---|---|---|---|
| phase misses | 4 | 630, 205 | 3,330, 1,218 | 39% |
| first solve fathoms | 176 | 5,560, 1,683 | 0, 0 | 46% |
| phase fathoms | 6 | 660, 197 | 1,040, 359 | 15% |

First solves 57% of the DNN CPU, phases 43%; DNN 99% of the run's CPU.

| arm | nodes | DNN iterations (first + phase) | CPU s | iterations / base | CPU / base | nodes both bounded: iterations / base |
|---|---|---|---|---|---|---|
| base | 186 | 11,220 (6,850 + 4,370) | 3,681 | 1.000 | 1.000 | 186 nodes: 1.000 |
| xl2 | 186 | 11,220 (6,850 + 4,370) | 3,594 | 1.000 | 0.976 | 186 nodes: 1.000 |

