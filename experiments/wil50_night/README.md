# wil50_night: screening settings on complete wil50 subtrees (Nov 24)

Nathan, Nov 24: "I want to gun for any remaining perf improvements - either in the solver itself, or via parameter
tuning." Each arm is `par/base.par` (a copy of cognac's `infra/phase1/recipes/wil50_d2r10.par`) plus the lines at
its end, run cold from `-F` on complete subtrees of stage B's depth-3 frontier nodes, one process a run, four at a
time (`run.sh`, `queue.sh`, `jobs.txt`), incumbent 48,818, solver as cognac 626230c8. `python3 analyze.py` prints
the tables at the end (`results.md`).

Subtrees: a220 `21:22,31:23,40:15` and b83 `8:9,21:11,45:40,34:0` (the bound lab's two,
`experiments/bound_lab/wil50_policy`), and for the branching arm c135 `21:20,31:24,45:22`, d178 `21:3,25:49,45:0`,
e212 `13:40,21:1,25:0` (stage B nodes 135, 178 and 212 short when converged).

**Read the iterations, not the CPU.** The container has four cores, all in use, and its CPU times moved by about
10% between runs that did the same work: `wcg80` and `wcg150` are the same run on both subtrees (the same nodes
and iterations) and took 4,063 and 3,660 s on a220. The container also runs the AMX kernel, which c7a does not.

## What came out

DNN iterations (first solves + cut phases) against base, summed over the subtrees an arm ran on:

| arm | lines added to base | subtrees | nodes | iterations / base |
|---|---|---|---|---|
| `xl2` | `branch.xbar_lines = 2` | all five | 631 for 767 | **0.80** |
| `g150x` | `cut.cp_wc_gap = 150`, `cut.cp_extrap = 100:1:30:2` | all five | 813 for 767 | **0.93** (a220 and b83: 0.91) |
| `g150` | `cut.cp_wc_gap = 150` | a220, b83 | the same | 0.93 |
| `wcg150`, `wcg80` | wc = 1.4 and the gap key at 150, 80 | a220, b83 | the same | 0.91 |
| `wc2g150` | wc = 2.0, gap 150 | a220, b83 | the same | 0.92 |
| `wcg250` | wc = 1.4, gap 250 | a220, b83 | the same | 0.93 |
| `extrap` | `cut.cp_extrap = 100:1:30:2` | a220, b83 | the same | 0.99 |
| `gw150` | `dnn.giveup_warm = 150` | a220, b83 | the same | 0.99 |
| `refresh80` | `dnn.refresh = 80` | a220, b83 | the same | 1.02 |
| `wc025` | wc = 0.25 (every phase at w) | a220, b83 | the same | 1.03 |

1. **`branch.xbar_lines = 2`: two of five subtrees shrink, none grows.** d178 140 -> 48 nodes (0.39 of the
   iterations: the subtree's own root takes c9 for c48), b83 181 -> 137 (0.73: two depth-5 nodes take c1 for c19
   and c20); a220 keeps 232 nodes with another line at one depth-4 node (0.99); c135 and e212 are the same runs
   as base. The five subtrees have 17 branching nodes under base and four of them change line, so this is four line
   choices, and the setting is
   global: it also widens rule 9's shortlist at depth 1, where it changes no line (wave 2, below). On tho40, 2 against
   4 was neutral on a solved sample of 150 nodes (`experiments/estimates_nov19`); 1 against 2 there was 723 -> 581
   nodes on three depth-3 subtrees. A solved sample is the test that would settle it for wil50.
2. **Near phases at w: `cut.cp_wc_gap = 150`, 0.93.** A cut phase that starts within the gap of the cutoff keeps
   the children's penalty factor (w = 0.25) for wc = 0.7. Every phase below a subtree's root starts within 80 of
   the cutoff, so 80 and 150 are the same run; the gain is 0.975 on a220 and 0.89 on b83, in the phases alone
   (first solves unchanged, the same nodes).
3. **Far phases: wc stays 0.7.** The only far phases here are the two subtree roots (254 and 202 short after the
   first solve). Shortfall by phase iteration:

   | root | wc | 200 | 400 | 600 | 800 | end |
   |---|---|---|---|---|---|---|
   | a220 | 0.7 (base) | 174 | 138 | 120 | 111 | 103 at 1,120 |
   | a220 | 1.4 | 173 | 132 | 114 | 104 | 98 at 1,050 |
   | a220 | 2.0 | 177 | 133 | 115 | 105 | 98 at 1,050 |
   | a220 | 0.25 | 190 | 163 | 147 | 137 | 117 at 1,490 |
   | a220 | 0.7 + extrap | 169 | 129 | 111 | 103 | 98 at 1,030 |
   | b83 | 0.7 (base) | 123 | 89 | 71 | 62 | 56 at 1,010 |
   | b83 | 1.4 | 134 | 100 | 82 | 71 | 62 at 1,090 |
   | b83 | 2.0 | 141 | 107 | 90 | 80 | 69 at 1,140 |
   | b83 | 0.25 | 113 | 87 | 75 | 68 | 66 at 910 |
   | b83 | 0.7 + extrap | 120 | 85 | 67 | 58 | 56 at 900 |

   wc = 1.4 is a little ahead on one root and behind on the other, 2.0 is no better, 0.25 stalls short on both
   (so the key matters: `wc025`, every phase at w, loses on a220 what it gains below). `wcg150`'s 0.93 on a220
   against `g150`'s 0.975 is its root phase ending 5 nearer and the tree below taking other lines (186 of the 232
   nodes in common), not the near phases.
4. **`cut.cp_extrap` helps the far phases** (the same end in 900 iterations for 1,010 on b83; 5 nearer on a220),
   which is 1% of a subtree with one far phase alone and 2% on top of the gap key (`g150x` 0.95 and 0.875).
   Every depth-2 and depth-3 node of the whole tree runs a far phase, so its share there is larger than here.
5. **Nothing in `dnn.giveup_warm = 150` or `dnn.refresh = 80`**: within 3% either way, opposite signs on the two
   subtrees.

Where base spends its iterations: first solves 55% and 52% of the DNN CPU on a220 and b83 (first-solve fathoms
47% and 41%), cut phases 45% and 48%; the DNN is all of the run's CPU.

## Wave 2

**The top of the tree does not move under `branch.xbar_lines = 2`.** Each of the root's 15 children, cold, under
base and under `xl2` (`top.sh`, `top/`): the same line and the same bounds at every one (`21:0`, `21:11`, `21:20`,
`21:21` r45; `21:22`, `21:23`, `21:24` r31; the other eight c0), and the root keeps r21. Depth 2 is rule 10, which
has no shortlist. So the setting acts at depth 3 and below, which is where a solved sample tests it.

**`g150x` on all five subtrees: 0.93 of the iterations, 813 nodes for 767.** By subtree 0.95 (a220), 0.875 (b83),
0.99 (c135), 0.97 (d178), 0.93 (e212). The extra nodes are one miss: d178's depth-4 node `24:48` starts 36 short;
base's phase grinds it to a fathom in 2,130 iterations, `g150x`'s at w = 0.25 stalls 2.4 short after 1,330 and
the node branches into 46 children that fathom in their first solves. On the 140 nodes both arms bound, 0.88.

**The two together (`g150xl2`) keep both gains** on the two subtrees `xl2` shrank: b83 137 nodes and 0.65 of
base's iterations (`xl2` alone 0.73), d178 48 nodes and 0.36 (0.39).

## The best-guess recipe (Nov 24)

Nathan, Nov 24: "I'd rather go aggressive here." `wil50-est-261124` runs under
`infra/phase1/recipes/wil50_d2r10_g150xl2.par` (cognac 0a4c0f83): `wil50_d2r10.par` with `cut.cp_wc_gap = 150`,
`cut.cp_extrap = 100:1:30:2` and `branch.xbar_lines = 2`, on four workers from a 1.5 h solve.

**Its path check: PASS** (`pathcheck/warm_wil50_d2r10_g150xl2.tsv`; `experiments/tho40b_analysis/path_check/dive.py
--warm`, from 48,818 toward `wil50.sln`): 48 nodes to depth 47, every one branched and kept the path child, the
leaf enumerated 48,816; the highest bound on the path 48,816.0000; 3,442 CPU s (`wil50_d2r10`: 3,389). The lines
are `wil50_d2r10`'s down to depth 4 (r21, c0, r31, c40, c49) and others from depth 5 (c9 for c48). A first
attempt stopped at depth 2 on a settings mismatch of my making (the recipe's comment was reworded while the
check ran, and a state carries the file's hash): the recipe file must not be edited while the estimate's
paused solve may be resumed or sampled.

## What follows from it

- `suggest/wil50-tune.spec` and `suggest/wil50_d2r10_{g150x,xl2,g150xl2}.par`: a tuning round on a solved sample of
  the paused frontier of `wil50-est-261124` (not on cognac, and it needs that estimate to have run): base, the
  two candidates and both together, as `tho40-tune-261119` was. Its sample is a placeholder until the estimate
  reads.
- `blas/phase0_blas.patch`: `infra/phase0`'s benchmark with BLAS variants (`--blas`); on cognac since 744952ac
  (and 1767cb10, the launcher under macOS's bash), launched Nov 24 as run 20261124-1425-eu-north-1.

## The tables (`analyze.py`)

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
| g150xl2 | 137 | 8,171 (5,392 + 2,779) | 2,224 | 0.648 | 0.737 | 47 nodes: 0.842 |
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
| g150x | 28 | 2,236 (1,306 + 930) | 778 | 0.988 | 1.011 | 28 nodes: 0.988 |
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
| g150x | 186 | 10,499 (6,407 + 4,092) | 3,418 | 0.971 | 0.927 | 140 nodes: 0.884 |
| g150xl2 | 48 | 3,917 (2,612 + 1,305) | 1,309 | 0.362 | 0.355 | 1 nodes: 0.934 |
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
| g150x | 186 | 10,419 (6,897 + 3,522) | 3,240 | 0.929 | 0.880 | 186 nodes: 0.929 |
| xl2 | 186 | 11,220 (6,850 + 4,370) | 3,594 | 1.000 | 0.976 | 186 nodes: 1.000 |

