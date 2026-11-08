# The Knuth estimator (`bnbdrv -e`): validation and the wil50 run (Nov 7, 2026)

`bnbdrv -e` (`latest/estimate.C`) expands the tree exactly to `start_dive_level`,
then makes random dives: at each level it bounds one node, branches, and follows
one child chosen with probability p, weighting what it sees by 1/p. Parameters are
in `./estimate.par` (format in RUNNING.md). Machine for everything below: a
2-core Xeon @ 2.1 GHz (AMX), one process per core, OpenBLAS (`scipy-openblas32`;
the reference LAPACK is about 7x slower per DNN node), known optimum as `-u`.

## How children are sampled

- **Default (QPB score):** p ∝ ((opt − predicted bound)/(opt − root bound))^e, the
  predicted bound being the QPB lookahead bound from `branch()`.
- **`QAP_EST_XBAR=a` (added Nov 7):** at a node with DNN marginals (rules 8 and 9),
  child (i,j) gets p = (1 − mix) x_ij^a / Σ + mix/kids, mix = `QAP_EST_XBAR_MIX`
  (default 0.1). The uniform part keeps every p positive, so the estimate stays
  unbiased. Nodes without marginals fall back to the QPB score.
- `QAP_EST_SEED` seeds `rand()` (unseeded runs repeat the same dives);
  `QAP_EST_HIST` names the per-dive history file (depth, node estimate, time estimate).

## Results

QPB only (`PARAM/<inst>.par`), QPB sampler, 2,000 dives from level 2:

| instance | actual nodes | estimate | actual CPU-s | estimate |
|---|---|---|---|---|
| nug12 | 592 | 597 ± 12 | 0.10 | 0.11 |
| nug15 | 4,786 | 4,610 ± 189 | 0.55 | 0.56 |
| nug18 | 240,357 | 219,331 ± 18,172 | 25.3 | 27.8 |
| nug20 | 656,094 | 586,210 ± 32,378 | 118.4 | 116.6 |

DNN recipe (`nug_deep_xbar_all.par -x 1:7:1:200:1`, the Nov 4 environment of
`../wil50_feasibility/tuning/README.md`), 200 dives from level 2:

| instance | sampler | nodes by level (2/3/4/5/6) | nodes | CPU-s |
|---|---|---|---|---|
| nug27 | **actual** | 83 / 189 / 48 / - / - | **331** | **271** |
| nug27 | QPB score | 81 / 167 / 332 / - / - | 593 ± 246 | 347 ± 117 |
| nug27 | x̄, a=1, mix 0.1 | 81 / 178 / 0 / - / - | 269 ± 85 | 206 ± 24 |
| nug27 | x̄, a=2, mix 0.1 | 82 / 133 / 2 / - / - | 227 ± 74 | 294 ± 49 |
| nug30 | **actual** | 225 / 1,880 / 1,242 / 338 / 125 | **3,820** | **3,773** |
| nug30 | QPB score | 225 / 1,764 / 756 / 0 / 0 | 2,755 ± 796 | 2,735 ± 481 |
| nug30 | **x̄, a=1, mix 0.1** | 225 / 1,217 / 1,924 / 227 / 63 | **3,665 ± 1,037** | **3,764 ± 635** |
| nug30 | x̄, a=1, mix 0.5 | 223 / 2,263 / 1,434 / 706 / 0 | 4,636 ± 1,055 | 3,609 ± 516 |

(The nug27 and nug30 full solves reproduce the M2 Pro trees exactly: 331 / 314 and
3,820 / 3,678 nodes / DNN fathoms; CPU here is about 2x the M2 Pro.)

- With the DNN recipe the QPB sampler is close to uniform (the QPB bounds are far
  from the cutoff) and rarely follows the children that go deep: on nug30 one dive
  in 200 reached level 4 and none reached 5-6, so it missed 463 nodes outright.
- The marginal sampler at a=1, mix 0.1 reached depth ≥ 4 in 13 dives of 200 and
  was within 4% on nodes and 0.3% on CPU on nug30. On nug27 its level-4 nodes sit
  under low-marginal children and it found none; the uniform share is the guard.
- The reported SE comes from the same heavy-tailed dives and understates the error
  when the deep part is undersampled (nug27, x̄ rows). Read a single run's SE as a
  lower bound; compare independent seeds.

## wil50

Production environment (`../wil50_feasibility/fullrun/settings.sh`), rule 9 through
depth 6 (`nug_deep_xbar9_all.par`), DNN at depths 1-7, dives from level 1 so each
process repeats only the root (`estimate.par`: `1 10000 1 5 1.5 1000`):

```sh
. experiments/wil50_feasibility/fullrun/settings.sh
cd <workdir with PARAM, SYM linked from latest/ and the estimate.par above>
QAP_EST_SEED=<distinct per process> QAP_EST_XBAR=1 QAP_EST_XBAR_MIX=0.1 \
QAP_EST_HIST=dives_<seed>.txt \
  latest/bnbdrv -e -p PARAM/nug_deep_xbar9_all.par -x 1:7:1:300:1 \
  -u $FULLRUN_UB -s $FULLRUN_SYM $FULLRUN_INSTANCE
```

Each process takes about 1.3 GB. Combine the processes' history files: the estimate
is the mean of column 2 (nodes) and column 3 (CPU-s) over all dives, and its SE the
standard deviation over sqrt(dives). The root has 15 depth-1 orbits.

QPB-sampler pilot (Nov 7, 4 dives, 9-11 min each on this Xeon): depths 3, 2, 2, 3;
node estimates 22,090 / 501 / 448 / 32,706; CPU estimates 4.4e5 / 1.3e5 / 3.1e5 /
3.9e5 s. All four dives landed on children that fathomed at depth 2-3, which is why
the marginal sampler was built; these dives are not used in the estimate.

### First marginal-sampler estimate (Nov 7, 15 dives)

Seeds 31 and 32, `QAP_EST_XBAR=1 QAP_EST_XBAR_MIX=0.1`, two processes for two hours
on the 2-core Xeon (8 + 7 dives). Per dive: depth, node estimate, CPU-s estimate:

```
4 971 224719          5 12536 3.64048e+06
8 1902617 1.94622e+07 3 1103 306837
7 3875416 1.24056e+08 6 58652 1.15428e+07
4 96605 1.60073e+07   4 4852 116864
2 697 250155          6 6592088 1.08021e+08
6 79952 2.45838e+07   6 561589 4.68152e+06
4 3203 795689         2 225 243826
4 258243 2.07253e+06
```

- **9.0e5 ± 4.7e5 DNN-era nodes, 2.1e7 ± 1.0e7 CPU-s on this Xeon** (about 6 hours
  on 1,000 such cores). The SE is a floor: two dives carry 73% of the CPU sum.
- Depths reached: 2 (2), 3 (1), 4 (5), 5 (1), 6 (4), 7 (1), 8 (1); the QPB sampler's
  four pilot dives all stopped at depth 2-3.
- Nodes by level, each process's estimate at its fifth dive: level 2 174 / 83,
  level 3 181 / 1,845, then 1e4-1e6 per level at depths 4-8. The rule-9 tree with
  the DNN from depth 1 keeps a depth-3 frontier of hundreds to a few thousand nodes,
  against the 29,580 orbits `../wil50_feasibility/roadmap/tree_size.md` assumed
  (its section 3 names this as the main upside risk). That is where most of the gap
  to the roadmap's 1e10-6e10 CPU-s comes from, and it is the cheapest thing to
  check directly: stage A of the sampling plan (depths 0-2 exact) gives the real
  depth-2 and depth-3 counts.

## To do: log per-level estimates per dive

The dive history file records only each dive's depth and its total node and CPU
estimates, so per-level variances can't be computed afterwards. bnbdrv does print
the per-level table (nodes, %fathom, time by level), but only cumulatively and only
every `print_everything_every` dives: the Nov 7 wil50 run printed it once per
process, after dive 5. Per-level estimates and their variances are valuable (where
the tree is, how sure we are of each level), so the next runs should log them:

- Simplest change: write each dive's per-level node and time contributions (its
  `path_stats` scaled by the dive's weight) as one line per dive to the history file,
  or to a second file next to it.
- Workaround with no code change: set `print_everything_every` to 1. The table is
  then printed after every dive, and dive k's contribution to a level is
  k·mean_k − (k−1)·mean_(k−1).
