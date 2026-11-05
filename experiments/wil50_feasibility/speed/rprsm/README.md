# The rPRSM pieces we never built (Nov 5 2026)

Item 1 of `../../roadmap/admm/alternatives.md` (row 2 of
`../../roadmap/admm/README.md`): the parts of Graham, Hu, Im, Li and
Wolkowicz's restricted Peaceman-Rachford method (rPRSM; INFORMS J. Comput. 34
(2022), arXiv 2006.01529) that `latest/dnn.C` did not have. [m] = measured
here, [e] = estimate. All runs on core 2 of the 4-core 2.1 GHz Xeon (no AMX),
one process at a time, one BLAS thread, with three other agents' jobs on the
other cores; iteration counts are deterministic, milliseconds carry that
noise.

## Verdict

**Kill: at the line cold, below it warm, and no CPU to speak of.**

- The **trace-constrained R-step** (`DNN_RPRSM=1`) does nothing on wil50:
  -0.9% iterations on the 8-node subset's fathoms, +10 at iteration 200 on
  its misses [m]; nothing on top of the DS piece [m]. With our
  Gilmore-Lawler certificate its multiplier `tau` lies along an exactly flat
  direction of the bound (below), so it cannot strengthen the certificate the
  way it strengthens Graham et al.'s box bound.
- The **doubly stochastic piece** (x-row and arrow of the Y-step projected
  onto the DS matrices, `DNN_RPRSM=10`) on all32, cold: the same 18 fathoms,
  **-10.2% iterations to fathom in total, -9.1% in the median (220 -> 200),
  -11.3% mean per node** [m]; every one of the 14 misses higher at iteration
  200 (median +33.5, range +23.7 to +66.7) [m]; but **-6.4% CPU on the
  fathoms and -4.1% over all 32 nodes** [m], since the iterations it saves
  are the early ones and the early iterations cost 15-20% more (equal after
  200) [m]. nug30 (m = 28): -11% iterations on 5 fathoms at equal CPU (+10%
  per iteration), +4 to +8 at 200 on 3 misses, no change of outcome [m].
  **Warm-started children**, as in the search: -6.5% iterations and -5.2% CPU
  on 3 fathoms, +1 to +4 at 200 on 3 misses, +3% CPU on them [m].
- The **restricted dual update** (`DNN_RPRSM=4`, the paper's actual
  "restricted" step, see below) is neutral alone on nug20 and cancels the DS
  gain when combined [m].

Against the kill rule (under 10% fewer iterations to fathom, a lower lb@200,
or any changed fathom): cold all32 passes by the total (10.2%) and fails by
the median (9.1%); lb@200 and the fathom set pass; warm children fail the
iteration test (6.5%). The CPU gain is 4-6% cold [m] and about 5% warm [m]
on fathoms, under the noise of the search's other decisions [e].

## What the paper does (read from the arXiv PDF)

The arXiv version (27 pages) was reachable this time by `curl` through the
proxy (WebFetch was refused for arxiv.org), so the memo's "unverified" items
were checked against Algorithm 3.1:

- **R-step**: `R = P_R(V'(Y + Z/beta)V)` over `R = {R psd, tr R = n+1}`
  (eq. 2.6, Section 3.2.1): the eigenvalues are projected onto the simplex
  `{lambda >= 0, sum = n+1}`, i.e. `R = (W - tau I)_+`. As the memo said.
- **The "restricted" dual update** (the memo's unverified detail): both dual
  steps are projected onto `Z_0 = {Z : Z_ii = 0, Z_0i = Z_i0 = 0, i >= 1}`,
  `Z^{k+1/2} = Z^k + gamma beta P_Z0(Y^k - V R^{k+1} V')`, the same after the
  Y-step, from a `Z^0` in `Z_A = {Z_ii = -(L)_ii, Z_0i = -(L)_0i}` (Theorem
  2.11: an optimal multiplier lies in `Z_A`). The arrow and diagonal of
  `L + Z` stay zero; in the Y-step those entries are `clip(VRV')`. That is
  the "restricted" of the name, a third piece we did not have. Their
  `gamma` is in (0, 1) (0.9 in the runs), `beta = n/3`, and the objective is
  rescaled to `P_V L P_V + sigma I` (eq. 3.10; not built: with our
  Gilmore-Lawler set the `sigma I` part is a no-op, see below).
- **The doubly stochastic projection is not in their algorithm.** Its Y-set
  is the box with the gangster zeros (eq. 2.7, Section 3.2.2). The DS
  constraints appear only in Remark 2.9 as redundant constraints one could
  add: `Mat(diag Y)` doubly stochastic, and each row/column a multiple of a
  vectorized DS matrix. Our "DS projection in the Y-step" is that remark,
  untested by the authors.
- Their bound is `g(Z) = min_{Y in box} <L + Z, Y> - (n+1) lambda_max(V'ZV)`
  at `Z^k`: our Proposition 2.2 with the box for `P`.

## What was built: `DNN_RPRSM` (bits, read in `dnn.C`; unset = off)

Not merged into the branch's `latest/dnn.C`: the idea is killed, so the
solver change is kept as `dnn_rprsm.patch` (the diff against a955e95; `git
apply` it there to rebuild the switch).

| bit | piece | how |
|---|---|---|
| 1 | trace-constrained R-step | `tau` from `trace(W)` and the tracked block's Ritz values (all K at a refresh); `R = W - tau I - sum_{theta<tau} (theta - tau) x x'`; the block holds the eigenvalues below `tau` plus the guard; `tau` capped at the block's top Ritz value if the block were too small (never: `capped 0` in every log) |
| 2 | x-row DS projection in the Y-step | the Y-step target's x-row `t`, projected onto `{x >= 0, row and column sums 1, eliminated x_ij = 0}` by exact block-coordinate ascent on the dual (`x = max(t - u_i - v_j, 0)`, simplex projections, warm duals; 1-3 sweeps per iteration), written to row 0 and column 0 |
| 8 | with 2: the arrow | also `Y_(ij),(ij) = x_ij` (Lemma 2.7: diag Y = first column), projecting `(2 t_0 + t_diag)/3` |
| 4 | restricted dual update | Z's arrow and diagonal set to `-L/scale` at the start (also over a warm start), reset after both dual steps |

The solver reads the same variable, so `DNN_RPRSM=10 bnbdrv -x ...` runs it in
the search. Prototype mode (`dnntest -P`) ignores it, since it reproduces
`dnn.py`'s iteration. The DS projection costs well under 1 ms per iteration at
m = 47 (2,209 entries).

## Why the certificate stays valid

Proposition 2.2 of `docs/dnn_bound/dnn_bound.tex` needs a symmetric multiplier
`Z`, a rigorous `lambda >= lambda_max(V'ZV)`, and a polyhedral set containing
the lifted permutations. Its proof uses `tr R = m+1` exactly, so
`<V'ZV, R> <= (m+1) lambda` holds for `lambda` of either sign.

- Bits 2, 4, 8 change only the iterates. The certificate multiplier is still
  `Z' = Z + beta (Y - VRV')` from the R-step, symmetric (the DS projection
  writes row 0 and column 0 from the same `x`, the restricted reset writes
  `Z_0i` and `Z_i0` from the same `L`). The Gilmore-Lawler part and Rump's
  test are untouched.
- Bit 1 changes the multiplier's spectrum: `V'Z'V = beta (W - R) =
  beta (tau I + sum_{theta<tau} (theta - tau) x x')`, so `V'Z'V - beta tau I`
  is negative semidefinite by construction, whatever the accuracy of the fp32
  eigenvectors (the default's negative-form argument, shifted by `tau`).
  Rump's shifts start at `delta = beta tau + rel ||S||_F` (`certify(...,
  lam0)`). Rump's test makes whatever `lambda` it returns rigorous, so `lam0`
  can only make the test fail (no certificate at that check), never make it
  wrong. With `tau < 0` the bound gains `(m+1) scale beta |tau|`, which the
  proposition allows. The double-precision refinement re-certifies its
  re-projected multiplier with `lam0 = 0`, and its trigger compares only the
  penalty beyond `(m+1) scale lam0`.

Two remarks from working it through:

- With the Gilmore-Lawler set, `Z -> Z + s I` is an exactly flat direction of
  the certificate: `M_00` and every `c_ij` gain `scale s`, so GL gains
  `(m+1) scale s`, and so does the penalty `(m+1) scale lambda`. `tau` moves
  along that direction (on nug20 `sdp_penalty` ranges from -61 to +166 over
  runs and iterations, while the bounds track the default's). So bit 1 cannot strengthen *our* bound as it strengthens Graham et
  al.'s box bound; it can only change the iterates. The same holds for the
  `sigma I` of their eq. 3.10.
- The DS projection leaves the row-plus-column part of Z's arrow
  unobserved by the Y-step (the projection is invariant to adding
  `a_i + b_j` to its input), so that part drifts. The GL certificate does
  not care, the **box certificate does**: `dnntest -b` on nug12 gives 567.35
  plain, -0.27 with `DNN_RPRSM=10`, 566.45 with `DNN_RPRSM=14` (the restricted
  update pins the arrow). Valid in every case, but the DS bits are for the GL
  certificate (the default) only. Margins stay at 4.4-4.5e-6 on wil50,
  against 4.9e-6 by default.

## nug20 root (`nug20_curves.txt`, `logs/nug20_*.log`)

Cold start, PRSM (1, 1), degree 3, `-L 2`, `-r 40`, one CholQR pass, no
cutoff, 400 iterations. First check at or above a level [m]:

| level | base | 1 trace | 2 x-row | 4 restr. | 10 x-row+arrow | 11 = 1+10 | 14 = 4+10 | 15 = all |
|---|---|---|---|---|---|---|---|---|
| 2500.0 | 130 | 120 | 120 | 130 | 110 | 110 | 130 | 130 |
| 2504.0 | 180 | 170 | 180 | 180 | 160 | 160 | 180 | 170 |
| 2505.0 | 220 | 210 | 220 | 220 | 200 | 200 | 220 | 220 |
| 2505.9 | 340 | 320 | 340 | 340 | 310 | 310 | 340 | 340 |
| ms/iteration | 4.9 | 5.2 | 4.8 | 4.7 | 5.1 | 5.3 | 4.7 | 5.7 |

The arrow is what matters in the DS piece (bit 2 alone is nearly neutral
after iteration 100); the restricted update is neutral alone and cancels the
DS gain when combined; the trace step is worth 0-20 iterations alone and
nothing on top of 10. A first version of bit 1 sized the tracked block by
all negative eigenvalues even when `tau < 0`: slightly stronger (2505.9 at
300 with 11) but 75% more per iteration on nug20; the block now holds the
eigenvalues below `tau` only.

## wil50, all32 (`run_all32.sh`, `summarize.py`, `runs/`)

The 32 depth-3 nodes of `../../tuning/giveup/all32.txt`, production-like
flags: degree 3, `-L 2`, `-g 1.0 -G 1.0`, `-r 40`, `-H 300`, `-i 2500 -u
48816`, `DNN_CHOLQR_PASSES=1`, one node at a time on core 2 [m]:

| config | fathomed | ms/it | median it | sum it (18 fathoms) | CPU fathoms | CPU misses | lb@200 on the 14 misses |
|---|---|---|---|---|---|---|---|
| base (unset) | 18 | 373 | 220 | 4,230 | 1,462 s | 1,867 s | - |
| `DNN_RPRSM=10` | 18 (same) | 389 | 200 | 3,800 (-10.2%) | 1,369 s (-6.4%) | 1,823 s (-2.4%) | median +33.5, min +23.7, max +66.7 |

Every fathom comes 10-40 iterations earlier (ratios 0.80-0.94, mean 0.887);
the give-up fires earlier on 5 of the 14 misses (e.g. node 16: 370 instead
of 470 iterations). Per-iteration cost by phase on the misses (cpu column of
the logs): iterations 1-100 342 -> 396 ms, 100-200 436 -> 450-470, 200-300
439 -> 438-450. The variant's iterate reaches the final rank structure
sooner (`npos` 1,480 instead of 1,680 at iteration 50, so a block of ~640
instead of ~440 tracked vectors), which is what makes its early iterations
dearer.

Development subset (nodes 2 9 14 18 19 21 24 32; `logs/subset8.out`), sum of
iterations on the 5 fathoms and lb@200 on the 3 misses against base [m]:

| config | sum it | lb@200, median |
|---|---|---|
| base | 1,150 | - |
| 1 (trace) | 1,140 (-0.9%) | +10.5 |
| 11 (trace + DS + arrow) | 1,060 (-7.8%) | +34.2 |
| 10 (DS + arrow) | 1,040 (-9.6%) | +34.0 |

## nug30, m = 28 (`run_nug30.sh`, `runs30/`, `logs/nug30.out`)

Eight depth-2 nodes (the eigensolver README's three and five from
`experiments/dnn_bound/predictor/nug30_depth2.csv`), the record setting's DNN
switches (PRSM (1, 1), degree 3, `-L 2`, `-B 0.21`, `-r 120`, one CholQR pass,
`-H 200`, `-u 6124`), cold [m]:

| | base | 10 | 11 |
|---|---|---|---|
| iterations on the 5 fathoms | 1,180 | 1,050 (-11.0%) | 1,050 |
| CPU on the 5 fathoms | 26.50 s | 26.39 s (-0.4%) | 26.08 s |
| ms/iteration | 21-24 | 23.6-25.7 | 22.8-25.8 |
| lb@200 on the 3 misses | - | +6.1, +7.6, +4.4 | +6.2, +7.8, +4.5 |

No fathom changes. At m = 28 the per-iteration premium (+10%) does not fade
by iteration 200 (`npos` stays 10-25 lower), so the iteration saving buys no
CPU.

## Warm-started children (`run_children.sh`, `runs_child/`, `logs/children_summary.txt`)

all32 is cold, but in the search most DNN nodes are children started from
their parent's Z and conditioned Y (`QAP_DNN_WARM_Y=2`), and the variant's
gain is in the early iterations. Two all32-type parents (depth 3, nodes 13 and
10 of `../../stronger_bounds/depth4/nodes.txt`, given up at 300) and six of
their children, `dnntest -K ... -Q b`, same flags; parent and children run
with the same switch [m]:

| | base | `DNN_RPRSM=10` |
|---|---|---|
| parents (cold misses), lb@200 | 48,503.1 / 48,501.5 | +29.1 / +35.5 |
| 3 fathomed children: iterations | 340, 340, 90 = 770 | 300, 330, 90 = 720 (-6.5%) |
| their CPU | 271 s | 257 s (-5.2%) |
| 3 missed children: lb@100 | - | +4.6, +1.3, +5.1 |
| their lb@200 | - | +2.9, +1.3, +3.8 |
| their CPU | 315 s | 325 s (+3%) |

No outcome changes. The children's base runs reproduce `../../warm_y/`'s
conditioned-Y counts (340, 300, 310, 90 for node 13's children). Warm, the
advantage at iteration 200 falls from +24 to +67 (cold misses) to +1 to +4,
and the iteration saving on fathoms to 6.5%: a warm start already gives the
iterate the doubly stochastic x-row and the arrow structure that the
projection enforces.

## Validity checks run

- `make dnncheck`: passes with the switch unset (before and after the last
  change), `DNN_RPRSM=10` and `DNN_RPRSM=15` (`logs/dnncheck_rprsm*.log`; part A max(lb - optimum) over
  2,592 node bounds -1.1e-13 and -1.8e-13). Before prototype mode ignored the
  switch, part B (agreement with `dnn.py`) failed with it set, as it must when
  the iteration differs; parts A, C, D, E passed.
- `python3 tests/brute_force.py --jobs 1 --solver-args "-x 0:20"` with
  `DNN_RPRSM=15` and `DNN_RPRSM=10`: 4,536 of 4,536 each
  (`logs/brute_force_rprsm*.log`).
- `validity_rprsm.py` (adapted from
  `roadmap/certificate/scripts/validity_bruteforce.py`): 12 random n = 8 QAPs
  (symmetric/asymmetric flows, some with a linear term), root, `-x 0:3`, and
  `-x 0:3` with three eliminated assignments; bits 0, 1, 2, 4, 10, 11, 15;
  three eigensolvers (subspace, ssyevd, dsyevd); a certificate at every
  iteration for 300 iterations: 756 runs, 215,918 certified values, every
  best value (compared at 17 digits from the RESULT line) below the
  enumerated optimum, worst -6.8e-10 (`logs/validity_rprsm.log`).
- `validity_bruteforce.py` itself (the `DNN_CERTX=47` variants at every
  check) with `DNN_RPRSM=11` and `=15`: PASS (`logs/validity_bruteforce_*.log`;
  note that its comparison is at the 4 decimals the variants print).

## Default path unchanged

With `DNN_RPRSM` unset, logs agree line for line (all bounds, iteration counts,
residuals, `npos`; timings stripped, `same.py`) between the binary built
before the change and after it (nug20 also with the final binary): nug20 root
400 iterations, wil50 nodes 19 and
11 (fathomed at 130 and 290), nug30 nodes `4:10,6:8` (fathomed at 330) and
`4:29,6:7` (given up at 200). Every changed expression reduces to the old one
for `tau = 0` (`0 - x`, `x - 0`, `theta < 0.0f`), and the `lam0` shift is
added only when non-zero.

## Not built

- The per-column DS structure of Remark 2.9 (each column block `x_ij` times a
  DS matrix) in the Y-step: symmetry couples the m^2 blocks, so it needs an
  inner Dykstra/ADMM over the whole matrix, several Y-step costs per
  iteration; the GL certificate already uses that structure.
- Their `gamma < 1` (0.9/0.9 was measured worse here before), `beta = n/3` and
  the rescaling (3.10).
