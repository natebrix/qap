# The eigensolver shortlist, measured (Oct 30 2026, Apple M1, Accelerate)

Question (experiments/wil50_feasibility/NOTES.md, section 6): the per-node
ideas for the DNN bound "compound to 3 to 5x at best". Which of them are
real at wil50's size (m = 47, K = 2117)?

## 1. The standalone studies did not run the solver's eigensolver

Every wil50, tho40 and nug30-predictor sample on this branch ran
`dnntest -f`, which is `-e 1`: a full single-precision `ssyevd` every
iteration. The solver (`bnbdrv -x`) uses `DNNParams` defaults, i.e. the
subspace tracker (`DNN_EIG_SUBSPACE`). On wil50 depth-3 node 2
(`18:43,11:29,46:44`), alone on this machine:

| eigensolver | ms per iteration | eigensolver share |
|---|---|---|
| `-f`, full ssyevd each iteration | 570 | 90% |
| default subspace tracker (degree 2, refresh 20) | 164 | 66% |

So the 250 CPU-s per node behind the wil50 estimate (1e10 CPU-s) was
measured in a mode 3.5x slower than the solver's. The estimate's CPU
figure, not its node count, is too high by about that factor.

## 2. But the default tracker is inexact at this size

Run to 2,500 iterations (primal and dual agree within 2 in every case),
shortfall below the cutoff 48,814.05:

| setting | at 800 | at 2,500 | ms/iter (alone) |
|---|---|---|---|
| exact (`-f`) | 223.9 | 202.3 | 570 |
| default: Chebyshev degree 2, refresh every 20 | 234.8 | 221.0 | 164 |
| degree 4, refresh 20 (`-d 4`) | 224.0 | **202.3** | 182 |
| degree 3, refresh 10 | 226.4 | 210.6 | |
| degree 2, refresh 5 | 250.9 | 233.2 | |

The degree-2 filter leaves the ADMM with a noise floor about 20 above the
exact value at m = 47 (the cpp report found no such effect at m <= 28). A
degree-4 filter reproduces the exact trajectory to the decimal for 11% more
per iteration. Refreshing more often does not help and refresh-every-5 is
worse, so inertia-gated refresh is at most a small time saver (the 120
scheduled refreshes are 52 s of 443). Logs: `node2_*.log`.

Where a degree-4 iteration goes at m = 47, late in a run: filter and
Rayleigh-Ritz 70 ms, form R (dsyrk, rank nneg ~570) 20 ms, scheduled
refresh 21 ms, forward + expand transforms 46 ms, Y/Z-step 13 ms,
certificates 6 ms. Side-switching does not apply: the negative side (~570
of 2117) is already the smaller one throughout.

## 3. Peaceman-Rachford: 1.5x fewer iterations, same bound

`DNNParams::gamma1` (dnntest `-G`, solver `QAP_DNN_GAMMA1`) adds a first
dual step `Z += gamma1 beta (Y - VRV')` between the R-step and the Y-step
(Graham, Hu, Im, Li, Wolkowicz 2022). The certificate still uses the
R-step multiplier `Z + beta (Y - VRV')` computed before that step, so
`V'Z'V = beta (W - R)` stays negative semidefinite by construction and
nothing in the proof changes. gamma1 + gamma must stay below 2:
(0.5, 1.618) diverges.

| node | ADMM (0, 1.618) | PRSM (1.0, 1.0) | PRSM (0.9, 0.9) |
|---|---|---|---|
| nug30 `4:10,6:8`, fathom | 500 iterations | 310-320 | 350 |
| nug30 `4:29,6:7`, miss | 6095.906 | 6095.904 | |
| nug30 `4:7,6:14`, miss | 6109.098 | 6109.085 | 6109.092 |
| wil50 node 11 `0:48,31:2,38:32`, fathom (degree 4) | 440 iterations, 80 s | 290, 50 s | |
| wil50 node 2, iterations to within 5 of the converged value (degree 4) | ~1,800 | ~1,000 | |

`nug30_prsm_iterations.txt`, `node2_cheb4_prsm.log`.

## 4. Block growth instead of forced refreshes

PRSM moves the spectrum of W more per iteration, and the tracker's
"block too small" rule (negative count grown into the guard) then forced
364 unscheduled full refreshes in 2,500 iterations (324 ms/iter). The
tracker now appends `guard` fresh pseudo-random vectors (deterministic
LCG) to the block instead, and the next filter, CholQR2 and Rayleigh-Ritz
absorb them; a full refresh runs only when the block cannot grow. With
PRSM on node 2: 120 scheduled refreshes, 96 growths, 177 ms/iter, same
converged value (202.1). This is on by default (it changes only the
inexact ADMM step, never a certificate); iteration counts on the nodes
tested are unchanged or one step lower.

## 5. In the solver: nug25, 8 workers, `nug_deep.par -x 2:5:1:300:1`

| | default | `QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4` |
|---|---|---|
| nodes | 728 | 696 |
| DNN nodes / fathomed | 721 / 688 | 689 / 658 |
| DNN CPU | 722 s | 509 s |
| CPU | 727 s | 514 s |
| wall | 104 s | 73 s |

`nug25_default.log`, `nug25_cheb4_prsm.log`. Both proved 3744. nug20
`-x 2:3:1:300` serial: DNN CPU 24.8 s -> 19.0 s, identical tree.

## 6. What this means for wil50

Per fathomable depth-3/4 node, against the `-f` mode the estimate used:
node 11 takes 50 s instead of 251 s (5.0x). Per clear miss stopped at 300
iterations: 53 s instead of 171 s (3.2x), and PRSM's bound at 300
iterations is 60 better (263 against 322 short), which sharpens the
give-up rule and any tree estimate built on it. So the "3 to 5x" is real
and is now in hand, but not from the listed items: 3.5x was already in the
solver's default mode (hidden by the `-f` flag in the standalone studies),
PRSM gives 1.5x on top, and degree 4 is what makes the fast mode exact at
m = 47. Inertia gating, side switching and transform fusion are each worth
a few percent at most now. None of this changes the fathoming depth, so the
wil50 verdict (a stronger relaxation is the lever) stands; the CPU side of
it is about 2e9 M1-CPU-s on the 4e7-node count rather than 1e10.

## 7. Validation

`tests/dnn_test.py` with a new part-A variant `prsm-cheb4` (`-d 4 -g 1.0
-G 1.0`): 1,728 node bounds, none above an optimum, no fathom at the
optimum. `tests/brute_force.py --solver-args "-x 0:20"` with the PRSM and
degree-4 switches set: 4,536 of 4,536. `make test` with default settings:
brute force, resume and parallel pass; the regression check fails only on
the known Apple-silicon root-bound difference (nug15 984.79 against 986.55,
as RUNNING.md records), with node counts within tolerance and no DNN bound
involved.

## Not done

- Making PRSM (1.0, 1.0) and degree 4 the solver defaults. That changes
  the tree (fewer DNN misses on nug25) and needs the standard A/B on the
  recorded nug27/nug28/nug30 runs; the warm start interacts with gamma1
  and was not retuned.
- gamma1/gamma other than (1.0, 1.0) and (0.9, 0.9); Graham et al.'s
  restricted variant with a different Y-step.
- The forward/expand transforms (46 ms, memory-bound) and the rank-nneg
  dsyrk (20 ms) are the remaining per-iteration targets, each worth
  perhaps 10%.
