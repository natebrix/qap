# Shipping the parent's eigenbasis to a warm child (Nov 6, cloud)

Question: the production candidate (`bestcore` in `../../pipeline_ab/run_pair.sh`)
runs its warm children at a low ADMM penalty, so they fathom in 20-40
iterations. Each warm solve starts with a full eigendecomposition of W:
about 0.7 s, or 0.5 s with `DNN_CORE_REFRESH_PARTIAL=1`, 6-9% of such a child
([m] in `../../tail_probe/README.md`, last paragraph). Does starting the child's
subspace tracker from the parent's tracked basis save that time, without
costing iterations or fathoms?

**[m]** = measured here (log named); **[e]** = estimate. Host: the 4-vCPU
2.1 GHz Xeon with AMX, one BLAS thread, every process pinned with `taskset`.
Code: commit `eb4828b` (`DNN_SHIP_BASIS=1`, off by default), plus the test fix
in section 8. **Killed, so the code is not on `cognac`:** it lives on the branch
`research/ship-basis` (`eb4828b`, test fix `64b73cb`); this directory's scripts
need a build of that branch (`git worktree add ../qap-basis research/ship-basis`).

## Verdict

- **Killed: the switch saves 0.3-0.8% of the warm children's CPU [m],
  against a kill line of 3%.** The start-up saving per warm solve is
  263-274 ms [m]: the first R-step drops from 530-565 ms (the full refresh)
  to 268-291 ms (the seed plus one filter step). That is 0.6% of the warm
  set's CPU [m], 3-4% of a 20-23-iteration child [m] (2.5% over the 8 short
  children run on one core), and 1.4-1.5% of the production DNN CPU [e].
  Even a free first R-step would give 1.7% [e].
- **No fathom lost, no fathom gained, the same 3,461 iterations [m].** One
  child fathoms one iteration later and one one iteration earlier
  (marginal crossings, rounding). The misses' lb agree with the switch-off
  run to -0.02..+0.03 at iteration 200 [m].
- **Valid by construction.** The certificate never reads the tracker
  (section 2). The default path is bit-identical to the old binary, old
  state files load, and a bad basis falls back to the full refresh with the
  switch-off run's trace (section 3).
- **Noise.** The base/base2 replicate differs by 13% between passes (the
  box's load); within a pass a total moves by about 0.5 points [m]. The
  saving is at that noise level, apart from the timed start-up.
- **bf16 basis**: the same iterations and fathoms as fp32 [m], half the
  bytes. **State files**: the basis adds 7.8 MB per wil50 parent at fp32
  (+30% on the 26 MB DNN state) or 3.9 MB at bf16 [e].
- **Gates**: all pass on `eb4828b` (section 8), after a fix to
  `tests/ship_test.py`. As committed, `make shipbasischeck` failed on a
  counting bug in the test, not in the solver.
- No `bestship` arm was added to the pipeline A/B.

## 1. What the switch does

`DNN_SHIP_BASIS=1` changes three places. Without it nothing runs.

- **Save.** `dnn_bound` keeps the tracker's last block in `DNNWarm`: the kb
  Ritz vectors of the last R-step (fp32, K x kb, K = (m-1)^2 + 1), and the
  tracked top eigenvector `vtop`. It saves nothing if the block is the whole
  space (kb = K) or the call failed.
- **Map.** `dnn_warm_child` maps each vector into the child's facial
  coordinates, with the same embedding that maps Z.
- **Seed.** A warm-started call outside a cut phase starts its first R-step
  from the mapped block (`seed_block`) instead of the full refresh.

### The mapping

Notation: V_p (n1_p x K_p) and V_c are the parent's and the child's facial
bases, with orthonormal columns. The parent's R-step works on
W_p = V_p' M V_p, where M = Y + Z/beta is the lifted matrix. Fixing facility
i to location j embeds the child's lifted space in the parent's by P: child
index 0 goes to parent 0 and (ij), and every other child index a goes to its
parent index map[a]. The warm start already uses it: Z_c = P' Z P.

- P keeps row and column sums, so the child's face lies inside the parent's:
  P V_c = V_p G, with G = V_p' P V_c.
- Hence V_c' (P' M P) V_c = G' (V_p' M V_p) G. At the start, the child's W
  is G' W_p G.
- So a parent vector x maps to x_c = G' x = V_c' P' V_p x: lift by V_p,
  restrict by P', project by V_c'. With the Householder form of V this costs
  O(m^2) per vector.
- The mapped vectors are not orthonormal: G'G = V_c' P'P V_c =
  I + V_c' e0 e0' V_c, because child index 0 collects two parent indices.
  The seed orthonormalizes them first.

`dnntest -mapcheck m:i:j` (`dnn_basis_map_check`) checks the identity
V_c' Z_c V_c = G' (V_p' Z V_p) G for random symmetric Z, to fp32 rounding.

### The seed and its fallback

`seed_block`, on the first R-step of the warm call:

1. CholQR of the mapped block (`DNN_CHOLQR_PASSES` passes, as the filter
   step uses).
2. Rayleigh-Ritz: H = Y'WY. Then the residual test, before the eigensolve.
   For orthonormal Y, |WY|_F^2 = |H|_F^2 + |WY - YH|_F^2, so
   res = sqrt(|WY|_F^2 - |H|_F^2) / |H|_F, which is the Ritz pairs'
   |WX - X Theta|_F relative to |H|_F. The seed passes if res <= tau
   (`DNN_SHIP_BASIS_TAU`, 0.25).
3. The Ritz values of H. The block is rejected if there is no negative one.
4. Default (`DNN_SHIP_BASIS_STEPS=1`): the R-step goes on with one Chebyshev
   filter step from Y. Its first product W Y is already there (T, from
   step 2). That step's own Rayleigh-Ritz trims or grows the block as after
   any filter step. With `DNN_SHIP_BASIS_STEPS=0`, the seed's Ritz pairs are
   this R-step's directly, with no filter step.

The block is rejected, and the R-step runs the full refresh as without the
switch, if any of these holds: the block is missing, has the wrong K, has
kb >= K, or has a non-finite entry; the Cholesky factorization fails; the
residual exceeds tau; no Ritz value is negative.

**`lam_hi` and `vtop`.** The filter needs an upper end lam_hi of W's
spectrum. After a full refresh it is the exact top eigenvalue. After a seed
that filters (the default), 6 shifted power steps (`DNN_SHIP_BASIS_POW`)
start from the shipped `vtop`, with the shift 1.05 |theta_min|. Without a
usable `vtop`, a fixed pseudo-random start is used. With `STEPS=0` the seeded
R-step does not filter, so it skips them. The next R-step's own power steps
start from the shipped `vtop`. The estimate is padded as the tracker's own
per-iteration estimate is: lam_hi = max(1.1 rho + 1e-3 |rho| + 1e-6,
theta_top + 1e-3 (|theta_top| + 1e-6)). The power-iterated `vtop` is then
carried, as after a refresh, by the tracker's 2 power steps per iteration
(`QAP_DNN_TIGHT_LAM=2`). An underestimate only weakens one filter step, and
the existing test theta_top >= lam_hi still forces a refresh.

**The refresh schedule** is unchanged. `since_refresh` is set so that it
reads 0 after the seeded R-step, as after the refresh it replaces. The next
full refresh comes 40 iterations later in both cases.

### State files

The basis also travels between processes. Version 6 is written byte for byte
as before unless both `DNN_SHIP_BASIS=1` and `QAP_DNN_SHIP` are set. Then the
file is version 7: each warm line gets one more field,
`bas<bits>:<K>:<kb>:<hastop>:<base64>`, or `-` for none. The basis is fp32
(exact, the tracker's own precision) by default, or bf16 with
`QAP_DNN_SHIP_BASIS_BITS=16`. Versions 1-7 are read. A reader older than
version 7 refuses a version 7 file instead of misreading it.

## 2. Validity

The bound does not depend on the tracker. A poor start can cost iterations,
never validity.

- The certificate reads only the multiplier Zc. `cert_lam` transforms Zc in
  double and runs Rump's Cholesky test on delta I - V'ZcV. `cert_poly` is
  the polyhedral (Gilmore-Lawler) part. Neither reads the tracker's state
  (`Xs`, `vtop`, the Ritz values, `lam_hi`, kb): `latest/dnn.C`, lines
  1305-1515.
- R is always formed in the negative form, R = W - sum over theta < 0 of
  theta x x'. So V'ZcV = beta (W - R) = -beta F F', which is negative
  semidefinite for any vectors x and any negative values theta. Rump's test
  then checks it in double on the actual Zc, with its margin.
- The double-precision refinement, when the penalty decides a fathom, uses
  `proj_psd` (a full eigendecomposition), not the tracker, and certifies its
  result again.
- Fathom decisions read certified values only. A seed changes the iterates,
  so it can change when a fathom is found, never whether a found one is
  valid.

## 3. Checks: identity, state files, fallback

All three ran on CPU 0 on the committed binaries (22:23-22:40 UTC). The logs
in `logs/` from before the restart came from an earlier build (its default
tau was 0.05), so they were rerun and overwritten.

**Identity** (`check_identity.sh`, `logs/ident/summary.txt`). Old binary:
`a4ad1db`, the branch point, built in scratch. New: `eb4828b`. Both ran with
the switch unset on nug20's root, on wil50 node 19 (cold, `w19`) and on f9's
parent with its child 37:21 through `-K ... -Q b` (`f9kid`). Each case ran
under the plain and the production environment, the old and new binaries
one after the other on the same core.

- All six cases: every `it` line and every RESULT identical, CPU fields
  stripped [m]. With the switch unset, the default path is bit-identical.
- With the switch on (`f9kid_prod_on.log`), the cold parent's trace is still
  identical to the old binary's [m]. A cold call is never seeded; only its
  saved state grows by the basis.
- `dnntest -mapcheck` on four (m, i, j) cases: `logs/ident/mapcheck.txt`
  for (6, 0, 0), (10, 3, 4), (12, 11, 5) and (16, 7, 15), the largest
  relative error is 2.3e-8 to 3.8e-8, fp32 rounding [m].

**State files** (`state_check.sh`, `logs/state/summary.txt`). The nug20 DNN
search of `tests/ship_test.py`, stopped after one 3-s slice with
`QAP_DNN_SHIP=32` [m]:

| file | written by | version | read by | result |
|---|---|---|---|---|
| v6 | new, switch off | 6 | old solver | resumes (exit 3, next slice) |
| v6 | new, switch off | 6 | new solver, switch on | resumes; 2 of 14 warm nodes seeded (parents bounded in the resumed process) |
| v7 | new, switch on | 7, 2 of 2 warm lines with a basis | coordinator `State` read and write | version, warm lines and nodes unchanged |
| v7 (rewritten by `State`) | | 7 | new solver | resumes; 16 of 17 warm nodes seeded |
| v7h (bf16 basis) | new, switch on, `QAP_DNN_SHIP_BASIS_BITS=16` | 7 | new solver | resumes; 16 of 17 seeded |
| v7 | | 7 | old solver | refuses: "not a version 1-3 state file" (exit 1) |

On nug20 the basis token is 221 kB per warm line at fp32 and 111 kB at bf16,
against 780 kB for the 32-bit DNN state [m].

**Fallback** (`fault_test.sh`, `logs/fault/summary.txt`). nug30's root (300
iterations) as the parent, two children, production environment [m]:

| block | seed | residual | children's trace |
|---|---|---|---|
| intact | used, both children | 0.013, 0.016 | (seeded) |
| pseudo-random vectors (`DNN_SHIP_BASIS_FAULT=1`) | rejected, both | 1.21, 1.22 | identical to the switch-off run |
| coordinates reversed (`=2`, a block in the wrong basis) | rejected, both | 0.74, 0.77 | identical to the switch-off run |

The fallback also runs in a real search. In the nug20 search of
`make shipbasischeck` (straight run, trace in scratch), 109 of 183 warm
children were seeded and 74 rejected (residual 0.28-0.30, above tau 0.25);
the rejected ones refreshed as without the switch [m]. On wil50 every one of
the 24 warm children was seeded, with residuals 0.011-0.035 [m].

## 4. The warm set (`run_all.sh`, `run_warmset.sh`, `summarize.py`, `warm/`)

Six wil50 depth-3 parents, each given up at 300 iterations as in the search,
and 24 of their depth-4 children. Each child starts from its parent's final
Z and conditioned Y (`dnntest -K ... -Q b`, the search's `QAP_DNN_WARM_Y=2`)
and, with the switch, from the parent's tracked basis. Production flags
(degree 3, `-L 2`, PRSM (1, 1), refresh 40, give-up after 300, one CholQR
pass, AMX), `bestcore`'s DNN environment, one BLAS thread, `-v -v` (the
profile lines; printing changes no number). Iteration counts are
deterministic; CPU is not.

Two passes of three configurations ran at the same time, one per core, so
that the CPU seconds of a pass are comparable:

- pass 1 (22:07-22:38 UTC): `base` (switch off, CPU 1), `ship` (the default,
  one filter step after the seed, CPU 2), `ship0` (`DNN_SHIP_BASIS_STEPS=0`,
  the seed's Ritz pairs as they are, CPU 3);
- pass 2 (22:38-23:05): `base2` (CPU 2) and `ship2` (CPU 1), the replicates
  with the cores swapped, and `ship16` (the parent's basis rounded to bf16,
  `dnntest -Wb 16`, as a state file with `QAP_DNN_SHIP_BASIS_BITS=16` would
  carry it; CPU 3).

The parents were identical in every comparison (iterations and lb) [m]: a
cold call is never seeded. Every one of the 24 children was seeded in every
switched configuration, with residuals 0.011-0.035 against tau 0.25 [m].
Of the 24 children, 15 fathom (20-145 iterations) and 9 are misses (300-330
iterations). The misses take 79% of the children's CPU [m].

Each configuration against the base of its own pass [m]
(`warm/summary_pass1.txt`, `warm/summary_pass2.txt`):

| configuration | children's CPU s (base / new) | saved | on the 15 fathoms | iterations | fathoms | moved by one iteration | first R-step ms, mean (base / new) | start-up saved per warm solve |
|---|---|---|---|---|---|---|---|---|
| ship, pass 1 | 1,181.10 / 1,171.82 | 0.79% | 2.92% | 3,461 / 3,461 | 15 / 15, the same set | f6/3 later (52 to 53), f6/2 earlier (24 to 23) | 565 / 291 | 274 ms (0.56% of the children's CPU) |
| ship0, pass 1 | 1,181.10 / 1,171.88 | 0.78% | 2.65% | 3,461 / 3,464 | 15 / 15, the same set | three later: f9/1 (29 to 30), f9/2 (36 to 37), f6/3 (52 to 53) | 565 / 155 | 410 ms (0.83%) |
| ship2, pass 2 | 1,026.41 / 1,023.47 | 0.29% | 1.24% | 3,461 / 3,461 | 15 / 15, the same set | as ship (the same iterations, deterministic) | 530 / 268 | 263 ms (0.61%) |
| ship16, pass 2 | 1,026.41 / 1,023.64 | 0.27% | 1.65% | 3,461 / 3,461 | 15 / 15, the same set | as ship | 530 / 280 | 251 ms (0.59%) |

The certified lb of the 9 misses, new minus base at iterations 50, 100 and
200 (the running best, smallest and largest difference) [m]:

| configuration | lb@50 | lb@100 | lb@200 |
|---|---|---|---|
| ship (and ship2, identical) | -0.66 to +0.06 | -0.09 to +0.11 | -0.02 to +0.03 |
| ship0 | -0.61 to +0.21 | -0.07 to +0.16 | -0.00 to +0.05 |
| ship16 | -0.64 to +0.06 | -0.10 to +0.11 | -0.02 to +0.05 |

The misses end where they ended without the switch, to a few hundredths.
The seed perturbs the start; by iteration 200 the trajectories have met.

**Where the start-up goes.** Without the switch, the first R-step is the
full refresh (ssytrd + sstedc + the partial back-transformation) plus form
R: 565 ms in pass 1 (396-720 per child) [m]. With the switch it is the
seed (85-121 ms: CholQR, the Rayleigh-Ritz product and eigensolve, 6 power
steps), one filter step and form R: 291 ms (245-357) [m]. A regular R-step
(the children's eigen part without the refreshes, per iteration) costs
228 ms in pass 1 and 200 ms in pass 2 [m]. The
net saving, 274 ms, is 0.8 of an iteration (341 ms per iteration for base's
children in pass 1) [m]. `ship0` drops the filter step too and saves 410 ms,
but its first iterate is rougher and three fathoms come one iteration later.
The map of the block in `dnn_warm_child` runs outside `dnn_bound`, so it is
not in these figures: about 690 vectors at O(m^2) each, some 15 ms per
child [e].

`summarize.py --profile` splits the children's CPU (`warm/summary_profile.txt`,
all 24 children) [m]. The switch removes exactly one full refresh per child
(96 to 72 refreshes, 0.52 s each in pass 1) and adds one filter step (3,365
to 3,389). Per filter step the cost is the same: base 183.2 and ship 182.6 ms
in pass 1, 163.6 and 163.5 in pass 2. Forward, expand, the Y/Z-step, the
certificates and the unprofiled rest each move by less than 0.4% of the
children's CPU.

## 5. Noise

- **Between passes (the base/base2 replicate).** The same configuration on
  identical work: base2's children took 0.869x base's CPU, 13.1% less, per
  child 0.66-1.03x (sd 9%); the parents 0.81-0.99x [m]. ship2 against ship:
  0.873x [m]. The box's load changed between the passes: pass 1 shared it
  with the gates and checks on CPU 0, pass 2 with the same-core runs. CPU is
  therefore compared only within a pass.
- **Within a pass.** Per child, ship/base has sd 4.6% and ship2/base2 2.0%
  [m]. The children's totals saved 0.79% and 0.29%. The two passes' figures
  differ by 0.5 points: that is the noise of one such total [m]. The saving
  is the size of that noise, and only the timed start-up saving (0.56% and
  0.61%) stands out of it.
- **Pass 1's extra on the short children does not reproduce.** In pass 1,
  the seeded runs (CPUs 2-3) spent 2-3% less per iteration than base (CPU 1)
  on the short children of f9 and f6, outside the first R-step [m]. With the
  cores swapped in pass 2 it was +0.3% on f9; on one core, +1.9% and -0.2% in
  the two pairings below [m]. In every comparison, iterations 2-10 ran
  0.1-0.2 s faster per child after a seed than after a refresh (24 children:
  0.18 s in pass 1, 0.13 s in pass 2). The slowdown
  is spread evenly over W*Y, CholQR and Rayleigh-Ritz, so it looks like a
  machine effect of the 0.5-s refresh, not fewer operations. Iterations 11-20
  and later gave part or all of it back [m] (`summarize.py --profile`, the
  last column). It is not counted as a saving of the switch.
- **Same core (`run_samecore.sh`, `samecore/`, CPU 0, 22:40-22:59).** f9 and
  f6 in the order base, ship, ship, base [m]:

  | pair | children's CPU (8 short children) | of which the first R-step |
  |---|---|---|
  | sB / sA (ship / base) | 0.957 (4.3% saved) | 2.51% |
  | sB2 / sA2 | 0.977 (2.3% saved) | 2.47% |
  | sA2 / sA (base replicate) | 0.981 | 0 |
  | sB2 / sB (ship replicate) | 1.001 | 0 |

  On short children the switch saves 2.5% in the first R-step [m]. The rest
  was 1.9% cheaper in the first pairing and 0.2% dearer in the second: not a
  saving that reproduces.

## 6. The bf16 basis

`ship16` seeds from the parent's basis rounded to bf16 (8 significant bits),
as a state file with `QAP_DNN_SHIP_BASIS_BITS=16` would carry it [m]:

- the same 3,461 iterations, the same 15 fathoms and the same two moves as
  the fp32 basis;
- residuals 0.0115-0.0351 against 0.0114-0.0350;
- the children's CPU 1.0002x ship2's (same pass);
- the misses' lb within -0.64 to +0.06 at iteration 50 and -0.02 to +0.05 at
  200, as with fp32.

On wil50 bf16 costs nothing measurable and halves the basis in the state
files. On nug20 (`make shipbasischeck`) it moved the search by 20 DNN
iterations of 8,710 [m].

## 7. State-file size per parent

A wil50 depth-3 parent has m = 47, so K = 46^2 + 1 = 2117. Its tracked
block at the end of the solve has kb = 678-694 vectors (the six parents,
`PARENT-BASIS` lines in `warm/ship_*.log`) [m]. A warm line then carries:

| field | fp32 / 32-bit | bf16 / 16-bit |
|---|---|---|
| DNN state (upper triangles of Z and Y, n1 = 2210), `QAP_DNN_SHIP` | 26.1 MB [e] | 13.0 MB [e] |
| basis, kb = 686 plus the top vector, `QAP_DNN_SHIP_BASIS_BITS` | 7.8 MB [e] (+30%) | 3.9 MB [e] (+15% of the 32-bit state) |

The [e] figures are the format's arithmetic: K (kb + 1) values at 4 or 2
bytes, times 4/3 for base64. The nug20 state files check the formula: 221 kB
measured per basis token at kb about 127, 247 kB from the formula with
kb = 127 [m]. Production ships the DNN state at 32 bits (`QAP_DNN_SHIP=32`
in `run_pair.sh`), so an fp32 basis would add about 30% to each parent's
warm line, and a bf16 basis about 15%.

## 8. Gates on the code commit

`gates.sh` (single-process set on CPU 0, 22:07-22:17 UTC, beside the
measurement; `logs/gates/summary.txt`), then `make test` (SET=multi, CPUs
0-3, after the measurement). Every result in the table is [m]. Production
environment = `bestcore`'s DNN
switches: `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2
DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1
DNN_R32_NONCERT=1`.

| gate | result |
|---|---|
| `make test` (brute force, resume, parallel, coord, regression; 23:05-23:07) | PASS (165 s): brute force 4,536 of 4,536, resume 5 of 5 identical, all parallel checks, all coordinator scenarios (the AWS one skipped: no boto3 on this box), regression exact |
| `make dnncheck` (default settings; parts A-E) | PASS (94 s) |
| `make shipcheck` (state v6) | PASS (160 s): the sliced run equals the straight run node for node |
| `make shipbasischeck` (state v7, the basis in the state files) | FAIL as committed, a test bug; PASS after the fix below (207 s): 109 of 183 warm nodes seeded, the sliced run with the fp32 basis equals the straight run node for node |
| `make continuationcheck` | PASS (7 s) |
| `python3 bench/bench.py check` | PASS (9 s): 7 instances match `bench/reference/` exactly (tolerance 0.0%) |
| `tests/dnn_test.py --parts ACDE` with `QAP_DNN_WARM_Y=2 DNN_SHIP_BASIS=1` + production | PASS (64 s) |
| `tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1`, same environment | PASS (120 s): 4,536 of 4,536 |
| the same brute force with `DNN_SHIP_BASIS=1` alone (default DNN settings) | PASS (127 s): 4,536 of 4,536 |

**The test fix** (`tests/ship_test.py`, committed with this study). As
committed in `eb4828b`, `make shipbasischeck` failed with "no node was seeded
from a basis in the straight run". The solver was right: its trace lines end
` seed 1 res 0.0779`, and the test counted lines ending in ` seed 1`, so it
found none. A second bug was waiting behind it: the bf16 case `ship64b16`
has `QAP_DNN_SHIP=64`, so the test held it to bit identity with the straight
run, although its docstring says a bf16 basis is only reported (it moves the
trajectory: 8,730 DNN iterations against 8,710). The fix matches
` seed 1 res R` and compares node for node only the exact case.

## 9. Decision: killed on the CPU line

The kill line: less than 3% CPU saved on the warm children, or any child
that fathoms later or loses a fathom beyond rounding.

- **CPU: killed.** On the 24 warm children the switch saved 0.79% (pass 1)
  and 0.29% (pass 2) of the CPU [m]. The noise of such a total within a
  pass is about 0.5 points; between passes it is 13% [m]. The part that
  stands out of the noise is the timed start-up: 263-274 ms per warm
  solve, 0.56-0.61% of the warm set [m]. Restricted to the 15 children
  that fathom, it is 2.9% and 1.2% [m]; on 8 short children on one core,
  2.3% and 4.3% [m]. Carried to `bestcore`'s node mix (276 warm first solves,
  median 20 iterations, 39% of the DNN iterations in cut rounds [m]), 263-274 ms
  per warm solve is 1.4-1.5% of the DNN CPU and 2.3-2.5% of the first
  solves' [e] (`est_production.py`). Every measure is under 3%.
- **Fathoms: none lost, none gained.** The same 15 fathoms and the same
  3,461 iterations in total. One child fathoms one iteration later (f6/3,
  52 to 53) and one fathoms one iteration earlier (f6/2, 24 to 23), the same
  in both passes and with bf16 [m]. Both are marginal crossings: f6/3's lb
  at its fathom is 0.18 above the cutoff. A start perturbed at the seed's
  residual (0.011-0.035) moves such a crossing by one iteration either way. I count
  this as rounding. `ship0` moves three fathoms later by one iteration
  each and none earlier: it fails this part too.
- **Why the saving is small.** The tail probe's 6-9% was the whole first
  refresh. A seeded first R-step still costs about 290 ms, the seed plus one
  filter step, so the switch removes about half of the 565 ms [m]. Even a
  free first R-step, costing no more than a regular one (about 230 ms, so
  335 ms saved), would be 1.7% of the DNN CPU and 2.8% of the first solves
  [e]. No
  start-up trick passes 3% on this node mix. Cut rounds are not seeded
  (`!par.cut_round`). `bestcore`'s trace has 23 rounds in 6 phases, so
  seeding them would add about 0.25% [e].

So `DNN_SHIP_BASIS` stays as committed: opt-in, off by default, with the
default path bit-identical and every gate passing. No `bestship` arm was
added to `../../pipeline_ab/run_pair.sh`. The state-file format (version 7)
stays readable; nothing writes it unless both switches are set.

## 10. Files

- `run_all.sh`: the warm-set measurement, two passes of three configurations
  at the same time on CPUs 1-3; it calls `run_warmset.sh` (the six families,
  the production flags and environment). Logs: `warm/<config>_<family>.log`,
  `warm/run_all.out`, and the summaries `warm/summary_pass1.txt`,
  `warm/summary_pass2.txt` (each configuration against the base of its own
  pass) and `warm/summary_noise.txt` (the replicates across passes).
- `run_samecore.sh`: the same-core replicate (f9 and f6, base ship ship base on
  CPU 0). Logs: `samecore/`, summary `samecore/summary.txt`.
- `summarize.py`: per child and in total against a base, with the kill line;
  `--profile [--short]`: where the children's CPU goes.
- `est_production.py`: the saving per warm solve carried to the production
  node mix (`../../pipeline_ab/runs_wfix/bestcore/trace.txt`).
- `check_identity.sh`, `state_check.sh`, `fault_test.sh`: section 3; logs in
  `logs/ident/`, `logs/state/`, `logs/fault/`, each with a `summary.txt`.
- `gates.sh`: the repository's gates (SET=single on one core, SET=multi for
  `make test`); logs in `logs/gates/`, one line per gate in
  `logs/gates/summary.txt`.
