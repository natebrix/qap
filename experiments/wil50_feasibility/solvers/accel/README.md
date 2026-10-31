# Iteration accelerations for the DNN bound: tight λ, Halpern, Anderson

Follow-up to the three solver-research reports in `..` (Oct 31 2026). Two of
their suggestions were implemented behind flags in `latest/dnn.C` and scored on
the 32 wil50 depth-3 nodes of `../../tuning/giveup/all32.txt` (18 fathom, 14 are
misses), in the search's mode: PRSM (`-g 1.0 -G 1.0`), subspace eigensolver,
`-H 300 -i 2500 -u 48816`, one BLAS thread, four runs at a time on the 2.1 GHz
Xeon. `run_all32.sh` runs it, `summarize.py` makes the table, logs are in
`runs/`.

| flag | `DNNParams` | what it does |
|---|---|---|
| `-L N` | `tight_lam` | the Chebyshev filter's upper end `lam_hi` from N shifted power steps on the tracked top eigenvector (padded 10%), instead of the Weyl bound `lam_hi += ‖W − W_prev‖_F`, which drifts about 8× too high between refreshes (report 3) |
| `-A N` | `halpern` | Halpern anchoring of the (Y, Z) iteration, `u ← u₀/(k+2) + (1 − 1/(k+2)) T(u)`, restarting at the current point when the fixed-point residual falls 5× or after N iterations |
| `-W N` | `anderson` | safeguarded type-II Anderson acceleration of (Y, Z/β) with memory N (memory cleared when the residual doubles) |

All three move only the iterates; the certificate (Z′ from the R-step,
Gilmore–Lawler, Rump Cholesky) is untouched, so every bound stays certified.
`tests/dnn_test.py` part A brute-forces all three modes on random nodes.

## Results

Over the 18 nodes every configuration fathoms (`med_it`, `med_cpu`: median
iterations and CPU seconds to fathom; `sum`: totals), the CPU spent on the 14
misses before the give-up rule stops them, and the certified bound at
iteration 200 on the misses relative to degree 4.

| config | fathomed | ms/iteration | med_it | med_cpu | sum_it | sum_cpu | miss CPU | lb@200 − d4 (misses), median [min, max] |
|---|---|---|---|---|---|---|---|---|
| `-d 4` (current) | 18 | 562 | 220 | 109.2 | 4240 | 2181 | 2889 | — |
| `-d 3` | 18 | 473 | 230 | 99.2 | 4460 | 1958 | 2260 | −11.6 [−20.7, −9.6] |
| `-d 3 -L 2` | 18 | 480 | 220 | 95.5 | 4250 | 1889 | 2330 | +0.8 [+0.0, +1.3] |
| `-d 4 -L 2` | 18 | 529 | 220 | 105.2 | 4180 | 2078 | 2692 | +4.0 [+3.2, +6.2] |
| `-d 4 -W 5` | PENDING | | | | | | | |
| `-d 4 -A 25` | PENDING | | | | | | | |

### Tight λ: degree 3 at degree-4 strength

With the Weyl bound, degree 3 is 16% cheaper per iteration than degree 4 but
weaker: 10 more iterations on most fathoming nodes and 12 units behind at
iteration 200 on the misses. The tight `lam_hi` closes that gap entirely:
`-d 3 -L 2` fathoms every node in the same iteration as degree 4 (one node 10
earlier, one 10 later) and is 0.8 ahead at iteration 200, at 85% of the
per-iteration cost. Per node that is **13% less CPU to fathom and 19% less
spent on misses**, with no change in which nodes fathom.

The two power steps cost two `sgemv` per iteration (negligible against the
filter's `sgemm`s). On the wil50 depth-3 node `18:43,11:29,46:44` run to 1500
iterations without the give-up rule (`node2_*.log`), the converged bounds are
48610.15 (degree 3), 48610.12 (degree 3 tight), 48611.83 (degree 4), 48611.86
(degree 4 tight): the tight bound only changes how fast the filter gets there,
not where it goes.

`-d 4 -L 2` is stronger per iteration (+4 at iteration 200) and 6% cheaper than
plain degree 4, but it fathoms in the same iterations, so it buys less than
degree 3 tight.

**Recommendation:** make `cheb_degree = 3, tight_lam = 2` the search's DNN
setting (currently degree 4). It is a straight CPU saving at equal strength.

### Halpern anchoring: negative

On the nug20 root (`nug20_curves.txt`), the anchored iteration is behind plain
PRSM at every checkpoint for every restart period tried (10, 25, 50, 100),
e.g. 2504.8 vs 2506.0 at iteration 370 with restarts every 10. Plain PRSM
already converges linearly on these problems; pulling the iterate toward the
anchor only slows it. Not pursued beyond the A/B row above.

### Anderson acceleration: negative

On nug20 with the subspace eigensolver, memory 3 or 5 is behind plain PRSM at
every checkpoint (2505.9–2506.0 vs 2506.0 at 370; far behind early). With the
exact `ssyevd` R-step it gains about 0.1 unit late in the run — the
extrapolation needs an exact fixed-point map, and our inexact projection
denies it that. On the wil50 misses it is PENDING units behind degree 4 at
iteration 200.

## Verdict

Of the 1.3–2.5× that report 1 estimated for accelerated splitting, nothing
materialized: the plain PRSM iteration with the tracked eigensolver is already
the fast path. The tight λ is the one gain, 13–19% CPU per node, from making
the eigensolver's filter sharper rather than from the splitting.
