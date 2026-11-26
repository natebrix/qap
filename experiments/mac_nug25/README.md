# nug25's root DNN solve ends `status error` on the laptop

Reported by the paper session (branch `claude/dnn-plus-paper`, `docs/dnn_plus/root_runs.md`), Nov 26 2026: on
an Apple M2 Pro (macOS, Apple clang, Accelerate) nug25's root DNN solve under `plain.par` (nug30's recipe, the
bound at the root only, no cuts, 4,000 iterations) stops at iteration 2,662 with `status error` and bound
1,707.3828, twice the same. On Linux with OpenBLAS it runs its 4,000 iterations (3,625.86 here). 41 other
problems ran normally on the laptop.

Three changes in `latest/dnn.C` on this branch (the file is otherwise cognac's), the third below:

- **The fix.** With `dnn.lazycert = 1` most checks skip the rigorous test and store the multiplier with the
  best estimate, to be certified when the solve ends. An error exit skipped that step and returned the last
  value certified in full, an early one. It now certifies the stored multiplier on an error exit too. A stored
  multiplier is a copy from an earlier check, so the bound is valid whatever went wrong afterwards.
- **The diagnostic.** When the R-step projection fails, a line on stderr names the LAPACK call, its `info`,
  why that iteration ran a full eigendecomposition (scheduled, forced by the eigenvalue bound, or after the
  subspace step failed), and how many entries of W and of its fp32 copy are not finite; a second line gives
  the iteration. At most 20 lines a process.

## Checked on Linux (cloud container, OpenBLAS), cognac 8d8f3876 with this file

| run | status | bound |
|---|---|---|
| cognac's own binary | `iter_limit`, 4,000 | 3,625.8573 |
| with the changes | `iter_limit`, 4,000 | 3,625.8573, no diagnostic line |
| projection made to fail at iteration 2,662, old exit | `error` | **1,707.3828** (the laptop's value to the digit) |
| the same, with the fix | `error` | **3,625.8540** |

`make regression dnncheck`: OK. The forced failure was a temporary edit, not committed. A NaN written into W
at that iteration did not make OpenBLAS's LAPACK return an error at all, so what fails on the laptop is not
known from here: that is what the diagnostic is for.

## What fails (the laptop's run of `run.sh`, Nov 26)

```
dnn: R-step projection failed in proj_psd_subspace: sstedc returned info 125714 (K = 577; a scheduled refresh);
     W has 0 non-finite entries of 332929, largest finite |W| 1.99664; its fp32 copy 0 non-finite, largest 1.99664
dnn: the solve stops with status error at iteration 2663 (m = 25, cold start)
first lb 3625.8517 ... status error iters 2662        (one BLAS thread; the fix's bound, where it was 1,707.38)
first lb 3625.8555 ... status iter_limit iters 4000   (Accelerate's default threading: no failure)
```

Accelerate's single-precision divide and conquer (`sstedc`) does not converge on a finite matrix: info
125714 = 217 x 578 + 288, the block of rows 217 to 288 of the tridiagonal. It happens with
`VECLIB_MAXIMUM_THREADS=1`, which `docs/dnn_plus/root_bound.sh` and `infra/laptop/solve.sh` both set, and not
with the default threading (the first version of `run.sh` left the variable out and saw no failure), nor
under OpenBLAS.

**The remedy** (third change in `latest/dnn.C`): when a divide and conquer returns info > 0, the same matrix
goes to the QR algorithm: `ssteqr` on the tridiagonal `ssytrd` left (saved before `sstedc`, which destroys it)
in the refresh, `ssyev` / `dsyev` in the full-eigendecomposition modes. A line on stderr each time (at most 20
a process). Checked on Linux by making the divide and conquer "fail" and return garbage, a temporary edit:

| forced failure | result |
|---|---|
| none | nug25 3,625.8573, as cognac's binary |
| `sstedc`, every third refresh (11 of them) | nug25 runs its 4,000 iterations, 3,625.8573 |
| `dsyevd`, every 7th call (nug15, 3,000 iterations) | 1,140.56099785 (unforced 1,140.56099785) |
| `ssyevd`, every 7th call | 1,140.56019 (unforced 1,140.56021) |

`make regression dnncheck sb`: OK.

**Confirmed on the laptop** (Nov 26, cognac 459f4089 with this branch's `latest/dnn.C` at 7a98bf2a):

```
=== one-thread
dnn: sstedc returned info 125714 (K = 577): ssteqr instead, info 0
first lb 3625.8559 primal 3625.7101 status iter_limit iters 4000 cpu 16.91
=== default-threads
first lb 3625.8555 primal 3625.7266 status iter_limit iters 4000 cpu 17.69
```

The one failure in 4,000 iterations, the solve to its limit, the bound within 0.0004 of the default-threading
run's. Not yet on cognac.

## On the laptop

From any checkout of the repository (it builds in a temporary worktree and touches nothing else; about a
minute and a half):

```
git fetch origin claude/cool-babbage-ilnogj && git show FETCH_HEAD:experiments/mac_nug25/run.sh | bash
```

It runs the solve on one BLAS thread, then with the default threading, and prints each one's stderr lines and
trace line. Expected with the remedy: the one-thread run reports `sstedc returned info 125714 ...: ssteqr
instead, info 0` and runs its 4,000 iterations.
