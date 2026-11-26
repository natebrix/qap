# nug25's root DNN solve ends `status error` on the laptop

Reported by the paper session (branch `claude/dnn-plus-paper`, `docs/dnn_plus/root_runs.md`), Nov 26 2026: on
an Apple M2 Pro (macOS, Apple clang, Accelerate) nug25's root DNN solve under `plain.par` (nug30's recipe, the
bound at the root only, no cuts, 4,000 iterations) stops at iteration 2,662 with `status error` and bound
1,707.3828, twice the same. On Linux with OpenBLAS it runs its 4,000 iterations (3,625.86 here). 41 other
problems ran normally on the laptop.

Two changes in `latest/dnn.C` on this branch (the file is otherwise cognac's):

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

## On the laptop

From any checkout of the repository (it builds in a temporary worktree and touches nothing else; about a
minute):

```
git fetch origin claude/cool-babbage-ilnogj && git show FETCH_HEAD:experiments/mac_nug25/run.sh | bash
```

It prints the diagnostic lines and the solve's trace line. Expected if the failure repeats: `status error` at
2,662 again, now with a bound near 3,625 and the name of the failing call.
