# A second review of b86352e and of pool_retry.patch (cognac session, Nov 16)

An independent read of `git diff e1596e2 b86352e -- latest` (8 files, +483/-13) and of
`../tai_solves/pool_retry.patch`, by a subagent of the cognac session (scratch files and test programs
outside the repository); two of its claims re-read here against the code. With the cause of the
`continuationcheck` failure that `README.md` here left open.

## Verdict

- **No bug that can make a bound invalid.** The first reviewer's assumption about the scaled transform
  (`dnn.dscale`), which they did not review, holds:
  - **The algebra checks.** The scaled basis spans the scaled range; the reflector ordering matches; the
    scaled Y/Z step is the plain step on the scaled problem; the certificate multiplier
    Zc = D(Zt + beta(Yt - VRVt))D gives V'ZcV negative semidefinite.
  - **The numbers agree:** a numerical test at m = 7 and 12 with random D, and nug15 and tai20b roots
    under all five scaled variants stay below the optimum.
  - **Every certificate is computed in original units.** `w.sc` is cleared before every call to
    `certify`, `cert_lam` and `cert_poly`; there is no `return` inside the loop; every `break` passes
    the back-conversion before `save` is written. Children read the parent's original-unit `save`;
    state files hold finished `save`s only.
  - **What a mistake there could still do:** cost iterations *and change the search* (branching on the
    solution's marginals, cut separation, child starts), not a bound.
- **The other new keys are sound.** Regional triangle cuts are valid for any facility subset, and pools
  match regions by membership. `cut.cp_extrap` keeps the cut multipliers non-negative. `cut.handover`,
  `cut.giveup_depth` and the `cut.cp_wc_*` keys only stop a phase or choose its penalty.
- **The 18 findings in `README.md`:** 17 agree, one (17) partly. In a rerun of the all-keys tai15b case
  every region was non-empty, so the `status error ... iters 0` node took the "no cut separated at round
  0" path (`dnncut.C:2289`), not the empty-region one (2263); inferred, not instrumented.

## New findings (none affects validity)

1. **`cut.cp_extrap` does almost nothing with `cut.cp_lazycert = 1`.** A lazy guard check has value -inf
   (`dnn.C:2286, 2327`), which reverts the jump and blocks every later one (1960, 2330-2333). Latent: no
   settings file sets both, and `cut.cp_lazycert` defaults to 0.
2. **The `-e -r` fix is partial.** `estimate.C` sets the root bound after step 1's `bnb()` has run
   (checked here: the assignment follows the call). Resumed nodes shallower than `start_dive_level` are
   bounded in step 1 with root bound 0, which affects `pick_param`'s level and the DNN max-gap filter.
   tai50b's dives start at level 1, so they are unaffected. Fix: set `top_stats->root_bound` before the call.
   **Fixed in cognac (Nov 16):** set before the call. It mattered more than "the level and the filter":
   on nug15 from a depth-0 state, all six depth-1 nodes of step 1 took strategy 2 (the deepest level's)
   where the run uses 0 and 1, with weaker bounds and other branching lines. `tests/estimate_test.py`
   now checks that step 1 bounds a state's nodes as a plain resume does.
3. **More silent off-paths of `dnn.dscale`** than finding 10 lists: m < 3, the prototype certificate, and
   `set_scaling` returning false (`dnn.C:559, 586`, a Householder vector formed without the sign choice
   that avoids cancellation). Only `-v` reports it.
4. `dnn.beta_sched = rb` with `dnn.dscale` mixes scaled norms and original-unit residuals (`dnn.C:2170`).
5. `cut.cp_extrap` keeps two to four extra (m^2+1)^2 arrays: 100-200 MB per process at m = 50 during a phase.
6. `cut.cp_wc_gap` is in absolute objective units (`dnncut.C:2318`; `cut.cp_wc_deficit` is relative): the
   300 tuned on tho40 means nothing on the tai-b instances.
7. Style: `cut.handover` ignores fields past the seventh; the scaled Y-step computes its `-v` diagnostics
   always.
8. Lab tool only: `dnncuts -phase -loadwarm` reports max(first_lb, res.lb) with `first_lb` from the file,
   so a state file from another node prints a best bound that is not this node's.

## pool_retry.patch

**Sound, and applies cleanly on b86352e** (no fuzz; builds). Every result kept is a complete
`dnn_bound` result, and `node.bound` takes the maximum. The parent's state is const and separate from
`save`; the failed pool is not passed on. With no retry firing, the all-keys nug15 DNN trace is identical
to unpatched b86352e's.

**How often it would fire.** Its trigger, replayed on existing DNN traces:

| traces | warm child solves | would retry |
|---|---|---|
| tho40 (`rule2_ab`, `tho40_search`) and wil50 (`wil50_policy`) | 6,545 (4,080 pooled) | 0 |
| tai35b | 300 | 32 |
| tai40b | 630 | 0 |

**Change before a long run:**
1. **A settings key, default off.** As written it adds none, so the settings record does not show it:
   patched and unpatched workers resume each other's state files without warning. Its second step (solve
   cold) fires on *every* warm-started solve that ends 1% under its inherited bound, pooled or not
   (checked here: the condition is `sp && ...`), so it reaches every run with `dnn.warm = 1`, tho40's
   and tai50b's included.
2. **Record each retry** in the DNN trace and the stats: which step, and the failed bound. Today a
   retried node shows `pool 0 kept 0`, and still `warm 1` after a cold re-solve; stderr stops after 20.
3. **Tie the trigger to the gate it protects.** The node's phase is decided on the solve's own bound
   (`bnb.C:884-885`). With `cut.zone = 1` (the tai solves) only the catastrophic failures (-3e11, -inf)
   lost their phase, and 1% catches those. With a zone under 1% (0.003-0.06 elsewhere), a solve ending
   0.5% under its inherited bound loses its phase with no retry. A trigger on "the inherited bound is
   inside the zone and the solve's is not" fits both.
4. Then the identity runs on a tho40 slice and a wil50 subtree (the replay predicts no change), and the
   tai35b case again.

**In cognac (Nov 16), for tho40-c:** the integration build's `latest/` (205bd9b, line for line b86352e)
and the retry (5f644a1), applied as a patch, not merged. Items 1 and 2 are done: `dnn.retry` (default 0)
switches both retries, and each one adds `retry S rlb X` to the node's DNN line (S: 1 without the pool,
2 cold, 3 both; X the first failed bound). Item 3 is not: at tho40's `cut.zone = 0.06` the 1% trigger lies
inside the zone. Checked: `make test`, `dnncheck`, `shipcheck`, `nodetrace`, `continuationcheck` pass; on
the tho40 subtree `26:10,8:0,30:32` under cut06's settings (38 nodes, 37 pooled first solves) the node and
DNN traces equal cognac aa5113d's but for pids and CPU, and `dnn.retry = 1` fired no retry there.

## `make continuationcheck`

**It has failed since 34e347c (Nov 8), on any machine; not a solver bug.** The check compares a DNN solve
split into warm-started chunks (10 + 15 + 20 iterations, `warm_primal = 1`) with one uninterrupted 45-
iteration solve, to 1e-12.
- **Bisected** (`git bisect` over the commits touching `latest/` since e6f4341, where the test was
  added): passes through 20fecaa, fails from 34e347c, "The nug30 record recipe is the solver's default".
- **The cause:** at 34e347c, reverting each new default alone in the environment, only `DNN_BETA_SCHED=`
  passes. The record recipe's `dnn.beta_sched = w=0.36` runs every warm-started call, these chunks
  included, at 0.36 x the base penalty by design, so the chunked iterates differ from the cold ones. Every
  certificate stayed safe throughout (648 safety checks in every run).
- **Fixed in the test** (cognac f07d83f): it pins the schedule off before solving. It passes on cognac
  and on b86352e.
