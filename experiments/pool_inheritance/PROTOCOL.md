# Parent cut-pool inheritance: predeclared protocol

Based on cognac `0be4e7d` (after `cdf8e18` and budget-study merge `5038718`).
Only `latest/dnncut.{h,C}` and `experiments/dnn_cuts/dnncuts.C` change solver
implementation. The new file options default off. No search integration,
production defaults, symmetry files, or branching rules change.

## Question and mathematical limits

Can the frozen parent pool buy the child's cut-phase bound for the cost of
one quiet solve? This tests roadmap/bound_strength.md item 1.

The exact child Lagrangian value is at least the parent's, with the following
precise qualification. Restrict the parent DNN matrix to the child face via
`Y_parent = T Y_child T'`. A retained triangle has exactly the same polynomial
after remapping its region and the opposite-side variables. When one variable
is fixed to 0 or 1, the four triangle residuals become, up to order,
`-W_jk`, `W_jk-p_j`, `W_jk-p_k`, and `p_j+p_k-W_jk-1`.
These are nonpositive by entrywise nonnegativity **and the assignment identities**
(including moments on the complementary region). Dropping a term `mu*g`,
with `mu>=0` and `g<=0`, increases the child objective pointwise. Thus the exact
kept-pool relaxation retains the parent lower floor even after these drops.

**A finite-iteration child certificate need not exceed the parent certificate.**
The exact infimum is monotone; lower approximations to different infima need
not be. A frozen certificate below the parent's matching certificate is logged
as a diagnostic requiring investigation, not automatically called a remap bug
and never silently replaced by a parent floor in the gain gate. Regression
checks establish identical round-trip objective coefficients/constant and DNN
values, and test the child lift in both orientations. Any unexpected cutoff
crossing on n10 or d2 is separately flagged as a validity question.

Similarly, the child's measured own-phase certificate is a **comparison target**,
not a proved upper bound on the inherited relaxation. Independently generated
finite pools need not contain each other, and a finite certificate is not an
upper bound on a pool's optimum. We report the requested parent/child comparisons
without presenting them as a rigorous two-sided numerical bracket.

## Nodes and fixed sampling rule

All use `qaplib/wil50.dat`, `-grid 5:10 -side 0`, cutoff **48,814.05** from
`-u 48816`. Node specifications are zero based. The parent is
`39:27,30:48,25:49`. Known child rows are checked against depth4/nodes.txt.

| Child | Additional fix | Independent historical comparison target |
|---|---|---:|
| d1 | `7:5` | Fathom: >48,814.05 |
| d2 | `26:5` | 48,789.9 (`dual/tri_d2.log`, rounded as specified) |
| d4 | `47:5` | Fathom: >48,814.05 |
| s5 | `5:5` | None: no independent own-phase run in this matrix |
| s10 | `10:5` | None: no independent own-phase run in this matrix |
| s11 | `11:5` | None: no independent own-phase run in this matrix |

The three additional facilities are `random.Random(20261103).sample(candidates,3)`,
where candidates are the ascending integers 0..49 excluding 39,30,25,7,26,47.
This yields **5,10,11**, selected before seeing any transfer result.

Indexing is explicit: regions use original indices of `side`; triangle members
use original indices of the **other** side. With side 0 and children `a:5`, all
six remove opposite-side member 5, while their regions lose different rows a.
Accordingly the kept/dropped counts must agree across these six children.

## Parent and three child arms

All processes are single-threaded (`-T 1` plus BLAS/OpenMP thread limits), use
the existing tracked eigensolver defaults, start cold in Y/Z, and leave `-warmY`
off. We skip optional parent-Z files: this first experiment isolates transfer
of the cut pool and multipliers. Within continued phases the usual Z-only
same-node warm start remains in use.

Common flags: `-grid 5:10 -side 0 -u 48816 -T 1 -v -delta 50 -nocap`.

| Arm | Purpose | Additional budget/options | Maximum iterations |
|---|---|---|---:|
| P | Parent n10 B phase; starts with empty pool | `-I 400 -r 12 -i 300 -k 1000 -freeze 0:1500 -savepool FILE` | 5,500 |
| A | Child uncut control | `-I 1500 -r 0 -i 300 -k 0` | 1,500 |
| B | Child frozen inherited pool | `-loadpool FILE -I 0 -r 0 -i 300 -k 0 -freeze 0:1500` | 1,500 |
| C | Child inherited pool, then Polyak updates/new cuts | `-loadpool FILE -I 400 -r 12 -i 300 -k 1000 -freeze 0:1500` | 5,500 |

There are **19 runs: one parent and 18 child runs**, once per arm. B performs
exactly one quiet solve, with no initial uncut solve, separation, or multiplier
update. C's initial 400 iterations already use the inherited pool. A's 1,500
iterations are a finite budget, not a proof of exact convergence. Standard
early fathoming stays enabled in every arm; any such value is a stopping
certificate, not a fully converged value. Gains in the gate use the actual
retained certificates, without imputing a higher value after a fathom.

The source pool is the multiplier vector actually used by the parent's final
solve. With the B protocol this is the best Polyak round restored for quiet
polishing. The file records that round, the quiet flag, and its matching
certificate, rather than assigning the overall best value to an unrelated
last multiplier vector. Zero-multiplier cuts are retained and counted.

## Predeclared gates and reporting

Let `g_i = B_i - A_i` for the six children. **KILL** if either:

1. `median(g_i) < 40`; or
2. `B_d2 < 48759.9` (48,789.9 minus 30, exactly as requested).

Otherwise **PASS**. Report every gain, the median, min/max and sorted spread.
Differences within 0.000002 of a gate are rounding-ambiguous instead of declared
passes from six-decimal output. Incomplete/failed runs cannot qualify. Validity
questions are reported separately from the numerical gate decision.

For each arm report the retained certificate, fathom outcome, cumulative DNN
iterations, and whole-process user+system CPU (including pool I/O, separation,
assembly, and quiet work). For B/C also report kept/dropped counts, source
pool provenance, and B minus the matching parent certificate. The requested
parent reference is approximately **48,743.88**; the actual freshly saved
matching certificate defines the diagnostic comparison.

For C on d1/d2/d4, report the first **logged round** reaching the historical
target, cumulative iterations, and a polled whole-process CPU observation at
that stage. Fathoms terminate at their first checked cutoff crossing (ten-iteration
check schedule); non-fathom target crossings are known only at round boundaries.
CPU observations may include work after the boundary because the runner polls;
only terminal whole-process CPU is exact process accounting. Preserve all rounds.

The three sampled siblings have no independent own-phase target in these 18
child runs. Report their full continuation progress and any fathom, without
inventing an own-phase optimum. This matrix measures time to known targets;
it does **not** provide a fresh paired child-only phase CPU baseline for all six.
Historical laptop CPU cannot be used as a paired cloud speedup measurement.

## Execution order and provenance

[plan.json](plan.json) contains the **exact default commands and queue**, generated
and committed before launch. The parent and six A runs are ready immediately.
B/C depend on successful completion of P and the saved pool, and cannot start
from a partial file. Queue shuffles use seed **20261104**, independent of the
facility sample. Within the shuffled dependent-node order, B then C are queued.
No arm, node, or parameter is selected from observed outcomes.

At most eight processes run simultaneously, pinned to distinct physical cores
and limited by the eight-CPU quota. The same 8 GiB container successfully ran
eight wil50 workers in the budget study. The initial uncut controls overlap
the parent; later B/C runs use free cores. The expected total is about 4–6 hours
on this container, depending on early fathoming; the runner stops at six hours
and preserves every started log if incomplete. Shared bandwidth, core assignment,
and tail concurrency limit interpretation of single-run CPU timings.

The manifest records the code commit, binary/input/plan/pool SHA-256 hashes,
commands, CPU placements, start times, wall time, whole-process CPU and maximum
RSS. Every started run keeps its raw log; an event stream records stage
observations. The parent pool is immutable throughout the child matrix.

The v1 file format accepts symbolic regional triangle types 0–3 only. Atomic,
five-point/local-7, nonzero `mu0`, and external-objective combinations fail closed
when pool I/O is requested. Unknown versions/types, malformed indices,
non-ancestor nodes, and invalid/nonfinite/excessive multipliers are rejected.
No arbitrary coefficient file can bypass the existing triangle validity argument.
The original library behavior is unchanged when both file paths are empty.

Before launch: run `make test`, the existing continuation regression, the new
`tests/dnn_pool.C` regression via `tests/dnn_pool.mk`, and driver smoke checks.
The new regression verifies exact round-trip Lagrangian coefficients and values,
72 child remaps with 8,640 permutation lifts across both orientations, and
one-solve frozen behavior. The default-off path is compared with the prior driver.

To reproduce, use the host-appropriate OpenBLAS linkage (exact cloud flags will
be retained under results/), then:

```bash
make -C latest dnncuts
make -C latest -f makefile -f ../tests/dnn_pool.mk poolcheck
python3 experiments/pool_inheritance/run.py --out /tmp/qap-pool-repeat --workers 8
python3 experiments/pool_inheritance/summarize.py /tmp/qap-pool-repeat
```

The output directory must be new. `--out` changes only the recorded pool/output
path, not the node/arm queue or parameters. The full results go in RESULTS.md
after completion; production integration and roadmap changes are left to Fable.

## Execution note: intermediate CPU telemetry

After the uncut controls, before any inherited solve, we found that this executor
does not expose the runner's child PIDs through `/proc`. The launched runner is
pinned at `ef25c61`; its intermediate CPU observations are therefore unavailable,
and its recorded round observations are collected at process exit. Terminal
`wait4` user+system accounting and all raw solver logs are unaffected. No node,
arm, parameter, stopping rule, or numerical gate changes.

For the first logged target round, the report instead bounds prefix process CPU
between the sum of the preceding DNN call times and terminal process CPU minus
the subsequent DNN call times. All unallocated CPU overhead can fall on either
side. Allow 0.05 seconds per printed round plus 0.01 seconds accounting slack.
This uses disjoint main-thread CPU intervals within total process CPU; it does
not assign a fictitious exact time to a target. Fathoming runs also have their
exact terminal process CPU. The analyzer is updated accordingly; the runner's
future-run telemetry uses the kernel process clock and still records round
progress if that clock is unavailable. The active run continues unchanged.

## Execution amendment: executor loss and four completion replays

At 2026-11-03 18:21:40 UTC the executor session was lost. The runner heartbeat
and all four active C logs stopped together; process session 16235 became
unavailable. Fifteen runs had completed successfully, including **all A and B
runs**. No solver error was logged. A full immutable copy of that attempt is
retained under `results/interrupted_original/`, including its manifest, event
stream, partial C logs, parent pool, and runner output.

Before any replay, [recovery_plan.json](recovery_plan.json) declares the four
unfinished runs: d4/C, s5/C, d2/C, d1/C. [recover.py](recover.py) repeats only
these runs from scratch, with the **same binary, pool, commands, and individual
CPU affinities**, now four concurrently. It preserves the original six-hour
deadline. Nodes, budgets, seeds, numerical gates, early stopping, and production
code are unchanged. Each replay is compared with every logged numerical round
of its interrupted prefix. The completed A/B evidence is never rerun or selected.

This is an explicit deviation from the planned single attempt per arm: a
completed study will have **19 completed runs and four interrupted attempts**.
Whole-process CPU is measured for each completed attempt. Terminal CPU for the
four interrupted attempts is unavailable; their logged DNN CPU supplies only
a lower bound on the additional work, reported separately. C timing after
recovery also has a different concurrency history. The frozen numerical gate
is unaffected; no whole-search speedup is inferred from these timings.

### Second executor loss and final detached retry

The first replay session stopped at 18:35:15 UTC, before any of its four C
solves completed a round. Resuming its session returned `executor key changed
during session recovery`. Its unchanged files are under
`results/interrupted_recovery1/`. No solver error or completed C value occurred.
Before relaunch, `recovery_plan3.json` pins a final retry of the same four
commands, binary, pool and cores, using a detached runner. The original six-hour
deadline remains; another infrastructure loss will be reported as incomplete,
without further retries. A completed matrix therefore has 19 completed runs
and **eight interrupted attempts (27 launches)**. The second interruption has
no completed DNN-call CPU observations, so its CPU overhead is unknown.

A recovery-code review also removed an unavailable Python process-clock call
before it was exercised. Round observations now explicitly have unknown CPU;
terminal wait4 CPU and the conservative prefix intervals remain the accounting
method. This changes no solver code, inputs, bound or gate.
