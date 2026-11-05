# Pipeline A/B: the Nov 4 candidates on the node-10 subtree

The convergence and certificate studies (`../roadmap/convergence2/REPORT.md`,
`../roadmap/certificate/REPORT.md`) left four candidates. Each is tested
end to end on the node-10 subtree through the coordinator, against the
production configuration on the same machine at the same time.

| arm | what it adds | judged against | kill line (DNN iterations) |
|---|---|---|---|
| `base` | production: rule 9, bundle D, parent-Y warm start, pool inheritance (`../coord_search/run_coord_search.sh`) | - | - |
| `zone` | `QAP_DNN_CUT_ZONE=0.005` (from 0.003): misses up to 244 short run a cut phase, node 10 (232 short) among them, so their children inherit a pool (the "equip phase") | base | 1.15x |
| `warmy` | zone + Y kept between cut rounds (`QAP_DNN_CUT_WARMY=1`), 12 rounds of 200 iterations, quiet round at most 1,000 (Astra's Y12s), first Polyak step damped (0.3) for phases that start from an inherited pool | zone | 1.05x, and the tree no larger |
| `certx` | the certificate polish within 10 units of the cutoff (`DNN_CERTX=32 DNN_CERTX_USE=1 ...`) | base | 1.03x |

Iterations are the deciding measure: they do not depend on the machine.
Running two arms at once on one machine also makes their CPU seconds
comparable. The smoke runs on nug20 (`runs_smoke/`, not committed) checked
the plumbing only.

## Running it (laptop)

```
git pull origin cognac
make solver                      # latest/bnbdrv; on macOS with Accelerate
experiments/wil50_feasibility/pipeline_ab/run_pair.sh base zone
```

- Each arm gets half the performance cores (`SLOTS=N` overrides).
- The script keeps the Mac awake (`caffeinate`) while it runs. Keep the lid
  open, or the power and an external display connected.
- `kill PID` (the script prints it) or ^C stops it cleanly. The workers
  drain, and the same command resumes the arms where they stopped. A
  finished arm is skipped.
- At the end it writes `summary.txt` and prints the `git add/commit/push`
  line for the results: `runs/ARM/` trace, results, logs and hosts. The
  stores, work directories and binaries stay out of git.
- The second pair, `run_pair.sh warmy certx`, needs the first-step damping,
  which comes in a later commit (`git pull` and `make solver` first; the
  script refuses `warmy` on a binary without it).

`summarize.py runs/` can be run at any time; unfinished arms are listed as
such.

## The damping setting (Nov 5, cloud box, `damp_test.sh`, logs in `damp/`)

Child 49:45 of node 10 with node 10's equip-phase pool
(`../roadmap/convergence2/runs/E_n10.pool`), the warmy arm's phase: Y kept
between rounds, rounds of 200 iterations, stall rule 3, uncapped target,
quiet round of 1,000 run regardless (in the search it runs only when the best
round's primal is above the cutoff, which none of these is).

| damping (first step : later steps) | rounds 1-4 | best round | after the quiet round |
|---|---|---|---|
| none | −102.4, −14.3, +13.3 (stall) | 48,744.8 | 48,760.9 (+29.5) |
| 0.3 : 1 | +16.2, +31.8, −227.6, −68.9 | **48,763.3** | **48,776.2 (+44.8)** |
| 0.3 : 0.5 | +16.2, +30.8, −72.2, −10.9 | 48,762.2 | 48,775.1 (+43.7) |
| 0.3 : 0.3 | +16.2, +30.0, +2.8, +30.2 | 48,761.7 | 48,771.3 (+39.8) |

- The first step is what matters: damped, the inherited multipliers improve
  at once instead of losing 102 units. The phase ends 18 units higher with
  one more round (800 iterations against 600), and the multipliers it hands
  to 49:45's children are those of a better round.
- Damping the later steps removes the round-3 overshoot but does not improve
  the best round. The warmy arm uses `QAP_DNN_CUT_DAMP1=0.3` (first step
  only).
- 49:45 stays about 50 units short either way, so it branches in every
  variant. Whether the better pool helps its children is for the search A/B
  to show.

## Results

(pending: `run_pair.sh base zone` on the laptop, then `warmy certx`)
