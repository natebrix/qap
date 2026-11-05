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

## Results

(pending)
