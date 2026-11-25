# tai80b's root with two-dimensional regions (Nov 25 2026)

Question: tai80b's root closes 46% of its published gap with the bound lab's regions (`comb7` =
`line:1,top:7,linegap:4:0.03`, 138 regions), where tai60b's closes 81%. tai80b's location clusters lie
in the plane (first MDS axis share 0.78; tai60b 0.92), and `comb7` cuts along one axis only. Does a
family cut along several directions close more? The bound lab named this as the next test and did not
run it (`experiments/bound_lab/tai60b/README.md` on `claude/cool-babbage-ilnogj`, the tai80b sections).
**[m]** measured here, **[e]** estimated.

## Answer

**No. With `line:4,top:7` (350 regions, every one of `comb7`'s among them) the root is 1.4M behind
`comb7` at phase iteration 500, at twice the CPU per iteration. tai80b stays out of reach.** [m]

| phase iterations | `comb7` (138 regions) | `line:4,top:7` (350) | difference | share of the published gap, `comb7` / `line:4,top:7` |
|---|---|---|---|---|
| 100 | 781,163,383 | 781,116,658 | -46,725 | |
| 250 | 787,033,428 | 786,291,076 | -742,353 | 0.4% / -1.9% |
| 500 | 793,523,957 | 792,127,493 | -1,396,463 | 21.0% / 16.6% |

- Behind at every check from iteration 200 on, by 1.1-1.4M from 300 to 500; the two gain about the same
  per 50 iterations by 500 (0.70M and 0.71M from 450 to 500), so the curve runs parallel below `comb7`'s
  and shows no sign of crossing it.
- 7.3 s of CPU per phase iteration against 3.4 (250 to 500). At equal phase CPU (3,637 s), `comb7` is
  at iteration 1,030 and 798.6M, 6.5M ahead.
- Stopped at iteration 510 by the bound lab's rule for this test (stop if not ahead of `comb7` at 500;
  `combh`, 200 regions, was stopped by it on Nov 16, 19k behind).
- What tai80b would need: tai60b's 81% is 812.4M, 11.0M above `comb7`'s plateau (801.4M at 2,320
  iterations). A family that is behind at 500 is not going to find that.

**Why a superset can be behind.** The CP-ADMM phase (`cut.cp = 100:1000`) adds, every 100 iterations,
the 1,000 most violated regional triangles over all regions, ranked by violation
(`latest/dnncut.C`, `separate`). With 350 regions the new directions' triangles take slots `comb7`'s
would have had, so this is not "`comb7` plus more cuts" but a different choice of 1,000. A run with
`cut.cp = 100:2500` would separate the region family from the cap; not run, since the verdict above
needs an 11M gain, not a recovered 1.4M.

## How it was run

- Solver: `cognac` at 70a14ffe with `diag.patch` (the bound lab's `live.patch`: `QAP_CP_LIVE=FILE`
  logs the phase's certified running best at every check; and its `threads.patch`:
  `QAP_BLAS_THREADS=N`, default 1 as before). Built against the `scipy-openblas32` wheel: the first
  attempt linked the system's reference BLAS (single-threaded, about 7x slower, `RUNNING.md`) and was
  stopped after 2 hours inside the first solve.
- `par/root80l4.par`: the bound lab's `par/root80.par` with `cut.regions` alone changed.
- `regions/tai80b.line4top7.txt`:
  `python3 -I experiments/bound_lab/regions/gen_regions.py qaplib/tai80b.dat line:4,top:7`. The same
  tool with `line:1,top:7,linegap:4:0.03` reproduces the bound lab's `tai80b.comb7.txt` byte for byte.
- From `latest/`, alone on a 4-core cloud container (15 GB; 8.0 GB resident at iteration 500):
  `OPENBLAS_NUM_THREADS=4 QAP_BLAS_THREADS=4 QAP_CP_LIVE=live.txt ./bnbdrv -p PAR -u 818415043 -d 0 -o root.state ../qaplib/tai80b.dat`
- `runs/l4t7/live.txt` this run; `runs/root80/live.txt` the bound lab's `comb7` run (its
  `runs/root80`, 4 threads), copied for the comparison. The first solve is the same in both (774,516,682.66
  at the first check); this container took 1,423 s of CPU for it against 1,796, so the CPU columns
  compare within a run better than across.
