# Same-node continuation: shorter cut-phase budget results

## Result

**Y12s is the only shortened budget that passes all three gates on every node.** It uses `-warmY`, 12 cut rounds of 200 iterations, and a quiet cap of 1,000, with the common 400-iteration uncut start. Every run begins with an empty cut pool. No solver code or production default changed.

Across the six paired nodes, Y12s used **19,902.171 versus 31,014.061 whole-process CPU seconds**: **64.17% of B's cost**, a **1.5583x phase-cost factor**. Total DNN iterations fell from **23,500 to 15,010** (36.13%), close to the 35.83% CPU reduction. The three near misses were fathomed earlier; all three expected short-of-cutoff nodes stayed below cutoff. **No validity question was triggered.**

All 30 predeclared runs completed successfully in **16,367.868 wall seconds (4 hours 32 minutes 48 seconds)**, consuming **115,232.458 CPU seconds (32.01 core-hours)** in total. At most eight single-thread processes ran concurrently, each on a distinct physical core. The matrix was committed before launch at `947801383c181a784b4ce2dea6f072f5bb9b72b1`, based on cognac `4095416756ef487251912396a727fac7ceacacaf`.

## Gates and roadmap interpretation

Gates compare each arm with its **fresh B run on the same node**: retained certificate at least B minus 2; the same fathom outcome, with strictly fewer iterations to fathom on 14, 6 and 9; and strictly lower whole-process CPU. First-crossing iterations include the initial 400 iterations. The cutoff is **48,814.05**. Quiet-only and retained values are separately reported below; an early fathom has no quiet solve. Full definitions and the original queue are in [PROTOCOL.md](PROTOCOL.md) and [plan.json](plan.json).

| Arm | Nodes passing all three gates | Baseline CPU covered by gates 1+2 | CPU/B on that passing subset | All-node phase factor |
|---|---|---:|---:|---:|
| Y12 | n14, n6, n9, n13, d2, n10 | 100.00% | 0.9050 | 1.1050x |
| Y6 | n14, n6, n9 | 29.68% | 0.5881 | Not qualified |
| Y8s | n14, n6, n9, n10 | 55.12% | 0.5941 | Not qualified |
| Y12s | n14, n6, n9, n13, d2, n10 | 100.00% | 0.6417 | 1.5583x |

Y12 is the full-budget control, not a shortened-budget winner. Y6 and Y8s save CPU but fail the bound gate on part of the sample. Their passing-subset ratios above cannot be applied to the whole run; the coverage column makes the omitted baseline cost explicit. All candidate/node pairs passed gates 2 and 3 in this batch. The failures were exclusively bound losses:

| Arm | Node | Bound minus fresh B | Failed gate |
|---|---|---:|---|
| Y6 | n13 | -16.196700 | Bound (minimum allowed difference: -2) |
| Y8s | n13 | -9.185342 | Bound (minimum allowed difference: -2) |
| Y6 | d2 | -10.963812 | Bound (minimum allowed difference: -2) |
| Y8s | d2 | -4.296116 | Bound (minimum allowed difference: -2) |
| Y6 | n10 | -14.672468 | Bound (minimum allowed difference: -2) |

For Y12s, the predefined near-miss and short-of-cutoff groups give different cost ratios:

| Y12s sample | B CPU s | Y12s CPU s | CPU/B | Phase-cost factor | Conditional factor at 46% |
|---|---:|---:|---:|---:|---:|
| All six nodes | 31,014.061 | 19,902.171 | 0.6417 | 1.5583x | 1.1973x |
| Short-of-cutoff nodes: 13, d2, 10 | 21,807.803 | 15,071.717 | 0.6911 | 1.4469x | 1.1656x |
| Near misses: 14, 6, 9 | 9,206.258 | 4,830.454 | 0.5247 | 1.9059x | 1.2798x |

The conditional factor is `1 / (0.54 + 0.46 * CPU_ratio)`. The all-six ratio suggests **about 1.20x**; using only the three short-of-cutoff nodes suggests **about 1.17x**. If those misses represent the roadmap's dominant class, the latter is the more relevant pilot estimate. Both are consistent with the proposed 1.1–1.3x item, provided the measured phase-cost ratio transfers to production and the 46% cost boundary matches this experiment. These are conditional projections, **not measured whole-search speedups**.

The evidence supports taking **Y12s into a production-policy/subtree A/B**. Preserving twelve cut/multiplier updates while shortening each inner solve worked across this sample; dropping rounds did not preserve the bound consistently. The policies also change subsequent cuts and multiplier trajectories, so this is not a clean causal decomposition of the mechanism. Keep production defaults unchanged pending the search comparison. This experiment does not establish a 4–12x reduction by itself or a result on sko49.

## Per-node and per-arm measurements

Best Polyak round is the earliest round attaining the highest printed cut-round bound, excluding the initial uncut solve and quiet solve. Round/quiet bounds have four printed decimals; final retained bounds have six. CPU is whole-process user plus system time, including separation and assembly. Eligible CPU/B is defined only where gates 1 and 2 pass. The first crossing is exact on the solver's ten-iteration certification schedule; no claim is made about unchecked intermediate iterations.

| Node | Arm | Best Polyak round (bound) | Quiet-only bound | Retained bound | Fathomed? | First crossing iteration | Total iterations | Process CPU s | CPU/B | Gates 1/2/3 | Eligible CPU/B |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|---|---:|
| n14 | B | 11 (48777.6488) | 48814.3434 | 48814.343442 | Yes | 4750 | 4750 | 6407.083 | 1.0000 | Reference | — |
| n14 | Y12 | 11 (48814.0938) | Not run | 48814.093776 | Yes | 3640 | 3640 | 5046.146 | 0.7876 | pass/pass/pass | 0.7876 |
| n14 | Y6 | 6 (48804.8201) | 48814.4214 | 48814.421355 | Yes | 2830 | 2830 | 3518.464 | 0.5492 | pass/pass/pass | 0.5492 |
| n14 | Y8s | 7 (48813.9654) | 48814.3354 | 48814.335355 | Yes | 2810 | 2810 | 3821.873 | 0.5965 | pass/pass/pass | 0.5965 |
| n14 | Y12s | 9 (48814.0656) | Not run | 48814.065613 | Yes | 2110 | 2110 | 2829.608 | 0.4416 | pass/pass/pass | 0.4416 |
| n6 | B | 3 (48815.3714) | Not run | 48815.371384 | Yes | 1260 | 1260 | 1524.070 | 1.0000 | Reference | — |
| n6 | Y12 | 2 (48814.3713) | Not run | 48814.371253 | Yes | 860 | 860 | 1071.117 | 0.7028 | pass/pass/pass | 0.7028 |
| n6 | Y6 | 2 (48814.3713) | Not run | 48814.371253 | Yes | 860 | 860 | 979.147 | 0.6425 | pass/pass/pass | 0.6425 |
| n6 | Y8s | 2 (48814.3713) | Not run | 48814.371253 | Yes | 860 | 860 | 1145.519 | 0.7516 | pass/pass/pass | 0.7516 |
| n6 | Y12s | 2 (48814.6134) | Not run | 48814.613353 | Yes | 780 | 780 | 950.680 | 0.6238 | pass/pass/pass | 0.6238 |
| n9 | B | 2 (48814.3058) | Not run | 48814.305754 | Yes | 990 | 990 | 1275.105 | 1.0000 | Reference | — |
| n9 | Y12 | 2 (48814.3695) | Not run | 48814.369455 | Yes | 810 | 810 | 1114.989 | 0.8744 | pass/pass/pass | 0.8744 |
| n9 | Y6 | 2 (48814.3695) | Not run | 48814.369455 | Yes | 810 | 810 | 916.781 | 0.7190 | pass/pass/pass | 0.7190 |
| n9 | Y8s | 2 (48814.3695) | Not run | 48814.369455 | Yes | 810 | 810 | 1015.247 | 0.7962 | pass/pass/pass | 0.7962 |
| n9 | Y12s | 2 (48816.2869) | Not run | 48816.286927 | Yes | 720 | 720 | 1050.166 | 0.8236 | pass/pass/pass | 0.8236 |
| n13 | B | 12 (48708.8084) | 48754.6194 | 48754.619389 | No | — | 5500 | 7480.778 | 1.0000 | Reference | — |
| n13 | Y12 | 12 (48757.0704) | 48762.5585 | 48762.558527 | No | — | 5500 | 7325.793 | 0.9793 | pass/pass/pass | 0.9793 |
| n13 | Y6 | 6 (48724.8136) | 48738.4227 | 48738.422689 | No | — | 3700 | 4306.408 | 0.5757 | FAIL/pass/pass | — |
| n13 | Y8s | 8 (48737.8887) | 48745.4340 | 48745.434047 | No | — | 3400 | 4495.178 | 0.6009 | FAIL/pass/pass | — |
| n13 | Y12s | 11 (48750.0484) | 48759.7927 | 48759.792738 | No | — | 3800 | 5424.863 | 0.7252 | pass/pass/pass | 0.7252 |
| d2 | B | 11 (48733.9886) | 48793.6819 | 48793.681902 | No | — | 5500 | 6436.866 | 1.0000 | Reference | — |
| d2 | Y12 | 11 (48791.8343) | 48799.3967 | 48799.396718 | No | — | 5500 | 6361.599 | 0.9883 | pass/pass/pass | 0.9883 |
| d2 | Y6 | 6 (48768.3173) | 48782.7181 | 48782.718090 | No | — | 3700 | 4461.581 | 0.6931 | FAIL/pass/pass | — |
| d2 | Y8s | 7 (48779.2169) | 48789.3858 | 48789.385786 | No | — | 3400 | 2850.807 | 0.4429 | FAIL/pass/pass | — |
| d2 | Y12s | 12 (48782.3703) | 48793.8459 | 48793.845910 | No | — | 3800 | 4490.928 | 0.6977 | pass/pass/pass | 0.6977 |
| n10 | B | 8 (48698.1175) | 48743.8798 | 48743.879795 | No | — | 5500 | 7890.159 | 1.0000 | Reference | — |
| n10 | Y12 | 8 (48736.1673) | 48746.3743 | 48746.374299 | No | — | 5500 | 7146.577 | 0.9058 | pass/pass/pass | 0.9058 |
| n10 | Y6 | 3 (48720.0758) | 48729.2073 | 48729.207327 | No | — | 3700 | 4564.515 | 0.5785 | FAIL/pass/pass | — |
| n10 | Y8s | 8 (48736.1673) | 48743.9844 | 48743.984403 | No | — | 3400 | 4174.485 | 0.5291 | pass/pass/pass | 0.5291 |
| n10 | Y12s | 12 (48740.7732) | 48748.3537 | 48748.353745 | No | — | 3800 | 5155.926 | 0.6535 | pass/pass/pass | 0.6535 |

## Scope and timing limits

Five nodes were held out from the preceding full-phase A/B; node 13 is a calibration repeat. All 30 runs were fresh, once per node and arm, under the predeclared rotating/shuffled queue. All five uncut controls match exactly at logged precision on each node. The earlier node-13 B and Y12 final bounds are reproduced exactly. Fresh B values, rather than the slightly different quoted step-7 archival values, define every gate here.

These are single observations on an eight-CPU, 8 GiB cloud container (Intel Xeon Platinum 8573C, OpenBLAS 0.3.30, one BLAS thread). Core assignment, shared bandwidth and decreasing tail concurrency affect timing. For example, n9/Y6, n9/Y8s and n9/Y12 all have identical logged bounds and 810 iterations, but cost 916.781, 1,015.247 and 1,114.989 CPU seconds: a 21.6% max/min spread. The small CPU gains of full-budget Y12 on nodes 13 and d2 should not be treated as established speedups. Y12s's aggregate iteration reduction closely tracks its CPU reduction, which is stronger evidence for the direction of its benefit, while the exact factor still needs replication in search.

The comparison charges the initial 400-iteration solve as requested. Production may enter with a different warm state or gate the quiet solve differently. A subtree A/B must retain the same certificate/fathom safeguards and charge the complete production policy, including any change in search work. No parameters were tuned after seeing this matrix, and no extra solver repetitions were run.

## Reproducibility and retained evidence

- [Predeclared protocol](PROTOCOL.md) and [exact command queue](plan.json).
- [Manifest](results/runs/manifest.json): source/binary/input/plan hashes, all commands, core placements, whole-process CPU, wall time, RSS, stage values and outcomes.
- [All 30 raw logs and event stream](results/runs/), [machine-readable summary](results/runs/summary.json), and [SHA-256 checksums](results/runs/SHA256SUMS).
- [Build settings](results/build.txt), [host metadata](results/host.json), [continuation regression output](results/preflight.log), [parser/gate checks](results/analysis_checks.log), and [final validation record](results/validation.json).
- [Runner](run.py), [gate evaluator](summarize.py), and [analysis regression checks](check_analysis.py).

Every raw log was reparsed against the final manifest; all 30 match and exit zero. Exact commands match the prelaunch queue, binary/plan/input hashes match, and the event stream confirms no concurrent reuse of a core. The continuation regression and parser/gate checks passed before launch. Solver source is unchanged from the cognac base.

To reproduce, build with OpenBLAS using the settings appropriate to the host (the exact cloud linkage is in `results/build.txt`), then use a fresh output directory:

```bash
make -C latest dnncuts
python3 experiments/same_node_continuation/budget_study/run.py --out /tmp/qap-budget-repeat --workers 8
python3 experiments/same_node_continuation/budget_study/summarize.py /tmp/qap-budget-repeat
```

The runner verifies the node specifications, predeclared queue, CPU quota and physical-core count, pins each process, sets BLAS thread limits, and enforces a six-hour wall limit. Lower `--workers` if fewer physical cores or insufficient memory are available; a timed-out partial matrix cannot qualify an all-node winner.
