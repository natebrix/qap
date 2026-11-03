# Continuation budget study: predeclared protocol

Based on cognac `4095416756ef487251912396a727fac7ceacacaf`. No solver implementation or production default changes. This protocol and the exact command queue in [plan.json](plan.json) are committed before launching the matrix.

## Question and sample

Can a shorter Y/Z cut phase preserve the baseline's certified value and fathoming outcomes while reducing total CPU? The roadmap assigns 46% of dominant-class CPU to cut phases.

| Node | Zero-based assignments | Expected baseline class |
|---|---|---|
| 14 | `47:23,10:44,40:34` | Near miss: fathomed |
| 6 | `20:12,40:28,17:25` | Near miss: fathomed |
| 9 | `26:27,8:38,34:28` | Near miss: fathomed |
| 13 | `21:28,8:41,28:19` | Leaves short; calibration repeat |
| d2 | `39:27,30:48,25:49,26:5` | Leaves short |
| 10 | `39:27,30:48,25:49` | Leaves short |

All32 line numbers and depth4 row 2 are checked against the supplied specifications before launch. Node 13 was used in the earlier experiment; it is a calibration repeat, not held out. The other five nodes were not used in that A/B. Every node/arm is run fresh once, including node 13, to avoid mixing historical CPU measurements into the gates.

## Arms

| Arm | Initialization | Cut rounds × iterations | Quiet cap | Total iteration cap |
|---|---|---:|---:|---:|
| B | Z-only | 12 × 300 | 1,500 | 5,500 |
| Y12 | Y/Z | 12 × 300 | 1,500 | 5,500 |
| Y6 | Y/Z | 6 × 300 | 1,500 | 3,700 |
| Y8s | Y/Z | 8 × 300 | 600 | 3,400 |
| Y12s | Y/Z | 12 × 200 | 1,000 | 3,800 |

All arms start from an empty cut pool and 400 cold uncut iterations. Common options are `-grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -delta 50 -nocap`; `-i`, `-r`, `-freeze 0:N` and `-warmY` vary only as declared above. Remaining driver defaults are unchanged. The cutoff is 48,814.05. The existing deterministic initialization and triangle separator are used; no extra randomized cuts or seed settings are added.

The library stops immediately at the first certified cutoff crossing. Cumulative iterations through that terminating call therefore give the exact first-crossing iteration on the solver's checkpoint schedule; certification checks occur every ten iterations by default. No unsampled intermediate crossing is claimed.

## Gates and reporting rules

For each candidate, compare to the fresh B run on the same node:

1. The retained certified bound at phase exit (`FINAL best_lb`) must be at least B minus 2. Both the quiet-only returned value and the retained bound are reported. If a run fathoms before the quiet phase, its quiet value is **not run**, and its terminal certified bound is used. This preserves already valid certificates instead of discarding them. The literal B-minus-2 threshold is used even if both arms fathom; bounds are not clipped at the cutoff. Differences within 0.000002 of the threshold are flagged as rounding-ambiguous rather than passed from six-decimal output.
2. The fathom outcome must match B. On nodes 14, 6 and 9, B must itself fathom and the candidate's first-crossing iteration must be strictly smaller. If an expected near-miss baseline does not fathom, that is a baseline outcome mismatch and the iteration gate cannot pass.
3. Whole-process user + system CPU must be strictly below B, including separation, assembly and all quiet work.

**Any certified cutoff crossing on node 13, d2 or 10 is a validity question, never an accepted fathom or a gate pass.** Raw logs remain unchanged. If B crosses on one of these nodes, no candidate comparison on that node can pass. No all-node winner may be declared from incomplete or failed runs.

Best Polyak round means the earliest cut round (1 onward) with the highest reported certified value, excluding the initial uncut solve and quiet phase. The quiet-only value, final retained bound, crossing iteration, total iterations and process CPU are retained separately.

The eligible CPU fraction is candidate CPU / B CPU only when gates 1 and 2 pass; otherwise it is unavailable. A partial passing subset is reported with its baseline-CPU coverage and is never presented as an all-node result. Only a budget passing all three gates on all six nodes qualifies as a universal winner here. Y12 is a control, not a shortened-budget winner.

For a qualifying arm, phase-cost factor = sum(B CPU) / sum(candidate CPU). If its phase ratio is r, the conditional whole-run factor is `1 / (0.54 + 0.46*r)`. This is an Amdahl projection assuming the 46% share and that this small sample's ratio generalizes; it is not a measured subtree speedup.

## Execution and provenance

Eight single-thread processes run at most one per physical core, pinned to eight distinct cores within the container's eight-CPU quota. The prior measured maximum RSS of about 859 MiB per wil50 worker fits the 8 GiB limit for eight workers with headroom for the controller. OpenBLAS and OpenMP thread limits are set to one, as is `-T 1`.

The queue uses five waves, each covering all six nodes with rotating arm assignments; within-wave order is shuffled with the predeclared scheduling seed 20261102. Thus B is not confined to the beginning and shortened arms to the end. Completion-driven scheduling can still change contention and core assignment; single-run CPU ratios remain subject to those limitations. No parameters or queue entries are selected from observed outcomes.

The previous full wil50 runs cost about 7,000 CPU seconds each on this container, substantially longer than the supplied 15–40 minute estimate. Early fathoming on nodes 6 and 9 should reduce total cost. The expected matrix time is roughly 4–6 hours with eight workers; the runner preserves partial logs and stops at six hours if necessary. A partial matrix cannot answer the universal-winner question.

The manifest records the source commit, binary hash, plan hash, commands, core assignments, start times, whole-process CPU, maximum RSS, stage observations and outcomes. A JSONL event stream provides coarse wall times of observed completed rounds. Those observation times are polling timestamps, not exact crossing wall times. Every started run retains its complete raw log.

Build with OpenBLAS and run:

```bash
make -C latest dnncuts
python3 experiments/same_node_continuation/budget_study/run.py --out /tmp/qap-budget-study --workers 8
python3 experiments/same_node_continuation/budget_study/summarize.py /tmp/qap-budget-study
```

Use a fresh output directory. The runner requires Linux affinity and resource-accounting APIs. The exact cloud build flags are retained in `results/build.txt`. The existing continuation regression is rerun before launch, with output in `results/preflight.log`.

## Results

Completed: all 30 runs finished successfully. **Y12s alone among the shortened budgets passed all three gates on all six nodes**, using 64.17% of baseline CPU (1.5583x phase-cost factor). The conditional 46%-share projection is 1.1973x across all six nodes, or 1.1656x on the three short-of-cutoff nodes. No unexpected cutoff crossing occurred. See [RESULTS.md](RESULTS.md) for the complete matrix, interpretation and evidence. The protocol above is unchanged from its prelaunch commit.
