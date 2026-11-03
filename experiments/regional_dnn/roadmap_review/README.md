# Can wil50 or sko49 fit into a 1,000-core week?

Research synthesis for Nathan Brixius, Fable, and Mythos, 2 November 2026.

**Source snapshot:** cognac commit `ca50467835d640a252f1f4d687618f593a395d7c`. Three Astra agents independently examined bound strength, tree/search policy, and numerical speed, then exchanged validity checks. This is a research review with small local experiments, not a production solver change or a demonstrated 4–12x improvement. No cloud fleet was launched.

## Conclusion and budget

There is a credible experimental path to the low end of the one-week target, but no measured basis yet to promise the high end. The best next work is to reuse expensive information and make better decisions with existing certificates before investing in a full new SDP implementation. Keep cross-region consistency as a genuinely different strengthening experiment.

The stated post-cut estimate of 2.2e9–7e9 reference-core-seconds remains a planning range, not a measured full-tree total. At 1,000 physical cores for seven days, capacity is 604.8 million core-seconds at perfect utilization, or **574.56 million at 95% utilization**. Thus the additional reduction required is **3.83–12.18x**. The already measured 1.81x benefit from cuts must not be counted again.

| Additional end-to-end reduction | Wall days at low estimate, 95% utilization | Wall days at high estimate, 95% utilization |
|---|---:|---:|
| 2x | 13.4 | 42.6 |
| 4x | 6.7 | 21.3 |
| 8x | 3.4 | 10.7 |
| 12.2x | 2.2 | 7.0 |

These are arithmetic scenarios, not predicted speedups. If an entire revised numerical implementation proves 2x faster on the target fleet, remaining algorithmic/search work must improve by approximately 1.91–6.09x. Measure the combined implementation: do not multiply gains that save the same iterations or descendants.

## What changed since the previous report

- The completed node-10 subtree A/B gives 139,149 -> 76,755 CPU seconds and 2,193 -> 365 nodes: **1.81x less CPU, 6x fewer nodes**. The gap between these factors is a central constraint on further work.
- Triple-consistency trials were negative. Five-point counts were mostly neutral. Slow30 dual optimization added only 5–8 units under the tested settings. Local-7 cuts helped at fixed multipliers on uncut iterates but lost 15–25 units when mixed into the tested triangle protocol.
- Those results justify lowering priorities, not declaring the exact larger relaxations exhausted. A lower-bound plateau is not an upper certificate on the optimized cut-pool value; even an exact primal solution of a frozen Lagrangian bounds that fixed-multiplier problem, not the best multiplier choice.
- Two existing roadmap memos already propose inheritance and CPU acceleration. The recommendations below refine those and add same-node continuation, full-certificate conditional child screening, tail-aware branching diagnostics, and cross-region categorical consistency.

## Ranked experiments

### 1. Reuse solver state within a node, and cut information between nodes

These are different experiments and should be tested separately first.

**Same-node continuation:** `dnn_bound` reloads Z but initializes Y afresh on each call; the cut rounds and quiet round call it repeatedly. Preserve compatible Y and Z and, in a later step, the spectral workspace across rounds. When the objective scale changes, rescale dual state correctly; preserve the best round's complete state for its quiet round. A negative test of a parent's Y projected into a child does not answer this same-node question. Small-instance evidence is in the speed memo; wil50/sko49 benefit remains unmeasured.

**Parent-to-child inheritance:** restrict the parent's selected cut pool and multipliers to each child, with exact fixed-variable substitution and constant accounting. A triangle touching the fixed membership reduces to a redundant nonpositive pair-mass penalty and may be dropped with a documented strengthening argument. Do not generalize that rule to arbitrary higher-order cuts. Warm-start the matching modified objective, and carry the required metadata across steal/checkpoint boundaries. The existing branch's state-version issues show why this is part of correctness, not just serialization convenience.

First gate: same-node paired continuation on wil50 nodes 13/d2 and sko49 children 11/13; separately compare inherited frozen pools against both child uncut and fresh child-cut controls on siblings of a failed phase. Move to complete subtrees only after logging time to certificate, additional fathoms, and failures. No inherited-gain factor is established yet.

### 2. Screen children with the full parent certificate

The current outer-assignment reduced-cost child bounds are not the only information available. Condition the full polyhedral part of the parent's certified objective on a candidate assignment, tightening the inner assignments as well. Reuse the parent's verified spectral allowance in the parent embedding. This potentially yields stronger child bounds without child SDP iterations.

Use the actual multiplier/slack that produced the best certificate, not merely the latest warm-start Z; preserve objective shifts, cut constants, and arithmetic allowances. The bound memo and tree memo give the proof obligations and implementation distinctions. Compare against existing child bounds by certified eliminations per second. A useful extension is to screen two candidate branching lines, then spend short DNN probes only on their unresolved expensive children.

### 3. Branch to reduce expensive surviving subtrees

Average child bounds need not predict total work: a few difficult children can dominate. Score candidate lines using calibrated survivor-cost estimates or a tail-sensitive criterion, then measure completed subtree CPU. Do not use a smooth formula such as 47^(gain/130) as an established conversion from bound units to runtime.

The cached child-bound analysis shows that five children carry 76–92% of the best-line surrogate cost. However, simple tail-QPB scores did not improve the cached selections; the useful new test is conditional certificates or limited actual DNN probes, not a scoring-only relabeling. These are surrogate analyses, not measured search savings. One sko49 sample contributes 80% of the +130 model total and two wil50 samples contribute 93%, making direct hard-class sampling essential. The `-F` subtree protocol also resets local depth and forces root rule 4; it does not exactly replay a full-search global-depth-three branching schedule. Keep the original A/B as valid evidence for its own protocol, and add a global-depth-faithful replay before extrapolating.

### 4. Remove numerical work, then benchmark accelerated kernels honestly

Prioritize low-risk skipped work and full same-node continuation before a large port. The speed audit identifies residual accumulation that can be restricted to check iterations when acceleration does not consume it, plus previously proposed redundant-copy/transform work. Microkernel speedups must be converted through whole-node timing.

AMX or GPU kernels remain possibilities, not established throughput gains. The current research environment does not expose AMX CPU flags, despite an older memo's environment report. Validate instruction availability on the intended fleet. Filter precision tests do not automatically validate all Rayleigh–Ritz or orthogonalization operations. Keep a certified numerical path and measure both whole-process CPU and wall time: the DNN's thread CPU timer alone omits helper-thread/device work.

In the measured subtree, roughly 54% of CPU is ordinary DNN bounding and 46% cut-phase/other work. Doubling only the latter buys about 1.30x overall; even deleting it entirely buys at most about 1.85x if all other work stays unchanged. GPU speedups change hardware-equivalent wall time; they are not automatically reductions in the reference-core algorithmic work or cost.

### 5. Test consistency across regions, not just within one region

For nested thresholds S subset T, each facility has three labels: S, T minus S, and outside T. For five or six facilities, enumerate the capacity-feasible label configurations (at most 243 or 729). Their unary and pair moments are linear aggregates of the existing Y. Separating their joint marginal polytope enforces a single common distribution across two threshold memberships, unlike independent within-region local-7 tests.

A row/column pair similarly gives four labels, with at most 1,024 configurations for five facilities. The existing small-polytope separator architecture is a starting point. Use exact or outward-rounded support maxima over admissible configurations for cut validity, consistent coefficient normalization, and the original node indexing/capacities.

Gate this on the triangle-capped quiet iterate. Freeze the existing triangle multipliers, introduce new cuts along a small nonnegative multiplier ray, and retain the old certificate as fallback. This avoids confusing damage to the old multipliers with the value of the new family. A negative measured result is a reason to stop; no positive gain is currently demonstrated.

### 6. Retain the global regional moment SDP as a gated research option

Section 10b of the main report still applies: separation first, independently verified cut second, incremental certified gain third. Full linked integration follows only if it beats repeated offline separation on useful work. Failure of selected local families does not prove the moment extension redundant, but neither does a violated moment extension establish a useful objective gain.

## Where sko49 enters

Today's sko49 depth-four cut logs at the inspected snapshot leave children 13, 9, 11, and 10 respectively 3.55, 20.27, 28.49, and 43.27 units short of cutoff 23,384.05. Child 13 is a sensitive diagnostic; child 11 or 10 is a stronger substantive target. A cheap gain below 20 can therefore be useful even though it would fail a blanket wil50 threshold. All four children share a parent, so expand beyond that parent before extrapolating.

The distance matrix is a 7x7 Manhattan grid, with locations in the first matrix. If using the transposed instance for the solver's symmetry machinery, change the location axis and assignment mapping consistently. Symmetry is already included in some baseline estimates; do not count the factor eight again.

## Evidence needed before committing a fleet

1. Pin instance, code, settings, global-depth semantics, and a reproducible snapshot of cut and solver state.
2. Run small correctness checks for objective transfer and new certificates, followed by paired node-level time-to-certificate tests.
3. Test complete representative subtrees: include the harder wil50 node-1/node-5 class, not only the successful node-10 class, and multiple sko49 parents.
4. Measure a combined candidate on a disjoint evaluation sample. Record bounds, fathoms, branches, CPU, wall time, memory, and restart overhead. Preserve failures and timeouts in the accounting.
5. Validate target-host throughput under realistic process counts and publish a sampling-based total-work estimate with uncertainty. At 95% utilization, the conservative work estimate must fit 574.56 million target-core-seconds for a one-week claim.

The promising path is **less repeated optimization + stronger cheap decisions + fewer expensive survivors + faster verified arithmetic**. The low-end target is plausible enough to justify these gates. The upper-end target remains an open research outcome, not a forecast supported by today's measurements.

## Detailed investigations

- [Bound strength, inheritance, and cross-region consistency](bound_research_findings.md)
- [Tree policy and conditional certificate screening](tree_research_findings.md)
- [Numerical speed and same-node continuation](speed_research_findings.md)

Each memo distinguishes repository measurements, new small experiments, mathematical arguments, and hypotheses. Supporting experiments are included where available; none is a production solver change.

Reproduction commands: [REPRODUCE.md](REPRODUCE.md).
