# Bounds experiment results ledger

The completed screens have not demonstrated a material whole-solver improvement. Coarse cuts lost at equal CPU on both measured tho40 nodes and tai60b's depth-2 branch. Regional equalities reached a native warm tho40 depth-5 cutoff 24 iterations sooner, while the warm depth-3 run reached the same stopping check with a gain of just 1.66 bound units. CP scaling reached the same certified target at the same check on its completed replays.

## Certified node results

All bounds and cutoffs are in the instance's raw objective units. CP CPU is the solver's single-thread CPU for the cut phase. Candidate overhead is added explicitly where measured. Each A/B pair used the same reduced instance and cutoff. Cold seeds, same-node restarts, and native inherited replays have different initial states, so compare only within each row.

| Experiment and protocol | Node, budget, initial bound | Certified base → candidate | CPU and cutoff result |
|---|---|---|---|
| Coarse, same-node continuation | tho40t `26:10,8:0,2:32,6:39,33:24`, 500 CP iterations after shared checkpoint 239,908.5305, 2,265 regional cuts | 240,346.0100 → 240,352.8080, **+6.80** at 500. **−29.70** at regional CPU budget | 57.90 → 66.30 s CP, plus 2.98 s separation and 0.55 s setup. Last candidate check in budget was 410 at 57.83 s including overhead. Cutoff 240,514.05, neither fathomed. |
| Coarse, same-node continuation | tho40t `26:10,8:0,2:32`, 500 CP iterations after shared checkpoint 237,829.1402, 2,400 regional cuts | 238,341.0981 → 238,350.4487, **+9.35** at 500. **−24.73** at regional CPU budget | 82.60 → 90.00 s CP, plus 3.64 s separation and 0.37 s setup. Last candidate check in budget was 430 at 81.34 s including overhead. Cutoff 240,514.05, neither fathomed. |
| Coarse order six, filtered same-node continuation | Same tho40t depth-5 checkpoint, 12 cuts retained after binary-union screening, 500 CP iterations | 240,346.0100 → 240,346.7646, **+0.75** at 500. **−60.81** at the reused regional CPU budget | 57.90 s regional baseline reused. Candidate spent 65.20 s CP, 10.19 s separation and screening, and 0.29 s setup. Last check in budget was 360. Neither fathomed. |
| Coarse, same-node continuation | tai60b `44:45,43:58`, 120 CP iterations after shared checkpoint 585,899,572.8148, 3,000 regional cuts | 594,804,614.4795 → 595,274,673.4528, **+470,058.97** at 120. **−297,998.81** at regional CPU budget | 188.90 → 200.20 s CP, plus 1.41 s separation and 1.63 s setup. Last candidate check in budget was 110 at 184.30 s including overhead. Cutoff 608,215,054.05, neither fathomed. |
| Coarse, late same-node continuation | Same tai60b checkpoint and 40 cuts, 600 CP iterations, no retuning | 603,671,585.5294 → 603,777,531.9516, **+105,946.42** at 600. **−51,004.53** at regional CPU budget | 987.10 → 1,032.20 s CP, plus 1.41 s separation and 2.16 s setup. Last candidate check in budget was 570 at 985.28 s including overhead. Cutoff 608,215,054.05, neither fathomed. |
| Regional equalities, cold fixed-budget diagnostic | Same tho40t depth-5 node, 200 cold DNN iterations then 400 CP, initial 237,392.2812, no inherited pool | 239,890.3934 → 239,917.5659, **+27.17** at 400 | 51.15 → 49.76 s CP in development timing. Both hit the cap below 240,514.05. The 239,800 target was first reached at check 360 → 350. Five earlier targets reached the same check. |
| Regional equalities, native warm production-recipe replay | Same tho40t depth-3 node, reconstructed parent state and 1,482 retained cuts, initial 237,871.8951 at 110 ordinary iterations | 238,577.8109 → 238,579.4667, **+1.66** at 690. At check 400 the gain was +1.50. | 109.99 → 114.56 s CP and 128.80 → 132.46 s process CPU. Both handed over at 690 and missed cutoff 240,514.05 by 1,936.24 and 1,934.58. |
| Regional equalities, native warm production-recipe replay | Same tho40t depth-5 node, reconstructed parent state and 1,718 retained cuts, initial 240,121.4553 at 110 ordinary iterations | Both fathomed cutoff 240,514.05. Final certificates were 240,514.1271 at check **687** and 240,514.2926 at check **663**. | 80.07 → 72.91 s CP and 91.74 → 85.80 s process CPU. Observed reductions were 3.49% in iterations, 8.94% in phase CPU and 6.48% in process CPU. |

The initial coarse arms added 40 validated cuts at zero multipliers while preserving the initial regional pool. Both arms refreshed regional triangles every 100 iterations. They ran in `dnncuts` using the production region and numerical settings, with a fresh same-node schedule, explicit stopping and no child inheritance. At tho40 depth five, 25 of the 40 candidates were implied within the local exact binary-union LP control. Four had positive numerical PSD witness violations in the stronger local control. In the order-six follow-up, 8 of 20 cuts were screened out and the numerical search found no positive PSD witness. The tai60b test had three nonempty labels with capacities 14, 8 and 36. It reached every distinct regional-arm certified target in more CPU, including setup and separation.

The long tai60b continuation used concurrent one-thread arms. Its fixed-iteration certificate gain peaked at 476,320.08 at check 130 and fell to 105,946.42 at check 600. The coarse arm reached the regional terminal target at check 580, twenty iterations earlier, taking 1,002.24 s through certification including overhead versus 986.79 s on the regional arm. Both stopped more than 4.4 million below the cutoff. Each executed 60 finite certificates with zero failures.

An affine-normalization audit on the same 40 depth-5 cuts found original-to-minimum lifted Frobenius norm ratios of 1.0165 median and 1.0471 maximum. The dyadic round trip and exact RHS recomputation gave nearly the same ratios. This rules out a large norm reduction from local affine reparameterization of the tested pool. [Normalization audit](verification/coarse_norm.json)

Each equality warm pair used the same executable, serialized input, parent warm blob and recipe. At depth three the enabled arm reached the off arm's terminal bound at check 680 instead of 690. The actual production stop still took 690 iterations on both arms, and the trace has only total phase CPU. At depth five the unchanged first solve varied from 11.45 to 12.74 s CPU. The reduction in cutoff iterations is a useful near-tight signal, with timing from one pair.

## CP scaling

The native node replay used `19:12,18:26`, one reconstructed parent warm state, 833 retained cuts, tai35b comb regions, and the tai50b recipe modified for this instance and a 200-iteration phase cap. The ordinary solve was 110 iterations with bound 278,672,093.3741 on every arm. Scaling ran throughout each CP call.

| CP exponent | Certificate at check 200 | Gain over off | CP CPU to common target | First and CP CPU |
|---:|---:|---:|---:|---:|
| Off | 280,065,070.2172 | — | 16.28 s | 25.62 s |
| .25 | 280,066,383.9761 | +1,313.76 | 16.04 s | 24.07 s |
| .50 | 280,067,658.0827 | +2,587.87 | 16.47 s | 25.01 s |

The common target was the off arm's check-200 certificate. All arms first reached it at check 200. The best candidate still missed the true cutoff of 283,315,445.05 by over 3.24 million. The small CPU differences are below the 15–20% promotion criterion. Even the unchanged ordinary solve varied from 8.03 to 9.34 s.

The native warm tho40 depth-3 scaling screen reused the equality pair's parent state with a 200-iteration cap. The stock arm reached 238,277.7695 at check 200 in 30.44 s of CP CPU. The .25 and .50 arms first reached that target at check 200. Their bounds were 238,282.4542 (+4.68) and 238,280.4580 (+2.69), and their CP CPU was 32.97 and 30.25 s. All three missed cutoff 240,514.05.

On tai60b `44:45,43:58`, a separate same-node scaling continuation started both arms from the coarse experiment's checkpoint at 585,899,572.8148 with 3,000 regional cuts. The .25 arm switched back to the ordinary metric at iteration 31. At check 60 the off certificate was 591,774,673.3183 and the on certificate was 591,901,546.0175 (+126,872.70). Both reached the off arm's check-60 target at check 60 in 90.54 and 88.95 s CPU. At the off arm's check-50 target, the scaled arm needed check 60. Both remained over sixteen million below the actual cutoff 608,215,054.05. This was a short fixed-pool diagnostic, before the first regional refresh at 100. The native cold tai60b baseline was stopped after about 905 seconds of wall time without a concluded node trace. The checkpoint pair supplies a comparison for the shared same-node start.

## Complete subtree and timing scope

Both equality arms closed the warm tho40t subtree `26:10,8:0,3:32` against incumbent 240,516 using the production recipe. The executable and recipe and the input state and warm blob hashes matched. The trace recorded identical ordered node assignments, branch decisions, and fathom causes. The comparison uses fresh trace counts and CPU, excluding the inherited history of one node and 496.267 s.

| Regional equality setting | New nodes | CP phases | Ordinary DNN CPU | CP CPU | New-node CPU | Process CPU |
|---|---:|---:|---:|---:|---:|---:|
| Off | 38 | 2 | 190.15 s | 80.99 s | 272.35 s | 272.65 s |
| On | 38 | 2 | 195.88 s | 80.62 s | 277.68 s | 277.86 s |

The subtree saved 15 CP iterations (541 → 526). Ordinary DNN work grew from 1,270 to 1,436 iterations as the states passed to descendants changed. Measured new-node CPU rose 1.96%. The identical tree and the small net CPU difference leave no evidence of material benefit on this subtree.

- The tai35b arms ran consecutively in single repetitions with machine isolation unverified. The cold equality pair overlapped development activity. The warm depth-3 equality enabled arm overlapped the coarse depth-3 gate. At depth five the equality arms overlapped the tai60b seed, and the off arm also overlapped the order-six coarse run. The two equality subtree arms ran consecutively, with each overlapping a coarse tai60b continuation. The measured CPU values describe development runs.
- The completed `nug12` controls proved 578 on stock, scaling and equality. Stock and scaling each evaluated 27 nodes, popped 10 inherited fathoms, and ran 24 DNN solves with nine CP phases. The equality arm evaluated 31 nodes, popped four inherited fathoms, and ran 28 DNN solves with nine CP phases. This is a small complete-search correctness control.

## Promotion decisions

- **Coarse labels.** The tested pools fail the equal-CPU gate on tho40 and tai60b at both continuation lengths. The tho40 fixed-iteration gains are far below the approximately +500–1,000 objective-unit screening target for costly nodes. The larger early fixed-iteration tai60b gain declined with additional iterations. The selected pools do not justify a subtree test.
- **Regional equalities.** The native warm depth-5 result is a near-tight signal. The slower far-miss phase and the unchanged 38-node subtree with no measured net saving argue against a larger tree trial on the current evidence.
- **CP scaling.** The completed tai35b and tho40 native warm replays and the tai60b same-node experiment showed no earlier crossing of their terminal stock targets or material CPU improvement. A critical-node native inherited replay remains unmeasured. Subtree validation would require a matched warm gain that clears the CPU-to-certificate gate.

Reproduction details are in the [native replay protocol](REPLAY.md) and the [independent certificate checks](verification/README.md). The numerical inputs include the coarse `summary.json` files, the `tai60b_d2_600` continuation traces, the `scaling-tai60b-checkpoint` driver traces, and the native run manifests and DNN traces for `equality-cold-{off,on}`, `equality-warm-{off,on}`, `equality-tree-{off,on}`, `scaling-tai35b-{stock,025,050}`, `scaling-tho40-{stock,025,050}`, and the `nug12` controls.
