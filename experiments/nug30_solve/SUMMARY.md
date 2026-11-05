# nug30 on the `cognac` branch: summary, Oct 25 – Nov 5, 2026

*Written Nov 5 by a subagent from the conversation transcript and the repository; numbers checked against the logs.*

**Scope.** This covers nug30, and the general solver work as it bears on nug30. The wil50 and sko49 studies are left out unless a finding was also measured on nug30.

**Conventions.**
- Numbers are **measured** unless tagged *(est.)*.
- Sources are files (paths relative to `/home/user/qap`) or commits on `cognac`.
- Prompt times are UTC. Laptop run dates follow the tuning README.
- AI systems are named by role:
  - "the assistant" is the cloud session;
  - "the laptop session" is an assistant session on Nathan's Apple M2 Pro that pushed results to `cognac`;
  - "subagent";
  - "an external AI reviewer" is a separate system Nathan consulted.

## Overview

- **What nug30 is.** A QAPLIB instance (Nugent, Vollmann and Ruml, 1968): n = 30, locations on a 5×6 grid, optimum **6124**.
- **The 2000 solve.** Anstreicher, Brixius, Goux and Linderoth proved the optimum with 1.19×10¹⁰ nodes. It took about 6.9 days on an average of 653 grid workstations, roughly 7 HP9000 C3000 CPU-years. Sources: `experiments/nug30_solve/nug30_2000_levels.txt` (from Table 5.12 of Nathan's thesis) and `docs/dnn_bound/dnn_bound.tex`.
- **The goal.** Nathan, Oct 28: "It's been, amazingly, 26 years since we solved nug30. I'd like to see what's possible now." The work started from Nathan's QPB branch-and-bound code and a list of ideas from Jean-Pierre Goux, Nathan's collaborator on the 2000 solve.
- **The method.**
  - QPB at every node.
  - A certified DNN bound on the nodes QPB leaves open at depths 1–7, used only to fathom.
  - Geometric branching rules, later branching on the DNN marginals.
  - Warm starts and a checkpointed parallel driver.
  - Every run starts from the known optimum (`--ub 6124`), as in 2000.
- **Best result (Nov 4).** Apple M2 Pro laptop, 8 workers: **1,838 CPU s, 322 s wall (5 min 22 s), 3,844 nodes**. 3,843 of those nodes got the DNN bound and 3,701 were fathomed by it.
  - Log: `experiments/wil50_feasibility/tuning/refine/nug30_best_nov4_ship32_j8_m2pro.log`, commit `ff7d517`.
  - Parameter file and range: `latest/PARAM/nug_deep_xbar_all.par`, `-x 1:7:1:200:1`.
  - Switches: `QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32`.
- **Best cloud result (Nov 5).** The record setting on a 4-core 2.1 GHz Xeon: 3,820 nodes, 3,270 CPU s, 14 min 9 s wall (`experiments/nug30_solve/nov5_cloud/`). Before it, the best cloud run was Nov 1's 4,528 nodes, 11,748 CPU s and 50.5 min (`experiments/wil50_feasibility/stronger_bounds/search/nug30_base.log`).
- **Laptop progress.** From Oct 29 to Nov 4, CPU went from 8,355 s to 1,838 s (4.5×) and wall time from 19 min 42 s to 5 min 22 s.

## Timeline of key prompts

| Date | Nathan asked | Outcome |
|---|---|---|
| Oct 25 17:36 | Look at Jean-Pierre's "nug30 Laboratory" page | The page's view: per-node speedups top out near 1.3×; a big gain needs a stronger bound |
| Oct 25 18:33 | "Implement your picks ... brute force and regression checks" | Brute force found 7 kinds of defect in the starting code (113 of 1,008 runs failed); all fixed. LAP warm start: 0.80× time |
| Oct 26 03:29 | nug24 with the latest code vs the baseline | 1.3× faster than the corrected original |
| Oct 27 03:52 | Have a subagent look for bigger ideas (paraphrase) | Frank–Wolfe warm start, an eccentric-location rule, and a DNN bound at the top, "possibly 5–20× at nug30" *(est.)* |
| Oct 27 10:02 | Try the eccentric idea ("net new strategies") | Rule 6; nug24 1.32× faster |
| Oct 27 22:15 | "Better branching up top? Bound improvements?" | Lookahead bug fixed; rule 7; a DNN feasibility study on nug20 |
| Oct 28 03:13 | "I would like to run the nug30 DNN test" | Python study: root DNN 5,949 vs QPB 5,332; a hybrid estimated at 10–100× *(est.)* |
| Oct 28 13:50, 17:39 | Build the C++ bound and nug parameters, then wire DNN into the search ("proceed as you suggest") | `nug_deep.par`, `latest/dnn.C`, DNN in `bnb.C`. **First nug30 solve finished Oct 29 00:16** |
| Oct 29 00:31 | A clean single-setting run, a write-up, an explainer page and verification by subagents ("a landmark problem and a new standard") | Clean run in about 80 min; reviews; `dnn_bound.pdf`; the "nug30, Again" page |
| Oct 29 05:05 | "Any opportunities to further improve the parameter settings for nug30?" | DNN to depth 7, plus DNN warm starts |
| Oct 29 05:49–05:53 | The external AI review; "a re-run of nug30 will be warranted" | Cancellation and warm-start fixes; 5,817 nodes; −17% CPU in the A/B |
| Oct 29 13:06–14:27 | Restructure the paper; thesis pages for the 2000 run | Table 6: nodes per level, 2000 vs 2026 |
| Oct 29 15:03 | Parameter intuitions for the laptop: "I'd like to get a fast run locally" | Notes in order: BLAS, branching quality, H = 200, more warm-start state |
| Oct 29 15:28 | Where do the speedup factors come from? | ~9,000× ≈ 20 (hardware) × 3 (QPB search) × 150 (DNN) *(est.)*. The 2000 run was ~7 CPU-years, not 1 |
| Oct 29 16:24 | "We now have a 20 minute run locally" | Paper Table 5 and the page updated |
| Oct 29 16:47–17:22 | Is DNN used only to fathom? Is QPB a predictor? Run the line-by-line test | Confirmed; AUC 0.994; line test (§4) |
| Oct 29 18:55 | "Are we tapped out on this bound?" | At the root, yes; triangle cuts give +0.2 and +0.4 |
| Oct 31 10:25 | Anything on branching choice itself, "look ahead vs eccentric"? | DNN marginals, which led to the laptop's rule 8 (−16% on nug30) |
| Oct 31 21:29 | Why ours is faster than Fujii et al.'s nug30 with the same relaxation | Certificate valid at every iterate, cheaper iterations, smaller tree; later measured at 230× less compute |
| Oct 31 23:50 | Jean-Pierre's runs: 2 h 35 → 1 h 30 → 44 min 39 s | His D6/D8 ported; his tunings tested |
| Nov 1 00:07, 04:14 | "Port D8 and D6 now"; other Goux ideas, "I approve" | `QAP_DNN_EXT` and `QAP_DNN_GU_ZONE`; β + refresh brought nug30 to 2,895 s |
| Nov 1 10:54 | Next steps "all the way through wiring cut phase" | Cut phase on nug30: +42% (cloud), +34% (laptop); closed |
| Nov 3 21:19 | "Pass the parent's primal solution to its children, not just the dual" | `QAP_DNN_WARM_Y`: 2,156 s |
| Nov 4 | (no prompt: the laptop session's own run) | `QAP_DNN_SHIP=32`: 1,838 s, 322 s |

## Findings

### 1. The DNN bound

- **What it is.** The facially reduced DNN relaxation (Zhao–Karisch–Rendl–Wolkowicz; Oliveira–Wolkowicz–Xu), solved by ADMM. Its certificate is valid at every iterate:
  - it takes the R-step multiplier, which is negative semidefinite on the face by construction;
  - it bounds λmax with a Cholesky test (Rump's criterion);
  - it evaluates the lifted cost by Gilmore–Lawler and subtracts floating-point margins.

  A node is fathomed when the certified value passes 6,122.05; nug30 costs are even. Sources: `experiments/dnn_bound/cpp_report.md`, `docs/dnn_bound/dnn_bound.pdf`.
- **Root.** DNN 5,949 (2.9% gap) against QPB 5,332.33 (12.9%). The root is converged: after 3,000 iterations the certificate is 5,949.4 and the primal value 5,948.9.
- **Python study** (`db8df94`, `experiments/dnn_bound/nug30/report.md`).
  - Depth 2: 34 of 36 sampled nodes fathomed, at about 130 s each. QPB separates the outcomes: misses fall below about 5,548.
  - Over all 225 depth-2 nodes the fathom rate was estimated at 82% *(est.)*.
  - Depth 3: 19 of 20 fathomed.
  - QPB alone was estimated at 30–120 CPU-days *(est.)*.
- **C++ bound** (`a57cc96`).
  - 43 ms against 339 ms per iteration at m = 28.
  - About 6× fewer iterations to certify, and 13–36× less CPU per node.
  - The same fathoming decision as Python on all 162 study nodes.
  - Sources of the gain: products with the face basis done with Householder reflections in O(m⁴), a certificate that needs no extra eigendecomposition, Gilmore–Lawler instead of a box bound, a single-precision subspace eigensolver, and β = 0.3·m/3.
- **Role in the search.** The DNN bound runs after QPB and before branching, and it only fathoms. Branching and child eliminations use QPB and its U matrix. It has two indirect effects:
  - the DNN value raises the bound that children inherit, which steers their strategy row;
  - its multiplier warm-starts the children.
- **Verification.**
  - Subagents reviewed the math and the code (2,592 bounds checked under sanitizers) and ran adversarial tests.
  - The external AI reviewer found two cancellation counterexamples. They are fixed by rounding the data to guaranteed lower bounds (`0a179fc`, `tests/dnn_test.py` part D).
  - None of this touches nug30: its data are integers and every sum stays below 2⁵³.
  - After the fixes, the clean settings reproduced the clean run's tree exactly.
- **QPB predicts the DNN outcome** (`c714294`, all 225 depth-2 nodes). As a fathom classifier, parametric QPB has AUC 0.994 (96.9% correct at its best threshold). The lookahead estimate scores 0.988, GLB 0.969 and PB 0.920.
- **Where the time goes** (Oct 29, `-x 2:7:1:300:1`, 2.8 GHz Xeon). The DNN bound is 99.7% of the CPU, and depths 3–4 account for 76% of it. Time per DNN node falls from 11.9 s at depth 2 to 2.2 s at depths 6–7.
- **2000 and 2026, level by level** (Table 6 of the report; 2026 counts from `experiments/nug30_solve/v4_final.state`):

| Level | 2000 nodes (fathomed) | 2026 nodes (fathomed) |
|---|---|---|
| 2 | 225 (0%) | 225 (73%) |
| 3 | 6,120 (0%) | 1,608 (94%) |
| 4 | 157,947 (1.7%) | 2,679 (99%) |
| 7 | 286,494,566 (61%) | 96 (100%) |
| 10 | 2,660,681,036 (85%) | none |

### 2. Eigensolver and ADMM speed

- **Cost per iteration at m = 28.** 30 ms on the M2 Pro with Accelerate, 40 ms there with scipy-openblas32, and 73 ms on the 2.8 GHz cloud Xeon (`RUNNING.md`). Per DNN node the laptop is about 3× faster than that Xeon.
- **Peaceman–Rachford step (PRSM, γ₁ = γ = 1).**
  - On nug30 node `4:10,6:8` it fathoms in 310–320 iterations instead of 500. Misses converge to the same values, and the certificate is unchanged (`experiments/eigensolver/README.md` §3).
  - With a degree-4 filter, in the full search (`cbc0a05`): 8,081 → 5,334 s (1.51×). DNN nodes fell 22%, and a 530,000-node QPB tail disappeared.
- **Degree 3 with tight λ:** −9% (4,008 → 3,653 s) at equal strength (`285ba15`).
- **Jean-Pierre's β × 0.7 plus a full eigendecomposition every 120 iterations:** −20% (3,632 → 2,895 s), same tree (`141d18d`). The two changes were not measured separately.
- **Smaller items.**
  - One-pass CholQR is part of the Nov 3 bundle and was never measured alone.
  - Block growth changed little on nug30.
  - The Intel AMX kernel cannot run on Apple silicon. On one nug30 node its bound was within 0.015 of the fp32 path.

### 3. Give-up rule and DNN depth range

- **The flag.** `-x D1:D2:gap:H:warm` runs the DNN bound at depths D1–D2. A node is given up after H iterations if its extrapolated bound stays below the cutoff.
- **Depth range.**
  - The first solve used depths 2–5, then 2–6.
  - The clean run used 2–6; the QPB search below the depth-6 misses still took 6.5M nodes.
  - With 2–7, all 96 depth-7 nodes were fathomed for 178 s, and the tree shrank to 5,817 nodes.
  - With 1–7 on the refined setting, the bound costs 36 s on the 9 depth-1 nodes. Its marginals improve rule 8's branching: depth-3 misses fell from 59 to 45.
- **The give-up threshold H.**
  - H = 300 was set from nug25 runs.
  - Under PRSM, H = 200 is clean on nug30: −8% CPU, no new misses (`e0b6baa`).
  - H = 150 is unsafe: extra misses and +4% CPU.
  - The gap-ratio rule at 0.65 gave about 40 false give-ups on nug30 and a 1,040,545-node QPB tail.
- **D8, `QAP_DNN_EXT`** (keep the DNN bound below a last-level miss).
  - On nug30 at H = 200 there is no last-level miss, so the tree is identical.
  - On nug28 it cut the tree from 3,500,289 to 1,634 nodes.
  - A bug dropped the flag when nodes passed through state files; it was fixed Nov 2 (`99d3b30`).
- **D6, `QAP_DNN_GU_ZONE=0.003`** (confirm a give-up near the cutoff). It changed no decision on nug30 and cost 4% CPU (`de2d47e`).

### 4. Branching

- **Rule 6** (the most eccentric free location; costs nothing to compute): nug20 2.1× faster, nug24 1.32×. It loses on the sparse-flow chr and esc instances.
- **Lookahead fix and rule 7.**
  - A bug made lookahead score essentially arbitrary lines. Fixing it gave 10–12% on nug20 and nug21.
  - Rule 7 runs lookahead on the 3 most eccentric locations plus 2 rows chosen by U. It became the default at depths 3–4 (nug24: 0.937× the time).
  - `nug_deep.par` adds rule 6 down to depth 6 (QPB-only nug25: 0.67× the CPU). It was the setting for every nug30 run through Oct 31.
- **Line-by-line test** (`f604947`; 5 depth-2 misses, 1,260 children).
  - At near misses the choice of line makes no difference.
  - At the far miss it does: 6 child misses on the best line against 17 on the worst. QPB ranks the lines with correlation 0.92–0.98.
  - The eccentric line was close to the best.
- **DNN marginals x̄** (`37f8b8f`; one nug30 node, 11 lines). x̄ predicts child gaps as well as child QPB does (+0.745 against −0.750), and it comes free. The eccentric and lookahead lines cost 1.8–2.2× the best line.
- **Rule 8** (`622af94`; branch on the line with the largest Σx̄²).
  - nug30: −16% CPU (4,904 → 4,133 s); depth-3 misses fell from 85 to 59.
  - nug28: **+37%**.
  - Rule 9 hedges between the two: nug30 −13%, nug28 +8%.

### 5. Warm starts

- **QPB Frank–Wolfe from the parent's solution:** nug24 1.16× faster.
- **DNN from the parent's multiplier.**
  - The first integration silently dropped the parent's state; fixed in `0a179fc`.
  - On nug30, iterations per fathom fall about 40%, but time per node only 13–16%.
  - Combined with depth 7: −17% CPU in a same-machine A/B (28,428 → 23,618 s, `f975248`).
- **Scheduling.** A node handed to another worker loses its parent's state, so reruns differ (5,817 against 5,819 nodes). This is the likely cause of the laptop's two depth-7 misses on Oct 29, though it was never tested.
- **The parent's primal Y as well** (`QAP_DNN_WARM_Y=2`). Together with one-pass CholQR and D6/D8: 2,870 → 2,156 s (−25%) (`4ae70b9`).
- **Shipped state** (`QAP_DNN_SHIP=32`, `8d59c8e`). The parent's state travels in the task files, so about 300 tasks start warm: −15% CPU on an identical tree. Parallel efficiency fell from 6.22× to 5.72×.

### 6. Cut phases

- **Atomic triangle cuts:** +0.22 and +0.41 on two depth-2 misses, nothing at the root. At the DNN solution they are violated by only 0.004–0.02 (`experiments/dnn_cuts/README.md`).
- **Regional triangles** on the 5×6 grid of the second matrix. Run standalone, they fathom miss `4:7,6:14` in 4 rounds and lift `4:29,6:7` by 25.9 of the 26.3 it needs.
- **Inside the search they lose on nug30.**
  - Cloud: +42% CPU (11,748 → 16,723 s), with 26% fewer nodes (`7b88bd0`).
  - Laptop: +34% with cuts at every depth, +1% with cuts at depth ≤ 2 only (`ea171cf`).
  - The reason: a nug30 depth-3 child costs about 0.7 s, so fathoming a depth-2 miss saves about 20 s, and a cut phase costs about 17 s.

### 7. Parallel runs

- **Driver.** `bench/parallel.py` builds a frontier (9 root children after symmetry), runs one process per worker with work stealing, and checkpoints every 600 s (300 s on the laptop). Runs survived 2-hour task limits, container restarts and changes of CPU.
- **Efficiency** (CPU time over wall time). Cloud, 4 workers: ×3.9. Laptop, 8 workers: ×7.1 at first, ×6.2–6.9 later, ×5.72 with ship.
- **Repeatability.**
  - Without warm starts, trees match across machines.
  - With warm starts a few marginal nodes move, and CPU varies about 2% (2,870, 2,895 and 2,943 s for the same setting).
  - The root QPB rounds differently on Apple hardware: 5,336.58 against 5,332.33 on the Xeon.
- **The coordinator** (`coord/`, Nov 4) was tested on a wil50 subtree, not on nug30. The ship feature came out of that work.
- **Comparisons.**
  - Jean-Pierre reports 2 h 35 min (5,550 nodes), then 1 h 30 min, then 44 min 39 s (5,570 nodes), on a 4-core 2.1 GHz Xeon. These are his figures.
  - Fujii et al. (2021, Table 4): 26,181 nodes and 5.4M core-seconds. Against our 2.8 GHz run's 23,600 core-seconds, ours used 230× less compute.

### 8. nug30 records

Cloud Xeon, 4 workers:

| Date | Setting | Nodes | DNN nodes / fathomed | CPU s | Wall | Source |
|---|---|---|---|---|---|---|
| Oct 29 | `nug_deep -x 2:5`, then `2:6`; mixed machines | 45,486,332 | 5,615 / 5,403 | 40,555 | ~3.0 h | `nug30_solve/run1–4.log`, `f8107dd` |
| Oct 29 | `-x 2:6:1:300`, clean, mostly 2.1 GHz | 6,532,419 | 5,764 / 5,551 | 18,614 | ~80 min | `clean_run*.log`, `708d348` |
| Oct 29 | same settings on 2.8 GHz (slower per iteration) | 6,532,419 | same | 28,428 | ~2.1 h | `ab_run*.log` |
| Oct 29 | `-x 2:7:1:300:1`, 2.8 GHz | 5,817 | 5,807 / 5,595 | 23,618 | 1 h 41 min | `v3_run1.log`, `fba0f53` |
| Nov 1 | + PRSM and degree 4, 2.1 GHz | 4,528 | 4,518 / 4,355 | 11,748 | 50.5 min | `search/nug30_base.log` |

Laptop (M2 Pro, 8 workers, Accelerate). Each row adds one change to the row above. Logs are in `experiments/wil50_feasibility/tuning/` unless noted.

| Date | Change | DNN nodes / fathomed | CPU s | Wall | Source |
|---|---|---|---|---|---|
| Oct 29 | `-x 2:7:1:300:1` (585,550 nodes) | 5,834 / 5,619 | 8,355 | 19:42 | `nug30_solve/laptop_run1.log` |
| Oct 30 | PRSM + degree 4 | 4,545 / 4,381 | 5,334 | 13:21 | `eigensolver/nug30_prsm_j8_m2pro.log` |
| Oct 31 | H = 200 | 4,519 / 4,356 | 4,904 | 11:17 | `giveup/nug30_ab/` |
| Oct 31 | rule 8 | 3,935 / 3,796 | 4,133 | 9:56 | `lines/nug30_xbar_H200_j8_m2pro.log` |
| Nov 1 | rule 8 at every node | 3,936 / 3,797 | 4,008 | 9:27 | `refine/nug30_xall_c4*` |
| Nov 1 | degree 3 + tight λ | 3,964 / 3,824 | 3,653 | 8:46 | `refine/nug30_xall_c3t2*` |
| Nov 1 | DNN from depth 1 | 3,842 / 3,700 | 3,632 | 8:43 | `285ba15` |
| Nov 1 | β 0.21 + refresh 120 | 3,820 / 3,679 | 2,895 | 6:59 | `141d18d` |
| Nov 3 | warm Y + one-pass CholQR + D6/D8 | 3,843 / 3,701 | 2,156 | 5:46 | `4ae70b9` |
| Nov 4 | shipped state | 3,843 / 3,701 | **1,838** | **5:22** | `ff7d517` |

On Oct 30, an Apple M1 desktop with the Oct 29 settings took 18,383 s and 42 min 54 s (`nug30_solve/m1_run1.log`).

## Tried and closed for nug30

| Idea | Why it was closed | Source |
|---|---|---|
| Shortest-path LAP | No gain: the existing LAP already ends with a shortest-path search | Oct 25 |
| Shorter Frank–Wolfe ("patience") | Neutral; switched off | Oct 26 |
| Per-node diagonal shift; face-restricted lookahead bounds | The shift cost more time than it saved on nug20; the face bounds were weaker | `experiments/diag_shift/`, `face_bound/` |
| DNN reduced-cost child bounds | 0 of 784 children fathomed; the certificate's outer assignment problem is flat | `0a179fc` |
| Atomic triangles; region branching | Triangles add only +0.2 to +0.4. A region node is no stronger than its weakest cell (530 s against 440 s) | `dnn_cuts/`, `region_branching/` |
| Level-2 RLT | Weaker: 5,770–5,846 against DNN's 5,944 (published values) | `wil50_feasibility/NOTES.md` |
| Momentum on top of PRSM | Worse at every cap (310 → 350 iterations), plus one false give-up | `4498650` |
| Gap-ratio 0.65; H = 150 | False give-ups | tuning README |
| Rule 9 | Not adopted for nug30: rule 8 does better there (−16% against −13%) | `b01f8f8` |
| Cut phase in the search | +42% (cloud) and +34% (laptop) CPU | §6 |

## Discrepancies and unverified points

- **RLT-2 on nug30.** On Oct 30 the conversation quoted about 5,934. The later, verified `NOTES.md` says 5,934 is an SDP value and gives RLT-2 as 5,770–5,846. This summary uses the later figure.
- **"Branching at the DNN depths is not a lever"** (`RUNNING.md`, based on a nug28 test). The later rule 8 result on nug30 (−16%) supersedes it there, though rule 8 lost on nug28.
- **DNN at depth 1.** The Oct 29 notes listed it under "Things that don't help". The Nov 1 measurement made it part of the best setting, through better branching.
- **Stale documents (fixed Nov 5).** When this summary was written, `experiments/nug30_solve/README.md`, `RUNNING.md`, `dnn_bound.pdf` and the "nug30, Again" page all stopped at Oct 29–30, and the later records were only in the eigensolver and tuning READMEs. All four were brought up to the Nov 4 record on Nov 5.
- **Fujii et al.'s node count.** 26,297 in their Table 5 ("generated" nodes, read from a screenshot) against 26,181 in their Table 4.
- **Estimates and reported figures.**
  - The 9,000× breakdown is an estimate.
  - So is the ~53 min the assistant estimated for our PRSM setting, for comparison with Jean-Pierre's run. His timings are his own reports.
  - On Oct 29 the assistant asked whether the laptop run used commit `2d7579c`; I found no confirmation.
- **Replication.** None of the laptop records was repeated, and single-run gains of a few percent are close to the ~2% run-to-run noise.

## Open items for nug30

- Re-run the recommended setting on nug28, where rule 8 lost 37%.
- Measure separately the changes that were only tested together: β against refresh, and warm Y against one-pass CholQR and D6/D8.
- Repeat the ship run.
- Re-fit the gap-ratio rule on nug30 traces.
- Test the certificate polish (`DNN_CERTX`) on nug30; it was only tested, and closed, on wil50.
- Two untried branching ideas: score lookahead lines by how many children clear a QPB threshold, and add rule 8's lines to lookahead's candidates at depth ≥ 3.
- ~~Bring the paper, `RUNNING.md`, `nug30_solve/README.md` and the "nug30, Again" page up to the Nov 4 record.~~ Done Nov 5.
- Jean-Pierre's next challenge, "Solving it in one minute", has no plan yet.

## Where to find things

- `experiments/nug30_solve/`: the README and logs for the Oct 28–30 runs, the 2000 per-level file, and `v4_final.state`.
- `experiments/wil50_feasibility/tuning/README.md`: the laptop nug30 records from Oct 31 to Nov 4 (logs in `giveup/nug30_ab/`, `lines/`, `refine/`, `goux/`).
- `experiments/eigensolver/README.md`: PRSM and the degree-4 filter (§3, §8).
- `experiments/wil50_feasibility/stronger_bounds/`: the cloud cut-phase A/B (`README.md`, `search/nug30_*.log`).
- `experiments/dnn_bound/`: `nug30/report.md`, `cpp_report.md`, `predictor/`.
- `docs/dnn_bound/dnn_bound.{tex,pdf}`: the technical report (Tables 4–6).
- `RUNNING.md`: the run table and the list of `QAP_DNN_*` switches.
- Code: `latest/dnn.C`, `bnb.C`, `branch.C`, `dnncut.{h,C}`, `latest/PARAM/nug_deep*.par`, `bench/parallel.py`.
- Pages: "nug30, Again" https://claude.ai/artifact/GS5uCJbhjczESdH5BYz7i6 and the DNN explainer https://claude.ai/artifact/RaC3rRG3VPkwihVXSgtmKz.
