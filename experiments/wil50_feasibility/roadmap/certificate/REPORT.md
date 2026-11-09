# Certificate study: more certified bound from the same iterates (Nov 4)

Question: can the certificate (the Gilmore-Lawler / Rump bound computed from a
DNN iterate) extract more lower bound from the iterates we already have?
Running notes are in `NOTES.md` and the proofs in `theory.md`.

Every number is marked **[m]** (measured; the log is named) or **[e]**
(estimate). Logs are in `logs/`. All runs used the production solver flags
(PRSM 1/1, degree 3, tight lambda 2, refresh 40, AMX bundle D, give-up 300) on
the 2.10 GHz AMX host, cores 2-3.

## Verdict

1. **Outside its polyhedral part, the certificate wastes almost nothing.**
   - The semidefinite term plus the rounding margins cost 3.6e-6 to 4.95e-6
     units on every wil50 log. No Cholesky test failed and none needed
     refinement [m].
   - The free dual variables are already handled (`theory.md` section 2):
     - the shift is exactly neutral;
     - the Y00 multiplier (the one Kim-Kojima-Toh bisect on) is blocked,
       because V'Z'V is exactly zero on R's range;
     - with exact LAPs, the assignment-constraint multipliers reduce to one
       remaining freedom, a skew-symmetric transfer.
2. **That skew polish is the one lever.**
   - It is the Hahn-Grant / Adams-Johnson RLT1 dual ascent applied to the
     Lagrangian matrix M = L + Z' instead of the QAP cost. It is valid for
     every lifted permutation (exact identities plus an explicit rounding
     allowance).
   - It raises the certificate at the same iterate on every node tried. **On
     19 of 19 fathoming nodes it crosses the cutoff 1 to 4 checks earlier
     [m].**
   - Iterations saved: 6.7% on cold first solves, 18% on warm-started
     children [m].
   - **In the search, the 46 pooled children of the phase miss [9:45] need
     11.2% fewer iterations and are 1.086x faster in CPU. The whole
     sub-subtree is 1.05x, with the same tree [m: E6].**
   - Estimated effect on the node-10 subtree: **about 1.04x [e]**. That is
     small, but the code is built, validated and opt-in.
3. **The cut-phase lag comes from the multiplier, not from how it is
   evaluated [m].**
   - At round ends the polish recovers only 10-17% of a 47-69-unit lag.
   - Of the round-end lag, about 84% is certificate deficiency and 11% is
     primal overshoot [m: E2b quiet round].
   - No certificate-side trick replaces the quiet round; with the polish it
     ends only about 6% sooner [m]. Closing the lag is a question of ADMM
     convergence.
4. **Free child bounds fathom almost nothing [m].** Complementary slackness
   makes them flat at convergence (proved).
   - At node 10's give-up (229 short): 0 of 47 children cross. The exact
     per-child Gilmore-Lawler bound is still 160-210 short, even though it
     adds a median +45.8 over the parent.
   - At a converged phase miss: 2 of 46 children cross, worth about 1% of
     that child class's CPU.
5. **These are dead [m]:**
   - a line search toward the post-step multiplier;
   - averages of recent multipliers;
   - the projected post-step multiplier (the old prototype certificate, here
     with the exact GL);
   - a skew matrix carried from one check to the next;
   - exact inner LAPs outside the AUTO window;
   - sharper eigenvalue bounds.

## Ranked table

| # | Idea | Effect on wil50 | Cost | Why it is valid | Cheapest test / kill criterion | Status |
|---|---|---|---|---|---|---|
| 1 | **Skew polish of M = L + Z'**, used as the certificate at checks within W ≈ 10 of the cutoff. Cheap version "HGR": 12 row-column-reduction passes plus 1 exact pass. Full version "HG": exact Hungarian passes. | **[m]** At equal iterates, 19/19 fathoming nodes cross earlier. Cold: 350 of 5,240 iterations saved (6.7%). Warm: 60 of 330 (18%). In a fathomed phase's quiet round: about 60 of 1,060 iterations (2% of the phase). Single nodes, CPU including the polish (E3): 7:45 1.10x, 47:45 0.97x, 49:45,0:12 1.05x. Search sub-subtree [9:45] (E6): 1.05x overall, pooled children 1.086x. **[e]** About 1.04x on the node-10 subtree. | HGR13: 0.5 s per evaluation at m = 46 (HG8 1.4 s, HGR49 2.1 s). About 2-5 evaluations per fathoming node, with an early exit once the cutoff is crossed, plus one per cut round. | Each pass is an exact identity for every lifted permutation: block-row LAP reductions, the outer LAP, spreading the leftover over the m-1 used entries, and averaging the (a,b)/(b,a) pairs (yy' is symmetric). So B plus a nonnegative residual equals <M, yy'>. Pass 1 is the exact GL, so certify's margin formula applies; later passes add an explicit 2u allowance. Every value at every check was tested against enumeration (n = 8) and five size-12 QAPLIB instances (tai12b is a sharp case: +615 over the converged certificate, still 4 below the optimum). The default path is unchanged: the plain E3 runs reproduce pool_search's trace exactly. | Node-10 subtree A/B, 4 cores, about 2 h, with `DNN_CERTX=32 DNN_CERTX_USE=1 DNN_CERTX_QUIET=1 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=10`. **Kill below 1.03x**, or on any change other than earlier fathoms. | Built, opt-in (env vars, which bnbdrv reads) |
| 2 | Cheaper LAP passes (warm-started JV or auction, or AMX/GPU), so that 32-50 passes are affordable | **[m]** On 37:45 at iteration 280, 2/4/8/32 exact passes give +4.3/+7.3/+9.4/+12.5. HG32 is 0.12 short of fathoming one check earlier still. **[e]** +1-2% on top of #1. | 1-2 days | Same as #1 | Microbenchmark of a pass; worth doing only if #1 is adopted | Idea |
| 3 | Free child bounds from the parent's certificate: today's reduced-cost `child_lb`, the outer LAP row fixed exactly, or the exact per-child GL (row and column fixed at both levels) | **[m]** At node 10's give-up: 0/47 on the branching line, 0/2,209 by reduced costs. At the phase miss [9:45] after its quiet round: median gain +0.58 over the parent (flat); 2/46 cross, and the search fathomed both in 10 and 30 iterations (13-15 s of 1,043-1,585 s). | About 0.1-0.3 s per child for the exact version | Same proposition as the parent's certificate: every child permutation is a parent permutation. | Done | **Closed** (`child_lb` stays unused by bnb.C) |
| 4 | A different multiplier from the same iterate: line search Zc + t(Z_{k+1} − Zc); averages of the last 4 Zc; the projected post-step Z; a warm-carried skew K | **[m]** All are below the standard certificate where fathoms are decided. At the fathom check: line search −0.7, average −1 to −6, projection −2 to −9. The warm-carried K is −24 on nug20. The average helps only at convergence (+2.4). | 0.3-2.4 s | Concavity, Rump's test and skewness respectively | Done | Closed |
| 5 | Sharper PSD/eigenvalue bounds (Rump/Jansson type); the shift, Y00 and assignment multipliers | **[m]** The margins are 4e-6 units. **Proved** (`theory.md`): the shift is neutral, Y00 is blocked, and the assignment multipliers reduce to the exact GL plus #1. | – | – | – | Closed |
| 6 | Exact inner LAPs at every check | **[m]** Identical to the standard value at every check in the fathom window: AUTO is already exact within 48.8 units. Outside the window they only feed the give-up rule and the Polyak targets. | 0.25 s per check | Valid | – | Closed |
| 7 | Closing the cut-round lag on the certificate side | **[m]** The polish recovers 10-17% of the lag; no round end becomes a fathom. | – | – | – | Closed. It is dual convergence after each Polyak update (the ADMM side). |

## Measurements

### E1-E5 and E1w: variants evaluated alongside the same iterates, without changing the run

Each entry is the first check above the cutoff: standard / HG8 / HGR13
(`logs/e1_table.md`, `logs/e1b_table.md`).

- **Depth-4 cold children of node 10:**
  - 47:45: 450 / 420 / 430. At the last check before the fathom, HG adds
    +3.2 and HGR +3.0 on a 6.7 lag.
  - 7:45: 310 / 280 / 280 (+9.3 / +8.7 on a 9.6 lag).
  - 37:45: 290 / 260 / 270.
  - 20:45: 130 / 120 / 120. 45:45: 150 / 140 / 140. 16:45: 160 / 150 /
    150. 26:45: 160 / 140 / 140.
  - 36:45: 170 / 150 / 160. 6:45: 180 / 160 / 160. 46:45: 200 / 180 / 180.
    27:45: 240 / 220 / 220. 17:45: 260 / 240 / 240.
- **Depth-5 cold children of [49:45]:**
  - 0:12: 650 / 610 / 620.
  - 0:22: 380 / 350 / 350. 0:17: 470 / 440 / 450. 0:18: 510 / 480 / 490.
    0:11: 530 / 500 / 490.
- **Children warm-started from node 10's Z and conditioned Y (E1w):** 7:45
  120 → 100; 47:45 210 → 170. Near the end the polish adds +4.5 to +13,
  which is 30-60% of the lag.
- **How much of the lag the polish recovers** near fathoms: 45-100% on
  depth-4 nodes; 15-20% on the slow depth-5 node (lag 13-15).
- **Cost per evaluation** at m = 45-46: exact GL 0.25 s, HGR13 0.5 s, HG8
  1.4 s, projection 1.5 s, line search 2.4 s. For comparison, ten iterations
  cost about 3.3-3.8 s.

### E2: the cut-phase lag (dnncuts, production phase protocol)

- **[48:45]** (`logs/e2_n48_45_m171.log`):
  - Rounds 1-5 end 44-135 short, with lags of 59-69; the polish adds +6.3
    to +9.8.
  - The phase stalls after round 5. The quiet round from round 5's
    multipliers fathoms at 1,060 iterations; the polished value is above the
    cutoff at the 1,000 check.
  - The gain shrinks during the quiet round: +31.5 at iteration 100, +2.4 at
    700, +1.2 at 1,000.
- **[9:45]** (`logs/e2_n9_45_m187.log`):
  - It reproduces the search's phase exactly: round 2's best is 48,750.28,
    63.77 short.
  - The quiet round converges to a certificate of 48,793.9 with primal
    48,796.3, still 20.1 short.
  - At round 2's end the certificate was 43.6 below the converged value and
    the primal 5.7 above it; the polish added +8.3.

### E3: plain run against the polished certificate, same core (`logs/e3_*`)

| Node | Plain | Polish (HGR13, W = 15) |
|---|---|---|
| 7:45 | 310 its, 96.3 s | 280 its, 87.4 s |
| 47:45 | 450 its, 145.9 s | 430 its, 151.0 s |
| 49:45,0:12 | 650 its, 188.9 s | 620 its, 179.1 s |

- On 47:45 the polish itself cost 5.9 s, and the time per iteration drifted
  8% between the two runs. Iteration counts are deterministic, CPU is not.
- **W = 15 is too wide for slow crawls; use about 10.**

### E4 and the Q3 numbers

- Node 10 at its 300-iteration give-up: `logs/e4_n10*.log`.
- Phase miss [9:45]: `logs/e2_n9_45_childgl_row19.log` and
  `logs/e2_n9_45_join_row19_*.log`.

### E6: search A/B on the [9:45] sub-subtree

- Setup: serial bnbdrv in `run_pool_search.sh`'s configuration; the second
  arm adds the polish within W = 10 (`logs/e6_*`, `logs/e6_summary.txt`).
- Same tree in both arms: 47 nodes, 46 fathomed, 3,492 FW iterations.
- **Total CPU 1,490.7 s → 1,418.7 s.**
- The 46 pooled children: 3,400 → 3,020 iterations, 943 → 868 s. 32
  children did not change (fast fathoms); 14 saved 10-60 iterations.

## Answers to the four questions

**Q1. Can the certificate extract more from a given iterate?**
- Only through the skew polish; every other route is closed above.
- Per check it adds +3 to +13 units near fathoms and +6 to +10 at cut-round
  ends.
- The cheap version (HGR13) costs 0.5 s per evaluation, about 15% of the
  ten iterations between checks, so it has to be confined to checks near the
  cutoff.

**Q2. Where does the lag come from, and can it be closed cheaply?**
- After each Polyak update, PRSM restarts from the old (Y, Z). The primal
  re-settles in about 100 iterations; the dual multiplier needs about 1,000.
- At [9:45] about 84% of the round-end lag is certificate deficiency. The
  polish recovers about a fifth of that; the rest is the multiplier itself.
- Nothing on the certificate side is cheaper than the quiet round. Ways to
  shorten it belong on the solver side: smaller multiplier changes per
  round, fewer and longer rounds, or a dual-only ascent at frozen multipliers.

**Q3. Can we fathom more without new solves?**
- No, and there is a proof why:
  - At a dual-optimal multiplier, the GL/RLT1 minimizing face contains the
    DNN optimum Y*.
  - So every child whose marginal x*_ab is positive gets exactly the
    parent's value.
  - Away from convergence the child bounds are less flat (node 10: median
    +46), but there the parent is far shorter still (229).
- Measured: 0/47 children at node 10's give-up and 2/46 at the converged
  phase miss, about 1% of that class.
- An exact per-child correction was built (`scripts/child_gl.py`),
  measured and closed.

**Q4. What does the SDP-certification literature offer?** (From memory; no
paper was re-read or downloaded.)
- Jansson, Chaykin and Keil's VSDP rigorous bounds have exactly our
  structure, with the trace bound m+1 known exactly.
- Rump's Cholesky test is already in use.
- Neumaier and Shcherbina's safe LP bounds correspond to the old box
  certificate, which GL dominates.
- Kim, Kojima and Toh's y0 bisection is blocked here (Q1).
- BiqCrunch and BiqBin derive child bounds from the parent's dual the same
  way, and theirs are just as flat.
- Exact rational rounding (Peyrl-Parrilo) does not raise values; it matters
  only for an independent audit of the final proof.
- The RLT1 dual ascent (Adams-Johnson; Hahn-Grant) applied to the Lagrangian
  residual is the one thing that adds value. It is item 1.

## Recommendation

- **Do not expect the certificate to raise the bound level.** The slack
  left at a given iterate is 2-13 units near fathoms.
- **If about 1.04x is worth two hours of four cores:** run the node-10
  subtree A/B with the opt-in environment above, and adopt the polish at
  1.03x or better.
- **A production version (about half a day)** would:
  - replace the environment variables with a `bnb.C` switch;
  - run HGR only when AUTO's exact value is within W of the cutoff, plus at
    round ends (already lean, with an early exit);
  - optionally add a faster LAP so that more passes are affordable (#2).
- **Close:**
  - all child-bound work;
  - the line-search, averaging and projection variants;
  - eigenvalue-bound refinements;
  - the roadmap item "skew-symmetric certificate polishing at non-converged
    iterates" as a 1.05-1.15x lever on the cut phase. It measured about 2%
    on the phase, and about 1.04x overall where it acts on fathoms.

## Files

- **`latest/dnn.C`:** an opt-in `DNN_CERTX` block.
  - Bits: 1 exact GL, 2 HG, 4 line search, 8 average, 16 dump of M, 32 HGR,
    64 warm K, 128 projected post-step multiplier.
  - Switches: `_USE`, `_WINDOW`, `_EVERY`, `_PASSES`, `_RPASSES`, `_QUIET`,
    `_CALLS`.
  - `CertInfo.lam` added. The default path is unchanged.
- **`latest/dnntest.C`:** `-K` accepts a `;`-separated list of children.
- **This directory:**
  - `NOTES.md`, `theory.md`;
  - `run_e1.sh`, `run_e1w.sh`, `run_e2.sh` … `run_e6.sh`, `queue.sh`;
  - `scripts/`: `validity_bruteforce.py`, `validity_qaplib.py`,
    `summarize.py`, `e1_table.py`, `child_gl.py`, `join_children.py`,
    `use_smoke.sh`;
  - `logs/`.
  - The 36 MB dumps of M were not committed.

## Checks after the merge into cognac (Nov 4, 23:40 UTC)

- **The default path is unchanged.** The plain E3 runs match
  `pool_search/trace.txt` to the printed digit: 7:45 has lb 48,815.8282 in
  310 iterations, and 47:45 has 48,814.2874 in 450. The `dnn.C` diff only
  adds lines, and every new path sits behind `DNN_CERTX`. `make solver
  dnntest` builds and `make dnncheck` passes after the merge.
- **The validity test was rerun independently.**
  `scripts/validity_bruteforce.py latest/dnntest DIR 8`: 8 random n = 8
  instances, root and one fixed node each, 519 values per run with bits
  1+2+4+8+32. None is above the enumerated optimum (worst slack 0; PASS).
  This catches errors in the dual-ascent identities, not rounding errors:
  the rounding allowances are about 1e-6 units, far below what
  enumeration resolves.
- **Correction to Q2.** In the production protocol PRSM does not restart
  from the old (Y, Z) after a Polyak update: it keeps Z and resets Y to
  e0e0' (`dnncut.C`: `par.warm_primal = cp.warm_primal && r > 0`, and
  `cp.warm_primal` was 0 in the search). The convergence study
  (`../convergence2/REPORT.md`, F1) shows that this reset causes most of
  the round-end lag. The measured split of the lag (84% certificate
  deficiency) stands.
- **Open before the polish goes into the proof run:** an independent read of
  the HG/HGR rounding allowance (`mg_hg`, `mg_hgr` and the per-pass `he` in
  `dnn.C`). The E3 run on 47:45 shows why it matters: the polished bound
  crossed the cutoff by 0.0026 units (48,814.0526 against 48,814.0500). That
  is three orders of magnitude above the allowance, but it is the kind of
  margin a final audit should recompute.

## Code status (Nov 9)

The `DNN_CERTX` experiment that produced these numbers (the variant
certificates, the Hahn-Grant polish, the line search, the averaged
multipliers) was removed from `latest/dnn.C` on Nov 9, ahead of the settings
registry (`docs/settings_plan.md`, step 1): no production recipe used it, and
its best variant, the skew polish (mode 32), was worth about 1.04x in
iterations to fathom. The code is in the history before that commit; the
logs in `logs/` and this report stand on their own.
