# Review of the DNN bound as implemented (Nov 10)

**Question.** Where, in the bound as it runs today (`fullrun/fullrun.par`, HEAD
`bfe38d4`), is there (a) faster ADMM convergence or (b) a stronger bound that
would matter in the wil50 run?

**How.** I read the solver (`latest/dnn.C`, `dnncut.C`, `bnb.C`) and the
record, checked every idea against the record before counting it, and ran
the cheap tests on the cloud box (4 Xeon cores with AMX, one process per
core). Sections 2 and 3 hold what was measured; section 4 ranks what is left.
`[m]` is measured here, `[r]` is from the record, `[e]` is my estimate.

**Status: in progress.** Sections 1, 2.1 (session 1), 2.2-2.4, 3.1, 3.2, 5 and
6 are final. Session 2 of section 2.1 (the cut phase's stall window and
sweep count) is running; the cut-phase findings, the bound-side findings and
the ranked list (sections 3.3, 3.4 and 4) follow when it and a wider-zone
test finish.

**So far.**
- The ADMM iteration and the typical leaf are tuned out: nothing tried is
  worth more than 2% of the tree (sections 2 and 3.1).
- Half of the leaf work is in the 14% of leaves that end within 3 units of
  the cutoff, so a few units of bound shorten the tree twice: fewer nodes
  and shorter marginal leaves (section 1).
- A warm solve that will miss can stop at iteration 110 instead of 300,
  worth up to 6% on the tail class, if cold solves keep their 300 (sections
  2.3 and 3.2).

## 1. Where the time goes

The node-10 subtree in the production configuration (`pipeline_ab/run_pair.sh
cpcutt`, wil50t, 94 nodes, baseline arm of section 2) [m]:

| work | count | DNN iterations | DNN CPU |
|---|---|---|---|
| first solves that fathom (the leaves) | 92 | 3,759 (58%) | 1,092 s (56%) |
| first solves that give up | 2 | 600 (9%) | 176 s (9%) |
| cut phases | 2 | 2,070 (32%) | 670 s (35%) |
| total | 94 nodes | 6,429 | 1,938 s |

The tail class has the same shape: the Nov 7 estimate for node 5 (316-329
short, `pipeline_ab/`) splits into about 50% leaves, 39% phases and 10%
give-ups (my split of its trace).

Two facts about the leaves decide what is worth doing there:

- **The median leaf is 20 iterations, and the mean is 41.** 13 of the 92
  leaves (14%) take 60 iterations or more and carry 49% of the leaf
  iterations; 3 leaves (220, 310 and 354 iterations) carry 24%. Every slow
  leaf stops within 3 units of the cutoff: its bound converges to a value
  barely above it. [m]
- **A leaf's fixed cost is 1.4 s, 18% of a 20-iteration leaf** (CPU = 1.40 s
  + 0.3154 s per iteration over stage B's 242 pooled leaves) [r].

So half of the leaf work is a few marginal leaves, which a slightly stronger
bound would shorten, and the other half is short solves that are already
near their floor.

## 2. What was measured

Everything here is at production settings on HEAD's binary. Iteration counts
are deterministic; CPU ratios are same-session (arms side by side, one per
core).

### 2.1 The node-10 subtree, one change per arm

Session 1 (`runs/s1_*`):

| arm | nodes / fathomed | phases (fathomed) | iterations (first + phase) | DNN CPU | against baseline |
|---|---|---|---|---|---|
| baseline | 94 / 92 | 2 (0) | 6,429 (4,359 + 2,070) | 1,938 s | |
| warm beta 0.15 instead of 0.25 (`w=0.15`) | 94 / 92 | 2 (0) | 5,985 (3,915 + 2,070) | 1,893 s | 1.074x iterations, 1.024x CPU |
| one extra certificate at iteration 5 of warm solves (`DNN_ADMM_CROSS_FIRST=5`) | 94 / 92 | 2 (0) | 6,267 (4,197 + 2,070) | 1,908 s | 1.026x, 1.016x |
| give-up floor 150 instead of 300, every call (`-x 0:4:1:150:1`) | 94 / 92 | 6 (5) | 8,900 (5,338 + 3,562) | 2,601 s | 0.72x, 0.75x |

- The baseline reproduces the Nov 7 record (6,426 iterations).
- `w=0.15` removes 12% of the leaf iterations but only 4.5% of the leaf CPU:
  the iterations get dearer. A wash.
- The floor of 150 loses because the cold root stops 318 short instead of
  214, misses the zone, and its phase is replaced by six phases one level
  down. A lower floor must not apply to cold solves.

Session 2 (`cut.stall` 2 and 5, `cut.cp_sweeps = 3:1`) is running.

### 2.2 A warm child's certificate, iteration by iteration

Eight depth-4 children of two depth-3 parents (`run_warmset.sh` families f9
and f6), a certificate at every iteration, at the production warm beta
(`w=0.25`). The record had this profile at `w=0.7` only. Values are relative
to the parent's final bound; the cutoff is at +65 (f9) and +86 (f6).

| child, by iteration | 1 | 2 | 5 | 10 | 15 | 20 | 25 | crossing |
|---|---|---|---|---|---|---|---|---|
| f9 37:21 (the fastest) | -448 | -251 | -128 | 0 | | | | 14 |
| f9 37:14 | -420 | -315 | -175 | -56 | -4 | +34 | +55 | 28 |
| f6 43:22 (the slowest) | -318 | -227 | -120 | -31 | +7 | +30 | +49 | 51 |

Over the eight children the start costs 320-480 units, the child is back at
its parent's level after 10-15 iterations, and it crosses at 14-51. The
multiplier's norm changes by under 0.1% over the whole solve.

Then the same eight children with the dual steps scaled down for the first
iterations (`DNN_ADMM_WARM_DAMP=D:K`: factor D for K iterations):

| arm | crossing iterations, f9 then f6 | sum |
|---|---|---|
| production | 28 33 14 15, 17 23 51 23 | 204 |
| no dual step for 1 iteration | 29 34 14 15, 17 24 52 23 | 208 |
| for 2 | 29 34 15 16, 18 24 52 24 | 212 |
| for 3 | 30 35 16 17, 19 25 53 25 | 220 |
| for 5 | 32 37 17 19, 20 27 54 26 | 232 |
| factor 0.3 for 3 | 30 35 15 16, 18 25 52 24 | 215 |

Holding the multiplier fixed makes the certificate recover faster (-133
against -175 at iteration 5 for 37:14), but it is heading for a plateau
about 95 below the parent, and it falls back as soon as the dual steps
resume. Every arm crosses later than production. Production passes that
plateau at iteration 8 of 28, so the other 20 iterations are dual progress,
and the first dual steps are not what delays a leaf. Closed.

### 2.3 The give-up rule on warm solves (replay, no new runs)

The rule cannot fire before iteration 110 and production's floor is 300 for
cold and warm solves alike (`dnn.giveup`). Replayed on the 24 warm children
of the Nov 6 warm set at `w=0.25` (`giveup_sim.py`):

| floor | iterations, 9 misses | fathoms lost of 15 | bound at the give-up against floor 300 |
|---|---|---|---|
| 300 | 2,740 | 0 | |
| 200 | 1,890 | 0 | 2-9 lower (8 of 9) |
| 150 | 1,490 | 0 | 5-18 lower (8 of 9) |
| 110 | 1,170 | 0 | 12-29 lower (8 of 9) |

The same floors on the six cold parents would stop them 150-200 units lower
at 110, which is the loss of section 2.1's last arm.

### 2.4 The zone boundary in stage B (tally, no new runs)

Stage B's depth-3 children by their depth-2 parent's first-solve shortfall;
parents up to 244 short ran a phase and passed their pool on, the others did
not (`zone_rd.py`):

| parent short | children | pooled | fathom in the first solve | mean iterations | mean CPU |
|---|---|---|---|---|---|
| 215-230 | 23 | yes | 0.87 | 82 | 26.9 s |
| 230-244 | 31 | yes | 0.84 | 89 | 30.3 s |
| 244-258 | 20 | no | 0.70 | 123 | 39.1 s |
| 258-272 | 10 | no | 0.70 | 150 | 46.5 s |
| 272-286 | 28 | no | 0.50 | 184 | 58.1 s |
| 286-300 | 14 | no | 0.14 | 265 | 82.8 s |

There is no cliff at 244: the step across the boundary (0.84 to 0.70) is
within the sampling error of 20-31 children per bin. The sample cannot say
whether a phase at a 244-300 parent would pay.

## 3. Findings

### 3.1 The leaves: the median leaf is at its floor, the cost is in the slow tail

- **Nothing I tried shortens a typical leaf by more than 2% of the tree.**
  A lower warm beta is a wash, an early certificate is 1.6%, and damping the
  first dual steps is negative (sections 2.1 and 2.2). The record closes the
  rest: the two-phase warm beta, the shipped eigenbasis, the polished
  certificate, child bounds (section 5).
- **The reason is in section 2.2's profile.** The start's loss is recovered
  in the first third of a leaf's iterations; the other two thirds are dual
  progress at the rate the iteration has, and no start or step-size change
  shortened them.
- **Half of the leaf work is in the 14% of leaves that take 60 iterations or
  more**, and every one of them ends within 3 units of the cutoff. Their
  cost is set by how far above the cutoff their bound converges, so they are
  the first to gain from a stronger pool: a few units shorten them by tens of
  iterations each. Section 3.3 measures this through the parent's phase.

### 3.2 Give-ups: a warm solve does not need 300 iterations to know

- **A warm first solve that will miss is recognised at iteration 110 in the
  replay**, with no fathom lost on 24 children (section 2.3). Production
  waits until 300 because `dnn.giveup` is one number for cold and warm
  solves, and a cold solve does need it (section 2.1, last arm).
- **Worth at most the 190 iterations per warm miss**: nothing on node 10
  (one warm miss), up to 6% on the tail class where warm misses are 10% of
  the CPU [e]. The 12-29 units the bound loses at the earlier stop must not
  move nodes out of the zone, so the zone test should then use the rule's
  predicted value.
- This needs a warm-only floor (a few lines and a settings key) and a
  tail-subtree A/B. Not built here.

### 3.3 The cut phase, 3.4 Bound strength, and 4. Ranked list

To come with session 2 and the zone test.

## 5. Checked against the record and closed

Ideas that look open from the code and are not. Each was tested and has a
kill on file; none should be proposed again without new evidence.

| idea | verdict in the record | where |
|---|---|---|
| two-phase warm beta (0.25, then 0.7 at iteration 60) | 1.002x per matched node against a 1.02x line | `pipeline_ab/` |
| damped dual steps at a warm start | slower at `w=0.7` (0.5 for 10 iterations); slower at `w=0.25` in six arms (section 2.2) | `speed/nov6_admm/` 5.3, here |
| restricted PRSM pieces, Anderson, Halpern, momentum, over-relaxation | killed | `speed/rprsm/`, `roadmap/admm/` |
| other solvers for the same SDP (bundle, Burer-Monteiro, SDPNAL+, CGAL) and other eigensolvers | closed | `roadmap/admm/` |
| the parent's eigenbasis shipped to a warm child | 1.4-1.5% against a 3% line | `speed/nov6_basis/` |
| fused transforms (`DNN_FUSE`) | 375 to 382 ms per iteration | `amx/README.md` |
| a different multiplier from the same iterate (line search, averages, the projected post-step multiplier) | below the standard certificate where fathoms are decided | `roadmap/certificate/` |
| certificate polish within 10 units of the cutoff (`certx`) | 0.983x in the tree, removed | `pipeline_ab/` |
| child bounds from the parent's certificate | 0 of 2,209 fathom; flat at convergence (proved) | `roadmap/certificate/` |
| more parent iterations (1,000 instead of 300) | break-even | `roadmap/convergence2/` |
| pooled children's first solve in CP form (`pool.cp = 1`) | 1.053x in DNN CPU, 1.004x in wall on node 10 (side 0): shorter phases, slower leaves; not adopted | `pipeline_ab/` |
| facility-side and location-side cut families, RLT2/RLT3, Lasserre 2, disjunctive cuts, parity, cutoff half-space | no family meets the bar; degree-2 ceiling | `new_cuts/SYNTHESIS.md` |
| large-Q box (count) cuts in the Polyak rounds | -32.3 and -22.6 on nodes 10 and 13 | `new_cuts/box_pool/` |
| region branching | a block node is about as strong as its weakest cell (one nug30 node, cold) | `experiments/region_branching/` |
| another line at depths 3-5 | rule 9 within about 5% of the best line | `tuning/lines/` |
| flow automorphisms | none | `roadmap/tree_size.md` |
| zone 0.010 at the 455-short tail node | a wash (Polyak-rounds phases, one node) | `tail_probe/` |

Two of these verdicts are narrower than their titles, and section 4 uses
that: the box cuts and the wide zone were tested only with Polyak-rounds
phases, and region branching only on nug30, cold, at a far miss.

## 6. Harness findings

1. **`pipeline_ab/run_pair.sh` does not pin `dnn.beta` and `dnn.gu_zone`.**
   Since 34e347c changed the registry defaults, an arm at HEAD runs 0.21 and
   0.003, not the full run's 0.3 and 0. Every arm in this review exported
   `QAP_DNN_BETA=0.3 QAP_DNN_GU_ZONE=0` first. Two exports in `arm_env`
   would fix it; not done here.
2. **A coordinator run with `-F` failed at HEAD.** The solver rejects `-F`
   with `-r` since f67e0d7 and the worker passed both. Fixed in bfe38d4
   (`coord/worker.py`, with a `fixes` scenario in `tests/coord_test.py`).
3. From the depth-1 line test: `fullrun_check` fails at HEAD because
   `settings.golden` predates `QAP_DNN_WARM_DIR`, and `sampling/common.py`'s
   `write_one` rejects QAP-STATE 9.
