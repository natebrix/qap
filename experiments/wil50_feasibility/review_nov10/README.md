# Review of the DNN bound as implemented (Nov 10)

**Question.** Where, in the bound as it runs today (`fullrun/fullrun.par`, HEAD
`bfe38d4`), is there (a) faster ADMM convergence or (b) a stronger bound that
would matter in the wil50 run?

**How.** I read the solver (`latest/dnn.C`, `dnncut.C`, `bnb.C`) and the
record, checked every idea against the record before counting it, and ran
the cheap tests on the cloud box (4 Xeon cores with AMX, one process per
core). Sections 2 and 3 hold what was measured; section 4 ranks what is left.
`[m]` is measured here, `[r]` is from the record, `[e]` is my estimate.

**Status: in progress.** Everything is final except the full-subtree zone
test (section 3.5) and the ranked list (section 4), which follow when that
test finishes.

**So far.**
- The ADMM iteration and the typical leaf are tuned out: nothing tried is
  worth more than 2% of the tree (sections 2 and 3.1).
- Half of the leaf work is in the 14% of leaves that end within 3 units of
  the cutoff, so a few units of bound shorten the tree twice: fewer nodes
  and shorter marginal leaves (section 1).
- A warm solve that will miss can stop at iteration 110 instead of 300,
  worth up to 6% on the tail class, if cold solves keep their 300 (sections
  2.3 and 3.2).
- In the cut phase only three Dykstra sweeps per iteration help: +7.7 units
  where the phase stops, 1.027x on node 10's subtree (section 3.3).
- One union node fathoms a line's strong children four times cheaper than
  one by one, but the strong leaves are the cheap ones: about 1.08x (section
  3.6).

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
  + 0.3154 s per iteration, my fit over stage B's 242 pooled leaves).

So half of the leaf work is a few marginal leaves, which a slightly stronger
bound would shorten, and the other half is short solves that are already
near their floor.

## 2. What was measured

Everything here is at production settings on the binary built at `bfe38d4`
(before 6335dea's change to rule 9). Iteration counts
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

Session 2, the cut phase (`runs/s2_*`; this session's CPU is 6% above session
1's for the same work, so each session has its own baseline):

| arm | nodes / fathomed | phases (fathomed) | iterations (first + phase) | DNN CPU | against baseline |
|---|---|---|---|---|---|
| baseline | 94 / 92 | 2 (0) | 6,429 (4,359 + 2,070) | 2,055 s | |
| three Dykstra sweeps per iteration (`cut.cp_sweeps = 3:1`) | 94 / 92 | 2 (0) | 6,094 (4,064 + 2,030) | 2,002 s | 1.055x iterations, 1.027x CPU |
| stall window 400 instead of 600 (`cut.stall = 2`) | 94 / 92 | 3 (1) | 6,618 (5,259 + 1,359) | 2,066 s | 0.971x, 0.995x |
| stall window 1,000 (`cut.stall = 5`) | 94 / 92 | 2 (0) | 7,200 (4,110 + 3,090) | 2,330 s | 0.893x, 0.882x |

The two phases of each arm, in units short of the cutoff:

| arm | root: iterations, CPU | root ends | pooled miss `45:49`: first solve, phase end (iterations) |
|---|---|---|---|
| stall window 400 | 780, 285 s | 83.6 | 44.5, 33.8 (400) |
| baseline | 1,110, 402 s | 73.8 | 39.5, 23.9 (960) |
| stall window 1,000 | 1,690, 607 s | 65.4 | 34.9, 20.7 (1,400) |
| three sweeps | 1,110, 429 s | 66.1 | 34.5, 20.7 (920) |

- Three sweeps reach the long window's bound in the baseline's iterations:
  +7.7 units for 7% more phase CPU.
- On the 92 leaves, three sweeps save 295 iterations, all of them in the 13
  slow leaves (1,827 to 1,515); the other 79 are unchanged (1,932 and 1,949).
  The long window saves 249 the same way and spends 1,020 more in the phases.
- With the short window one marginal leaf (`0:49`, 354 iterations in the
  baseline) misses by 2.3 and needs a phase of its own.

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
at 110, which is the loss of section 2.1's last arm. With the rule's window
cut from 50 to 30 iterations it fires at 70 and still loses no fathom here
(790 iterations on the 9 misses), but this set has only three slow fathoms
(120-150 iterations), too few to trust a shorter window.

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

### 3.3 The cut phase: a flat optimum on node 10, and one cheap gain

The phase became one CP-ADMM solve on Nov 7. Its stall rule and its beta
factor (`wc=0.7`) came over from the Polyak rounds and were not varied in
the search since; three sweeps were measured on node 13's frozen pool (+6.2
at 1,000 iterations) and never run in the search (`roadmap/admm/cpadmm/`).

- **Three Dykstra sweeps per iteration: +7.7 units where the root's phase
  stops, 1.027x on the subtree** (section 2.1, session 2). No node changes
  outcome on this subtree, so what shows is only the 13 slow leaves getting
  17% shorter. The gain thins out below the root: +5.0 at the pooled miss's
  first solve and +3.2 after its phase. Under the tree model 3-8 units are
  worth 1.05-1.25x where nodes do change outcome [e]; one 94-node subtree
  cannot show that.
- **The stall window is at a flat optimum.** A shorter window saves phase
  iterations and pays them back in the leaves; a longer one costs more than
  it returns.
- **A pooled miss's own phase buys little.** `45:49` arrives 39.5 short with
  the root's pool, and its phase ends 23.9 short after 960 iterations, 15%
  of the subtree. Without the root's pool the same node needs 1,540
  iterations to end 26.2 short (session 1's last arm), so the cuts reach the
  same place either way and a pool only moves the work. Cutting that phase
  to 400 iterations moves the cost to its children about one for one (the
  short-window arm).

### 3.4 Bound strength: no cut family is left at degree 2

- **The triangle pool's reach ends below the cutoff at depth 3**: node 13
  certifies 48,771.7 after 3,000 CP iterations, 42 short, nearly flat [r].
  The record's synthesis finds no further family at degree 2 (section 5),
  and I found none.
- **Box cuts inside the CP phase are the one untried combination, and I rate
  it unlikely.** The pool format carries boxes and the CP phase has no
  multiplier step for them to disturb, which is what killed them in the
  Polyak rounds. But the +17 they showed was measured on top of an
  unconverged triangle iterate (48,714 certified, where the triangles alone
  later reach 48,772 in CP form), and by round 12 of the pool test the
  Lagrangian with boxes was no higher than without [r]. The cheap first
  step is to separate boxes at a CP phase's final iterate on nodes 10 and
  13: if the violations are small there, the family is closed for good.
- So what a node gets from the cuts is set by two things that can still
  move: how converged its phase is when it stops (three sweeps, section
  3.3) and which nodes run a phase at all (the zone, below).

### 3.5 The cut zone at far nodes: one phase at the parent for many at the children

Production runs a phase at a node within 244 units of the cutoff after its
first solve (`cut.zone = 0.005`). A node farther out branches without one,
and each child that then lands inside the zone runs its own phase from
scratch. The record's wide-zone test was at one node 455 short with
Polyak-rounds phases, a wash (`tail_probe/`); the band between 244 and about
400 was not tested, and not in CP form.

First-level test (`zone/first_level/`, wil50t, the node-10 arm's settings):
two depth-3 nodes from stage B's frontier, bounded cold, then their children's
first solves, with the parent's phase (zone 0.009) and without (0.005).
Stopped after 16 minutes for the full-subtree test below, so the same 13 and
11 children are in both arms:

| parent | its phase | children that fathom in the first solve, without / with | first-solve iterations, without / with | where both miss, the pool is worth |
|---|---|---|---|---|
| `13:1,21:0,45:2`, 325 short | +124 in 750 iterations, 277 s | 5 / 10 of 13 | 2,839 / 1,424 | +57 to +71 |
| `0:0,21:1,45:30`, 362 short | +132 in 730 iterations, 274 s | 2 / 7 of 11 | 2,884 / 1,841 | +53 to +67 |

- **Without the parent's phase most children miss by 15-150 and each needs
  a phase of its own**; with it most of them fathom in the first solve and
  the rest miss by 8-94.
- **One 740-iteration phase stands in for about twenty.** Section 3.3 found
  that a pool moves the work without extending the reach, so the node count
  below should not change; what changes is that the cuts are found once at
  the parent instead of once per child.

(The complete subtree of the 325-short node under zone 0.005 and 0.0075 is running.)

### 3.6 A search-side candidate: one node for a line's strong children

A union node is the parent with all but a chosen set of one line's children
forbidden (`DNNParams::elim`), warm-started from the parent's own state. It
is a solve of the parent's size, and if it crosses the cutoff every child in
the set is fathomed. The record closed this as region branching on one
nug30 node, solved cold at a far miss (`experiments/region_branching/`), and
noted that the wil50 check was never run. This is that check, warm, at a
parent whose children all fathom.

Probe (`union/`, a scratch pass `-G` in the test driver, the solver
untouched): parent `26:27,8:38,34:28` on `wil50.dat`, 72 short at its
give-up, line row 37. Its 47 children all fathom alone, in 1,439 iterations
and 474 s. "Strong" is by how the children did alone:

| set | as one union node | members alone | CPU, union / alone |
|---|---|---|---|
| 2 strongest | fathoms, 22 it, 8.5 s | 40 it, 13.7 s | 0.62 |
| 4 strongest | fathoms, 29 it, 10.9 s | 80 it, 27.2 s | 0.40 |
| 8 strongest | fathoms, 39 it, 13.4 s | 160 it, 55.1 s | 0.24 |
| 16 strongest | gives up 1.5 short at 330 it, 118 s | 325 it, 112 s | fails |
| 4 marginal (alone: 25-28 it, under 3 above the cutoff) | fathoms, 121 it, 42.8 s | 105 it, 35.5 s | 1.21 |
| 2 slowest | fathoms, 140 it, 50.2 s | 129 it, 39.7 s | 1.27 |
| 8 slowest | gives up 28.5 short at 300 it, 108 s | 442 it, 139 s | fails |

- **A union is about as strong as its weakest member, a little less.** That
  is the record's finding for blocks, and it holds here.
- **Among strong members it converges almost as fast as one child**: 22, 29
  and 39 iterations for 2, 4 and 8 of them.
- **So it pays for exactly the children that fathom with room to spare**:
  here the 8 that are 19 or more above the cutoff at iteration 20, four
  times cheaper as one node. One marginal member is enough to sink a set.
- **The prize is small, because the strong leaves are the cheap ones.** On
  node 10's two lines 16 and 21 of the 46 leaves stop 23 or more above the
  cutoff. Those 37 leaves are 19% of the leaf CPU and 11% of the subtree
  [m], so union nodes are worth about 1.08x there once the unions' own cost
  is paid [e]. The marginal leaves that carry the cost cannot be grouped.
- **What it needs.** A way to pick the strong children before bounding them
  (the parent's marginals, or the union's own first 20 iterations; being checked), nodes that carry eliminations in `bnb.C` and in
  the state file, and a test with a pool, which the test driver cannot do.
  A few days of work.

## 4. What is left, ranked

By expected effect on the run's CPU. "Subtree" is node 10's.

| # | what | evidence | expected | cost | next step and kill line |
|---|---|---|---|---|---|
| 1 | a wider cut zone, so that far nodes that will branch run a phase and pass their pool on | first level: 10 of 13 and 7 of 11 children fathom in the first solve against 5 and 2 [m]; the full subtree is being measured (section 3.5) | | a settings line | |
| 2 | three Dykstra sweeps in the cut phase, `cut.cp_sweeps = 3:1` | +7.7 units where the root's phase stops for 7% more phase CPU, 1.027x on the subtree with no outcome changed [m]; +6.2 on node 13 [r] | 1.03x; up to 1.25x where nodes change outcome [e] | a settings line | node 5's subtree with and without; adopt at 1.03x in DNN CPU on matched nodes with no fathom lost |
| 3 | a give-up floor for warm solves only, 110-150 instead of 300 | replay: 9 warm misses fall from 2,740 to 1,170 iterations, no fathom lost of 15 [m] | up to 1.06x on the tail class, nothing near the cutoff [e] | a settings key and a few lines; the zone test on the rule's predicted value | tail subtree A/B; kill under 1.02x |
| 4 | one union node for a line's strong children | 8 of 47 children in 39 iterations against 160 [m] | about 1.08x [e] | a few days: eliminations in the node and the state file, a predictor | a pooled parent first; kill if the strong set is under a quarter of the leaves |
| 5 | the early certificate at iteration 5 of warm solves | 1.016x on the subtree [m] | 1.016x | promote `DNN_ADMM_CROSS_FIRST` from an experiment name to a settings key | none; take it with the next settings change |
| 6 | box cuts at the CP phase's final iterate | untried, and the evidence is weak (section 3.4) | 0 to +10 units [e] | an hour with `dnncuts` | close if the largest violation is under 0.05 |
| 7 | branching marginals from the phase's final iterate (`bnb.C` takes `dnn_xbar` from the first solve before the phase's state replaces it) | not measured | unknown | two lines and a key | node 10: if the line changes, run both subtrees |

Not worth doing, measured here: a lower warm beta, damped dual steps at a
warm start, a shorter or longer stall window (sections 2.1 and 2.2).

**Where the next factor is.** Nothing above is worth 1.5x. Each item moves
a few percent, and together they are perhaps 1.2-1.3x if they stack. The
bound's reach is set by the degree-2 relaxation and the triangles, and the
ADMM that computes it is near its floor. The line choice at depths 1-2 is
the more likely place for a larger factor: the tho30 and tho40 line-truth
studies find the best line at 0.6-0.8x of rule 9's nodes, and at 21:1 row
45 scores 1.6-1.9x on the modelled subtree
(`sampling/results/depth1_lines.md`).

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

1. **`pipeline_ab/run_pair.sh` did not pin `dnn.beta` and `dnn.gu_zone`.**
   Since 34e347c changed the registry defaults, an arm at HEAD ran 0.21 and
   0.003, not the full run's 0.3 and 0. Every arm in this review exported
   `QAP_DNN_BETA=0.3 QAP_DNN_GU_ZONE=0` first. Fixed on Nov 10: `arm_env`
   now exports both, and the settings file it builds for `cpcutt` equals the
   one this review's baseline ran.
2. **A coordinator run with `-F` failed at HEAD.** The solver rejects `-F`
   with `-r` since f67e0d7 and the worker passed both. Fixed in bfe38d4
   (`coord/worker.py`, with a `fixes` scenario in `tests/coord_test.py`).
3. From the depth-1 line test: `fullrun_check` fails at HEAD because
   `settings.golden` predates `QAP_DNN_WARM_DIR`, and `sampling/common.py`'s
   `write_one` rejects QAP-STATE 9.
