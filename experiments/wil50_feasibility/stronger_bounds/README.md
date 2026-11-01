# Stronger DNN bounds for wil50: regional triangle cuts

Evaluation of Astra's research memo `research_memo_astra.md` (Oct 31 2026) and
its experiment bundle `astra/` (as received; see `astra/README.md`). The memo
proposes strengthening the facially reduced DNN bound with constraints on
*regions*: for a set S of locations, the membership variables
z_i = Σ_{a∈S} x_{ai} ∈ {0,1} have first and second moments

    p_i = Σ_{a∈S} Y_(ai),(ai),    W_ij = Σ_{a,b∈S} Y_(ai),(bj)

that the DNN relaxation already contains, and the Boolean triangle
inequalities on (p, W)

    W_ij + W_ik − W_jk − p_i ≤ 0   (three apexes),
    p_i + p_j + p_k − W_ij − W_ik − W_jk ≤ 1

are valid for every permutation but **not implied** by the DNN relaxation
with the atomic triangles (`../../dnn_cuts/`): the memo gives an exact
rational witness (K₃,₃ as both matrices, optimum 10, DNN bound 9, regional
cut raises it to 10). The regional version matters for wil50 because its
first matrix is the Manhattan distance on a 5 × 10 grid, which is exactly the
sum of 13 threshold indicators (row prefixes t = 1..4, column prefixes
t = 1..9): every term of the objective is a regional count, so the cuts speak
the objective's language. The memo's further steps (mixed-region and box
inequalities, selective triple consistency, a regional second-order moment
SDP block) build on the same regions.

Astra's pilot on wil50's hardest depth-3 node `34:40,28:22,14:29` (m = 47,
give-up bound 48359, cutoff 48816): 1500 uncut PRSM iterations give
48409.2; adding the 1000 most violated regional triangles at a uniform
multiplier 1 for 700 more iterations gives **48449.1 (+39.8)**; the same at
multiplier 3 gives 48238.1, at multiplier 2 with 100 cuts 48412.5. An uncut
PRSM continuation for 700 iterations stays at 48408.2.

## Implementation: `dnncuts -grid R:C`

`experiments/dnn_cuts/dnncuts.C` (the triangle-cut prototype used for nug30)
gained the regional family. `-grid 5:10` builds wil50's 13 grid regions in
original location indices and maps them to the free locations of the node;
`-prefix` builds the prefixes of the free locations for validation on
arbitrary instances. The cuts are folded into the objective in Lagrangian
form, `L + Σ μ_t G_t` with the shift lowered by Σ μ_t h_t, so the existing
certificate applies unchanged and **every printed value is a valid bound**.
Each round solves the DNN relaxation (warm-started), separates the most
violated regional triangles in the last primal iterate (`-k` per round), and
updates the multipliers by a Polyak step toward `min(cutoff, best + delta)`,
halving the step when the bound drops; `-mu0 M -t 0` reproduces a fixed
uniform multiplier.

### Validation

* K₃,₃ witness (`astra/toy6.dat`, ADMM mode `-e 0 -g 1.618 -G 0`, prefix
  regions): plain DNN 9.0000, two violated regional triangles (violation 1/8,
  as the memo predicts); fixed μ = 1 gives 9.25, the Polyak step with
  `-delta 0.5` reaches 9.9993 in two rounds (optimum 10).
* `validate_cuts.py`: 12 random n = 7 instances from the `tests/brute_force.py`
  families, root and one- and two-fix nodes, prefix and atomic cuts, every
  round's bound checked against enumeration: **84 certified values, 0 above
  the optimum**. Of 36 nodes, 6 have a DNN gap; the prefix cuts close 6.4% of
  it on average, the atomic triangles 10%.

## Node 1 reproduction (`run_node1.sh`, logs in `node1/`)

Four runs, one core each, on the 2.1 GHz Xeon:

| run | protocol |
|---|---|
| `astra` | 1500 PRSM iterations uncut, then the 1000 most violated regional triangles at multiplier 1 for 700 iterations (Astra's +39.8 trial) |
| `polyak` | same start, 8 rounds of 700 iterations, Polyak step toward best + 100, new cuts entering at 0 |
| `polyak1` | as `polyak`, new cuts entering at multiplier 1 |
| `control` | uncut PRSM for 7100 iterations (`dnntest -H 0`), the equal-budget reference |

Certified bound after each round (round 0 = the 1500 uncut iterations,
48409.27 in every run; `control` at the same iteration counts):

| round | iterations | `control` (uncut) | `astra` (μ = 1 fixed) | `polyak` | `polyak1` |
|---|---|---|---|---|---|
| 0 | 1500 | 48409.27 | 48409.27 | 48409.27 | 48409.27 |
| 1 | 2200 | 48409.41 | **48449.17 (+39.9)** | 48438.24 | 48297.01 |
| 2 | 2900 | 48409.45 | | 48472.31 | 48408.84 |
| 3 | 3600 | 48409.47 | | 48369.02 | 48458.21 |
| 4 | 4300 | 48409.48 | | 48452.81 | 48492.62 |
| 5 | 5000 | 48409.49 | | 48495.47 | 48511.31 |
| 6 | 5700 | 48409.49 | | 48518.33 | 48494.24 |
| 7 | 6400 | 48409.49 | | **48530.25 (+121.0)** | 48422.87 |
| 8 | 7100 | 48409.49 | | 48406.74 | 48416.05 |

* **Astra's trial reproduces exactly**: +39.89 against their +39.83, same
  1000 cuts at μ = 1, primal 48450.8 against the certificate 48449.2 (the
  modified problem is solved, so that is the value of this Lagrangian).
* **The uncut bound is converged at 1500 iterations**: 7100 iterations buy
  0.22. The whole gain is from the cuts; the equal-budget control is flat.
* **Optimizing the multipliers triples the gain**: the Polyak run reaches
  **48530.25, +121.0**, in round 7 with 5311 active cuts (a new 1000 added per
  round from the 7000–11000 violated), still rising 12 per round before the
  round-8 overshoot (the step halves after a drop, so a longer run would
  recover). Entering new cuts at μ = 1 (`polyak1`) overshoots in round 1
  (−112) and ends at +102; new cuts should enter at 0.
* Cost: each 700-iteration round is 410–430 s CPU, 2% more per iteration than
  uncut (the `L_add` term); the +121 cost 8 rounds, 3400 s on top of the 844 s
  uncut solve — about 5× the cost of the converged uncut bound and about 20×
  the search's per-node spend (the give-up rule stops this node at 300
  iterations).

Node 1 needs 48814 to fathom; the cuts closed 30% of its 405-unit gap. The
bound was still moving, and the multiplier update is the crudest possible
(uniform Polyak step, no cut deletion, 1000 cuts per round); neither the
ceiling of the family nor its price is established yet.

## Step 2: the 14 misses (`run_misses.sh`, logs in `misses/`, table by `summarize_misses.py`)

The `polyak` protocol of node 1 on each of the 14 misses of
`../tuning/giveup/all32.txt` (ids = line numbers): 1,500 uncut PRSM
iterations, then up to 8 rounds of 700 iterations, 1,000 new regional
triangles per round entering at multiplier 0, Polyak step toward best + 100;
the run stops when the certified bound clears the cutoff 48,814.05. One core
per node, four at a time, Xeon 2.1 GHz. (Three container restarts killed runs
in flight; the reruns reproduced the killed runs' rounds to the decimal, the
partial logs are in `misses/partial/`.)

| node | uncut, 1,500 it. | short | best certified | gain | round | fathomed | CPU s |
|---|---|---|---|---|---|---|---|
| 1 | 48,409.3 | 404.8 | 48,530.2 | +121.0 | 7 | no | 4,220 |
| 2 | 48,611.8 | 202.2 | 48,700.5 | +88.7 | 7 | no | 4,167 |
| 3 | 48,573.3 | 240.7 | 48,642.1 | +68.8 | 8 | no | 4,024 |
| 4 | 48,586.1 | 228.0 | 48,689.0 | +102.9 | 8 | no | 4,236 |
| 5 | 48,498.4 | 315.7 | 48,594.6 | +96.2 | 7 | no | 4,137 |
| 6 | 48,758.0 | 56.0 | **48,814.06** | +56.0 | 5 | **yes** | 3,011 |
| 7 | 48,611.4 | 202.6 | 48,712.2 | +100.7 | 7 | no | 3,998 |
| 8 | 48,598.1 | 216.0 | 48,674.8 | +76.7 | 8 | no | 3,965 |
| 9 | 48,770.4 | 43.6 | **48,814.07** | +43.6 | 7 | **yes** | 3,736 |
| 10 | 48,655.9 | 158.2 | 48,730.5 | +74.6 | 7 | no | 4,045 |
| 13 | 48,632.4 | 181.7 | 48,727.7 | +95.3 | 7 | no | 3,944 |
| 14 | 48,729.9 | 84.1 | 48,812.3 | +82.4 | 8 | no (1.75 short) | 3,976 |
| 15 | 48,753.5 | 60.6 | **48,814.06** | +60.6 | 5 | **yes** | 2,850 |
| 16 | 48,758.5 | 55.6 | **48,814.05** | +55.6 | 7 | **yes** | 3,681 |

* **Near misses (43–84 short at convergence): 4 of 5 fathomed**, in rounds 5–7
  (5,000–6,400 iterations). Node 14, 84 short, ends 1.75 short after round 8,
  still rising 3 per round.
* **Clear misses (158–405 short): +69 to +121, median +95**; none fathoms. The
  gain does not shrink with the gap (node 1, the worst, gained the most).
* The uncut 1,500-iteration bound is the converged value (node 1: 7,100
  iterations buy 0.2), so `first_lb` is a fair control; every gain is the cuts.
* **Cost**: 4,000–4,200 CPU-s per clear miss, 2,850–3,740 per fathom, against
  100–175 s for the same node under the search's give-up rule: 20–40×. The
  Polyak step overshoots on some nodes (node 8: −341 in round 3 and three
  rounds to recover; node 3: −294; node 7: −171), and 700 iterations per
  round were chosen for safety, not tuned.

## Assessment

The regional triangle cuts are the first thing that moves the *strength* of
the wil50 bound rather than its speed, and they move it on every node tried:
+95 at the median, three quarters of the ~130 units that one level of the
tree is worth (`../NOTES.md` section 5). Where the uncut bound is within
~85 of the cutoff the node fathoms outright: these depth-3 near misses are
what the depth-4 children of the clear misses look like, so the cuts should
move a large fraction of the depth-4 frontier up one level and, on the worst
tenth (node-1-like, 400 short), save the depth-6 tail. Modelled on the
`NOTES.md` tree that is 2–3× fewer node-equivalents, not the 47× a full
level gives: +95 against +130 sits on the threshold, and the outcome depends
on the gain at depths 4 and 5, which has not been measured.

At today's cost the cuts do not pay: 20–40× per node for 2–3× fewer nodes.
Two things decide whether they will, in this order:

1. **Cost.** Round 2 already gives +50 at 2,900 iterations; the warm-started
   rounds probably need far fewer than 700 iterations, and a damped step
   (smaller delta, or a step on the violated cuts only) would avoid the
   overshoot rounds. A protocol at ~3× the give-up cost on the marginal
   nodes looks reachable and is the next experiment: on nodes 6, 9, 15, 16,
   the same gain with `-i 200`–`300` and a smaller delta, measured in CPU.
2. **Selective use.** Run the cuts only where the uncut bound lands within
   ~100 of the cutoff (the converged near misses, and depth-4 children of
   clear misses after their give-up); branch the far ones as now. In the
   search that is a second phase after the give-up rule fires: if the gap at
   the give-up is below a threshold, add cuts instead of branching.

Then the memo's next inequalities on the same regions (mixed-region and box
inequalities, selective triple consistency, the regional second-order moment
block), each worth testing on the same 14 nodes: +30–50 more at the median
would make the full level. The exact fit of wil50's data to these regions
(Manhattan distance = 13 threshold indicators) is why this family works
where the atomic triangles did nothing (+0.2 on nug30); sko49 and tho40
share the grid structure, tai instances do not.

## Step 3: the cost of a fathom (`cost/run_cost.sh`, table by `cost/summarize_cost.py`)

The four near misses the cuts fathom, under cheaper protocols than step 2's
(1,500 uncut iterations, then 700-iteration rounds). Iterations and CPU to
the first certified bound above the cutoff; one core, four at a time:

| protocol | 6 | 9 | 15 | 16 | fathomed | fathoms: mean iterations, CPU |
|---|---|---|---|---|---|---|
| step 2: 1,500 uncut, 700-iteration rounds, delta 100 | 5,000 | 6,400 | 5,000 | 6,400 | 4 of 4 | 5,700, 3,320 s |
| 1,500 uncut, 300-iteration rounds, delta 100 | 3,900 | 3,000 | 3,900 | 3,600 | 4 of 4 | 3,600, 2,060 s |
| the same with a damped step, delta 30 | 4,200 | 3,600 | 0.003 short at 4,500 | 1.3 short | 2 of 4 | 3,900, 2,310 s |
| **from the give-up point**: 400 uncut, 300-iteration rounds, delta 50 | **2,500** | **1,900** | **3,400** | 0.003 short at 4,000 | 3 of 4 | **2,600, 1,510 s** |
| from the give-up point, 150-iteration rounds | no (48,795) | no (48,813) | no | no | 0 of 4 | – |

* Shorter rounds pay down to about 300 iterations: the warm-started PRSM
  does not need 700 per round (−37% CPU at 4 of 4). At 150 the certificate
  cannot follow the cuts: the primal sits 80–100 above the cutoff while the
  certified bound stalls below it, so nothing fathoms.
* Damping the step avoids the overshoot rounds but costs more rounds than it
  saves; delta 50–100 is right.
* Starting from the give-up point, where the search would hand the node
  over, is the realistic protocol: **1,900–3,400 iterations in all for a
  fathom, 1,100–2,000 CPU-s**, against about 250 s for the give-up itself
  (360–480 iterations under H = 300): **4–8× the give-up budget**, down from
  20–40× in step 2. Node 16 ends 0.003 below the cutoff twice (the same
  margin it fathomed by in step 2); it would fathom with one more round.

What this buys in the search, on the `NOTES.md` tree: the dominant work is
the ~8e5 depth-4 children of the clear misses, each a near miss that today
costs its give-up plus 47 depth-5 children (about 4,800 CPU-s). A cut phase
at 1,500 s that fathoms three of four of them replaces that with about
2,800 s per node; four of four, 1,600 s. That is **1.7–3× on the dominant
class** at the cost and gain measured here, provided the depth-4 gains match
the depth-3 ones, which is the next thing to measure (the 16 depth-4
children of nodes 6 and 9 in `../tuning/giveup/all32.txt`, ids 17–32, all of
which fathom uncut in 170–690 iterations, are not the right sample: a depth-4
sample under the clear misses is).

## Step 5: the cut phase in the search (`latest/dnncut.h`, `bnb.C`, `QAP_DNN_CUT_ZONE`)

The cut rounds are now a library (`latest/dnncut.{h,C}`, regions on either
side; the `dnncuts` driver runs over it and reproduces every earlier log to
the decimal), and `bnb.C` runs a **cut phase** after a DNN miss whose
certified bound ends within `QAP_DNN_CUT_ZONE` × |incumbent| of the cutoff:
from the node's own iterate, `QAP_DNN_CUT_ROUNDS` (12) rounds of
`QAP_DNN_CUT_ITERS` (300) iterations with `QAP_DNN_CUT_K` (1,000) new regional
triangles each and a Polyak step toward best + `QAP_DNN_CUT_DELTA` (50),
stopping at the cutoff; the regions are the prefixes of the
`QAP_DNN_CUT_GRID=R:C` grid over the original indices of side
`QAP_DNN_CUT_SIDE` (wil50: `5:10`, side 0; nug30: `5:6`, side 1, the second
matrix). The children keep the uncut multiplier; the results block counts
`cut phase: nodes, fathomed, CPU`. Off by default.

First runs:

* nug20, `default.par -u 2570 -x 2:2:1:300:1`, PRSM + degree 4,
  `QAP_DNN_CUT_ZONE=0.02 QAP_DNN_CUT_GRID=4:5 QAP_DNN_CUT_SIDE=1`: the two
  depth-2 misses get a cut phase, **one fathoms** (15 s for both phases),
  84,592 → 54,718 nodes.
* **nug30's two depth-2 misses, standalone** (`dnncuts -grid 5:6 -side 1`,
  800 uncut iterations then 300-iteration rounds, delta 20): `4:7,6:14`
  (15.1 short) **fathoms in round 4** (2,000 iterations, 86 CPU-s);
  `4:29,6:7` (26.3 short) reaches 6,121.59 in round 16 (+25.9, 0.46 short,
  primal 6,124.5, still rising). The atomic triangles had given +0.2 and +0.4
  on nug30 depth-2 misses (`../../dnn_cuts/README.md`); the regional family on
  the grid of the second matrix is a different thing. On nug30 the depth-2
  survivors carry the tree (58 of 225 under `nug_deep.par`, 1,552 depth-3
  children), so a cut phase at depth 2 is the next nug30 A/B.

## Step 4: depth 4 (`depth4/`, table by `depth4/summarize_depth4.py`)

Sixteen depth-4 children of four clear misses (nodes 10, 13, 2, 7 of
all32: 158, 182, 202, 203 short at convergence; one random free facility and
four random free locations per parent, seed 20261101), each run uncut to
1,500 iterations (`dnntest -v`, cold start) and with cuts from the give-up
point (400 uncut iterations, then 300-iteration rounds, delta 50, up to 12
rounds). Cutoff 48,814.05.

| | nodes | what happens |
|---|---|---|
| fathom uncut | **7 of 16** | 4 within 140–240 iterations, 3 slow ones at 510–630 (0.1–0.5 above the cutoff at convergence); with cuts the slow ones fathom in 4 rounds (1,600 iterations) |
| near misses, 1–22 short at convergence | 2 | **both fathom with cuts** (rounds 9, 3,100 iterations, 1,750 CPU-s) |
| near misses, 62–109 short | 7 | none fathoms: gains +64 to +105 over the give-up-point bound, ending 25–61 short after 12 rounds (4,000 iterations, 2,200 CPU-s each) |

Two things the sample settles. The level rule holds: parents 158–203 short
give children −6 to +109 short, about 130 units closer, and 44% of them
fathom uncut, which the `NOTES.md` model (children "170–370 short") was
pessimistic about. And the cut phase at depth 4, from a cold give-up point,
fathoms only the children within about 25 of the cutoff at convergence
(about 60 at the give-up point): the gains are the same as at depth 3 (+64 to
+105), but a cold 400-iteration bound at m = 46 sits 48–55 below its
converged value (26–36 at depth 3), so the 62–109-short children end 25–61
short. In the search the children are warm-started from the parent and the
stall rule stops a non-fathoming phase after three or four rounds, which
changes both numbers in the right direction, but the depth-4 picture is
**a 10–25% saving on the dominant class, not the 1.7–3× of the step-3
model**: on these 16 children, the cuts remove two of nine near-miss
subtrees (2 × 47 depth-5 children) for the price of nine cut phases.

What would change it: a larger gain per node (the memo's next inequalities;
+130 at the median would fathom the 62–66-short children), or cheaper
rounds. The cut phase from a *converged* bound (step 2's protocol) fathomed
depth-3 nodes 43–84 short at 2,850–3,980 CPU-s; on the 62–66-short depth-4
children that is about even with their uncut subtrees (47 depth-5 children at
~100 s), so a two-stage policy (give-up → converge → cuts, for nodes within
~70 of the cutoff at convergence) is the best use of the family at depth 4
as it stands.

### The cut phase in the search: A/B (`search/run_ab.sh`, cloud Xeon 2.1 GHz, 4 workers, PRSM + degree 4)

| run | nodes | DNN nodes / fathomed | total CPU |
|---|---|---|---|
| wil50, subtree below near miss 6 (`-F 20:12,40:28,17:25 -x 0:4:1:300:1`), base | 48 | 48 / 47 | 3,565 s |
| the same, `QAP_DNN_CUT_ZONE=0.002 QAP_DNN_CUT_GRID=5:10` | **1** | 1 / 0, the phase fathoms the node | **1,576 s** |
| nug30 `nug_deep.par -x 2:7:1:300:1`, base | 4,528 | 4,518 / 4,355 | 11,748 s |
| the same, `QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_GRID=5:6 QAP_DNN_CUT_SIDE=1 QAP_DNN_CUT_DELTA=20` | 3,342 | 3,332 / 3,169 (+46 by the phase) | 16,723 s (+42%) |

(`bench/parallel.py` merges statistics through the state files, which do not
carry the cut-phase counters, so the merged block has no "cut phase" line;
phase fathoms = fathomed − DNN fathomed.) On the wil50 near-miss subtree the
phase fathoms the node at depth 3 instead of bounding its 47 children:
**2.3×**. On nug30 the phase removes 26% of the nodes and costs 42% more
CPU; Fable's laptop A/B on the refined setting says the same (+34% at every
depth, break-even restricted to depth 2 with `QAP_DNN_CUT_MAXDEPTH`,
`../tuning/README.md`). The difference is per-node cost: a nug30 depth-3
child costs under a second, a wil50 depth-4 child about 75 s. The cut phase
is a wil50 tool, and the research memo (`research_bound_quality.md`) says
where its next gains are: the dual method and the certificate's lag, before
new inequalities.

## Step 6: the dual method's cheap fixes (`dual/run_dual.sh`, `dual/run_dual2.sh`)

`research_bound_quality.md` read the step-2/3/4 logs and found three losses
in the multiplier method, each with a cheap remedy: the Polyak target
`min(cutoff, best + delta)` makes the step vanish within 0.01 of the cutoff
(fathoming nodes crawled 4-6 rounds there); the certificate trails the
Lagrangian primal by 25-90 at 300-iteration rounds (87-108 at 150); and the
step overshoots (-100 to -390 rounds). `dnncut.C` gained `-nocap` (target =
best + delta), `-freeze R:N` (a quiet round of N iterations with the
multipliers frozen and no new cuts, after round R or at the end) and the
driver `-t theta`; all protocols below start from the give-up point (400
uncut iterations) unless stated.

| experiment | nodes | result |
|---|---|---|
| **uncapped target**, 300-iteration rounds, delta 50 | near misses 6, 9, 15, 16 | **4 of 4 fathom in 1,000-2,800 iterations** (step 3: 1,900-3,400 and one miss at 4,000); the bounds cross the cutoff cleanly (48,814.3-48,815.5) instead of crawling |
| | node 14 (84 short) | +84, ends 37.6 short (step 3: +77, 61 short) |
| **quiet round** (600 it.) after the 150-iteration protocol that fathomed nothing in step 3 | 6, 9, 15, 16 | **4 of 4 fathom**: the certificate goes from 48,785-48,813 to 48,814.5-48,817.3 in the one frozen round |
| quiet round after 12 capped 300-iteration rounds | depth-4 children 1, 2, 4, 11 (62, 109, 66, 66 short) | +35, +35, +17, 0 over the capped run; child 1 ends **11.5 short**, child 2 47, child 4 40, child 11 60 (its last round was an overshoot, so the frozen multipliers were poor); none fathoms |
| 5,000 cuts per round instead of 1,000 (step-2 protocol) | clear misses 2, 7 | +65 and +98 against +89 and +101: no gain, the Polyak step gets coarser |
| theta 0.4 instead of 1 (step-2 protocol) | overshooters 3, 8 | +75 and +74 against +69 and +77: the -337 and -389 rounds disappear, the eight-round value barely moves |
| both fixes together (uncapped + quiet round) | depth-4 children 1, 4, 11 and node 14 | identical to the quiet-round-only runs on the children (the cap is inactive while best + delta is below the cutoff, i.e. more than 50 short): child 1 **11.5 short**, child 4 40, child 11 60; node 14 goes from 37.6 short (uncapped alone) to **8.2 short** with the quiet round (+113.8 over the give-up point). None fathoms. The rounds still lose 2-3 overshoots each (-69 to -138), which is what theta 0.4 removes; theta 0.4 + quiet round on child 1 and node 14 is the next run |
| theta 0.4 + quiet round (`dual/run_dual3.sh`) | child 1, node 14 | child 1 **48,803.0, 11.1 short** (theta 1: 11.5); node 14 **48,805.3, 8.7 short** (theta 1: 8.2). No overshoot rounds, and the destination barely moves: the present cuts with Polyak rounds reach about 48,803 on the 62-short child whatever the step, so the last 9-12 units are this multiplier method's ceiling, which the memo's section-4 protocol measures and the 5-point count inequalities are meant to pass |

The two fixes are now the search's cut-phase defaults (`bnb.C`:
`cap_at_cutoff = 0`, a 600-iteration quiet round when the phase ends
unfathomed with its primal above the cutoff; `QAP_DNN_CUT_NOCAP=0` and
`QAP_DNN_CUT_FREEZE=0` restore the old behaviour). What they change: a
near-miss fathom from the give-up point now costs **2.5-3x the give-up
budget** (1,000-1,300 iterations on three of the four), down from 4-8x in
step 3 and 20-40x in step 2; and the inner lag is real and recoverable, which
moves the depth-4 62-66-short children to within 12-40 of the cutoff with
the present cuts. The pool size and the step factor are not the limit; the
memo's next items (a bundle method on the rounds, then the 5-point count
inequalities) are.

**Laptop, same day (`laptop/README.md`):** with the quiet round started from
the *best* round's multipliers (`DNNCutParams::freeze_best`, now the default;
from the last round's an overshoot freezes a bad multiplier) and 900-1,600
quiet iterations, **four of four depth-4 children 62-99 short at convergence
fathom** (d1, d4, d6, d11; +117 to +134 each, a full level) at 690-1,050
laptop CPU-s, about half of their 46 uncut depth-5 children: **about 2x on
the dominant class**. In the search, the cut phase uncapped with the quiet
round is 4.5-7x on the depth-3 near misses 6, 9, 15 (node 16 did not fathom
uncapped and is the one to re-test with `freeze_best`). The zone has to be
about 0.003 to reach children 117-134 short at the give-up point, and the
stall rule carries the failures in that band (600-700 CPU-s each). Next: a
clear-miss subtree in the search with zone 0.003 (`search/run_ab2.sh`), the
end-to-end test on the dominant class.
