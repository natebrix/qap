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
