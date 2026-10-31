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

RESULTS PENDING

## Step 2: the 14 misses (`run_misses.sh`, logs in `misses/`)

PENDING

## Assessment

PENDING
