# Nonlinear valid constraints for the wil50 DNN bound: nothing survives (Nov 6)

Question (roadmap, Nov 6): are there valid constraints that are *not* linear
in the lifted matrix Y, or linear only in a lift we do not form, that would
tighten the facially reduced DNN bound at wil50's depth 3-6 nodes, and that a
PRSM or the Lagrangian cut phase (`latest/dnncut.C`) could impose cheaply?
The residual gap after the triangle phase is 117-143 units at depth 3
(`../grid/REPORT.md`), and every linear family on Y tried so far is closed
(`../SYNTHESIS.md`, `../box_pool/README.md`).

**Answer: no.** Every nonlinear constraint on Y that is valid for all
permutation lifts and convex is implied by the three ingredients the DNN
already has (Y positive semidefinite, Y nonnegative, the lifted row and column
sums), and we can say why in one paragraph (section 1). The content that is
*not* implied lives at degree 3 and higher in x (triples of facilities), which
is the RLT2 / Lasserre-2 / disjunctive territory the Nov 4 study closed. One
linear family that had not been looked at, the k-partition pigeonhole
("clique") inequalities on the grid's rows, is violated at every saved iterate
(section 3), but it decomposes into box inequalities on single rows and
columns whose violation mass equals that of the threshold boxes that failed
their pool gate on Nov 4, and whose strongest members *are* those threshold
boxes; a fixed-multiplier screen in a bound confirms it (section 4). No
switch was added to the solver: nothing earned a prototype under the rules
set in advance.

Tags: **[m]** measured, with the log named; **[e]** estimated. Node specs
are `cell:facility` on `qaplib/wil50.dat` (first index the 5 x 10 grid, second
the flows). Iterates: node 13 = `21:28,8:41,28:19`, node 10 =
`39:27,30:48,25:49`, d2 = node 10 + `26:5` (m = 46), all dumped by the Nov 4
agents (`../grid/runs/*.Y.bin`, `../facility/runs/*/quiet.Y.bin`,
`../facility/runs/n10c/parent.Y.bin`; the binary dumps are not in git and
live in the main checkout). Scripts here: `nlc_analyze.py` (sections 1-3,
every candidate at one iterate), `line_boxes.py` and `interior_mass.py`
(section 3, violation mass), `make_line_ladd.py` and `run_lines_fixedmu.sh`
(section 4). Logs in `logs/`. Nothing was run on the laptop or on AWS.

## 1. Why degree-2 nonlinear constraints are implied (the argument)

Write y = (1, vec X), Y = E[y y'] under any probability measure on
permutations (the convex hull of the lifts is what a bound can see). Six
facts about the DNN set D = {Y ⪰ 0, Y ≥ 0, lifted row/column sums, gangster}:

1. **Every function of one facility's location is linear in x.** For facility
   i and any g: g(loc_i) = Σ_a g(a) x_ai. So for any vector of such functions
   u = (1, g_1(loc_1), ..., g_m(loc_m)) = Φ y, the moment matrix E[u u'] =
   Φ Y Φ' is a *compression* of Y: PSD whenever Y is. Every PSD, SOC or
   Jensen-type constraint on the first and second moments of single-facility
   quantities (regional memberships z_i^S, coordinates r_i, c_i, their powers,
   counts q = Σ_Q z_i) is therefore implied. (This is the memo's remark that
   [1 p'; p W] ⪰ 0 is redundant, in general form.) What is *not* implied about
   these moments is the integrality of u, and a polytope's integrality is
   expressed by linear facets: triangles, 5-point, local-7, boxes, all
   tested. The "integer variance" bound Var(q) ≥ θ(1-θ), θ the fractional
   part of E q, looks nonlinear but is exactly the piecewise-linear box family
   (q-b)(q-b-1) ≥ 0.
2. **Every pair block is an exact joint distribution.** For facilities i ≠ j
   the m x m block Y_(ai),(bj) is nonnegative, zero on the diagonal (gangster)
   and has marginals x_·i and x_·j (lifted sums); any such matrix is a convex
   combination of point masses (a, b), a ≠ b, each realised by some
   permutation. So every constraint that is valid for all permutations and
   involves one pair block (Jensen on distances, (E d)^2 ≤ E d^2 with
   d^2 = (r_i - r_j)^2 both linear in the block, moment bounds on the flow
   between two cells, ...) is implied.
3. **Spectral constraints on conjugates are majorization, and majorization
   is implied.** For permutations X B X' has B's spectrum, so for the
   relaxation M(Y) = E[X B X'] (linear in Y: M_ab = Σ_ij B_ij Y_(ai),(bj)) the
   valid convex statement is λ(M) ≺ λ(B) (Ky Fan: top-k sums of eigenvalues
   bounded by B's, trace equal). Proof that D implies it: with B = Σ_l λ_l
   u_l u_l', M = Σ_l λ_l P_l where P_l = E[(X u_l)(X u_l)'] ⪰ 0 is a
   compression of Y, tr P_l = E||X u_l||^2 = u_l' E[X'X] u_l = 1 and Σ_l P_l =
   E[X X'] = I, both from the gangster entries and the lifted sums; for any
   k-dimensional projection Π, tr(Π M) = Σ_l λ_l tr(Π P_l) with tr(Π P_l) ∈
   [0, 1] summing to k, hence ≤ the k largest λ_l. The same argument gives
   λ(E[X' A X]) ≺ λ(A) and, for any aggregation T of cells into super-cells,
   λ(E[X' T'T X]) ≺ λ(T'T): the eigenvalue bounds on the coarse problems are
   implied too. Trace-power identities tr((XBX')^k) = tr(B^k) are degree 2k in
   x; their only convex degree-2 consequence, ||E[XBX']||_F ≤ ||B||_F, is a
   case of majorization. The Anstreicher-Brixius QPB construction (B⊗A =
   I⊗S + T⊗I + Q, Q ⪰ 0, objective ≥ tr S + tr T + x'Qx) is implied by Y ⪰ 0
   through the Schur complement Y_xx - xx' ⪰ 0 and the gangster/sum identities
   <I⊗S, Y_xx> = tr S; the SDP relaxation dominating QPB is the classical
   result.
4. **The partial transpose and the realignment.** Two reshapings of the
   m^2 x m^2 block Y_xx that are *not* compressions: the partial transpose
   (Y^{T_B})_(ai),(bj) = Y_(aj),(bi) and the realignment R_(ab),(ij) =
   Y_(ai),(bj) = E[X ⊗ X]. For a permutation lift, Y^{T_B} is a symmetric
   orthogonal involution (eigenvalues ±1) and R is a permutation matrix, so
   -I ⪯ Y^{T_B} ⪯ I and ||R||_op ≤ 1 are valid, genuinely nonlinear (LMI)
   constraints that I had not seen used for QAP. Both are implied: Y^{T_B} is
   symmetric, nonnegative (its entries are entries of Y) and has unit row sums
   (Σ_bj Y_(aj),(bi) = Σ_j x_aj = 1 by the lifted sums), so it is doubly
   stochastic and Perron-Frobenius puts its spectrum in [-1, 1]; R is doubly
   stochastic for the same reason. Checked at the iterates (section 2):
   λ_max(Y^{T_B}) = 1.00001, the all-ones eigenvector, second eigenvalue
   0.83-0.84; σ_max(R) = 1.00001.
5. **Objective-level constraints.** The cutoff half-space <L, Y> ≤ c changes
   nothing: for any cut family with multipliers μ and a multiplier ν ≥ 0 on
   the half-space, the Lagrangian dual is max_μ L(μ) + max_ν ν (L(μ) - c),
   which equals max_μ L(μ) whenever every L(μ) < c, i.e., whenever the node is
   not already fathomed. Reduced costs (certificate child bounds) are a dual
   object the Nov 4 study measured at 0 of 2,209 eliminations and are the
   sibling agent's topic. Parity: wil50's objective is even and d(a, b) ≡
   parity(a) + parity(b) mod 2 on the grid, so f/2 ≡ (number of odd-degree
   facilities on odd cells) mod 2; with 24 odd-degree facilities among 50 both
   residues of f mod 4 occur (2,000 random permutations: {0, 2}), so there is
   no free shift of the cutoff from 48,814.05 to 48,812.05.
6. **Copositive and rank-one reasoning.** The exact cone is completely
   positive (Povh-Rendl); for any compression E[u u'] with binary u the CP
   slice with M_0i = M_ii is exactly the Boolean quadric polytope (Burer 2009),
   i.e. linear facets again; the K_r inner approximations of the copositive
   cone need N^3 = 10^10 variables at N = 2,210. Rank-one identities on Y
   (vanishing 2 x 2 minors Y_ab Y_cd = Y_ad Y_cb) are degree 4 in x and their
   linearisation is RLT; the lift-and-project closure of a line is its worst
   child (facility agent, +42 to +59 at node 10); a disjunctive cut from a
   linear split v'y ≤ c or ≥ c needs a gap in the achievable values of v'y,
   which only integer counts have (count splits: +2.4 to +2.9, closed), not
   the eigenvector directions where the relaxation's second moment deviates
   most from rank one.

So at degree 2 the DNN is "complete" for nonlinear convex constraints, and the
missing strength is a hierarchy gap (degree ≥ 3), consistent with the three
Nov 4 reports. (Had a convex nonlinear constraint g(Y) ≤ 0 survived, imposing
it would have needed no new cone: its supporting hyperplane at the current
iterate, g(Y_k) + <∇g(Y_k), Y - Y_k> ≤ 0, is a linear cut on Y, and the pool
with its Polyak multipliers and certificate takes it as it takes a triangle;
an LMI like -I ⪯ Y^{T_B} ⪯ I would enter as the rank-one cut v' Y^{T_B} v ≤
1 of its top eigenvector. The cost would have been that of the dense box
terms, which is what sank the boxes in the pool.) Literature to Oct 2026 (`logs/arxiv_qap_2023plus.txt`,
export.arxiv.org and optimization-online.org searched Nov 6): the 2024-2026
QAP relaxation papers are exactness conditions for the SDP (arXiv
2409.08802, 2408.05942), error bounds for rank-one DNN reformulations with
exact penalties (2403.11210), graph-alignment phase transitions (2606.15581)
and approximation results (2604.12584); no relaxation stronger than the DNN
at n ≥ 40, no new valid constraint class for the QAP lift. Partial-transpose
and realignment criteria appear only in the entanglement literature.

## 2. Numerical check of the implied families at the saved iterates [m]

`nlc_analyze.py` on six iterates (`logs/nlc_*.log`; node 13: tri5 = the grid
agent's 5-round triangle phase, quiet = the facility agent's 8 rounds + 1,500
quiet, uncut = 1,500 uncut iterations; node 10: tri5 and uncut converged; d2:
tri5). Residuals of the iterates: lifted sums off by 3e-5, row sums of x by
8e-4, entries ≥ 0, gangster exact.

| check | node 13 tri5 | node 13 quiet | node 13 uncut | node 10 tri5 | node 10 uncut | d2 tri5 |
|---|---|---|---|---|---|---|
| majorization E[XBX'] ≺ B, largest top-k excess | -3e-13 | - | - | 1e-13 | - | 1e-13 |
| majorization E[X'AX] ≺ A | 1e-13 | - | - | 2e-13 | - | -2e-13 |
| ‖E[XBX']‖_F against ‖B‖_F | 222.4 / 251.7 | | | 219.9 / 249.4 | | 214.7 / 243.8 |
| λ_max(Y^{T_B}), second eigenvalue | 1.00001, 0.834 | - | - | 1.00001, 0.844 | - | 1.00001, 0.843 |
| λ_min(Y^{T_B}) | -0.834 | | | -0.844 | | -0.843 |
| σ_max, σ_2 of the realignment | 1.00001, 0.99988 | | | 1.00002, 0.99982 | | 1.00002, 0.99986 |
| spread, k = 4: Dmax(4) - max over 4-sets of Σ E d | +5.3 | +5.3 | +5.5 | +4.8 | +5.2 | +4.8 |
| spread, k = 3, 5 (Manhattan) | +1.1, +4.0 | | | | | |
| rows-only / cols-only spread, k = 3..5 | +0.4..+1.6 / +0.4..+2.6 | | | | | |

(The last rows: for any k cells the total pairwise L1 distance is at most
Dmax(k), computed exactly by enumeration over the node's free cells; its
expectation is linear in Y. Over-spreading of light-flow sets is the
direction in which the DNN cheats on average, facility report A.3, but no
k-set is over-spread beyond the feasible maximum; the Nov 4 reports covered
the opposite direction, over-clustering, with the same result.)

The spectrum of Y at node 13's tri5 iterate: top eigenvalue 5.24, 567
eigenvalues above 1e-2, trace 48: the relaxation's optimum is far from rank
one, so rank-based reasoning has much to cut but no convex handle (section
1.6).

## 3. The one violated family: k-partition pigeonhole inequalities [m]

For a partition of the free cells into parts with capacities (the 5 rows, the
10 columns, the 12 mid5 super-cells), P_ij = probability that facilities i
and j share a part (= Σ_parts W^{part}_ij, linear in Y, a cross-threshold
object: same row = equal membership in all four row thresholds). Any t
facilities with t > #parts have at least f(t) same-part pairs, f(t) = min over
integer loads n_p ≤ cap_p, Σ n_p = t of Σ_p C(n_p, 2) (Chopra-Rao clique
inequalities of the k-partition polytope, with capacities). The 2-part case
with t = 3 is the sum triangle. Minimum over t-sets by greedy plus 1-swap
local search (an upper bound on the true minimum, so violations may be
larger):

| partition, t | f(t) | node 13 tri5 | node 13 quiet | node 13 uncut | node 10 tri5 | node 10 uncut | d2 tri5 |
|---|---|---|---|---|---|---|---|
| rows, t = 6 | 1 | 1.100 | 1.129 | 1.134 | 1.163 | 1.152 | 1.156 |
| rows, t = 7 | 2 | **1.894 (-0.106)** | 1.977 (-0.023) | 1.924 (-0.076) | 1.997 (-0.003) | 1.929 (-0.071) | 2.009 |
| rows, t = 8 | 3 | **2.973 (-0.027)** | 2.972 (-0.028) | 2.966 (-0.034) | 2.954 (-0.046) | 2.972 (-0.028) | 2.956 (-0.044) |
| rows, t = 9..11 | 4, 5, 7 | +0.15, +0.53, +0.10 | | | +0.22, +0.62, +0.20 | | +0.22, +0.65, +0.30 |
| columns, t = 11..16 | 1..6 | +0.79 .. +0.12 | | | +0.71 .. +0.03 | | +0.76 .. +0.13 |
| mid5 super-cells (12 parts; 11 at nodes 10, d2), t = parts+1..+6 | 1..6 | +2.3 .. +2.7 | | | +1.7 .. +2.1 | | +1.9 .. +2.3 |

So the exact-row partition is violated at every iterate, by 0.003-0.106 pairs
(the row problem's integrality: seven facilities spread over five rows with
too few same-row pairs), the columns and super-cells never. Not implied by
PSD alone: Jensen on the loads gives only Σ_p (μ_p^2 - μ_p)/2 ≥ 0.6 at t = 7
against the integer 2.

**It decomposes into box inequalities on single rows.** Σ_pairs P_ij = Σ_p
E[C(n_p, 2)], and the box inequality (n_p - b)(n_p - b - 1) ≥ 0 on the exact
row p with the same facility set T and b = ⌊E n_p⌋ gives E[C(n_p, 2)] ≥ g(E
n_p), g the piecewise-linear interpolation of C(n, 2) at the integers; a
separable convex piecewise-linear function with integer breakpoints attains
its minimum over {Σ μ_p = t} at an integer point, so Σ_p g(μ_p) ≥ f(t). Hence
the clique inequality is a sum of five "line boxes", and the line-box family
on the 15 exact rows and columns is at least as strong. `line_boxes.py` runs
the pool's box separation (`dnncut.C` type 8: add/drop local search, 40
starts per region and b) on the lines and on the 13 thresholds at the same
iterate (`logs/linebox_*.log`):

| iterate | line boxes: distinct violated, top-200 sum, max | threshold boxes (killed Nov 4): same | threshold triangles: top-200 sum, max, top-1000 sum |
|---|---|---|---|
| node 13 tri5 | 663, **8.23**, 0.075 | 580, **8.35**, 0.080 | 5.49, 0.063, 17.8 |
| node 13 quiet | 682, 7.98, 0.076 | 596, 8.40, 0.074 | 5.28, 0.054, 15.0 |
| node 10 tri5 | 577, 6.91, 0.081 | 540, 8.11, 0.081 | 6.82, 0.062, 20.7 |

Per line at node 13 tri5: rows 0 and 4 and columns 0 and 9 carry 371 of the
663 violated boxes and the four largest violations (0.066, 0.072, 0.075,
0.063); the interior rows 1-3 have 5-48 violated boxes each with maxima
0.027-0.047 and the interior columns 10-44 each with maxima 0.021-0.047. The
edge lines are the edge thresholds or their complements (row 0 = row<1,
column 0 = col<1, row 4 = the complement of row<4, column 9 of col<9), and a
box on a region and on its complement are the same inequality (q ↦ |Q| - q,
b ↦ |Q| - b - 1). So the strongest members of the family are exactly the
threshold boxes that lost -32.3 and -22.6 in the pool gate
(`../box_pool/README.md`), and the whole family's violation mass is that
family's (8.2 against 8.4 on the top 200). The interior lines alone
(`interior_mass.py`, `logs/linebox_n13_tri5_interior.log`): 268 violated
boxes, top-200 sum **4.1** against the threshold boxes' 8.4, median of the top
200 0.020 against 0.039, maximum 0.050 against 0.080. The predeclared
criterion for a prototype in `dnncut.C` was a violation mass clearly above the
killed family's on the genuinely new members: they have half of it. Killed on
that count; section 4 closes it in a bound as well.

## 4. Fixed-multiplier screen of the interior-line boxes, node 13 [m]

The literature agent's on-top screen (`../literature/wil50/run_box_ontop.sh`):
the grid agent's tri5 protocol (400 uncut iterations, 5 Polyak rounds of 300
with 1,000 triangles a round, target best + 50 uncapped, quiet round of 600)
with an external term (`dnncuts -Ladd`) of the 200 most violated interior-line
boxes separated at the tri5 iterate, μ = 0.3 in the doubled form (G = 2b on
the diagonal, -1 off it, h = b(b+1), the same convention as
`make_box_ladd.py`), every cut checked on 300 random permutations and one per
count before it is written; a matched control with the same binary and
bundle-D settings (`run_lines_fixedmu.sh`, host in `logs/host.txt`).
**Predeclared rule: KILL unless the quiet-round certified gain over the control
is at least +14**, the threshold boxes' gain in this very screen (48,728.10
against 48,713.85, `../literature/wil50/runs/n13_ontop_mu0.3_k200.log`),
which then failed the pool gate.

Result (`logs/n13_tri5_control.log`, `logs/n13_lines_mu0.3_k200.log`, term
construction in `logs/n13_lines_mu0.3_k200.mk`; 2.10 GHz Xeon with AMX, one
core, 0.34 s per iteration):

| arm | round 0 (400 it.) | best Polyak round | quiet round, certified | primal at the quiet round | process CPU s |
|---|---|---|---|---|---|
| control: tri5 protocol, this binary | 48,597.17 | 48,686.38 (round 3) | **48,713.95** | 48,735.35 | 895.7 |
| + 200 interior-line boxes at μ 0.3 (446,577 nonzeros, const 261.6) | 48,595.30 | 48,698.31 (round 5) | **48,724.58 (+10.6)** | 48,746.91 (+11.6) | 900.6 (+0.5%) |
| for comparison, Nov 4: + 200 threshold boxes at μ 0.3 (856,065 nonzeros), against its control 48,713.85 | 48,598.85 | 48,700.36 (round 5) | 48,728.10 (+14.3) | 48,752.51 | - |

**KILL: +10.6 < +14.** The interior-line boxes behave as a weaker copy of the
threshold boxes: the same shape of run (the control overshoots in round 4,
-17, and freezes round 3's 2,862 multipliers for the quiet round; the term
arm does not overshoot and freezes round 5's 4,680, so part of the +10.6 is
the Polyak path rather than the family, the confound the literature agent
noted on Nov 4), the same order of fixed-multiplier gain, and a family of the
same type then lost -32.3 and -22.6 when its multipliers joined the Polyak
rounds (`../box_pool/README.md`). Today's control reproduces the Nov 4 control
to 0.1 (48,713.95 against 48,713.85 on a host without AMX), so the screen is
stable. No pool gate was run: a pass here was the condition for one.

## 5. What is closed now, and why

- **Spectral / majorization / Ky Fan / trace-power / PB-QPB-type constraints
  on E[XBX'], E[X'AX] and on any coarse aggregation:** implied by PSD +
  gangster + lifted sums (section 1.3, proof; numerically to 1e-13).
- **Partial transpose and realignment LMIs** (-I ⪯ Y^{T_B} ⪯ I, ||E[X⊗X]||
  ≤ 1): implied by nonnegativity + lifted sums (doubly stochastic); the only
  eigenvalue at 1 is the all-ones vector's, the next is 0.83.
- **SOC / small PSD / Jensen constraints on regional or coordinate moments,
  on counts, on single pair blocks:** compressions of Y or statements about
  one exact distribution: implied (section 1.1-1.2). The regional M2 block
  (degree 4) was closed Nov 4 by saturation.
- **Integer-variance bounds on counts:** the box family in disguise; boxes
  on thresholds killed Nov 4, on interior lines killed here (sections 3-4).
- **Objective-level:** the cutoff half-space adds nothing (dual argument);
  no parity shift of the cutoff (f mod 4 takes both values); reduced-cost
  eliminations closed Nov 4.
- **Copositive inner approximations, rank-one / intersection / disjunctive
  cuts on spectral sets:** linear facets (BQP) or RLT at degree 4, or capped
  by the worst child; closed with the Nov 4 arguments (section 1.6).
- **k-partition clique (pigeonhole) inequalities on rows, columns and
  super-cells:** rows violated by 0.003-0.106 pairs at every iterate, columns
  and super-cells never; the family is a sum of line boxes, whose violation
  mass and strongest members are the killed threshold boxes'; the genuinely
  new members (interior lines) give +10.6 certified at a fixed multiplier on
  node 13 against the threshold boxes' +14.3 in the same screen (section 4):
  killed by the predeclared rule.
- **Spread (Dmax) and perimeter inequalities:** never violated (slack 1-5).
- **Literature 2023-Oct 2026:** nothing stronger than the DNN at n ≥ 40.

What this means for the plan: the bound side of wil50 at depths 3-6 is at the
degree-2 ceiling; the remaining 117-143 units need degree-3 information
(coarse RLT2 on super-cells, 1-2 weeks, 0 to +30 estimated by the literature
agent) or branching, which the measured levers already exploit. The next
bound-side test worth its cost is not a cut family: it is the sampling plan's
tail nodes, where the per-level gains of the triangle phase (+117 to +134 at
depth 4) decide the tree.
