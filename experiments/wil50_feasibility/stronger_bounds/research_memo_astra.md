# Stronger bounds for wil50: regional moment SDPs and selective consistency cuts

**Research memo for Nathan Brixius and Fable — 31 October 2026**

**Repository examined:** [natebrix/qap, cognac](https://github.com/natebrix/qap/tree/cd43c7645d8af01a2ec22aa0ae9bba9fafea7022), commit `cd43c7645d8af01a2ec22aa0ae9bba9fafea7022`.

## 1. Assessment

There is a mathematically sound route to a stronger bound that preserves the large DNN block and adds substantially smaller objects: **enforce higher-order consistency of membership in sets of locations**.

For wil50, the natural sets are the 13 threshold regions that exactly decompose its 5×10 Manhattan distance matrix. I recommend a progression:

1. Add regional Boolean-quadric, box, and hypermetric inequalities through the existing `DNNParams::L_add` interface.
2. Add selected constraints requiring three facilities' pairwise assignment marginals to admit a common joint distribution. This captures part of RLT2 without allocating the full relaxation.
3. Add second-order moment SDP blocks for selected regions, and possibly small blocks involving several regions.

The third item is the substantive stronger SDP formulation. Its root-level regional matrices have dimension **1,276**, versus **2,402** for the existing facially reduced DNN matrix and approximately **three million** for a direct second-order assignment moment matrix before further linear reductions. The additional regional blocks have a known trace, supporting spectral corrections to approximate dual certificates.

I established strict strengthening with an exact rational example and a six-facility QAP with nonnegative data: ordinary DNN has value **9**, while one regional cut raises the bound to the true optimum **10**. The example satisfies every individual-assignment Boolean triangle. Thus the disappointing atomic-triangle experiment in the repository does not settle the regional question.

The wil50 pilot is also positive: on the difficult depth-three node `34:40,28:22,14:29`, a pool of 1,000 regional triangles with fixed nonnegative multipliers produced a certified lower bound of **48,449.0658**, versus the strongest computed uncut certificate of **48,409.2349**. The gain is **39.8309**. This is a bound on that conditional subproblem, **not a new root lower bound for wil50**. It does not fathom the node, and is still below the estimated 130–260-unit target.

This is not a claim that wil50 is solved, that a new best-known wil50 root bound has been established, or that the full regional moment formulation has been benchmarked. The distinction between a proved stronger relaxation, a computed improvement, and an effective branch-and-bound improvement is essential.

## 2. Baseline and literature

The relevant target is improvement at difficult nodes, not merely an attractive root value. The repository's `experiments/wil50_feasibility/NOTES.md` and subsequent convergence logs suggest that roughly **130–260 additional bound units** at depths 3–5 could materially change the search. Those numbers are extrapolations from samples and tree shape, not a guaranteed conversion from bound gain to runtime.

The literature baseline also needs care:

| Source or method | Relevant observation | Interpretation |
|---|---|---|
| QAPLIB, Montréal mirror [1] | wil50 incumbent 48,816; listed lower bound 48,121 | A reference point, not proof that no later bound exists |
| Graham et al., rPRSM manuscript, Table 4.3 [2] | Published numerical lower bound 48,126 | Beating the older 48,121 listing by a few units would not establish a new record |
| Fujii et al., Table 5 [3] | Estimated wil50 trees of 20,523,849–39,330,160 nodes under their branching rules | Supports the scale of the problem; these are estimates, not completed wil50 proofs |
| Mittelmann–Peng, Table 3 [4] | Manhattan-specific small SDPs give 46,526 and 46,467 | Useful structural ideas, but these computed bounds are weaker than modern DNN on wil50 |
| Date–Nagi, Table 5 [5] | RLT2 dual ascent gives nug30 5,846; nug20 2,520 | Performance depends on the instance and the accuracy of the relaxation solve; do not assert a general ordering between full RLT2 and DNN |
| Hahn et al. [6] | Their RLT3 coefficient storage is about 10.1 TB at order 50 | Full RLT3 is an unattractive first implementation target |

Changing ADMM to PRSM, Halpern iteration, Anderson acceleration, or another convergent algorithm can improve the bound reached within a time budget. It does not by itself increase the optimum of an unchanged relaxation. Likewise, a model described as a “new SDP” need not dominate the existing DNN model: several classic QAP models are equivalent [7], while some scalable sparse-graph formulations explicitly relax the larger model further [8].

The closest precedents for the proposed work are Jünger–Kaibel's QAP box inequalities [9], the hypermetric cuts and ADMM implementation in MADAM [10], and partial higher-order moment relaxations [11]. The application and selection strategy proposed here should be investigated as a synthesis of these ideas. I have not established novelty of this precise combination.

## 3. Exploit wil50's actual geometry

Direct inspection of `qaplib/wil50.dat` confirms

\[
D_{ab}=|r(a)-r(b)|+|c(a)-c(b)|,
\qquad r(a)=\lfloor a/10\rfloor,\quad c(a)=a\bmod10.
\]

**Indexing warning:** in this file, the first matrix is the grid distance matrix and the second is the flow matrix. In the repository's native convention, use \(X_{ai}\), with grid location \(a\) first and flow facility \(i\) second. Do not assume the first matrix is the flow matrix merely because a comment uses that interpretation.

Let \(\mathcal S\) contain the four row prefixes and nine column prefixes:

\[
S^r_t=\{a:r(a)<t\},\ t=1,\ldots,4;
\quad S^c_t=\{a:c(a)<t\},\ t=1,\ldots,9.
\]

Then, exactly,

\[
D_{ab}=\sum_{S\in\mathcal S}
\left|\mathbf 1_{a\in S}-\mathbf 1_{b\in S}\right|.
\]

Their capacities are 10, 20, 30, 40 and 5, 10, …, 45. Define a binary membership variable for each facility:

\[
z_i^S=\sum_{a\in S}X_{ai}\in\{0,1\},
\qquad \sum_i z_i^S=|S|.
\]

Writing the flow matrix as \(F\), the ordered-pair QAP objective becomes

\[
f(X)=\sum_{S\in\mathcal S}\sum_{i,j}F_{ij}
\bigl(z_i^S+z_j^S-2z_i^Sz_j^S\bigr).
\]

This identity explains why regional consistency might matter: it constrains variables that directly describe every term of the distance objective. It is not evidence of a particular numerical gain.

At a branch-and-bound node, keep the original grid coordinates of the remaining locations. For free facilities, \(z_i^S\) sums only over free locations in \(S\), and its capacity is the number of those locations. Fixed interactions remain in the reduced linear term and shift. Renumbering the reduced matrix as a new rectangular grid would be incorrect.

Retain the full DNN model. Summing bounds from 13 independently optimized balanced-cut problems can lose the coupling imposed by a single common assignment and need not beat the current bound.

## 4. Regional cuts are genuinely additional constraints

Let \(Y\) be the existing lifted assignment matrix, and define its regional moments

\[
p_i^S=\sum_{a\in S}Y_{0,(a,i)},\qquad
W_{ij}^S=\sum_{a,b\in S}Y_{(a,i),(b,j)}.
\]

On the DNN feasible set, diagonal assignment entries can replace the first-row entries in the definition of \(p\). For an approximate ADMM iterate, those quantities can differ; residuals must be reported when separating cuts.

The regional first-order moment matrix

\[
\begin{bmatrix}1&p^T\\p&W\end{bmatrix}
\]

is a linear compression of \(Y\), and is already PSD. **Adding only this PSD condition is redundant.** The new content comes from integrality conditions on the regional variables.

A related trap concerns Manhattan distance-matrix SDPs. If \(C=\alpha J-D\succeq0\), the constraint \(\alpha J-\widetilde{\mathbb E}[X^TDX]\succeq0\) is already implied by the full assignment moment PSD condition and lifted assignment equations: factor \(C=\sum_\ell v_\ell v_\ell^T\) and write the left side as \(\sum_\ell\widetilde{\mathbb E}[(X^Tv_\ell)(X^Tv_\ell)^T]\). Each summand is a PSD compression of the assignment moment matrix. This explains why simply appending that particular small geometric SDP constraint cannot supply the missing strength.

### Boolean triangles and hypermetrics

For distinct facilities \(i,j,k\), valid inequalities include

\[
W_{ij}+W_{ik}-W_{jk}\le p_i
\]

and its two other choices of apex, together with

\[
p_i+p_j+p_k-W_{ij}-W_{ik}-W_{jk}\le1.
\]

Each is linear in the existing \(Y\); no additional PSD block is needed. Each of the three observables may also use a different region. This permits cuts coupling different thresholds.

For stronger inequalities, let \(s_i=1-2z_i\) and

\[
Q_{ij}=1-2p_i-2p_j+4W_{ij}.
\]

For any integer vector \(b\) whose entries sum to an odd number,

\[
b^TQb\ge1
\]

because \(b^Ts\) is an odd integer. Supports of size 3, 5, and 7 yield the familiar triangle, pentagonal, and heptagonal families. MADAM provides an ADMM precedent for adding these dynamically [10]. This does not mean MADAM's implementation can be transplanted into the QAP code without adapting its linear operators.

### Integer counts and box inequalities

For facility set \(T\), region \(S\), and integer \(r\), the count

\[
q=\sum_{i\in T}z_i^S
\]

is an integer. Thus

\[
(q-r)(q-r-1)\ge0.
\]

Linearizing its expectation gives

\[
2\sum_{i<j\in T}W_{ij}^S
-2r\sum_{i\in T}p_i^S+r(r+1)\ge0.
\]

These inequalities prohibit an implausibly small variance around a fractional count. They connect directly to the QAP box-inequality literature [9]. Signed counts and mixed regions supply further valid quadratic inequalities. Separation is the difficult part, and geometric regions are a practical starting set rather than an exhaustive search.

There are \(13\cdot4\binom{50}{3}=1,019,200\) ordinary same-region triangles at the root, small enough for exhaustive separation after aggregation. Exhaustively testing pentagons and heptagons is much less attractive; use randomized or local-search separation and retain a limited active set.

## 5. Exact strictness example, including a genuine QAP gap

Use six locations split into two groups of three. Split six facilities into two groups of three as well. Let **both** QAP data matrices be the adjacency matrix of \(K_{3,3}\): entries are 1 across the groups and 0 within a group.

If \(a\) facilities from the first facility group occupy the first location group, the ordered-pair objective is

\[
f=2\bigl(a^2+(3-a)^2\bigr)
=9+4(a-3/2)^2.
\]

Since \(a\in\{0,1,2,3\}\), the true optimum is 10. A DNN solution can behave as though \(a=3/2\) with zero variance, attaining 9.

Here is an explicit rational witness. Define a 6×6 correlation matrix \(Q\) with diagonal 1, off-diagonal entries \(-1/2\) within either facility group, and cross-group entries 0. Set

\[
p_i=1/2,\qquad W=(J+Q)/4.
\]

Set \(Y_{00}=1\), \(Y_{0,(a,i)}=Y_{(a,i),(a,i)}=1/6\). Conflicting assignment pairs have zero entries. For nonconflicting pairs,

\[
Y_{(a,i),(b,j)}=
\begin{cases}
W_{ij}/6,&a,b\text{ in the same location group},\\
(1/2-W_{ij})/9,&a,b\text{ in different location groups}.
\end{cases}
\]

Exact rational verification gives the spectrum

\[
2\ (1\text{ time}),\quad 1/4\ (8),\quad
3/16\ (16),\quad0\ (12).
\]

It satisfies the lifted row/column equalities, diagonal links, gangster zeros, and nonnegativity. All atomic Boolean triangles hold; their maximum violation is \(-1/12\). However, for the first three facilities and first location group,

\[
p_1+p_2+p_3-W_{12}-W_{13}-W_{23}-1=1/8.
\]

Writing the slack of this regional inequality as

\[
g=1-z_1-z_2-z_3+z_1z_2+z_1z_3+z_2z_3,
\]

the QAP objective is \(10+8g\). The cut forces \(g\ge0\), so it closes the gap exactly. A Lagrangian multiplier of 8 produces the constant objective 10 on permutations.

The native C++ certificate reproduced this behavior:

| Run | Certified bound | Iterations |
|---|---:|---:|
| Ordinary DNN | 8.9999998013 | 175 |
| DNN plus this regional cut | 9.99999999998 | 150 |
| Enumeration of all 720 permutations | 10 exactly | — |

The supplied validation also checks all four regional triangle types on all 720 permutations: 2,880 checks. This example proves that the proposed cut family can strengthen the existing model on nonnegative QAP data. It does not estimate the gain on wil50.

## 6. The stronger SDP: second-order regional moments

For each selected region \(S\), introduce moments \(u_I^S\) for facility subsets \(I\) of size at most four. The empty-set moment equals one. Index a matrix by

\[
\mathcal I_2=\{I\subseteq[m]:|I|\le2\},
\]

and define

\[
M_2^S(u)_{I,J}=u_{I\cup J}^S,
\qquad M_2^S(u)\succeq0.
\]

Using set union incorporates \(z_i^2=z_i\). Link the low-order moments to the original DNN matrix:

\[
u_\emptyset^S=1,\qquad u_{\{i\}}^S=p_i^S,
\qquad u_{\{i,j\}}^S=W_{ij}^S.
\]

For residual region capacity \(k_S\), impose the cardinality identities

\[
\sum_{i=1}^m u_{I\cup\{i\}}^S=k_Su_I^S,
\qquad |I|\le3.
\]

The proposed bound minimizes the existing QAP objective over **the original DNN constraints plus these linked regional constraints**. Every permutation extends to a feasible moment sequence by setting \(u_I=\prod_{i\in I}z_i\). Therefore it is a valid relaxation, and its optimum is at least the original DNN optimum. Adding regions or additional valid constraints cannot weaken the exact optimum.

This is genuinely higher order. Merely introducing free matrices called “moments” without the repeated-entry identities and links would not provide the stated strength.

The regional triangle slack \(g\) in Section 5 satisfies \(g^2=g\) modulo the Boolean identities. Since \(g\) has degree two, \(M_2\succeq0\) implies

\[
\widetilde{\mathbb E}[g]
=\widetilde{\mathbb E}[g^2]\ge0.
\]

Thus this SDP removes the strictness witness too. This argument applies to the other Boolean triangle types. I am not claiming that second-order moments automatically impose every pentagonal or heptagonal inequality; those can be added explicitly.

### Size and trace

At the root, \(m=50\):

| Object | Size before additional facial reduction |
|---|---:|
| Current reduced DNN PSD block | 2,402×2,402 |
| One regional \(M_2\) block | 1,276×1,276 |
| Distinct moments of degree at most four per region | 251,176 |
| One dense regional matrix in doubles | 13.03 MB |
| Thirteen dense regional matrices, one copy each | 169.33 MB |
| Full assignment order-two monomial basis, using compatible pairs | 3,003,751 entries |

The last row would require about 72 TB for a single dense square matrix before more reductions. It is not the formulation recommended here. Storage for a solver is larger than storage for one PSD matrix: account for primal/dual copies, eigenvectors, linear operators, and factorization workspaces.

The crude sum of dense cubic eigendecomposition costs for all 13 regional blocks is about 1.95 times one full eigendecomposition of the current block. **This is not a runtime prediction.** The existing subspace tracker, extra consensus equations, conditioning, and convergence rate can change the comparison substantially. Start with one or two selected regions, not all 13 by default.

There is an especially useful exact identity:

\[
\operatorname{tr}M_2^S=1+k_S+\binom{k_S}{2}.
\]

It follows from the Boolean and cardinality equations. With an orthonormal facial reduction, the reduced PSD block has the same trace. If \(k_S>m/2\), representing the complementary region replaces \(z\) by \(1-z\); the complete degree-two moment formulation is equivalent, but its monomial basis has a smaller trace and may be better conditioned.

For additional control over cost, introduce moment blocks only on selected groups of regional observables. A degree-two block on seven binary variables has dimension 29. Alternatively, enforce the exact local binary marginal polytope for those seven variables using 128 nonnegative configuration weights. Link their first and second moments to \(Y\). Such blocks can include different regions and therefore constrain their coupling. This is a selective relaxation, not a claim of equivalence to the full higher-order hierarchy.

## 7. Solving and certifying the strengthened bounds

### Fastest experiment: keep the existing certificate

Write any valid cut as \(\langle G_t,Y\rangle\le h_t\). For nonnegative multipliers \(\mu_t\), define

\[
L_\mu=L+\sum_t\mu_tG_t,\qquad
s_\mu=s-\sum_t\mu_th_t.
\]

For every feasible permutation,

\[
s_\mu+\langle L_\mu,Y\rangle\le s+\langle L,Y\rangle.
\]

Consequently, any valid lower certificate for the modified objective is a valid bound for the original node. Reuse `L_add`, the shift adjustment, objective-scaling bookkeeping, and the current spectral/GL certificate. This requires no new cone in the first prototype.

It remains valid if multiplier optimization stops early or produces poor multipliers. Those failures affect strength, not validity. Retain the maximum certified value across trials, including the unmodified bound. Count all coefficient contributions conservatively in the rounding allowance, especially if constructing the modified objective involves cancellation.

A stabilized bundle method is a better eventual outer optimizer than an aggressive Polyak step aimed at an unattainable cutoff. Approximate, infeasible ADMM primal iterates are useful for proposing cuts but are not automatically exact supergradients of the dual function. Separate that heuristic from the lower certificate, and tighten inner solves before drawing conclusions about a plateau.

### Integrated regional SDP

There is an intermediate decomposition that is useful before building a joint consensus solver. Let \(\mathcal D\) be the original DNN feasible set, \(\mathcal M_S\) the regional order-two moment feasible set with its own normalization and cardinality constraints, and write their linking equations as \(A_SY=B_Su^S\). Dualizing just those links gives

\[
q(\lambda)=s+
\inf_{Y\in\mathcal D}\left\langle L+\sum_S A_S^*\lambda_S,Y\right\rangle
+\sum_S\inf_{u^S\in\mathcal M_S}\langle-B_S^*\lambda_S,u^S\rangle.
\]

Every \(\lambda\) gives a valid lower bound. The first subproblem uses the existing generalized-objective DNN solver; the regional subproblems are independent smaller SDPs. Certified lower estimates of **all** subproblem infima can be summed safely. This separates the difficult software changes and permits parallel regional solves, at the cost of an outer dual optimization loop. The regional subsolver and its certificate remain new work; the current `L_add` cut prototype does not implement them.

The larger model can use a two-block consensus splitting: one block contains the product of PSD cones and simple bound cones, and the other contains the linear moment identities, links, and objective. PSD projections remain independent eigendecompositions. The linear projection is no longer the current entrywise clipping step; it needs a structured solve, potentially iterative with a reusable preconditioner.

Do not append several Gauss–Seidel updates to ADMM and assume convergence. A sound two-block formulation or an established convergent multi-block method is needed. Partial higher-order relaxations are well precedented [11], but this QAP implementation still needs numerical development.

The trace identities provide a certification route. If an approximate dual slack \(S_b\) for block \(b\) satisfies \(S_b\succeq-\epsilon_b I\), and every feasible block has trace \(\tau_b\), then

\[
\langle S_b,R_b\rangle\ge-\epsilon_b\tau_b.
\]

Subtract the sum of these corrections from the candidate dual objective, alongside correctly handled polyhedral terms and arithmetic error. This is a route to a certificate, not a complete implementation: all equality multipliers, signs, scaling, and trace conditions must be included and checked. In particular, a primal SDP objective is never a fathoming certificate by itself.

## 8. Selective triple consistency: a practical RLT2 supplement

For three free facilities \(i,j,k\), let their pairwise placement marginals be

\[
P^{ij}_{ab}=Y_{(a,i),(b,j)},\quad
P^{ik}_{ac},\quad P^{jk}_{bc}.
\]

Require a nonnegative tensor \(t_{abc}\), supported only on distinct locations, satisfying

\[
\sum_c t_{abc}=P^{ij}_{ab},\quad
\sum_b t_{abc}=P^{ik}_{ac},\quad
\sum_a t_{abc}=P^{jk}_{bc}.
\]

Every permutation satisfies these equations. They test whether the three pairwise marginals can coexist, a condition not guaranteed merely by their separate assignment constraints. Full RLT2 has further shared consistency equations across all selected indices; this proposed selective construction should not be mislabeled full RLT2.

One facility triple at order 50 has 117,600 possible injective placements, versus approximately 2.305 billion degree-three assignment monomials for all triples. Adding a limited number of triples avoids the latter allocation.

An alternative is to separate projected cuts. If potentials satisfy

\[
\alpha_{ab}+\beta_{ac}+\gamma_{bc}\ge0
\quad\text{for every distinct }a,b,c,
\]

then

\[
\sum_{ab}\alpha_{ab}P^{ij}_{ab}
+\sum_{ac}\beta_{ac}P^{ik}_{ac}
+\sum_{bc}\gamma_{bc}P^{jk}_{bc}\ge0
\]

is valid. A bounded LP can find violated potentials, and the resulting inequality goes through the same `L_add` machinery. In the supplied separator, potentials are rounded to multiples of 1/1024 and repaired by checking every injective placement in integer arithmetic. Thus the validity of the rounded cut does not rest on the LP solver's feasibility tolerance.

## 9. Computational investigation

### Experimental scope

I built a separate driver against the pinned branch, using its existing `dnn_bound` and generalized objective hook. The principal test was the first difficult depth-three sample, `34:40,28:22,14:29`, with 47 free rows and columns. The original run used 1,500 PRSM iterations. Cut trials and controls restarted from that run's saved multiplier and used 700 iterations. The strongest earlier certificate was retained even when a continuation returned a weaker one.

| Run | Iterations | Returned lower certificate | Difference from strongest uncut certificate |
|---|---:|---:|---:|
| Initial uncut PRSM | 1,500 | 48,409.2349 | +0.0000 |
| Uncut ADMM continuation | 700 | 48,399.0700 | -10.1649 |
| 100 cuts, multiplier 2; ADMM | 700 | 48,412.5387 | +3.3038 |
| 100 cuts, multiplier 10; ADMM | 700 | 48,403.7761 | -5.4588 |
| Uncut PRSM continuation | 700 | 48,408.2001 | -1.0348 |
| 1,000 cuts, multiplier 1; PRSM | 700 | 48,449.0658 | +39.8309 |
| 1,000 cuts, multiplier 3; PRSM | 700 | 48,238.1095 | -171.1253 |


The strongest uncut certificate was 48,409.2348707950. The best cut trial improved it by **39.8309330592**, and improved the matched 700-iteration PRSM control by **40.8657009016**. Relative to incumbent 48,816, it closes approximately **9.8% of the previous conditional node gap**. The repository's cutoff is 48,814.05, so the node remains unfathomed by roughly 365 units.

These are **computed lower certificates**, not exact optima of the relaxations. The 1,000 cuts were the most violated same-region triangles in the original iterate, each with multiplier 1. They were not dynamically updated or optimized by a bundle method. Raising every multiplier together to 3 was unhelpful. This illustrates why a large violation count or a larger penalty does not guarantee a stronger finite dual bound; it does not establish a ceiling on the fully optimized cut relaxation.

The controls also matter: the warm-start interface restores Z but initializes Y afresh. Consequently a continuation is not an uninterrupted run, and some 700-iteration certificates were weaker than the original 1,500-iteration one. Comparing only against those weaker controls would exaggerate the gain. The 39.8309 figure uses the strongest uncut certificate obtained here.

### Separation diagnostics

| Separation set on the original iterate | Regions/observables | Violations above 1e-5 | Maximum violation |
|---|---:|---:|---:|
| Same-region threshold triangles | 13 regions | 11,206 | 0.0910535 |
| Same-region rectangle triangles, after duplicate/complement removal | 770 regions | 21,705 | 0.0910535 |
| Mixed-region triangles | 611 observables; 151,320,260 inequalities | 389,765 | 0.0910535 |

The original PRSM iterate had appreciable feasibility error. An uncut 700-iteration ADMM refinement reduced the maximum lifted row/column moment discrepancy to approximately **3.89e-5**, with maximum diagonal-link discrepancy **3.27e-5**. The largest regional triangle violation remained **0.0905820**. This is evidence that the diagnostic is not simply disappearing as the primal residual shrinks, although the refined iterate is still not an exact feasible DNN point.

Ten selected facility-triple separator LPs all returned valid violated cuts. Each rounded cut was checked on all **97,290** injective placements of its three facilities among 47 locations. Their initial violations were approximately 0.250–0.281, under a different normalization from the triangles. The first triple, local indices `4,25,46` (original flow indices `4,26,49`), still had a violation of **0.2724130** after the ADMM refinement. Those ten separation LPs took roughly 1.5–7.7 seconds each in this environment. Their cuts were **not** used in the reported strengthened-bound trials.

The exact witness checks passed, including 2,880 permutation checks and 2,000 additional checks of mixed-region coefficient construction, overlapping supports, and repeated facilities. The supplied makefile fragment was built successfully, and the delivered driver reproduced the toy cut certificate.

### What the pilot does and does not establish

It establishes a meaningful certified improvement on one difficult wil50 subproblem, using a valid cut family that differs from the repository's atomic-triangle experiment. The earlier atomic experiments were on nug30 and other instances, so their absolute gains should not be compared directly to wil50's different objective scale; the exact witness is the clean demonstration that the families differ.

It does not establish an improved wil50 root record, an optimality proof, a uniform gain across the 16 sampled nodes, a wall-time advantage for a complete search, or the performance of the full regional order-two SDP. The larger SDP and joint triple-consistency bound remain proposed implementations. The scripts are research prototypes, not production fathoming code.

The `cpu` values in the logs follow the repository's timer convention; runs used different BLAS thread counts and some ran concurrently. They should not be presented as controlled wall-time comparisons. For a production decision, rerun promising settings under a fixed resource budget and record elapsed time explicitly.

## 10. Recommended development sequence

1. **Use the existing 16 depth-three wil50 nodes as the evaluation set.** Include the clear misses, not just nodes already close to fathoming. Maintain a well-converged uncut reference and an equal-budget uncut control.
2. **Try regional and mixed-region cuts through `L_add`.** Begin with threshold sets, then rectangles or heuristically selected sets. Record certified gain, marginal residuals, active cuts, and time. Cut violation alone is not an objective improvement, and its magnitude depends on normalization.
3. **Add a small number of projected triple-consistency cuts.** Select facilities from violated regional cuts or high-cost conflicting marginals. Check rounded cut validity independently.
4. **If the required gain remains missing, prototype one regional \(M_2\) block.** A middle column prefix is a natural first test, followed by regions selected using actual violations. Implement cardinality facial reduction and trace-aware certification before using it to fathom.
5. **Increase strength selectively.** Compare extra regional blocks, mixed-region small blocks, and larger integer-count inequalities. Avoid committing to all regions at every node until the cost/benefit is measured.
6. **Measure the tree effect.** Compare additional fathoms, residual child counts, and total work. A 130-unit improvement at one node is not sufficient evidence for a factor-of-47 reduction across the search; heterogeneous gains and expensive bounding can change the result.

If a cut family adds only fractions of a bound unit after converged controls, retire it for this target. If it consistently supplies 100–300 units on difficult nodes at a modest multiple of the DNN cost, it deserves integration and a broader tree experiment. These are proposed research gates, not observed results for the full formulation.

There are useful fallback options. A disjunctive bound obtained by taking the minimum of all child bounds for one complete branching partition is valid; combining partitions can strengthen the relaxation further. Its cost can approach performing the branching work itself. Symmetry reduction can still accelerate a root solve, but a typical fixed node destroys much of the grid symmetry, and speed improvements must be distinguished from relaxation strength.

## 11. References and reproducibility

1. **QAPLIB**, [Montréal instance table](https://qaplib.mgi.polymtl.ca/). The [Lehigh mirror](https://coral.ise.lehigh.edu/data-sets/qaplib/qaplib-problem-instances-and-solutions/) retains an older wil50 lower-bound entry, so the mirrors should not be treated as interchangeable evidence of a current record.
2. **Graham, Hu, Im, Li, Wolkowicz**, *A Restricted Dual Peaceman–Rachford Splitting Method for a Strengthened DNN Relaxation for QAP*. [Author-hosted 2021 revision](https://www.math.uwaterloo.ca/~hwolkowi/henry/reports/RevisionApr2021newPRSMforQAP.pdf), especially Table 4.3, wil50 row.
3. **Fujii, Ito, Kim, Kojima, Shinano, Toh** (2021), *Solving Challenging Large Scale QAPs*. [arXiv:2101.09629](https://arxiv.org/abs/2101.09629), Table 5.
4. **Mittelmann and Peng** (2010), *Estimating Bounds for Quadratic Assignment Problems Associated with Hamming and Manhattan Distance Matrices Based on Semidefinite Programming*. SIAM Journal on Optimization 20, 3408–3426. [Author manuscript](https://optimization-online.org/wp-content/uploads/2008/05/1980.pdf), Table 3. Also relevant: [Gueye–Michelon's distance-variable formulation](https://www.sciencedirect.com/science/article/abs/pii/S219244062100040X); its structural motivation should not be confused with dominance over modern DNN.
5. **Date and Nagi**, *RLT2-based Parallel Algorithms for Solving Large Quadratic Assignment Problems on Graphics Processing Unit Clusters*. [arXiv:1710.03732](https://arxiv.org/abs/1710.03732), Table 5; journal version INFORMS Journal on Computing 31(4), 2019.
6. **Hahn, Roth, Saltzman, Guignard** (2013), *Memory-Aware Parallelized RLT3 for Solving Quadratic Assignment Problems*. [Manuscript](https://optimization-online.org/wp-content/uploads/2013/12/4144.pdf), Section 4. The underlying RLT3 bound is described by Hahn, Zhu, Guignard, Hightower, Saltzman, [INFORMS Journal on Computing 24, 202–209](https://doi.org/10.1287/ijoc.1110.0450).
7. **Povh and Rendl** (2009), *Copositive and Semidefinite Relaxations of the Quadratic Assignment Problem*. Discrete Optimization 6, 231–241. [Author summary and manuscript](https://optimization-online.org/2006/10/1502/).
8. **Bravo Ferreira, Khoo, Singer**, *Semidefinite Programming Approach for the Quadratic Assignment Problem with a Sparse Graph*. [arXiv:1703.09339](https://arxiv.org/abs/1703.09339). Relevant as a scalability comparison, not a proven strengthening of the full DNN model.
9. **Jünger and Kaibel** (2001), *Box-Inequalities for Quadratic Assignment Polytopes*. Mathematical Programming 91, 175–197. [University repository](https://kups.ub.uni-koeln.de/54792/), [1997 submitted manuscript](https://kups.ub.uni-koeln.de/54792/1/zpr97-285.pdf).
10. **Hrga and Povh** (2021), *MADAM: a Parallel Exact Solver for Max-Cut Based on Semidefinite Programming and ADMM*. Computational Optimization and Applications 80, 347–375. [Open-access article](https://doi.org/10.1007/s10589-021-00310-6), [implementation](https://github.com/HrgaT/MADAM).
11. **Campos, Misener, Parpas** (2023; online 2022), *Partial Lasserre Relaxation for Sparse Max-Cut*. Optimization and Engineering 24, 1983–2004. [Open-access article](https://doi.org/10.1007/s11081-022-09763-y), [implementation](https://github.com/cog-imperial/Partial-Lasserre-relaxation-sparse-maxcut). This supports the selective higher-order strategy, not a demonstrated wil50 result.
12. **Kim and Kojima** (2021), *Generating Cutting Inequalities Successively for Quadratic Optimization Problems in Binary Variables*. [arXiv:2107.08665](https://arxiv.org/abs/2107.08665). An additional route uses objective-cutoff-dependent valid inequalities; those must be scoped to the incumbent and node assumptions used to certify them.

Repository evidence is pinned to the commit above, particularly `experiments/dnn_cuts/README.md`, `experiments/wil50_feasibility/NOTES.md`, the `depth3_converged` logs, and `latest/dnn.h`. The research work was performed in a detached local worktree; no changes were pushed to `cognac`.

The experiment bundle contains the driver, separators, exact witness verification, toy QAP, validation output, and scalar logs. Large binary ADMM checkpoints are omitted; commands to regenerate them are included. The instance SHA-256 is `2ee770305d830283a74f7950543584035e4cfa4967338666b8c0d0a75d550a1b`.
