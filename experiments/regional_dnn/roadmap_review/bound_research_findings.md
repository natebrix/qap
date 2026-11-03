# Bound audit and additional opportunities, 2026-11-02 (America/Los_Angeles)

Repository audited: `cognac checkout at ca50467`, supplied revision `ca50467`. No repository files changed and no large QAP solves run. Tests are in `bound_audit_checks.py`.

## Findings that change the interpretation of the existing evidence

1. **The triangle-pool optimum has not been upper-bounded.** `dualtest/README.md` establishes achieved certified values, not a certified distance to the pool optimum. Thirty rounds gaining only 5–8 establishes slow progress under those settings. It does not establish “within 10,” “not 50 above,” or exhaustion. `local7/README.md` explicitly shows that enlarging a pool lowered the *achieved* bound by 15–25; this is optimization failure, since zero multipliers on the new cuts recover the old feasible dual solution. A second slow30 failure would not close this logical gap.

   Write `v_pool = min { f(Y): Y in D, g(Y)<=0 }` and `q(mu)=min_{Y in D}[f(Y)+mu'g(Y)]`. The code certifies a lower bound on `q(mu)`, hence on the QAP optimum. Its printed `primal` is an approximate, generally infeasible Y value for one modified objective. Even an **exact D-feasible** minimizer of one frozen Lagrangian only supplies an upper bound on that particular `q(mu)`, not on `sup_mu q(mu)` or `v_pool`. The required upper witness for the pool is a D-feasible, all-pool-cut-feasible Y evaluated using **f**, without the multiplier penalties.

   A practical cap diagnostic is convex primal recovery from several near-stationary iterates, followed by verified affine/PSD/nonnegative/pool feasibility. A second route is a master LP of upper planes from *verified D-feasible* witnesses: maximize `t` subject to `t <= f(Yk)+mu'g(Yk)`, `mu>=0`. If finite, its objective upper-bounds `sup_mu q(mu)`; its LP dual constructs a convex combination satisfying pool inequalities. Approximate infeasible iterates do not furnish these upper planes without error corrections. A loose or unbounded master is inconclusive.

2. **Literature comparisons do not prove dominance.** The original `research_memo_astra.md` and `NOTES.md` correctly distinguish computed RLT2 dual-ascent values from the RLT2 optimum and distinguish selected independent triple marginal tests from full RLT2's shared consistency. `roadmap/bound_strength.md` compresses this into an unsupported exhaustion claim. A weaker achieved bound in another relaxation is not a theorem that its intersection with DNN adds nothing. No recommendation to build full RLT2/RLT3 follows: their costs still make them poor first experiments. Likewise vertex fixing is not a reason to reject all grouped disjunctions or cut-generating lift-and-project methods.

3. **Additive inherited gain is unproved.** A child inherits at least the parent's *certified scalar bound*, and the exact conditioned fixed-multiplier Lagrangian has monotonicity. Neither says the child's uncut one-level improvement adds to the parent's cut gain. Nor is the child's own achieved cut bound an upper bound on the value of inherited multipliers or the child pool optimum. Thus the roadmap's 1.5–4x inheritance estimate needs an experiment; it is not justified by its purported bracket.

## Cut inheritance: what is true, and the constants

For one fixed assignment r, form the embedding T from child lifted coordinates to parent lifted coordinates. Parent coordinate 0 and r both equal child coordinate 0; deleted incompatible assignments are zero; the other coordinates retain their original IDs. For any symmetric parent matrix M:

`M_child = T' M T`,

`M_child[0,0] = M[0,0] + 2 M[0,r] + M[r,r]`,

`M_child[0,a] = M[0,p(a)] + M[r,p(a)]`.

Then `<M,T Yc T'>=<M_child,Yc>` exactly, including for fractional child matrices. Every child DNN-feasible Yc embeds into the parent's DNN feasible set. With `M=L+sum mu G` and constant `c=shift-sum mu h`, `q_child(mu)>=q_parent(mu)` follows. This applies without assuming any empirical convergence. The code's `dnn_warm_child` already implements this transformation for Z, including the constant entry. It stores only Z; this is an initialization, not automatically a transferred certificate.

Implementation choices:

- First transfer the **whole weighted objective** by T for a correctness reference; include its scalar shift once. Alternatively reconstruct tuples using original region and facility IDs and symbolically substitute fixed memberships.
- If reducing the original QAP separately, add only the transferred cut objective and cut constants to the already reduced QAP; do not add the original objective's fixed-assignment shift twice.
- `dnncut.C` currently carries `shift - hsum` separately from `Ladd`; higher-family RHS terms may be nonzero and must travel with the cut.
- Parent and child objective scales can differ. Preserve the scale metadata and the code's existing rescaling of Z.
- Do not assume a last quiet iterate corresponds to the exported best multipliers: export the exact pool, mu, Z, and objective snapshot together. The current internal `best_mu` / `best_save` pairing is the right starting point.

### Triangle cuts touching the fixed facility really can be dropped

For the other two facilities, define the four regional pair masses:

`q11=Wjk`, `q10=pj-Wjk`, `q01=pk-Wjk`, `q00=1-pj-pk+Wjk`.

They are nonnegative in the child DNN: they are sums of entries of Y over S×S, S×Sc, Sc×S, and Sc×Sc, using the assignment moment equalities. Condition each of the four Boolean triangles on `zf=0` or `zf=1`; its entire left-minus-right expression becomes one of `-q00,-q01,-q10,-q11`. It is therefore DNN-redundant.

Dropping its **entire** weighted expression, including its RHS contribution, removes a nonpositive function from the minimization objective and thus can only increase the child's fixed-multiplier relaxation value. Consequently dropping is safe and monotonicity survives, although exact objective equality with the conditioned parent is lost. Keeping the matrix term while forgetting the RHS is not this operation and may invalidate a bound. Keeping the RHS while deleting the matrix generally weakens it. Preserve both or remove both.

The roadmap's branching suggestion “avoid facilities with much touching triangle multiplier mass” does not follow: those touching triangles disappear as redundant nonpositive penalties; losing that mass is not inherently losing bound strength.

### This rule is false for higher-order cuts: exact small counterexample

The supplied check constructs an n=8, region-capacity-3, full 65×65 **rational QAP-DNN matrix**. Facilities are partitioned into two groups of four. Its regional moments are:

- `p_i = 3/8`;
- `W_ij = 1/16` within a group, `9/64` across groups, and `W_ii=p_i`.

Uniformly spread assignment moments within the three region locations and five complement locations: distinct-location, distinct-facility entries are respectively `W_ij/[3*2]`, `(p_i-W_ij)/[3*5]`, or `(1-p_i-p_j+W_ij)/[5*4]`. Diagonals/arrow entries are `p_i/3` or `(1-p_i)/5`; gangster entries are zero.

The script verifies PSD by **exact Fraction Gaussian/LDL elimination**, including every zero pivot's remaining row, and verifies all assignment-moment equalities in exact arithmetic. All regional triangles hold. Yet for either group of four:

`sum_i p_i - sum_{i<j} W_ij - 1 = 1/8 > 0`.

This is precisely the residual of a valid five-point count cut after fixing its fifth membership to zero. Hence a general implementation must substitute and retain residual 5-point/local7 cuts, not discard every cut touching a fixed facility. The example proves nonredundancy against full DNN; it does not claim the existing wil50 pool contains this particular cut.

Raw test output (0.63 s wall time):

```
Exact full QAP-DNN witness: PSD positive/zero pivots (49, 16) ; regional triangle max -1/16 ; conditioned five-point violation 1/8
Exact integer parent-child objective identities: 600
```

## Ranked next tests

### 1. Fixed-parent-certificate child screening (tree agent owns implementation)

The tree agent's genuinely additional idea passes the proof audit. Let `M=L+Zc` be the *actual certified* parent residual objective, including cut contributions if present, and suppose the existing certificate bounds `<Zc,Yp> <= tau_parent*delta`, with `tau_parent=m_parent+1`. For a child:

`b_child = c_parent + GL_child(T'MT) - tau_parent*delta - rounding_margin`

is valid. It requires no new eigensolve. Retain the parent's trace factor, not the child's, because the spectral statement is about embedded Yp. Retain the parent shift/cut RHS once; or perform an algebraically equivalent split into reduced-QAP and transformed-dual pieces. Each inner LAP, as well as the outer LAP, is now conditioned, unlike the existing `child_lb` outer reduced-cost shortcut.

Cheapest test: export a certified snapshot at failed triangle phases on wil50 d2 and sko49 child13 (only 3.5 short in the historical run), score all children of one candidate facility, and compare with existing child_lb and later actual child certificates. Measure *screened children per CPU second*, not bound units alone. The script checks all 600 parent-child objective identities for n=5 across every fixing and child permutation. It does not yet exercise the production GL or rounding code.

This overlaps inherited-pool tree savings; do not multiply independently estimated subtree reductions. It can cheaply decide which inherited children actually need DNN solves.

### 2. Inherited triangle pool, with honest gate and separated costs

Reuse the parent's exported best cut objective/pool and Z on six children of node10, then on the full node10 subtree if the gate passes. Include d1,d2,d4 plus three fixed random siblings, as already proposed. Record (a) equal-iteration certificates, (b) time to fathom, (c) avoided new cut phases, and (d) bound relative to the child's own independently optimized pool. Do not sell inherited gain as parent gain plus the uncut level gain.

Useful engineering precedents exist for strengthened SDP branch-and-bound and bundle methods, but the present mathematical proof is self-contained. This is primarily an amortization opportunity; no measured post-cut speedup is supplied here. Same-node persistence of Y and spectral work is a separate, speed-agent-owned optimization and must not be conflated with the old negative test of *child-conditioned* Y initialization.

### 3. Cross-region categorical local consistency: a genuinely untested cheap strength gate

Select two nested thresholds S⊂T and partition locations into `S`, `T\S`, and `Sc∩Tc = complement(T)`. For r facilities use one common distribution on their three-valued labels, respecting the **remaining** capacities of all three cells. All required unary and pairwise label moments are already sums of entries of full Y.

For r=5, enumerate at most `3^5=243` configurations; r=6 uses729. Eliminating the final label gives `2r+4*C(r,2)` nonconstant coordinates: 50 and72. A row/column split instead gives four cells and `4^5=1024` configurations for five facilities. The existing away-step/LP local7 separator generalizes directly: for any direction a, compute `beta=max_config a'v` by finite enumeration, yielding an independently valid linear cut.

Why this is not the retired experiment: binary local7 enforces separate distributions for each region; this asks for **one coupled distribution across regions**. The old full-location triple test has only three facilities; this has five or six, with coarser labels. Neither old negative test proves this constraint family redundant. Pairwise cross-region identities, nested same-facility products, and the small first-moment PSD compression are already supplied by full Y; those alone are not new cuts. A nonredundancy witness against every previously tested family is not established here; the diagnostic is meant to determine whether actual capped iterates expose a useful gap.

Cheapest decisive test: dump the quiet triangle-capped iterate on wil50 n13/d2 and sko49 child13/9. Separate 100–300 subsets per partition, concentrating on facilities with high conflicting pair moments and high objective weights, and compare normalized violation with a binary/mixed-triangle baseline on those **same moments**. Use mid partitions first. Then freeze the old triangle multipliers and line-search a single aggregate ray of the top20–40 new cuts, retaining the alpha=0 certified baseline. Advance only for reproducible +10–20 extra certified units at plausible phase cost; sko49 child13's 3.5-unit and child9's20.3-unit historical deficits make particularly informative endpoints. A negative gate rejects this tested configuration, not all linear inequalities.

No bound-unit or runtime gain is measured. This is a lower-cost alternative gate before implementing the proposed ~1,129×1,129 global regional M2 blocks, not an assertion it dominates them.

### 4. Repair heterogeneous dual optimization without risking the old bound

Perform the missing fixed-best-triangle-multiplier tests first. Preserve `(mu_tri,Z,best_certificate)` and add new-family multipliers along one normalized ray; include alpha=0, test a short increasing ladder, then optimize at most a few accepted coordinates/families. A joint optimization may always retain the old certified bound as its incumbent, even if its current iterates are worse. Scale by induced lifted-objective norm or measured perturbation/response, rather than counting cuts or using the same raw multiplier step across families.

This is a safer *diagnostic* than interpreting another slow30 plateau as a ceiling proof. If fixed-ray additivity is positive, investigate proximal/serious-step dual updates; if it is negative, the tested ray is unhelpful. Small frozen-Lagrangian primal-certificate gaps still do not certify the outer dual gap. No automatic +40 or full-level gain follows.

### 5. Grouped regional probes, behind the cheaper child screen

For a near-fathomed node choose `z_f^S=0 or1`. The current `DNNParams.elim` mask can express either domain by forbidding assignments of f inside or outside S, avoiding a new general conic constraint. Two certified solves give a parent lower bound equal to their minimum; one child above cutoff permanently eliminates that half-domain. For counts of a facility subset, `q<=a or q>=a+1` is a richer integer split but requires new constraint support.

This is valid and differs from single-assignment probing, but full-sized solves could cost more than the avoided subtree. Start only after fixed-parent GL screening, on nodes3–25 short, and prefer masks with high fractional uncertainty and asymmetric predicted certificate gains. Kill if neither half rises enough after a short equal-budget warm probe. It is a reserve idea, not a predicted week-saving factor.

## Budget and attribution

The planning baseline remains the user's estimated 2.2e9–7e9 post-cut core-seconds; it is not a measured full-instance total. At 95% useful utilization, a 1,000-core week supplies 574,560,000 useful core-seconds; the required improvement is **3.83–12.18x**. This audit provides no measured new factor toward that requirement. Inheritance, child screening, better certificates, and new cuts all target overlapping descendants. Combine them only through a paired subtree run's total CPU and generated/solved nodes, then validate on multiple parent classes on both instances. Do not exponentiate a median bound gain by `47^(gain/130)` as a forecast: threshold crossings, node weights, and the strong observed tail determine the actual work.

## Primary-source checks

- BiqBin, Gusmeroli et al., describes strengthened SDP relaxations with hypermetric inequalities, bundle bounding, and parallel branch-and-bound. It supports those methods as precedents, **not** a QAP speedup prediction or an unverified claim about its exact cut-inheritance implementation: https://arxiv.org/abs/2009.06240
- Date and Nagi explicitly describe a GPU **Lagrangian dual-ascent** algorithm for RLT2; these are achieved lower bounds, not certified upper witnesses for the RLT2 optimum: https://arxiv.org/abs/1710.03732
- de Klerk and Sotirov's vertex-fixing/symmetry work provides stronger SDP bounds for instances with suitable transitive symmetry. It does not imply all disjunctive relaxations equal ordinary branching on arbitrary depth3–5 nodes: https://optimization-online.org/2009/09/2399/
- Graham et al.'s strengthened restricted dual PRSM paper is directly relevant to distinctions between redundant subproblem constraints, multipliers and convergence: https://doi.org/10.1287/ijoc.2022.1161 . No new quantitative prediction here is attributed to that paper.
