# CP-ADMM for the cut phase: the frozen pool inside the Y-step (Nov 7, 2026)

The question is item 2 of `../alternatives.md` section 4: can the cut phase
reach its plateau in far fewer iterations if the cut pool goes inside the
ADMM? Today the phase re-solves the DNN with Lagrangian cut multipliers and
moves them by Polyak steps between rounds. CP-ADMM (de Meijer and Sotirov)
puts the cuts into the PRSM's Y-step as constraints and projects onto them by
Dykstra sweeps. This study certifies every check with the sweeps'
multipliers. The memo's kill line: not within 10 of 48,763 by 2,000
iterations, or more than 30% overhead per iteration over the plain PRSM
iteration.

Labels: [m] is measured here, in the log named; [e] is an estimate. All runs
are on this box (4-core Xeon with AMX, OpenBLAS), one core each, four at a
time. Iterations t are counted from the give-up point, for both protocols.

## Short answer

**It passes, by a wide margin.** On node 13 with its own round-12 pool
(13,000 regional triangles, frozen), from the give-up point:

- With one Dykstra sweep per iteration, CP-ADMM certifies 48,753 at
  t = 1,040 [m]. The Polyak protocol gets there at t = 4,860 [m], in its
  quiet round. That is 4.7x fewer iterations and 4.9x less CPU (547 s
  against 2,659 s) [m]. Three sweeps per iteration reach it at t = 850
  (5.7x, 5.8x in CPU) [m].
- CP-ADMM passes the Polyak protocol's final value (48,755.4 at t = 5,100)
  at t = 1,120 [m]. It has 48,767.1 by t = 2,000 and 48,770.5 by t = 3,000,
  still rising by about 2 per 500 iterations [m]. The three-sweep run
  certifies 48,771.7 at t = 3,000 [m]. So this pool's dual optimum is at
  least 48,771.7, not about 48,763.
- The overhead per iteration is +3.8% in the dual test's settings [m]. The
  CP work itself is +6.7%, and the eigensolver is cheaper on the CP
  iterates. In the production settings the overhead is +13.5%, because the
  plain run skips most Rump tests (DNN_CORE_LAZYCERT) and the CP run
  certifies every check [m]. The projection costs 18.4 ms per iteration:
  two dense P-steps of 8.7 ms and one pass over the 13,000 cuts of
  0.56 ms [m].
- The production settings give the same picture. CP-ADMM reaches 48,753 at
  t = 1,170 [m]. The Polyak run under the same settings reaches it at
  t = 4,290, in its quiet round, and ends at 48,762.1 [m].
- Every value is certified through the existing L_add certificate path.
  - Validity: 77 brute-forced small nodes, 12,656 certified values, none
    above the optimum [m].
  - The cut extraction matches `assemble` to 1e-15 relative [m].
  - The standard path, given CP-ADMM's multipliers, certifies 48,757.0,
    above anything the Polyak protocol certified [m].
  - With the switch off the solver is trace-identical to the parent
    commit [m].
- Two more nodes agree. Node 14 fathoms at t = 770 against the Polyak
  protocol's 4,390 (5.7x). Node d2 reaches the Polyak final value at
  t = 1,050 against 5,090 (4.8x) [m].
- The box-cone paper's cheap projection does not apply: none of our
  families has the bound-factor form. Its Appendix-A grouping of triangles
  would fit the regional triangles, but the cones would be nearly scalar
  [m].

**Two more nodes**, with the same settings and protocol and each node's own
round-12 pool (`run_node.sh`).

*Node 14* (`47:23,10:44,40:34`). Here the Polyak protocol fathoms, in its
quiet round at t = 4,390 (2,412 s). CP-ADMM s1p1 fathoms at t = 770 (430 s):
5.7x fewer iterations and 5.6x less CPU [m] (`polyak_new_14.log`,
`cp_s1p1_14.log`).

*Node d2* (`39:27,30:48,25:49,26:5`, the dual test's other node). Polyak:
48,790.8 after the quiet round, 23.2 short of the cutoff [m] (the laptop log
has 48,789.9). The dual test's slow ascent put d2's pool optimum near
48,795. CP-ADMM s1p1 reaches 48,785 (within 10 of that) at t = 880 against
4,670, 5.3x [m]. It passes Polyak's final value at t = 1,050 against 5,090,
4.8x [m]. It has 48,802.5 at t = 2,000, 11.5 short of the cutoff and still
rising [m]. With three sweeps it reaches 48,785 at t = 750 (6.2x) and
48,805 at t = 1,890. It ends at 48,807.9 after 4,500 iterations, nearly flat
(+0.7 between t = 3,000 and 4,000) and 6.1 short of the cutoff [m]
(`cp_s3p3_d2.log`). So d2 does not fathom on its frozen pool. But the pool's
dual optimum is at least 48,807.9, about 48,808-48,810 from the flat curve
[e], not about 48,795. d2 is 4-6 short at the cap, not 18. On d2 the
per-check certificate also moves more: from t = 1,500 on it lies up to 6
below the running best [m].

## Method

**The test.** Node 13 (`21:28,8:41,28:19`, m = 47, cutoff 48,814.05) with the
dual test's settings (`dnncuts` defaults: degree 4, refresh 20, gamma 1/1,
beta 0.3; `run_ref.sh`). Both protocols start at the give-up point: round 0
is 400 uncut iterations from cold, and the next solve starts warm with Z kept
and Y = 0. The Polyak protocol then runs 12 rounds of 300 iterations with
1,000 new triangles a round, steps toward best + 50 uncapped, and ends with
the quiet round of 1,500 (`stronger_bounds/dual/run_five.sh`). The frozen
pool is that run's own pool after round 12's new cuts: 13,000 triangles,
zero multipliers included (`-savepoolall 12:...`). CP-ADMM runs one solve
from the give-up point with that pool in the Y-step and nu = 0 at the start.

**The Y-step.** In scaled units (s = max |L|, L/s the objective), with
gamma1 = gamma = 1, the PRSM iteration is unchanged except for the Y-step:

    R = P_psd(V'(Y + Z/beta)V),  Zc = Z + beta (Y - VRV'),  Z <- Zc
    Y = argmin { ||Y - X|| : Y in P, <G_t, Y> <= h_t for all t },  X = VRV' - (L/s + Z)/beta
    Z <- Z + beta (Y - VRV')

P is today's Y-step set: 0 <= Y <= 1, Y_00 = 1, the gangster zeros and any
eliminated assignments. Without cuts the argmin is today's clip.

**Dykstra.** The sets are P and the half-spaces H_t, in cyclic order with P
first. Dykstra keeps one correction per set, and X = y + p_P + sum_t p_t at
all times. A half-space's correction is a multiple of its normal,
p_t = nu_t G_t with nu_t >= 0. Because P comes first in a sweep, its step is

    y = P(X - sum_t nu_t G_t)

whatever p_P was. So the warm start is the vector nu alone. The half-space
steps are Hildreth's:

    nu_t' = max(0, nu_t + (<G_t, y> - h_t) / ||G_t||^2),   y <- y - (nu_t' - nu_t) G_t.

One sweep is a P-step and K passes over the cuts. Each Y-step ends with a
P-step, so the returned Y is always in P. nu carries over between
iterations. When beta changes, nu is rescaled so that mu below stays put.

**Cheap sweeps.** A regional cut is linear in the moments of its region S:
p_i = sum_{a in S} Y_(ai),(ai) and W_ij = sum_{a,b in S} Y_(ai),(bj). So the
half-space steps run on the 13 regions' moments, not on Y. A step
y <- y - d G_t changes region S' by -d c_i k on p_i and by -d c_ij k(k-1)/2
on W_ij, with k = |S ∩ S'|. That is about 80 flops per cut. The P-step is one
pass over Y: it forms X - G'nu, clips it, and sums the (a, b) blocks of the
result for the moments. The blocks are grouped by the set of regions that
hold both a and b. The 5 x 10 grid's prefix regions give 48 such sets at
node 13 [m].

**Certificate.** For any symmetric Zc and any mu >= 0, every Y feasible for
the relaxation has <G_t, Y> <= h_t and Y = VRV' with trace R = m + 1, so

    <L, Y> >= <L + G'mu + s Zc, Y> - h'mu - s <Zc, VRV'>
           >= GL(L + G'mu + s Zc) - h'mu - (m+1) s lambda_max(V'Zc V).

This is the LP duality of `research_bound_quality.md` section 3.3, and the
inequality the cut rounds already certify, with L_add = G'mu. The KKT
conditions of the projection read
beta (Y - VRV') + L/s + Z + G'(beta nu) in -beta N_P(Y), so mu_D = s beta nu
are the cut multipliers in the original units. At every check the hook
writes L_add(mu_D) and h'mu_D. dnn.C certifies L + L_add(mu_D) + s Zc through
the existing path: GL with the margin over sum |L| plus the terms'
magnitudes, and Rump's Cholesky on V'Zc V. Then h'mu_D, rounded up, is
subtracted and the result rounded down. The dense term is summed in long
double, one rounding per entry. nu >= 0 by construction, so mu_D >= 0.

**Code (opt-in).** `DNNParams::ycuts` (dnn.h) is a hook. When it is NULL,
which is the default, nothing changes. When it is set, dnn.C forms X, calls
the hook's projection and does the Z-step in a separate pass. At every check
it certifies with the hook's L_add(mu_D) and h'mu_D. DNN_FUSE and
DNN_CORE_LAZYCERT are ignored there, and the accelerations, the prototype
certificate, DNN_CERTX and the child bounds are refused. `dnncut.C` holds
the regional-cut projector (`CutProjector`), its self-test and
`dnn_cp_bound`. The `dnncuts` driver gains these options:
- `-savepoolall R:F`: the whole pool after round R.
- `-cpadmm F`: the give-up point, then the CP solve, with `-cpiters`,
  `-cpsweeps`, `-cppasses`, `-cpdiag`, `-cpcheck`, `-cpquiet`, `-cpomega` and
  `-cpmu`.
- `-cpplain`: the same solve without the cuts.

## The box-cone paper (Charbonneau, Hu, Sotirov, arXiv 2610.41051)

Read Nov 7 from export.arxiv.org. The cheap projection needs inequalities
that share one slack. A row a'x <= b times every x_j and every 1 - x_j gives
t = b - a'x >= 0 and 0 <= u_j <= t. That is one box cone per row, projected by
a sort and a clip (their Theorem 1 and Proposition 1). None of our cut
families has this bound-factor form, so the guess is confirmed. The QAP's
rows are equalities, and facial reduction already removes their products
(their F Y = 0 is our V). The products of the bounds 0 <= x <= 1 are the
McCormick inequalities (their section 5.1). In the QAP DNN they follow from
Y >= 0 and the lifted row sums: Y_(ij),(kl) <= sum_l' Y_(ij),(kl') = x_ij, and
x_ij - Y_(ij),(kl) = sum_{l' != l} Y_(ij),(kl') <= sum_{l' != l} x_kl' = 1 - x_kl.
(The row sums hold because every column of a feasible Y lies in range(V).)
The 5-point and box (count) cuts are products of two slacks,
(q - b)(q - b - 1) >= 0, which the paper's section 5.2 places outside the
box-cone form. The local-7 cuts are facets of a 7-point hull, with no shared
slack.

The paper's Appendix A is closer to us. In correlation coordinates the four
triangles of a triple read |Z_ik + Z_jk| <= 1 + Z_ij and
|Z_ik - Z_jk| <= 1 - Z_ij, so triples with a common base pair share a slack.
Our regional triangles are that family on each region's moments
(Z_ij = 1 - 2p_i - 2p_j + 4W_ij), so they could be grouped. It would not
pay, for three reasons:
- The cones would be nearly scalar. In the round-12 pool, 11,101 of 12,027
  (region, triple)s hold one of their four inequalities [m]. Greedy base-pair
  cones hold 4.6 third vertices on average, median 3, and 852 of 2,601 hold
  one [m] (`pool_stats.py`).
- The projection lives in an image variable, and its coupling to Y is a
  solve with I + H'H. That solve is diagonal for Max-Cut, where each image
  entry is one matrix entry. For us each moment sums |S|^2 entries of Y and
  the 13 regions overlap.
- Our half-space steps on the moments already cost 0.56 ms per pass over
  13,000 cuts [m].

The paper's bounds are not certified ("should not be interpreted as
certified bounds"). Its MADAM runs recover scalar multipliers from the cone
multipliers for MADAM's own bound correction; our mu_D certificate is that
step.

## Validity

1. **Cut extraction** (`-cpcheck`, `cp_selftest` in dnncut.C), on node
   13's round-12 pool (13,000 cuts, `logs/cp_selftest_pool12_13.log`; the
   round-2 pool's run is `logs/cp_selftest_pool2_13.log`):
   - (a) L_add from the projector's terms against `assemble` at the same
     random multipliers. The max |diff| is 2.3e-12 on entries up to 2,352,
     1e-15 relative. Both have 3,816,180 nonzeros, and the magnitude sums and
     h'mu agree to 10 digits [m].
   - (b) 60 sampled cuts with their explicit sparse rows, read off the
     matrix `add_rcut` writes (441-10,374 nonzeros each). ||G||^2 is exact,
     and <G, Y> - h matches the moment formula and `rviolation` to 7e-15 [m].
   - (c) The P-step's grouped moments match `aggregate()` to 3e-15 relative
     [m].
   - (d) After a half-space pass (8,232 multipliers moved), the incrementally
     kept moments match `aggregate()` on y_P - G'(nu_new - nu_old) to 8e-15
     relative [m].
2. **Brute force** (`validate_cp.py`, `logs/validate_summary.txt`). The set:
   - 24 random instances at n = 7, 8, 9 in the four families of
     `tests/brute_force.py`, each at its root and at two random nodes, with
     prefix regions;
   - nug12 with its matrices swapped so that the 3 x 4 grid is the first
     matrix (grid 3:4 regions): its root (optimum 578) and four nodes with
     three facilities fixed (enumerated).

   Each node gets a 3-round Polyak pool, then three CP arms: 1 sweep; 3
   sweeps of 5 passes; and 1 sweep with the cutoff at optimum - 1.

   Results over 77 nodes and 231 runs [m]:
   - 12,656 certified values. max(certified - optimum) = -1.2e-9, and no run
     fathomed with the cutoff at the optimum.
   - All 75 fathoms at optimum - 1 have lb > cutoff.
   - The 3x5 arm closes a median 84% of the uncut gap on the 28 nodes with a
     gap, so the cuts bind.
   - On nug12's root CP-ADMM certifies 574.0 against the optimum 578. The
     uncut value is 567.9, and the Polyak rounds there drop to 485.7 [m].
   - A rerun with the final binary gives the same summary [m].
3. **The production path is unchanged** (`check_identity.sh`,
   `logs/ident_summary.txt`). The binaries of this tree, with the switch off,
   print the same traces as scratch builds of the parent commit 57adbe6 (every
   "it", "round", "quiet" and RESULT line, CPU stripped) [m]:
   - dnntest on nug20's root and wil50 nodes 19 and 10, in the plain
     settings, the production settings, and the production settings with the
     Nov 6 bundle (DNN_BETA_SCHED, LAZYCERT, NOSYMV, NOSYMW,
     REFRESH_PARTIAL): 9 of 9 identical;
   - the whole 5,500-iteration cut phase on node 13 (564 lines, `run_ref.sh`):
     identical;
   - with the final binaries again, a 4-round cut phase with a quiet round of
     300 on node 13 (196 lines, `CUT4=1`): identical.

   The CP code path with an empty pool reproduces the plain solve bit for
   bit (nug12 and node 13), so the restructured Y/Z step is exact [m]. The
   dnntest traces, the 4-round cut phase and the empty-pool checks were run
   with the final binaries (`logs/ident_summary.txt`).
4. **The large values, cross-checked through the standard path**
   (`logs/cp_s1p1_q_13.log`). After 1,500 CP iterations (certified
   48,762.50), the final mu_D were frozen as L_add, and the plain solve ran
   1,000 more iterations from the CP state (Z warm, Y = 0 as in every
   Polyak round). That is the cut rounds' own certificate path, with no
   hook. It certifies 48,756.97, above anything the Polyak protocol
   certified (48,755.40) [m]. Its Lagrangian primal settles at
   48,767.3-48,769.1, which is about the dual value L(mu_D). The CP
   certificate with the same mu_D (48,762.50) lies below that, as a lower
   bound on L(mu_D) must [m]. The run also repeats the s1p1 trace exactly up
   to t = 1,500, so the runs are deterministic [m].

5. **Gates** on the final code (`logs/gates.txt`) [m]:
   - `make test`: brute force 4,536 of 4,536, resume 5 of 5, parallel and
     coordinator all passed (the AWS scenario skipped, no boto3 here), and
     the regression check exact;
   - `make dnncheck`: OK;
   - `make shipcheck`: ok, ship64 identical to the straight run;
   - `make continuationcheck`: PASS;
   - `python3 bench/bench.py check`: 7 of 7 match the reference exactly.

## Results

### The curve (node 13, the dual test's settings)

Best certified value so far against iterations after the give-up point. The
"Polyak" column is the run on this box (`logs/polyak_ref_13.log`). The laptop
log (`stronger_bounds/dual/tri_13.log`) has round ends only; it gives 48,686.4
at t = 900, 48,708.97 at 3,600 and 48,755.42 after the quiet round, within
0.04 of this box's run. "uncut" continues the give-up point without cuts
(`-cpplain`).

| t | Polyak, best so far | Polyak, at the check | CP 1 sweep x 1 pass (s1p1) | CP 1 x 5 (s1p5) | CP 3 x 3 (s3p3) | uncut |
|---|---|---|---|---|---|---|
| 0 | 48,596.8 (give-up point) | | | | | |
| 100 | 48,596.8 | 48,559.4 | 48,579.9 | 48,575.6 | 48,591.3 | 48,524.1 |
| 200 | 48,614.2 | 48,614.2 | 48,651.6 | 48,649.2 | 48,664.2 | 48,580.9 |
| 300 | 48,636.5 | 48,636.5 | 48,685.2 | 48,683.6 | 48,697.5 | 48,603.1 |
| 500 | 48,643.9 | 48,643.9 | 48,719.9 | 48,718.7 | 48,730.4 | 48,621.3 |
| 750 | 48,666.1 | 48,646.0 | 48,741.1 | 48,739.9 | 48,748.8 | 48,628.8 |
| 1,000 | 48,686.4 | 48,595.6 | 48,752.0 | 48,751.0 | 48,758.2 | 48,631.3 |
| 1,500 | 48,686.4 | 48,681.8 | 48,762.5 | 48,761.6 | 48,766.9 | 48,632.6 |
| 2,000 | 48,692.5 | 48,679.1 | **48,767.1** | 48,766.4 | 48,769.8 | 48,632.7 |
| 2,500 | 48,702.3 | 48,638.7 | 48,768.9 | 48,768.3 | 48,770.9 | |
| 3,000 | 48,703.3 | 48,703.3 | 48,770.5 | 48,770.0 | 48,771.7 | |
| 3,600 (round 12 ends) | 48,709.0 | 48,709.0 | | | | |
| 5,100 (quiet round ends) | **48,755.4** | 48,755.4 | | | | |

All values [m]; `analyze.py` prints this table from the logs. The Polyak
column at the check falls after overshoots: rounds 4 and 9 lowered the
bound, and each round restarts Y at 0.

Iterations and CPU after the give-up point to first certify a level [m]:

| level | Polyak | CP s1p1 | CP s3p3 | CP s1p5 |
|---|---|---|---|---|
| 48,713 (give-up + 116) | 3,930 it, 2,166 s | 450 it, 241 s | 380 it, 208 s | 460 it, 248 s |
| 48,743 (kill line - 20) | 4,370 it, 2,401 s | 790 it, 417 s | 650 it, 354 s | 810 it, 431 s |
| **48,753 (kill line)** | **4,860 it, 2,659 s** | **1,040 it, 547 s** | **850 it, 461 s** | **1,070 it, 565 s** |
| 48,755.4 (Polyak's final) | 5,090 it, 2,779 s | 1,120 it, 588 s | 920 it, 498 s | 1,170 it, 617 s |
| 48,763 | not reached | 1,550 it, 810 s | 1,230 it, 665 s | 1,630 it, 856 s |

All [m]. CPU is thread CPU from the give-up point, setup included. To the
kill-line level, s1p1 takes 4.7x fewer iterations and 4.9x less CPU than
Polyak, and s3p3 5.7x and 5.8x. More half-space passes per sweep do not help
(s1p5 runs a little behind s1p1). More sweeps per iteration do help (s3p3).
The extra sweeps in s3p3 cost 19 ms per iteration [m], and it still needs
16% less CPU than s1p1 to reach the kill-line level.

### Overhead per iteration

Milliseconds per iteration, CPU (the dnn.C profile at verbose 2 over each
solve; `logs/cp_*_13.log`) [m]:

| part | plain (`-cpplain`) | CP s1p1 | difference |
|---|---|---|---|
| forward transform | 38.7 | 37.8 | -0.9 |
| eigensolver (R-step) | 403.8 | 389.2 | -14.6 |
| expand | 26.9 | 27.3 | +0.4 |
| Y/Z-step | 27.2 | 50.3 | +23.0 |
|   of it, the projection (2 P-steps of 8.7 ms, 1 pass of 0.56 ms) | - | 18.4 (19.5 with the residual diagnostic) | +18.4 |
| certificate, per iteration (every 10th) | 12.9 | 24.2 | +11.3 |
| **total (solve CPU / iterations)** | **510.2** | **529.5** | **+19.3 (+3.8%)** |

The CP-specific work is +34.3 ms (+6.7%) [m]: +23.0 in the Y/Z-step and
+11.3 at the checks. Of the Y/Z-step's part, 18.4 ms is the projection. The
rest is the two passes that the fused plain step does not need (forming X and
the separate Z-step). A CP check costs 242 ms against 129 ms [m]. It forms
the cut term, which is the dense L_add plus three passes over the n1^2
matrix. And from about t = 1,800 the value sits within 48.8 of the cutoff,
where GL switches to exact inner LAPs (`DNN_LAP_AUTO`). The split between the
two was not measured. The eigensolver is 14.6 ms cheaper on the CP iterates:
they track a smaller block (665 vectors on average against 690, 560
negative eigenvalues against 580) [m]. The net is +3.8%. s3p3 costs
550.1 ms (+7.8%) [m]. The kill line is 30%, so this passes.

### How far the projection is from converged

The Y-step's output after one sweep against a 30-sweep projection from the
same start, every 200 iterations (`-cpdiag 200:30`), and the cut residuals of
the output [m]:

| t | arm | norm(Y - Y_30) / norm(Y_30) | max violation / norm(G), output | the same after 30 sweeps | cuts violated by > 1e-4 | active cuts (nu > 0) |
|---|---|---|---|---|---|---|
| 200 | s1p1 | 3.8e-3 | 1.0e-3 | 1.5e-4 | 2,675 | 2,228 |
| 1,000 | s1p1 | 1.5e-3 | 4.1e-4 | 5.6e-5 | 1,619 | 2,530 |
| 2,000 | s1p1 | 8.1e-4 | 2.7e-4 | 4.9e-5 | 1,186 | 2,687 |
| 3,000 | s1p1 | 4.8e-4 | 1.7e-4 | 3.9e-5 | 921 | 2,740 |
| 1,000 | s3p3 | 2.3e-4 | 1.6e-4 | 6.8e-6 | 1,820 | 2,517 |
| 3,000 | s3p3 | 5.4e-5 | 3.8e-5 | 1.5e-6 | 778 | 2,632 |

One sweep is not a converged projection. Its output is 0.1-0.4% away from
the 30-sweep projection, and 30 sweeps leave the max violation 5-8x lower.
That does not matter here, because the ADMM itself is far from converged
(primal residual 2e-2 at t = 1,000) and every iteration warm-starts the
next. The inexactness shrinks as the iterates settle. About 2,500 of the
13,000 cuts carry a positive multiplier [m]. The certificate does not depend
on how well the projection converged: any nu >= 0 is valid.

### Production settings

The production DNN settings (`DNN_AMX=1 DNN_CHOLQR_PASSES=1`, degree 3,
tight lambda 2, refresh 40, gamma 1/1, and the bundle `DNN_R32_NONCERT=1
DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_CORE_LAZYCERT=1
DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1`), with the same
protocol, node and kind of pool. The Polyak run writes its own round-12 pool
(`logs/pool12_prod_13.txt`). The lazy certificate certifies only the decisive
checks, so the Polyak curve comes from a second run without it
(`polyak_prodnl_13.log`). Its iterates are the same: round 12 ends at
48,729.17 and the quiet round at 48,762.09 in both [m].

| t | Polyak, best so far | CP s1p1 | CP s3p3 |
|---|---|---|---|
| 300 | 48,651.7 | 48,688.9 | 48,702.4 |
| 500 | 48,663.0 | 48,719.6 | 48,732.1 |
| 1,000 | 48,702.7 | 48,748.8 | 48,757.0 |
| 1,500 | 48,702.7 | 48,758.7 | 48,764.5 |
| 2,000 | 48,711.5 | **48,763.8** | 48,768.3 |
| 3,000 | 48,724.8 | 48,768.8 | |
| 3,600 (round 12 ends) | 48,729.2 | | |
| 5,100 (quiet round ends) | **48,762.1** | | |

| level | Polyak | CP s1p1 | CP s3p3 |
|---|---|---|---|
| 48,713 | 2,070 it, 656 s | 450 it, 149 s | 360 it, 126 s |
| **48,753** | **4,290 it, 1,348 s** | **1,170 it, 380 s** | **860 it, 296 s** |
| 48,763 | not reached | 1,920 it, 623 s | 1,350 it, 461 s |

All [m] (`cp_s1p1_prod_13.log`, `cp_s3p3_prod_13.log`). To the kill-line
level, s1p1 takes 3.7x fewer iterations and 3.5x less CPU, and s3p3 5.0x and
4.6x. The Polyak protocol does better here than in the dual test's
settings: round 12 ends at 48,729.2 (against 48,709.0), and the quiet round
passes 48,753 after 690 of its iterations (against 1,260). CP-ADMM still
gets to 48,753 in 27% of the iterations.

Per iteration, against three baselines [m]:
- the plain solve with the lazy certificate (`cp_plain_prod_13.log`):
  0.2891 s. CP s1p1 costs 0.3282 s, **+13.5%**; s3p3 0.3439 s, +19.0%.
- the Polyak rounds with the lazy certificate: 0.3000 s. CP s1p1 is +9.4%.
- the Polyak rounds certifying every check, as CP-ADMM does: 0.3154 s. CP
  s1p1 is +4.1%.

The Y/Z-step adds 23.0 ms (the projection is 18.3) and the checks add 19.6
ms [m]. All three readings are inside the kill line.

### Warm multipliers

`-cpmu` starts nu from the pool's own multipliers, the Polyak round-12 mu
(`logs/cp_s1p1_mu_13.log`). This is the situation of a pooled child whose
parent hands down its multipliers. Z is still the uncut give-up point's, and
Y = 0 [m]:

| t | 100 | 300 | 500 | 1,000 | 1,500 |
|---|---|---|---|---|---|
| s1p1 from nu = 0 | 48,579.9 | 48,685.2 | 48,719.9 | 48,752.0 | 48,762.5 |
| s1p1 from Polyak's mu | 48,592.4 | 48,695.1 | 48,726.9 | 48,755.6 | 48,764.6 |

The head start is +10 to +13 early and +2 to +4 from t = 1,000 on.
It reaches 48,753 at t = 920 against 1,040 [m]. The multipliers alone carry
little because they are not consistent with the uncut Z. The first checks
are lower (48,479 against 48,488 at t = 50). A pooled child would also bring
the parent's Z.

### Variants

All on node 13 in the dual test's settings, from the give-up point [m]:

| arm | 48,753 at | t = 1,000 | t = 2,000 | per iteration | note |
|---|---|---|---|---|---|
| s1p1: 1 sweep, 1 pass | 1,040 | 48,752.0 | 48,767.1 | 0.530 s | the memo's design |
| s1p5: 1 sweep, 5 passes | 1,070 | 48,751.0 | 48,766.4 | 0.531 s | passes do not help |
| s3p3: 3 sweeps, 3 passes | 850 | 48,758.2 | 48,769.8 | 0.550 s | best per CPU |
| s1p1, omega = 1.5 | not reached | 48,746.0 | | 0.526 s | over-relaxed steps: ahead to t = 850, then oscillates |
| s1p1, round-2 pool (3,000 cuts) | 1,150 | 48,749.1 | 48,763.2 | 0.523 s | the first 3,000 cuts carry most of it |
| s1p1 from Polyak's mu | 920 | 48,755.6 | | 0.535 s | see above |

The over-relaxed run's half-space step grows from 5e-4 to 8e-2 after
t = 850, and its certificate falls back to 48,640-48,660
(`cp_s1p1_w15_13.log`). The round-2 pool is the pool after round 2's cuts
in the same Polyak run (`logs/pool2_13.txt`, `polyak_pool2_13.log`). It ends
4 below the full pool at t = 2,000. A CP phase would therefore need only a
few separations, early.

## Decision

Against the memo's kill line:
- **Within 10 of 48,763 by 2,000 iterations: passes.** s1p1 certifies 48,753
  at t = 1,040 and 48,767.1 by t = 2,000 [m]. The production settings give
  1,170 and 48,763.8 [m]. On d2 the same level (within 10 of the dual test's
  estimate) comes at t = 880 against 4,670 [m]. Node 14 fathoms at t = 770
  against 4,390 [m].
- **Per-iteration overhead at most 30%: passes.** It is +3.8% in total and
  +6.7% for the CP work itself, in the dual test's settings [m]. In the
  production settings it is +13.5%, of which 6.8 points come from
  certifying every check [m].

**Decision: pass.** The memo expected 1.3-2x fewer phase iterations [e]. The
measured gain to the Polyak protocol's plateau on node 13 is 4.7x in
iterations and 4.9x in CPU with one sweep per iteration, and 5.7x and 5.8x
with three [m]. Nodes 14 and d2 give 5.7x (to the fathom) and 4.8x (to
Polyak's final value) with one sweep [m].

Why it works. The Polyak phase moves the multipliers once per 300
iterations, along the subgradient at a 300-iteration iterate. Two of the 12
steps overshoot: rounds 4 and 9 lower the bound [m]. Each round starts Y at
0, and the certificate ends each round 16-54 below the Lagrangian primal;
the quiet round recovers 46 [m]. CP-ADMM has no rounds. The
multipliers are the projection's own KKT multipliers, they move a little at
every iteration, and they converge together with Z to the saddle point of
the cut-constrained DNN. The memo's worry was the fast60 proxy, where
frequent Polyak moves did worse. That proxy took noisy subgradient steps
from inexact solves. Dykstra is a projection, and its multipliers are
consistent with the current iterate at every step.

The value is also higher. CP-ADMM certifies 48,771.7 at t = 3,000 and is
still rising [m]. The Polyak protocol's 48,755.4 and the dual test's slow
30-round ascent (48,763 on the laptop's pool) were both below this pool's
dual optimum. The regional triangles are still capped well short of the
cutoff on this node: 42 short instead of 59 [m]. The dual test's conclusion
(the family is capped) stands. But the cap on this pool is at least 48,771.7.
That is 16 above the Polyak protocol's value and 9 above the slow ascent's
48,763 on the laptop's pool. On d2 the cap is at least 48,807.9 [m], 6 short
of the cutoff, against the slow ascent's 48,795. With fresh cuts separated
during the solve, d2 is the node to watch.

## What a production version would need

Not built tonight. The pieces:

1. **Pool changes between rounds.** Today each round adds 1,000 cuts,
   separated on the round's last Y. A CP phase is one solve. It would
   separate on its own iterate every k iterations and append the new cuts
   with nu = 0. The node's regions do not change, so the projector's
   grouping stays. A pass costs 0.56 ms per 13,000 cuts [m], so the pool's
   size is not a cost. Cuts whose nu stays at 0 can be dropped. Not
   measured: whether separating inside the solve keeps the speed. That is
   the first test, against this study's frozen-pool runs, on node 13 and on
   d2 (6 short at its frozen pool's cap). The round-2 pool's run suggests a
   few early separations carry most of the gain [m].
2. **Pooled children.** Remap the parent's pool as today (`remap_pool`, with
   zero multipliers kept) and start nu = mu / (s beta) from the parent's
   mu_D (`-cpmu` does this). From the Polyak round-12 multipliers it gains
   +10 at t = 200-300 and +3.6 at t = 1,000 over nu = 0 [m]. The larger step
   would be the child's own first solve in CP form. Today that solve freezes
   the inherited multipliers (`dnn_pool_bound`), and the first solves are
   61% of the iterations in `bestcore` [m, `pipeline_ab/README.md`].
3. **The quiet round.** Not needed. A quiet round from the CP state at
   t = 1,500 (mu_D frozen, Y restarted) reaches 48,757.0 after 1,000
   iterations [m]. Continuing the CP solve for those 1,000 reaches 48,768.9
   [m]. The quiet round exists because each Polyak round restarts Y and its
   certificate trails the Lagrangian. A CP phase is one solve with no
   restarts.
4. **The warm start.** The CP state is (Y, Z, nu). Within a phase there is
   one solve, so there is no restart and QAP_DNN_CUT_WARMY is moot. For the
   children, hand on Z (and Y under QAP_DNN_WARM_Y) as today, plus mu_D as
   the pool's multipliers (`export_pool` with mu = mu_D). For a stop and
   resume, nu goes into the shipped state (QAP_DNN_SHIP).
5. **Stopping rules.** The stall rule and the primal gate count rounds.
   A CP phase needs per-check versions: the gap-ratio rule on the running
   certified best, and the gate on the CP primal, re-tuned.
6. **Cheaper checks.** Certifying every check is 6.8 of the 13.5 points of
   production overhead [m]. The lazy rule carries over: the polyhedral value
   with the current mu_D is an upper bound on the certified one.
7. **Scope.** The projector handles regional cuts on side 0 (wil50). nug30's
   side-1 grid and atomic cuts need a second P-step loop or a transposed
   layout. The P-step implements eliminated assignments as dnn.C does, but no
   run here exercised them.

**Expected effect on the cut phase's share.** In the production `bestcore`
arm on the node-10 subtree, the phases are 6,200 of 15,963 DNN iterations
(39%) [m, `pipeline_ab/README.md` run 1]. They average about 1,030
iterations, and none fathoms: they are stall-stopped misses that hand their
pools down [m, run 2's `wlowc`, 6 phases with 0 fathoms]. Two readings:
- **Same bound, fewer iterations.** A 1,000-iteration Polyak phase gets
  node 13 to 48,686.4 (round 3, t = 900). CP-ADMM is there by t = 310 [m]:
  2.9x fewer iterations, at +14% per iteration. The phases would cost about
  0.39 of today, and the subtree's DNN work would fall to
  0.61 + 0.39 x 0.39 = 0.76 of today, about 1.3x [e]. To round 12's value
  (48,709.0 at t = 3,600) it is 8.4x fewer iterations (t = 430) [m].
- **Same budget, higher bound.** In a 1,000-iteration budget CP-ADMM
  certifies 48,752 against 48,686, +66 [m]. Phases that end within about 60
  of the cutoff today would fathom instead of handing a pool down, and the
  pooled children would inherit better multipliers [e]. The effect on the
  tree needs a search run to estimate.

Either way the next step is the search, not more single-node runs: the
node-10 A/B of `pipeline_ab` with a CP phase (items 1, 2, 5) against
`bestcore`, scored on DNN CPU and fathoms.

## Separation inside the solve (follow-up, Nov 7)

**Why.** The frozen round-12 pool hands CP-ADMM all 13,000 cuts at t = 0.
The Polyak run had them only by t = 3,300, separated on its own later
iterates. The honest comparison separates during the CP solve, at the
Polyak protocol's cadence.

**Method** (`-cpsep k:K[:CAP[:DROP]]`, in dnncut.C and dnncuts.C; dnn.C is
untouched).
- The start pool is round 0's own separation at the give-up point: 1,000
  triangles separated on the give-up point's Y, exactly as round 0 of the
  cut rounds does. That is the same as starting from an empty pool and
  separating at t = 0 on the give-up point's Y. It is also the
  coordinator's pool-0 file (9 regions, 1,000 cuts).
- Every k iterations, the solve separates on that iteration's Y-step
  output: the projected Y (in P), which is the iterate the ADMM carries.
  It is the same object as a Polyak round's last Y, the output of the
  round's last Y-step.
- The separator and ranking are dnn_cut_rounds': the four triangle types on
  every i < j < k of every region, violation above 1e-4, sorted by
  violation, the K most violated that are not already in the pool (the same
  key). The moments come from the P-step that produced the Y, which the
  self-test checks against `aggregate()`.
- New cuts enter with nu = 0 at the end of the cyclic order. The region
  grouping does not change (a node's regions are fixed), and the sums are
  rebuilt from nu at every Y-step.
- Options: CAP bounds the pool (zero-multiplier cuts leave first, the
  longest at zero first); DROP drops a cut whose nu has stayed 0 for that
  many iterations (checked at each separation). Neither was needed for the
  runs below (both 0).

**Check.** Through t = 300 the cadence-300 run's trace equals the
coordinator's pool-0 run (`cp_s1p1_pool0_13.log`) bit for bit. Through
t = 100 so does the cadence-100 run's [m]. So the start pool is the same,
and the separations are the only difference.

### Node 13, the dual test's settings

The no-lookahead baseline is the coordinator's run (`cp_s1p1_pool0_13.log`,
commands in `polyak_r1_13.log` and `run_pool0.sh`'s step 2): CP s1p1 on round
0's 1,000 cuts, frozen, 2,000 iterations. The separating runs are
`cp_s1sep300_13.log` and `cp_s1sep100_13.log`, s1p1 with
`-cpsep 300:1000` and `-cpsep 100:1000`, 3,000 iterations (`run_sep.sh`).
Best certified value so far [m]:

| t | Polyak | frozen round-12 pool | pool-0, frozen | separating, every 300 | separating, every 100 |
|---|---|---|---|---|---|
| 100 | 48,596.8 | 48,579.9 | 48,568.3 | 48,568.3 | 48,568.3 |
| 300 | 48,636.5 | 48,685.2 | 48,664.7 | 48,664.7 | 48,680.1 |
| 500 | 48,643.9 | 48,719.9 | 48,693.1 | 48,707.0 | 48,716.8 |
| 1,000 | 48,686.4 | 48,752.0 | 48,715.0 | 48,747.7 | 48,750.7 |
| 1,500 | 48,686.4 | 48,762.5 | 48,720.6 | 48,760.7 | 48,762.0 |
| 2,000 | 48,692.5 | 48,767.1 | 48,722.4 | 48,766.3 | 48,766.8 |
| 3,000 | 48,703.3 | 48,770.5 | | 48,770.3 | 48,770.4 |
| pool at the end | 13,000 | 13,000 | 1,000 | 4,025 | 4,594 |

First iteration (and CPU since the give-up point) to certify a level [m]:

| level | Polyak | frozen round-12 pool | pool-0, frozen | every 300 | every 100 |
|---|---|---|---|---|---|
| 48,686.4 | 1,770, 970 s (round 3 ends at 48,686.36, t = 900) | 310, 167 s | 440, 248 s | 390, 221 s | 330, 187 s |
| 48,709.0 (round 12) | 3,600, 1,984 s | 430, 231 s | 770, 431 s | 520, 292 s | 450, 252 s |
| 48,713 | 3,930, 2,166 s | 450, 241 s | 910, 510 s | 550, 308 s | 480, 268 s |
| 48,743 | 4,370, 2,401 s | 790, 417 s | not reached | 900, 499 s | 830, 459 s |
| **48,753** | **4,860, 2,659 s** | **1,040, 547 s** | not reached | **1,150, 633 s** | **1,070, 589 s** |
| 48,763 | not reached | 1,550, 810 s | not reached | 1,670, 916 s | 1,590, 874 s |

**In-solve separation keeps the frozen pool's speed.**
- With a separation every 100 iterations, CP-ADMM reaches 48,753 at
  t = 1,070. The frozen round-12 pool reaches it at 1,040 and the Polyak
  protocol at 4,860. That is 4.5x fewer iterations and 4.5x less CPU than
  Polyak [m].
- At the Polyak cadence (300) it is t = 1,150, 4.2x [m].
- At t = 3,000 both separating runs are within 0.2 of the frozen pool [m].
- The pool-0 baseline stops at 48,722.4, its 1,000-cut cap [m]. The cuts
  separated later are what carry the rest, and the solve finds them on its
  own iterates.

The pool stays small. Separations at t = 600 already find only 1,900-3,000
violated triangles, most of them in the pool, and the pool ends at
4,000-4,600 cuts, not 13,000 [m]. The separation itself costs 3 ms per call,
0.1 s over 30 calls (0.03 ms per iteration) [m].

Per iteration, solve CPU over iterations [m]:

| run | s per iteration |
|---|---|
| plain | 0.510 |
| frozen round-12 pool | 0.530 |
| every 300 | 0.548 |
| every 100 | 0.549 |
| pool-0 (the coordinator's run, under other load) | 0.560 |

Against the frozen run, the Y/Z-step (50.8 against 50.3 ms) and the checks
(24.9 against 24.2 ms) are the same. The difference is in the eigensolver
(403.5 against 389.2 ms) and forward (40.1 against 37.8 ms): the iterates
and the load differ, not the separation. Against the plain iteration it is
+7.5%, well inside the kill line's 30% [m].

SEP_MORE

## Files

- `run_ref.sh`: the Polyak protocol on node 13, with the old binary and with
  this tree's, which also writes the round-12 pool. Modes `prod` and
  `prodnl` run the production settings, with and without the lazy
  certificate.
- `run_cp.sh`: the CP arms (s1p1, s1p5, s3p3, plain, `*_prod`; QUIET,
  EXTRA and TAG for the variants).
- `run_node.sh`: the same runs on another node (d2, 14): the Polyak
  reference with its pool, CP s1p1, and CP s3p3 (`cp3`).
- `validate_cp.py`: the brute-force validity check.
- `check_identity.sh`: the switch-off identity check. The "old" binaries in
  it and in `run_ref.sh` are scratch builds of the parent commit 57adbe6
  (`latest/` and `experiments/dnn_cuts/dnncuts.C` from `git archive 57adbe6`,
  then `make dnntest dnncuts`). Set OLD, OLD_DNNTEST and OLD_DNNCUTS to use
  others.
- `analyze.py`: the curve, level, overhead and residual tables from the logs.
  For example: `python3 -I analyze.py s1p1 s1p5 s3p3 plain` (node 13);
  `--ref polyak_prodnl_13.log s1p1_prod s3p3_prod`; or
  `--ref polyak_new_d2.log --node d2 --levels 48785 48790.8 -- s1p1 s3p3`.
- `pool_stats.py`: the pool's structure, for the box-cone question.
- `logs/`:
  - `polyak_{ref,new,prod,prodnl}_13.log`, `cp_*_13.log`, and the d2 and 14
    logs (`polyak_new_{d2,14}.log`, `cp_s1p1_{d2,14}.log`, `cp_s3p3_d2.log`);
  - `pool12_13.txt`, `pool12_prod_13.txt`, `pool12_d2.txt`, `pool12_14.txt`,
    `pool2_13.txt`: the frozen pools;
  - `cp_selftest_pool{2,12}_13.log`, plus `polyak_pool2_13.log`, the run
    that made the round-2 pool;
  - `validate/`: samples (nug12's root and one random instance) and the
    swapped nug12 instance; `validate_summary.txt`;
  - `ident/`, `ident_summary.txt`, `gates.txt`.
- Code: `latest/dnn.h` and `dnn.C` (`DNNYCuts`, `DNNParams::ycuts`),
  `latest/dnncut.h` and `dnncut.C` (`CutProjector`, `cp_selftest`,
  `dnn_cp_bound`, `DNNCPParams`, `-savepoolall`'s `pool_all_round`), and
  `experiments/dnn_cuts/dnncuts.C` (the options).
