# Facility side and disjunctive strengthening: working notes (agent 3, core 3)

Started Nov 4 2026, 05:27. Append-only log of what was read, found, run.

## 05:30 Read
- latest/dnn.C header, latest/dnncut.h, experiments/dnn_cuts/dnncuts.C header
- stronger_bounds/README.md steps 1-13, research_bound_quality.md (1.3, 1.4, 3.4-3.6, 5),
  roadmap/bound_strength.md, NOTES.md sections 1-5, triples/ and local7/ READMEs
- Key facts: regional triangles on the 13 grid thresholds (side 0 = first matrix =
  Manhattan grid) give +158.6 on node 13, +131 on node 10, +133 on d2 after the quiet
  round (laptop dual/tri_*.log); 5-point, triples, local-7 all washes or negative.
- research_bound_quality.md section 5 lists "facility-side regions on wil50" as
  do-not-bother (flow matrix random 0-9, no decomposition identity). To be tested.
- Machine: Xeon @2.8 GHz Cascade Lake (avx512, NO AMX), 4 cores, 15 GB.

## 05:35-05:45 Setup
- Conventions checked: first index = rows of A = wil50's 5x10 Manhattan grid (index r*10+c),
  second = rows of B = flows. Here "location" = first index, "facility" = second index;
  a fix a:b = facility b on location a; lifted index of x_(a,j) = 1 + a*m + j (local).
  scripts/ynode.py reproduces the logged primal from a dumped Y to 1e-6 (plumbing node).
- Speed here (no AMX on this Cascade Lake): dnntest node 13, -d 4: 477 ms/it with DNN_AMX=1;
  -d 3 -L 2 -r 40: 389 ms/it with DNN_AMX=1, 346 ms/it with DNN_AMX=0. With the grid agent's
  run on another core, 0.56 s/it measured in the phase. Using DNN_AMX=0 -d 3 -L 2 -r 40.
- Built src/fdrv (src/build.sh): links latest/*.o read-only + src/dnncut_f.C, a renamed copy of
  latest/dnncut.C made by src/make_dnncut_f.py (adds an export hook only). Modes: -phase
  (triangle protocol, stops before the quiet round, saves the pool's Lagrangian term + ADMM
  state), -cont (continue from a saved state with extra fixed terms), -pool (Polyak on an
  external sparse cut list on top of the frozen triangle term). Plumbing (m=10 node):
  fdrv -phase + -cont reproduces dnncuts' quiet round to all printed digits (44974.9659).
- Flow matrix (logs/flow_analysis.log): statistically a uniform random symmetric 0-9 matrix.
  Spectrum top 223.8 then 34.4 (random matrix: 227.4, 31.7); centered top 39.3 vs 36.9 random
  (semicircle edge 41.3); triangles with all flows >= 9: 32 (random expectation 29), >= 8: 168
  (165); densest 4/6/10/15-subsets as dense as in a random 0-9 matrix (mean flow 9.0/8.13/7.0/6.3
  vs 8.83/8.0/6.93/6.16). No clusters, no block structure: facility "clusters" are noise.
- Literature fetches: arxiv.org and kups.ub.uni-koeln.de do not resolve from this session
  (WebFetch blocked); WebSearch snippets work. Literature is cited from search results + memory.
- 05:44 launched runs/n13 (run_phase.sh, 8 rounds; 12-round run killed: at 0.56 s/it it would
  have taken ~37 min, above the 30-min limit).

## 05:50-06:00 Diagnostics at node 13's uncut converged iterate
Iterate: the grid agent's dump new_cuts/grid/runs/n13_uncut.Y.bin (dnntest 1,500 it, -d 3 -L 2
-r 40, lb 48,627.29, primal 48,628.38), copied to runs/n13_uncut1500.Y.bin (read-only use).
scripts/ydiag.py (logs/ydiag_n13_uncut1500.log), box_sep.py (logs/box_n13_uncut1500.log),
sep_dist.py --fam hull3 (logs/hull3_n13_uncut1500.log), gapdecomp.py (logs/gapdecomp_n13_uncut1500.log):
- Facility triples: min over all 16,215 triples of Delta_jk + Delta_jl + Delta_kl = 4.79 (valid >= 4,
  grid bipartite + injective): parity cut never violated. Checkerboard region: no triangle violated
  (max -0.0067). Facility regions (random facility sets, s over locations): no triangle and no box
  inequality violated (max -0.03..-0.10).
- Threshold regions: max triangle violation 0.034-0.080, 41-2,530 violated per region (as in the record).
  General box inequalities (any support, local search) on the thresholds: best normalized violation
  (half the box measure = triangle units) 0.03-0.064 with supports 3-19, i.e. no larger than the
  triangles. Parity regions: checker 0, rowpar 0, colpar 0.023 (support 3).
- Triple distance hull (LP over the 520 achievable distance triples): violated only by metric
  triangle facets H = (1,-1,1) (max 0.49, 20 of the 2,000 smallest-sum triples), which are sums of
  the 13 thresholds' apex triangles on the same facility triple: implied by the regional family.
- Gap vs a 49,036 local-search solution (scripts/node_ls.py): the gap (407.6) is a small difference
  of +6,424 / -5,807 over pairs; heavy pairs (B 8-9) sit at mean expected distance 3.46 in the DNN vs
  3.62 in the solution, light pairs (0-2) at 6.16 vs 6.01. Per-facility shares are uncorrelated with
  flow degree (corr 0.00): no facility group carries the gap; it is the diffuse bisection smoothing
  that the threshold cuts address.
- Lines study (tuning/lines/work, 300-it cold children): node 10's best lines row 49 and col 45 share
  their lowest child 49:45 (48,625.0 at 300 it); min child per line is parent -37..+57 at 300 it.

## 06:00 More diagnostics at node 13 uncut converged (same iterate)
- Facility-region triangles with T OPTIMIZED (scripts/facreg_search.py: for each location triple
  the most violated T is a 47-variable QUBO, multi-start local search; 60 location triples x 4
  types, half of them compact): max violation 0.0062 (logs/facreg_n13_uncut1500.log). So no
  facility region at all, including any region of a cut decomposition of B + c(J - I) (which
  would make the objective linear in facility-regional moments, as the 13 thresholds do on the
  location side), has a violated triangle.
- Dispersion (scripts/disp_check.py, logs/disp_n13_uncut1500.log): smallest DNN sum of expected
  pairwise distances over k facilities vs the grid minimum sigma_k: k=4 10.74 vs 8, k=5 20.62 vs 16,
  k=6 33.33 vs 25, k=7 50.05 vs 38. Slack 2.7-12: the DNN never places a group near compactness.
  "Facility groups placed too compactly" is refuted at this iterate.
- Lines study x-weighted (logs/lines_xweighted.log): on every line of the three parents, the
  x-weighted mean of the children's 300-iteration bounds is above the parent's converged bound on
  most lines (parent 1: +66..+124; parent 3: -21..+56; parent 10: -12..+34), so the parent's
  optimum is outside the hull of the line's children (if Y* = sum x_c Y_c then V >= sum x_c V_c):
  one level of lift-and-project on a line does cut it off, but the closure's value is the MIN
  child, which is the child the DNN favours most (x = 0.30-0.47 on the best lines; node 10's best
  lines row 49 and col 45 share the bottleneck child 49:45, x = 0.47).
- Box inequalities with BOTH sides free (scripts/box2_search.py, logs/box2_n13_uncut1500.log):
  count c = |pi(Q) cap P| for any location set P and facility set Q (Juenger-Kaibel boxes,
  unsigned), alternating QUBO search from 13 threshold starts + 17 random starts: violated only
  when P is a threshold (max 0.031 triangle units, |Q| = 15 and 12); every random-start search
  ends at zero violation. Together with facreg_search: at this iterate the DNN is consistent on
  every count except those of the 13 threshold regions, i.e. it "cheats" only in the variables
  the objective is written in.
- Cross-check: the grid agent (new_cuts/grid/NOTES.md, 05:30-05:55) finds on the same iterate the
  same triple minimum 4.79, zero violated parity-region triangles, col-parity max 0.023, and
  flow-weighted sub-QAP slack 140-1,700: independent agreement.
- Facility 45, the xbar line at all three lines-study parents (tuning/README.md calls it
  "location 45"; in the lines code "col" is the second index, i.e. a FACILITY), is the LIGHTEST
  facility of wil50 (flow row sum 178, rank 50 of 50; mean 222). The DNN sends it to the
  periphery: at node 10 its largest marginal is location 49 (the corner (4,9)), x = 0.47, which is
  also the bottleneck child of the best lines. Facility-side reading of the branching data: the
  lightest facility's line is the most decided one; rule 8/9 (xbar) already finds it.
- Node 13 uncut: corr(flow degree, DNN expected eccentricity of the facility's location) = -0.90,
  corr(flow degree, sum_a x_aj^2) = -0.48: heavy facilities central, light ones peripheral and more
  decided. This first-order structure is the GL/star part, already exact in the DNN.

## 06:07 Note on the deliverable
- Writing REPORT.md with the Write tool is refused by the harness for subagents ("return findings
  as text"). The report will be delivered as the text of the final hand-back to the coordinator,
  who can commit it as REPORT.md. NOTES.md (this crash-recovery log) is kept as asked.

## 06:12 Queue (runq.sh, jobs.txt) and plan
- runs/n13: triangle phase (8 rounds, -d 3 -L 2 -r 40; rounds track dual/tri_13.log to +-1:
  r1 +39.4, r2 +68.9, r3 +89.2, r6 +96.9, r7 +103.6), then quiet 1,500 with -C (child bounds
  from the cut certificate), dumping quiet.Y.bin.
- runs/valid: validate_fdrv.py (extra-region path, n=8 enumeration).
- runs/n10c: the cell disjunction at node 10's bottleneck cell 49:45 (x = 0.47): child 49:45,
  parent with 49:45 eliminated (-E), parent itself (-C 2 -X), 1,500 iterations each. One level of
  Sherali-Adams/N+ on the variable x_(49,45) has value min(child, eliminated); every disjunction on
  row 49 or column 45 (incl. the full line closures) is bounded by the child.
- next: count-split disjunction on node 13 on top of the triangle state (scripts/pick_split.py,
  count_term.py: both sides by Lagrangian terms at a few multipliers from the quiet state);
  runs/d2 phase + quiet for a second node's triangle-phase iterate.
- Star inequalities (scripts/star_check.py, logs/star_n13_uncut1500.log): hub j with weights
  B_j.: DNN value minus the conditional GL bound sum_a x_aj LAP_a is 110-264 per hub (minus
  min_a LAP_a: 144+). Implied (proof in the report) and far from binding.

## 06:13 n13 phase done; give-up iterate diagnostics
- runs/n13/phase.log: 8 rounds, best round 7 48,700.75 (+103.6 over round 0 48,597.17); laptop
  reference dual/tri_13.log r7 48,699.80 (+103.0). cpu 1,614 s. State saved (tri.*, 6,199 cuts).
- Give-up iterate (round 0, 400 it; logs/*_n13_giveup400.log): same picture as converged uncut
  (facility regions max 0.008, triple min 4.81, checker none, dispersion slack 2.8-12.3), BUT the
  iterate is not PSD (min eig of Y -0.0105; regional covariances down to -0.44 vs -0.0003 at 1,500 it),
  so large-support box "violations" (up to 12 box units at support 43-46, 2.9 at support 15) are
  primal infeasibility, not cuts: the count inequalities of near-full support are the cardinality
  equality sum_j z_j = |S| that the facial reduction already enforces. Diagnostics are meaningful
  at converged iterates only; box_sep.py now caps the support at 15.

## 06:18 Draft conclusions so far (to be confirmed at the triangle-phase iterate)
1. No facility structure: B is a random 0-9 matrix (spectrum, dense subsets, heavy triangles all
   at random-matrix values). The only facility-side structure is first order (flow degree 168-267,
   corr -0.90 with the DNN's expected eccentricity), which the DNN already prices exactly (stars).
2. The objective is <B, Delta> + linear, Delta_jl = sum over the 13 thresholds of P(j, l split by
   the threshold). The regional triangles ARE bisection inequalities of the flow graph B (one
   cardinality-constrained bisection per threshold), i.e. the facility side already enters every
   family tried, through B's weights. Facility regions T make the objective linear in their
   moments only through a cut decomposition of B + c(J - I) (O(n^2) regions, real weights).
3. At converged iterates the DNN is consistent on every facility-side count/triangle searched
   (QUBO-optimized facility regions max 0.006-0.008, both-sides-free boxes only on thresholds,
   parity/dispersion/star slack): it cheats only on the threshold memberships.
4. Disjunctions: one level of lift-and-project on a line = min over its children; bounded by the
   bottleneck child (the largest-x cell); N+ on a single variable = min(child, anti-child);
   dominated by branching. Measurements pending (node 10 cell 49:45; count split on n13).
- Facility-triple consistency LMP3_fac (scripts/lmp3.py: LP over the 97k injective placements of
  a facility triple, min L1 deviation of the three DNN pair blocks; logs/lmp3_n13_uncut1500.log),
  uncut converged node 13: heaviest-flow triples (flow sum 27) 0.0268, 0.0199, 0.0026; random
  triples 0.0023, 0.0002, 0.0013. The DNN's pair blocks ARE least consistent on the heaviest
  triples (about 10x random), but at the 1-3% mass level: each block has mass 1, so ~0.01 of
  probability would have to move. Order of magnitude per heavy triple in objective units:
  0.01 x (1-2 distance) x flows 9 x 3 pairs ~ 0.3-0.5.

## 06:30 n13 after the triangle phase (quiet round 1,500 it from round 7's multipliers)
- runs/n13/quiet.log: lb 48,745.63 (+148.5 over the give-up 48,597.17; laptop 12-round reference
  48,755.42), primal 48,749.34; 68.4 short of the cutoff.
- Certificate child bounds from the cut certificate (fdrv -C): 0 of 2,209 children fathomed;
  excess over lb median 0.054, 90% 0.22, 99% 19.4, max 43.5 (68.4 needed); 2,043 within 1 of lb.
  The cut certificate is as flat as the uncut one (tuning/README.md section 2).
- Diagnostics at the quiet iterate (logs/*_n13_tri.log): threshold triangles max 0.015-0.054
  (counts 45-438, down from 0.034-0.080); facility regions (QUBO over T, 120 location triples)
  max 0.0001; checker none; triple min 4.87; dispersion slack 2.9-11.0; hull3 only 3 metric facets
  (0.08-0.14, implied by threshold triangles). Signed boxes on thresholds with support 5-15:
  0.02-0.07 triangle units (the 5-point family was a wash in the pool). box2 found |P|=46,|Q|=43
  "violations" of 0.055: cancellation artifacts of near-full counts (E c^2 - (E c)^2 with
  E c ~ 42 and column-sum residuals 5e-4), not cuts.
- LMP3_fac at the quiet iterate: heaviest triples 0.023, 0.022, 0.011, 0.009 (uncut: 0.027, 0.020,
  0.003): the triangle phase does not repair heavy-triple consistency. The one facility-side
  family still violated -> pool test (scripts/lmp3_sep.py --k 60 on the quiet iterate, then
  fdrv -pool on top of the frozen triangle term, against an equal-iteration control).
- Count split picked at the quiet iterate (runs/n13/split_pick.log): col<2, Q+ {34,38},
  Q- 10 facilities, beta 8, E c 8.576, Var 0.104 (f(1-f) 0.244); jobs queued (mu 10/40/160 per side).

## 06:40 LMP3 on heavy facility triples and compact location triples (n13 quiet iterate)
- scripts/lmp3_sep.py --k 60 (logs/lmp3sep_n13_tri.log, runs/n13/lmp3_heavy60.cb): 60 cuts from the
  60 heaviest-flow facility triples (flow sum 27): violations min 0.0006, median 0.016, max 0.114,
  sum 1.46; the largest on triples inside the densest 4-set {5, 9, 20, 27} (0.065-0.114). Exact RHS
  (rounded potentials, max over all placements), perm check passes. ~6 s per LP.
- LMP3_loc (scripts/lmp3.py --loc, logs/lmp3loc_n13_tri.log): compact location triples (L-shapes,
  perimeter 4) 0.059, 0.053; random location triples 0.006, 0.0015. Cross-mandate observation (the
  location side): the 3-local consistency of compact location triples is violated ~10x more than
  random ones, in directions no facility-region triangle sees (facreg max 0.0001 at this iterate).
  scripts/lmp3_sep_loc.py writes those cuts (runs/n13/lmp3loc_compact60.cb, running).
- Queued: runs/n13/lmp3pool.log (fdrv -pool on the 60 heavy-triple cuts, 200 + 4 x 300 + quiet 800,
  delta 5, theta 0.5, warm Y) and runs/n13/lmp3ctl.log (equal-iteration control, 2,200 it).
- LMP3_loc cuts (runs/n13/lmp3loc_compact60.cb, logs/lmp3locsep_n13_tri.log): 60 L-shaped location
  triples, violations min 0.007, median 0.044, max 0.104, sum 2.76 (twice the heavy facility triples'
  1.46); perm check passes. Pool test queued (runs/n13/lmp3locpool.log).

## 06:48 Node 10's bottleneck cell 49:45 (the cell disjunction)
- runs/n10c/child4945.log (dnntest -x 39:27,30:48,25:49,49:45 -u 48816 -T 1 -d 3 -L 2 -r 40 -i 1500
  -H 0 -C 1): certified 48,690.99 at 1,500 it, primal 48,710.09 (m = 46, not fully converged; the
  child's DNN value is in [48,691, ~48,710]). Parent node 10 converged uncut 48,655.9 (misses table).
  So every disjunction on row 49 or on column 45 (the full line closures included; one level of
  lift-and-project on either line) is worth at most +35 certified / ~+54 by the primal estimate at
  node 10, which is 158 short. The child's certificate fathoms 0 of its 2,116 grandchildren.
- Lines study, all 30 lines of the 3 parents: the line's minimum child is its largest-x cell on 22
  of 30 lines, the second-largest on 7, the 5th on 1 (logs/lines_xweighted.log). The cap on any
  line disjunction is the child of the cell the DNN favours most.
- runs/n10c/elim4945.log (node 10 with x_(49,45) = 0 imposed, dnntest -E 49:45, 1,500 it):
  certified 48,672.21 (primal 48,670.03: converged). Anti-child +16.3 over the parent's 48,655.9.
  So one level of Sherali-Adams / N+ on the variable x_(49,45) (the hull of child and anti-child)
  is worth min(+35..+54, +16.3) = +16 at node 10 (158 short), for two DNN solves.
- The densest flow cluster {5, 9, 20, 27} (all six flows 9) at node 13: expected pairwise distance
  sum 13.9 (uncut) / 14.0 (after triangles) against the grid minimum 8; expected number of adjacent
  pairs 2.27 / 2.26 (at most 4 possible). Not placed too compactly; its LMP3 inconsistency
  (0.065-0.114 on its triples) is in the joint of the pair blocks, not in the distances.

## 07:18 Node 10 parent with the same settings; the cell disjunction in one table
- runs/n10c/parent.log (dnntest -x 39:27,30:48,25:49 -u 48816 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0
  -C 2 -X): certified 48,648.55, primal 48,651.14 (the -d 4 misses run gave 48,655.9).
  Same settings: anti-child (x_(49,45) = 0) 48,672.21 / 48,670.03 -> +21..+24;
  child 49:45 48,690.99 / 48,710.09 -> +42 certified, ~+59 by the primal.
  One level of N+ / Sherali-Adams on x_(49,45) (x = 0.47): +21..+24 for 2 DNN solves.
  Any disjunction on row 49 or column 45: <= +42..+59. Node 10 is 166 short (this run).
- Certificate child bounds at node 10 (logs/child_n10_parent.log): 0 of 2,209 fathomed; excess
  over lb median 0.16, 99% 0.48, max 7.6; 2,195 within 1. On the 446 lines-study children the
  certificate excess vs the child's own 300-it gain (median +117): Spearman 0.21; x vs gain -0.72.
- 07:22 split control with a Z-only restart from the converged quiet state (runs/n13/zonly/
  split_ctl_zonly.log, fdrv -cont quiet -iters 400): certified 48,713.37 against 48,745.63 at the end
  of the quiet round, primal 48,752.09: restarting Y from scratch costs ~32 certified units that 400
  iterations do not recover. All split jobs switched to -wy (Y and Z from the saved state) and rerun.

## 07:36 Count split on node 13 (outer region col<2, from the quiet state, Y+Z warm, 400 it)
- control (fdrv -cont quiet -wy -iters 400): 48,746.78 (primal 48,749.60), +1.15 over the quiet end.
- side c <= 8 (E c = 8.576): mu 10 -> 48,751.48 (+4.7), mu 40 -> 48,744.32 (primal 48,751.51; the
  certificate lags at the larger multiplier). The side's value is about +5: the DNN moves the
  fractional count mass 0.58 down for ~5 units. Dropped the mu = 160 runs (killed one; log in
  runs/n13/zonly/). Side c >= 9 running; then the mid-region split (col<4).
