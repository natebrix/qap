# grid/ working notes (location side: wil50's 5 x 10 Manhattan grid)

Agent: grid research agent, Nov 4 2026, core 2 only. Append-only log.

## 05:27 Read

- latest/dnn.C header (certificate: any multiplier, Cholesky + GL polyhedral part),
  latest/dnncut.h (regional triangles, 5-point, local-7; L_ext hook),
  experiments/dnn_cuts/dnncuts.C (driver; -Ladd file[:const], -dumpY prefix).
- stronger_bounds/README.md steps 1-13, research_bound_quality.md (1.3, 3.4-3.6),
  research_memo_astra.md 1-4, roadmap/bound_strength.md, NOTES.md 1-5,
  sym/README.md, explainer section 3 ledger, triples/README.md, local7/README.md.
- Index convention (checked in local7_separate.py and sym/README.md): in a
  spec `i:j` on qaplib/wil50.dat, i is the LOCATION (first matrix A = grid,
  cell i = row i//10, col i%10) and j the FACILITY (second matrix B = flow).
  Lifted index of (local location a, local facility i) is 1 + a*m + i.
- Triangle-phase reference protocol (dual/run_five.sh):
  `dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500`
  laptop: node 13 48,755.42 (58.6 short), d2 48,789.87 (24.2), node 10 48,744 (70.1).

## 05:30-05:55 Machine, tools, first measurements (node 13 = 21:28,8:41,28:19)

- This box: Xeon 2.8 GHz, NO AMX flags in /proc/cpuinfo; dnntest at m = 47
  runs 0.55 s/iteration on core 2 (DNN_AMX=1 has nothing to use). A 1,500-it.
  solve is ~14 min; the full step-7 triangle protocol (5,500 it.) ~50 min, too
  long for the 30-min window, so the triangle iterates here are from shorter
  protocols (stated per run).
- Tools written: gridlib.py (node reduction, row/col split of <L,Y>),
  qapsa.c -> bin/qapsa (SA heuristic with linear term; heuristic UBs only),
  decomp_ub.py, dispersion.py (exact Dmin(k) by margins + max-flow, checked
  against brute force), analyze_Y.py, threshold_split.py, subqap_viol.py,
  region_search.py.
- runs/n13_uncut.log: dnntest -x 21:28,8:41,28:19 -u 48816 -T 1 -d 3 -L 2 -r 40
  -i 1500 -H 0 -v -dumpY runs/n13_uncut  -> lb 48,627.29, primal 48,628.38, 829 s.
- Decomposition (decomp_ub.py, SA 40 runs x 1000 sweeps):
  root: full 48,816 (= optimum; row 16,304 + col 32,512), row-only 15,766,
  col-only 32,418, sum 48,184 (DNN root 48,135).
  node 13: full 49,022 (row 16,248 + col 32,774), row-only 15,902, col-only
  32,682, sum 48,584 < uncut DNN 48,627: the decoupled bound is useless.
- Node 13 uncut iterate (analyze_Y.py): row part 15,999.6 (98 ABOVE the
  row-only optimum), col part 32,631.7 (50 below the col-only optimum).
  Parity-region triangles (checkerboard, row parity, 2x2 checker): ZERO
  violated; col parity max 0.023. Clique dispersion: never violated (most
  clustered 3-set sums 4.79 vs Dmin 4; 9-set 92.0 vs 72). Single interior
  rows/cols: max 0.017-0.027 (prefix mid regions 0.073-0.080).
- Flow-weighted sub-QAP inequalities (subqap_viol.py, greedy high-flow cliques
  k = 4..12, Q by SA = optimistic): slack by 140-1,700 units. Local packing
  is not where the DNN cheats.
- threshold_split.py: per threshold, the DNN iterate's cut term vs the SA
  min-bisection of that term alone: mid regions BELOW their bisection optimum
  (row<3 -12.4, col<5 -11.1, col<6 -16.7, col<4 +0.5); outer regions 50-118
  above. Sum of SA minima 48,048 (fully decoupled bound).
  Accounting at node 13: true coupling value (OPT_node - OPT_row - OPT_col)
  >= 438 (SA), the DNN captures 47 of it: the gap is the row-column coupling.

## 05:55-06:30 Coupled threshold subsets; the triangle iterate

- region_search.py (free-region 1-flip search for the most violated
  triangle, from each named region, top 20 triples): the prefix regions are
  local optima (col<5 0.0801 stays 0.0801; col<6 0.070 -> 0.079 by moving
  toward col<5); from the parity regions the search climbs to 0.019-0.058.
- xor_regions.py: corner rectangles (row<t & col<s) and quadrant
  checkerboards (row<t xor col<s): max 0.067, sums 10-30x below the mid
  prefix regions. Mixed regions do not carry the coupling either.
- coupled_split.py / root_vs_node.py (SA optima of the K-part of the
  objective alone, K = a set of thresholds, coupled through the joint cell
  capacities = a coarse "super-cell" QAP). Node 13, UNCUT iterate:
  K = mid5 = {row<2,row<3,col<4,col<5,col<6} (12 super-cells): iterate 24,035.6,
  coupled optimum 24,250 (robust: 2 seeds x 30 runs x 1000 sweeps), decoupled
  24,052. Deficit -214 (node), -76 against the ROOT coupled optimum 24,112.
  Pairs row x col: -12..-62 (mid ones). K = 7 mid: -282 / -150 (root).
  K = 9 (all rows + col<3..7): -303 / -151. The best known node solution
  (SA 49,022) has mid5 part 24,294: the mid5 coupling carries 258 of the
  391 gap at the uncut iterate.
- runs/n13_tri5.log: dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300
  -r 5 -delta 50 -nocap -freeze 0:600 -x 21:28,8:41,28:19 -dumpY runs/n13_tri5
  -> 48,713.85 certified (+117 over the 400-it. give-up point; the laptop's
  12-round protocol reached 48,755.4). 1,725 s CPU.
- At the triangle iterate (analyze, threshold_split, coupled_split,
  root_vs_node): every single threshold is now AT or ABOVE its own
  bisection optimum (min margin col<6 +0.3, row<3 +5.3); the row problem
  +125 above its optimum, the column problem +6. Parity still 0 violated,
  dispersion still slack. Coupled row x col sets still below: mid5 -117
  (node) / +21 (root), K=7 -181 / -49, K=9 -213 / -61.
  => after the triangles, all 1-D (single threshold, row-only, col-only)
  hulls are satisfied in the objective direction; what remains is the
  row-column coupling, and it is not local in thresholds (pairs -16..-25,
  5 thresholds -117, 9 thresholds -213) nor in facilities (k-sets slack).

## 06:20-06:40 Objective splitting (Lagrangian decomposition along a coarse coupled problem)

- Idea: for a threshold set K, <L_K, Y> >= OPT_K(node) is a valid linear
  cut whenever OPT_K(node) is a proven lower bound of the K-part; with
  multiplier alpha its Lagrangian term is -alpha L_K, i.e. the instance
  A' = s D - s alpha D_K, and bound = DNN(A')/s + alpha OPT_K. (make_inst.py ->
  inst/wil50_mid5_a{k}of{s}.dat, inst/wil50_{row,col,mid5only}.dat.)
  Superadditivity: the gain is at most alpha (OPT_K - DNN(L_K)): only the
  K-problem's OWN relaxation gap can be harvested, so OPT_K must come from
  something stronger than the DNN (exact solution).
- runs/n13_mid5a1of2.log: dnntest -x 21:28,8:41,28:19 -c 1e12 -T 1 -d 3 -L 2 -r 40
  -i 1500 -H 0 -v -dumpY runs/n13_mid5a1of2 inst/wil50_mid5_a1of2.dat
  -> lb' 73,194.81 (903 s). bound = 36,597.41 + 0.5 OPT_K:
  with the node's SA optimum 24,250: 48,722.4 = +95.1 over uncut 48,627.3;
  with the root SA optimum 24,112: 48,653.4 = +26.1. (SA values: potential
  only, not certified.) The iterate's mid5 part rises 24,035.6 -> 24,094.8,
  rest falls 24,595.6 -> 24,552.8. => the mid5 coarse problem's own DNN gap
  at node 13 is at least 2 x 95 = 190.
- Started runs/n13_tri5_mid5a1of2.log (same tri5 protocol on the alpha = 0.5
  instance, -delta 100 for the doubled scale, -c 1e12): on top of triangles.

## 06:40-07:10 On top of the triangles; validation

- runs/n13_tri5_mid5a1of2.log: dnncuts -grid 5:10 -c 1e12 -T 1 -k 1000 -v -I 400
  -i 300 -r 5 -delta 100 -nocap -freeze 0:600 -x 21:28,8:41,28:19
  -dumpY runs/n13_tri5_mid5a1of2 inst/wil50_mid5_a1of2.dat  (2,040 s CPU)
  -> lb' 73,254.29 (round 0 at 400 it.: 73,153.17). Wil50 bound = lb'/2 +
  0.5 LB_K: with LB_K = 24,250 (node SA optimum) 48,752.1 = +38.3 over the
  plain tri5 run (48,713.85, same protocol); with the root SA optimum 24,112
  it is 48,683.1 (-30.7). Break-even LB_K = 24,173.4. At the 400-it. give-up
  point: 48,701.6 vs 48,596.8 (+104.8). The triangles gain less on top
  (+101 scaled = +50.6 original, vs +117 plain): the two overlap.
- The coarse problem's own DNN at node 13 is <= 2 (48,627.3 - 36,597.4) =
  24,059.8 (superadditivity at the uncut level), so certifying LB_K >= 24,173
  (break-even) means closing >= 113 of its >= 190-unit gap, and +38 needs the
  coarse optimum itself.
- validate_objpart.py -N 12 --seed 3 (2x3, 2x4, 3x3 grids, random flows,
  random K, 0-2 fixes, alpha in {0.25, 0.5, 0.75}, exact OPT_K by
  enumeration): 48 certified values, 0 above the node optimum
  (runs/validate_objpart.log).
- Started runs/n13_tri5_mid5only.log: the coarse problem (A = D_mid5) at
  node 13 with the same tri5 protocol (-delta 25): its own DNN + triangle bound.
- Coarse (mid5) SA optima at node 13's ancestors: depth 2: 24,198 (21:28,8:41),
  24,192 (21:28,28:19), 24,160 (8:41,28:19); depth 1: 24,130-24,164; root
  24,112; the node 24,250. Against the break-even 24,173 on top of the
  triangles, a depth-2 key gives +0..+12.5, the node's own key +38: the cut
  must come from the node's own coarse problem (or its own coarse subtree).

## 07:10-07:35 The coarse problem's own bound (prices the certification)

- runs/n13_tri5_mid5only.log: dnncuts -grid 5:10 -c 1e12 -T 1 -k 1000 -v -I 400 -i 300
  -r 5 -delta 25 -nocap -freeze 0:600 -x 21:28,8:41,28:19 inst/wil50_mid5only.dat
  -> round 0 (plain DNN, 400 it.) 23,927.5; after 5 rounds + quiet 23,988.9
  (primal 24,007.7). Coarse SA optimum 24,250 (d2: 24,230, 5 seeds agree).
  The coarse problem's relaxation gap is >= 242 (1.0%), larger in relative
  terms than the fine node's (node 13: 48,714 tri5 vs 49,022 SA = 0.63%).
  Certifying the break-even LB_K = 24,173 needs a coarse B&B closing 184 of
  it; the +38 needs the full 261. 12-way branching instead of 47-way, keys
  shared by ~70-230 depth-3 nodes, but 3-5 coarse levels: not a cheap cut.
  Coarse level-1 experiment (coarse_bb.py) dropped: one level cannot close 184.
- Started runs/d2_tri5.log (d2 = 39:27,30:48,25:49,26:5, tri5 protocol, dump).
- Coarse (mid5) SA optima: node 10 = 24,230 (best of 5 seeds; 24,236 on 4),
  d2 = 24,230 (5 seeds agree): the extra fix 26:5 is coarse-consistent with
  an optimal coarse solution of its parent, so d2 gains no coarse strength.
- 07:3x: the Write tool refuses REPORT.md in this subagent session ("return
  findings as text"); the report goes into the hand-back message for the
  coordinator to save as REPORT.md.

## 07:31-08:05 d2

- runs/d2_tri5.log (tri5 protocol, d2): 48,749.62 certified (round 0 at 400 it.
  48,656.94; laptop 12-round reference 48,789.87). 1,700 s CPU.
- Its iterate's ORIGINAL objective <L,Y> is 48,872.7, 123 above the certified
  value and 95 above the Lagrangian primal (48,778.2): the quiet round's
  multipliers over-satisfy the cuts (round 4 overshot), so every part
  measured at this iterate is inflated by ~100; node 13's iterate did not
  have this (<L,Y> 48,715.2 vs 48,713.9). Readings at d2 (analyze_Y,
  threshold_split, coupled_split, root_vs_node): parity 0 violated,
  dispersion slack (most clustered 3-set 5.19 vs 4), all thresholds +25..+123
  above their bisection optima, row and column problems +187/+192 above,
  mid5 -3 (vs SA 24,242 in this run; 24,230 best), K=7 -49, K=9 -141.
- Started runs/d2_tri5_mid5a1of2.log (objective splitting at alpha 0.5 on
  top of the triangles; bound = lb'/2 + 0.5 x 24,230).
- runs/d2_tri5_mid5a1of2.log (tri5 protocol on inst/wil50_mid5_a1of2.dat, -c 1e12
  -delta 100, d2): lb' 73,333.45 -> 36,666.73 + 0.5 x 24,230 = 48,781.7 =
  **+32.1** over plain tri5 (48,749.62); at the 400-it. point +71.1
  (48,728.0 vs 48,656.9). Break-even LB_K at d2 = 24,165.8.
- 08:30 coarse level-1 B&B on node 13's mid5 problem (facility 40, 12
  children, 400 it. each, plain DNN) started: runs/n13_coarse_level1_f40.log.
- runs/n13_coarse_level1_f40.log (coarse_bb.py 21:28,8:41,28:19 40 --iters 400;
  12 children, 120-130 s CPU each, ~0.31 s/it): children 23,968.1 (super-cell
  row 2 x col 4, the single cell 24) .. 24,043.5; LB_K after one coarse level
  = 23,968.1, +40.6 over the parent's 400-it. plain value 23,927.5 and below
  its tri5 value 23,988.9. At this rate the break-even 24,173 is 4-6 coarse
  levels away. Converged check queued: runs/coarse_conv.sh (parent and the two
  lowest children at 1,500 it.).
- runs/coarseconv_*.log (dnntest -c 1e12 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0 on
  inst/wil50_mid5only.dat): node 13's coarse problem, plain DNN converged
  23,949.9 (primal 23,941.6; coarse gap to 24,250 = 300, 1.24%); its lowest
  coarse child (facility 40 -> cell 24) 23,989.8 = +39.9: the converged
  per-level rise of the coarse LB is ~40, the same as at 400 it. Break-even
  (24,173) is 223 above the plain coarse bound: ~5-6 coarse levels.
  => the certification of the objective-splitting cut is a coarse B&B of
  10^4-10^5 nodes per key: item 1 is killed by its own criterion (< 60/level).
- 09:15 queued node 10 tri5 (dump) then node 10 objective splitting alpha 0.5.
- second-lowest coarse child (facility 40 -> cell 34) converged 24,003.6
  (400 it.: 23,977.4). Convergence adds ~22-26 to parent and children alike;
  the min-child rise per coarse level stays ~40.

## 09:22-09:55 node 10

- runs/n10_tri5.log (tri5, node 10 = 39:27,30:48,25:49): 48,711.63 certified
  (round 0: 48,612.89; laptop 12-round reference 48,743.94). Iterate <L,Y>
  48,688.7 (cuts net violated, no inflation).
- At the tri5 iterate: parity 0 violated; dispersion slack (3-set 5.15 vs 4);
  every threshold at or above its bisection optimum (min col<4 +2.9); row
  problem +137, column problem +77 above their optima; coupled mid5 -155
  (vs 24,242 in this SA run; -143 vs the best 24,230), K=7 -190, K=9 -260;
  mid5 against the root optimum -25. The node-13 pattern replicates.
- runs/n10_tri5_mid5a1of2.log started (objective splitting, alpha 0.5).
- runs/n10_tri5_mid5a1of2.log: lb' 73,264.53 -> 36,632.26 + 0.5 x 24,230 =
  48,747.3 = **+35.6** over plain tri5 (48,711.63); 400-it. point +87.0.
  Break-even LB_K 24,158.7. On top of tri5, three test nodes: +38.3 (13),
  +35.6 (10), +32.1 (d2); median +35.6, conditional on proving the SA coarse
  optima (24,250 / 24,230 / 24,230).
- validate_objpart.py -N 8 --seed 11 --gap-only (only instances whose plain
  DNN is below the optimum; 2x3, 2x4, 3x3): 32 certified values, 0 above the
  node optimum (runs/validate_objpart_gap.log). With K = all thresholds
  (OPT_K = OPT) the split bound tends to OPT as alpha -> 1 (inst 51: 575.3 at
  alpha 0.75, optimum 576), as it must. Total: 80 values, 0 invalid.

## 10:30 Done. Conclusions (the hand-back message carries the full report)

1. No location-side inequality family beyond the prefix triangles is worth
   carrying: parity, dispersion, sub-QAP packing, free/rectangle/XOR regions,
   l1/hypermetric on E[d] are unviolated or dominated at 4 iterates.
2. After the triangles, every single threshold and both 1-D problems are at
   or above their own optima in the objective direction (nodes 13, 10; d2
   inflated); the remaining gap is the row x column coupling of the middle
   thresholds, global in thresholds and facilities.
3. Objective splitting along the 12-super-cell (mid5) problem: +32..+38 on
   top of tri5 on 13/10/d2 IF the coarse optimum is proven; proving it needs
   a coarse B&B ~5 levels deep (coarse relaxation gap 1.0-1.2%, min-child
   rise ~40/level): killed by its pre-declared criterion.
