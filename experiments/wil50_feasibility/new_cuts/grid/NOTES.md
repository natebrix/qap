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
