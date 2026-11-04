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
