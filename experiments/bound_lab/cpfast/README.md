# cpfast: what limits the tail of the CP cut phase on tho40 (bound_lab agent `cpfast`, Nov 16 2026)

DRAFT 2 (01:40), being filled in as runs finish. [m] measured (log named), [e] estimated. Iterations
are the metric; the box was shared (load 5-8 on 4 cores).

## Answer (so far)

1. **The tail is the ADMM penalty of the phase.** tho40.par's `dnn.beta_sched = w=0.36` was tuned for
   warm children that fathom in 20-40 iterations; with no `wc` the cut phase runs at the same factor
   (beta = 0.36 x 0.21 x m/3). There the phase's primal residual is 15-19 times its dual residual from
   the first check to the last (`rhat/shat`), and the cut multipliers, which move by
   (penalty x violation / ||G||^2) per iteration, creep: h'mu is 4,473 at iteration 300, 9,058 at
   1,000, 11,021 at 2,000 [m `logs/A_base.log`]. At 8x the penalty it is 11,109 at iteration 300
   [m `logs/A_wc288.log`].
2. **With `dnn.beta_sched = w=0.36,wc=2.88` (an existing key) the phase is 2.4x to 4.6x shorter at the
   same cost per iteration** on off-path nodes of depths 3-5 [m]:
   - `26:10,8:0,38:32` (depth 3, 2,427 short), in the search, warm: fathoms after 801 phase iterations
     instead of 3,341 (`search/w38_wc2.88`, against `../tho40_search/runs/w38_z030`);
   - the warm subtree `26:10,8:0,18:32` end to end: **38 nodes and 5,891 DNN iterations (625 CPU s)
     against 74 nodes and 13,035 iterations (1,437 CPU s)**: 2.2x fewer iterations.
3. Not the bottleneck: the pool's growth, the accuracy of the Dykstra projection, (pending: the
   eigen-tracker).
4. **Caveats found**: (a) the change of penalty at the start of the phase costs its first 30-150
   iterations (the running best stays at the first solve's value), which is most of a near miss's
   phase; (b) on nodes on the path to the best known solution (depths 8-10, where the cut bound comes
   within 100 of the optimum) `wc=2.88` is slower than today's settings. Being measured: a rule
   that sets the phase's penalty from its own residual ratio (`wca`), and a rule that leaves near
   misses at `w` (`cut.cp_wc_gap`).

(tables follow)
