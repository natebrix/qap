# cpfast: what limits the tail of the CP cut phase on tho40 (bound_lab agent `cpfast`, Nov 16 2026)

DRAFT, being filled in as runs finish. [m] measured (log named), [e] estimated. Iterations are the
metric; the box was shared (load 5-8 on 4 cores).

## Answer (so far)

1. **The tail is the ADMM penalty of the phase, not the cuts.** tho40.par's `dnn.beta_sched = w=0.36`
   was tuned for warm children that fathom in 20-40 iterations; with no `wc` the cut phase runs at the
   same factor (beta = 0.36 x 0.21 x m/3). At that penalty the phase's primal residual is 15-19 times
   its dual residual from the first check to the last (`rhat/shat` in the trace), and the cut
   multipliers, which grow by (penalty x violation / ||G||^2) per iteration, creep: h'mu reaches 11,000
   after 2,000 iterations. With `wc=2.88` (8x the penalty) it is there after 300.
2. **`dnn.beta_sched = w=0.36,wc=2` to `2.88` makes the phase 2.4x to 4.6x shorter at the same cost per
   iteration** (table below; an existing key, no code).
3. The pool's growth, the Dykstra projection's accuracy, the eigen-tracker: not the bottleneck (below).

(tables and the rest follow)
