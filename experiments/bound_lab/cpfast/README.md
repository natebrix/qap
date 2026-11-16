# cpfast: what limits the tail of the CP cut phase on tho40 (bound_lab agent `cpfast`, Nov 16 2026)

INTERIM (being written as results come in). [m] measured (log named), [e] estimated.

## Answer so far

- **The phase runs at a penalty that is about 4x too low.** `dnn.beta_sched = w=0.36` (tuned so that warm
  children fathom in 20-40 iterations) also sets the phase's beta, because `wc` is unset: beta = 0.36 x 0.21 x m/3.
  At that beta the primal residual is 16-20x the dual one for the whole phase (`rhat/shat` in the trace).
  With `dnn.beta_sched = w=0.36,wc=1.44` (phase beta x4, nothing else changed) on the depth-3 node
  `26:10,8:0,38:32` [m, `logs/A_base.log`, `logs/A_wc144.log`]: 50% of the long run's gain at 170 iterations
  instead of 340, 75% at 350 instead of 910, 90% at 620 instead of 1,830; 10.6 short of the cutoff at 1,000
  where today's settings need about 3,250. Same cost per iteration (0.12 s).
- Not the bottleneck: the pool's growth (saturated by iteration 300), the accuracy of the Dykstra projection
  (exact sweeps give about 1.3x on nug12/nug15, over-relaxing the half-space steps diverges on tho40).

(details below as they are measured)
