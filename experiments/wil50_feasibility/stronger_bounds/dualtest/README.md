# The dual-method test: is the Polyak protocol at the pool's dual optimum?

Memo candidate 1 (research_bound_quality.md section 2.1 and the section-4
protocol B). The regional triangles with 12 Polyak rounds and a quiet round
leave node 13 58.6 short and child d2 24.2 short (`../dual/tri_13.log`,
`../dual/tri_d2.log`, +158.6 and +133.0 over the give-up point). A bundle
method, or multipliers moved inside the ADMM, can only recover the distance
between those values and the dual optimum of the same cut pool. This test
measures that distance with the present code: a slow, careful subgradient
ascent (30 rounds, theta 0.3, target best + 30) with and without new cuts
after round 12 (`dnncuts -addto 12`, new option), and a frequent-update
variant (60 rounds of 100 iterations) as a proxy for in-ADMM multipliers.
`run_dualtest.sh` runs the six jobs two at a time, about 45-60 laptop
minutes each.

Decision rule, on the certified value after the quiet round against the
12-round reference:
- `slow30_pool12` gains 30 or more on either node: the multipliers are the
  lever; build the bundle step (its job is to reach that value in 5 rounds
  instead of 30), then consider in-ADMM updates.
- `slow30_pool12` under 10 and `slow30` under 20: the Polyak protocol is
  already near the pool's dual optimum and the regional triangles are capped
  where they stand; what remains is higher-order consistency (the regional
  moment block, candidate 6, or Astra's section-10b moment separation).
- `fast60` better than `slow30` per CPU second: frequent small updates pay,
  which is the case for moving the multipliers inside the ADMM.

## Results (laptop, Accelerate, two jobs at a time, Nov 2; logs `<proto>_<id>.log`)

Reference: the step-7 12-round runs `../dual/tri_13.log` (best Polyak round
48,708.97, quiet 48,755.42, 5,500 iterations, 1,067 CPU-s) and
`../dual/tri_d2.log` (48,730.62, 48,789.87, 5,500 it., 1,021 CPU-s). Cutoff
48,814.05. "Drops" = rounds whose certified value fell below the previous
round's; none fell by more than 20 in the slow protocols.

| node | protocol | best Polyak round | after quiet | gain vs 12-round ref | short | drops | iterations | CPU-s |
|---|---|---|---|---|---|---|---|---|
| 13 | slow30 | 48,715.86 (r30) | 48,763.05 | +7.6 | 51.0 | 3 | 9,400 + 1,500 | 1,650 |
| 13 | slow30_pool12 | 48,715.28 (r27) | 48,762.39 | +7.0 | 51.7 | 3 | 9,400 + 1,500 | 1,649 |
| 13 | fast60 | 48,635.50 (r60) | 48,753.88 | −1.5 | 60.2 | 7 | 6,400 + 1,500 | 1,215 |
| d2 | slow30 | 48,735.74 (r27) | 48,795.28 | +5.4 | 18.8 | 3 | 9,400 + 1,500 | 1,539 |
| d2 | slow30_pool12 | 48,736.04 (r30) | 48,795.75 | +5.9 | 18.3 | 2 | 9,400 + 1,500 | 1,462 |
| d2 | fast60 | 48,656.91 (r0) | 48,786.01 | −3.9 | 28.0 | 17 | 6,400 + 1,500 | 1,008 |

Reading against the decision rule: the second case. `slow30_pool12` gains 7
and 6, `slow30` 8 and 5, at 2.4x the reference's Polyak budget and 1.5x its
total CPU. The 12-round Polyak protocol plus the quiet round is within about
10 units of the pool's dual optimum on both nodes, and new cuts after round
12 are worth under 1 unit (slow30 vs slow30_pool12). The slow ascent is
still creeping at round 30 (+0.3 to +0.5 per round on node 13), so the dual
optimum of the pool is a little above 48,763, but not 50 above. A bundle
step would reach the same value in fewer rounds; it would not change the
value. `fast60` is worse on both nodes and per CPU second: 100-iteration
rounds leave the certificate 60-120 below the slow protocol's Polyak values
(the inner lag at 100 iterations is larger than at 300), the step sees a
noisier subgradient (7 and 17 drops), and the quiet round recovers most but
not all of it. Frequent small multiplier updates do not pay with this
certificate; multipliers moved inside the ADMM would need a different
measure of progress than the certified value to be worth building.

Conclusion: the regional triangles are capped where they stand, about +165
on node 13 (51 short) and +139 on d2 (18 short). Candidate 1 is closed
alongside 2 and 5. What remains of the memo is higher-order consistency: the
regional moment block (candidate 6) or Astra's moment-separation route, whose
first question is whether the separated cuts' violations at these iterates
are large enough to matter (the triple cuts' were not).
