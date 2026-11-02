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
