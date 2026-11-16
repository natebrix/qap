**The cuts are valid and they do cut into the hull, but in the DNN relaxation they are worth
too little: closed as a material lever.** wil50 fails the gate outright; tho40 sits just under it.

| | measured [m] | the gate |
|---|---|---|
| wil50 node 13, on top of the triangle rounds, equal iterations | **+1.4** (scatter 5) | +10 |
| wil50 node 10, the same (3 rounds) | @@N10GAIN@@ | +10 |
| tho40 far-miss nodes d4_d, d5_d, on top of the triangle rounds, equal iterations | **+70 to +290 a round, median +135**; d4_d fathoms at 1,350 iterations where the triangles alone end 196 short | +150 |
| tho40 gate node, exchange cuts alone against the plain solve | +150 (the triangle phase: +1,327 and a fathom) | |

1. **Valid** [m]. The formula matches brute force exactly (asymmetric matrices, diagonals);
   2,248 certified values on n = 12 nodes never pass the enumerated optimum over
   exchange-optimal completions (or the optimum, on path nodes), and no value on the tho40 path
   to the best known solution passes 240,516 (stages 1, 3e). The best known solutions of tho40,
   wil50, sko42, sko49, tai50b and nug30 are exchange-optimal; wil50's has one exchange of gain
   0, i.e. a second permutation of value 48,816.
2. **At the level of the integer hull the idea is very strong** [m]: on n = 12 instances 71 to
   90% of off-path nodes at depths 3-5 lose their best completion and 29 to 77% have no
   exchange-optimal completion within 5% of the optimum (stage 2). The linear inequalities on Y
   keep most of that over the exact lifted hull: the same nodes are affected, 9 to 56% of them
   are excluded outright, and 55 to 70% of the finite rise survives (stage 2b).
3. **The DNN relaxation realises about a tenth of the violation it shows** [m]. 2 to 4% of the
   d m^2 cuts are violated at converged plain iterates, concentrated on one fixed facility
   (88 to 98% of the mass on facility 8 at corner 0 in the tho40 subtree `26:10,8:0`). The sum
   of violations is 810 to 1,413 cost units at tho40 depths 4-5 and the certified gain +110 to
   +217; it is 15 units at wil50 node 13 (gain +1.4) and 68 at node 10. wil50's dense, random
   flows leave almost nothing for an exchange to gain in expectation; tho40's sparse,
   heavy-tailed flows leave something, an order of magnitude less than the regional triangles
   give there (+1,300 to +1,900).
4. **What it would be worth** [e]. wil50 (and, unmeasured, sko42/sko49, whose flows are of the
   same kind): nothing. tho40: on top of a triangle phase, about one 150-iteration round in the
   climb and two or more at the plateau, at far-miss give-up nodes only; nothing at the near
   misses the phase fathoms in about 100 iterations. Of the order of 10% of the cut-phase
   iterations at far-miss nodes, less on the tree: not material. The gain is a level shift,
   so it needs multipliers that cost no iterations of their own (stepped inside the phase, as
   in the joint rounds here, or inherited).
5. **Not done**: stage 4 (children; the gate did not pass), sko49, the cuts as half-spaces
   inside the CP form's Y-step, a full-protocol (5,500-iteration) wil50 pair, the depth-4 node
   d4_e, the rotation and fixed-pair families on the large nodes (built and validated on small
   ones only).@@TAINOTE@@

**What I would do next**, in order: (a) nothing more on wil50/sko with this family; (b) if the
triangle phase goes into tho40's settings (`../tho40_nodes/`), look first at the side
observation of 3c (the rounds form with a target of best + 400 is ahead of the CP form at
equal iterations on the two far-miss nodes measured), which is worth more than these cuts;
(c) only after that, the exchange cuts as an add-on there: row half-spaces in the CP
projector plus inheritance (a child's cut is the parent's with the newly fixed facility's
term moved into c0), expected +135 a node at far-miss nodes; (d) the logical form at the
bottom of the tree, where the LAP elimination test starts to fire (16 of 64 assignments at
nug12 depth 4, 27 of 36 at depth 6) and nodes are cheap anyway.
