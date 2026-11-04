# Box (count) inequalities in the cut pool: killed

The literature agent's one remaining recommendation (`../literature/REPORT.md`,
recommendation 2; `../SYNTHESIS.md`): the Jünger-Kaibel 1-box inequalities
with large facility sets as an in-pool cut type, killed if the median gain
over the triangle phase on nodes 10, 13 and d2 is below +15.

## Implementation

`latest/dnncut.{h,C}`, cut type 8 (`dnncuts -box K[:QMAX[:STARTS]]`, search
`QAP_DNN_CUT_KBOX`): for a region's membership vector z and any facility set
Q and integer b >= 1, (q - b)(q - b - 1) >= 0 for the count q = sum_Q z_i,
which in the regional moments reads b sum_Q p_i - sum_{i<j in Q} W_ij <= b(b+1)/2
(the sum triangle is |Q| = 3, b = 1). Separated per region and round by
`box_gate.py`'s add/drop local search (b = 1..7, |Q| <= 24, 40 starts),
same Polyak step and certificate as the triangles; the pool's box terms are
summed per region before they are written into the objective.

Validity: `stronger_bounds/validate_cuts.py --families box` (n = 7 and 8, 36
and 24 nodes, `logs/validate_n{7,8}.log`): 108 certified values, none above
the enumerated optimum (the small cases have almost no DNN gap). nug12 root
with the cutoff at the optimum 578 (`logs/nug12_*.log`, grid 3:4 on the
distance side): triangles 573.58, triangles + boxes 576.98, boxes alone
576.48, never above 578. So on nug12 the family is valid and stronger than
the triangles, as the literature agent found with SCS.

## Protocol (predeclared in `run_box_pool.sh`)

Step 7's (`stronger_bounds/dual/run_five.sh`) with the search's per-node
settings (bundle D): from the give-up point (400 iterations), 12 Polyak rounds
of 300 iterations with 1,000 triangles a round, target best + 50 uncapped,
quiet round of 1,500 from the best round's multipliers; the box arm adds 200
box cuts a round. Each node's two arms ran side by side on the same host
(2.10 GHz Xeon with AMX). Gate: KILL if the median over 10, 13, d2 of (box -
tri) retained certificates is below +15.

## Result: KILL, decided by the first wave

| node | arm | round 0 | best Polyak round | quiet round = retained | Lagrangian primal, round 12 | process CPU s |
|---|---|---|---|---|---|---|
| 10 | triangles | 48,613.01 | 48,702.24 (round 12) | **48,749.32** | 48,763 | 2,881 |
| 10 | + boxes | 48,613.01 | 48,673.10 (round 5) | **48,717.03** (**-32.3**) | 48,763 | 2,945 |
| 13 | triangles | 48,597.17 | 48,708.56 (round 12) | **48,754.88** | 48,767 | 2,838 |
| 13 | + boxes | 48,597.17 | 48,679.07 (round 9) | **48,732.28** (**-22.6**) | 48,754 | 2,939 |

With two of the three gate nodes at -32.3 and -22.6, the median is at most
-22.6 whatever d2 gives, so the second wave (d2 and 14) was stopped after its
first round (`runs/*.stopped`) to free the cores for the pool-inheritance
measurement (`../../pool_search/`). `summary.txt` is `summarize.py`'s table.

## Reading

- **The boxes do not raise the Lagrangian value.** By round 12 the box
  arm's Lagrangian primal equals the triangle arm's on node 10 (48,763) and
  is 13 lower on node 13, with 2,374-2,381 box cuts added and box violations
  still up to 0.04-0.05 at round 12. The +17 the literature agent measured at
  a fixed multiplier on top of a converged triangle iterate does not appear
  when the boxes join the Polyak rounds from the give-up point.
- **They cost the certificate.** At the give-up point the boxes are violated
  five times more than the triangles (0.40 / 0.28 against 0.074), the first
  Polyak step is dominated by them, and round 1 falls 84-107 below round 0
  (the triangles gain 38-40). The certificate then trails the Lagrangian by
  83 and 66 per round on average over rounds 2-12 (up to 112) against 54 and
  51 for the triangles (the dense box term; the literature agent saw the same
  lag on d2), so the best certified round comes earlier (rounds 5 and 9
  against 12) and the quiet round starts from weaker multipliers: its
  converged value (primal 48,720 / 48,739) sits 34 and 20 below the triangle
  arm's (48,754 / 48,759).
- **Per-iteration cost +2-4%** (2,945 / 2,939 against 2,881 / 2,838 CPU-s).
- Like the 5-point family (step 7), the boxes are valid, not implied by the
  triangles, and strong on a small grid, but on wil50 at depths 3-4 the
  multiplier method cannot use them. With this the linear families on Y are
  closed, as `../SYNTHESIS.md` recommended after this test.
