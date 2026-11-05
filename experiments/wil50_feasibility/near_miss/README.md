# Near misses: cut phases that end a few units short (Nov 5)

From the node traces of five runs of the node-10 subtree (`../coord_search/`,
`../pool_search/`, `../pipeline_ab/runs/`): nine nodes ended their cut phase
within 25 units of the cutoff and branched. Every one of them ended at the
end of its 1,500-iteration quiet round (status iter_limit after 3-4 rounds
and the quiet round, 2,400 iterations in all), not on the stall rule.

| missed by | phases | their subtrees (DNN iterations) |
|---|---|---|
| < 2 | 2 (0.27, 0.68) | 10,700 and 7,150 |
| < 10 | 6 | 2,100-13,410, mean 6,300 |
| < 25 | 9 | mean 7,300 |

So a phase that ends within a few units, with its Lagrangian value above the
cutoff, might pay for a longer quiet round: another 1,000-1,500 iterations
against a subtree of about 6,000. Two cautions:

- The warmy arm (`../pipeline_ab/`, Y kept between rounds) attacks the same
  lag: with Y kept the certificate trails the Lagrangian value by 4-10 units
  at a round's end instead of 44-67 (`../roadmap/convergence2/REPORT.md`, F1).
  So the warmy run should be looked at first.
- `run_quiet.sh` tried to reproduce the 0.27 miss with `dnncuts` and failed
  as a test. Started cold, the node's uncut bound converged to 48,790.9 (23
  short) after 4,000 iterations. The cut rounds then lost 41-99 units, so
  the quiet round would restart from the uncut multipliers. In the search
  the node started from its parent's state, and that state is not
  available offline. The log is stopped in its quiet round. A faithful
  test needs the search's own phase with a longer quiet round: an
  environment switch in `bnb.C` and an arm.
