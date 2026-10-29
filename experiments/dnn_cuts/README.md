# Triangle inequalities on top of the DNN bound (prototype, Oct 29 2026)

Question: can the certified DNN bound be made strong enough to fathom much
higher in the tree (for wil50, sko49) by adding valid inequalities?

`dnncuts.C` (`make -C latest dnncuts`) folds the Boolean-quadric triangle
inequalities on the lifted matrix Y into the objective with Lagrange
multipliers mu >= 0, for lifted variables a, b, c with distinct facilities
and locations:

    (T1)  Y_ab + Y_ac - Y_bc - Y_aa <= 0
    (T2)  Y_aa + Y_bb + Y_cc - Y_ab - Y_ac - Y_bc <= 1

<L,Y> >= <L + sum mu G, Y> - sum mu h for every permutation, so the certified
DNN bound of the modified objective (shift lowered by sum mu h) is a valid
bound; it goes through the new `DNNParams::L_add` hook of dnn.C. Each round
solves the relaxation (warm-started), separates the most violated triangles
in the last primal iterate, and takes a Polyak subgradient step.

Validity: 76 certified values on random n = 7, 8 instances with optima by
enumeration, none above the optimum.

## Results (nug30, `-I 2500 -i 1500 -k 3000`, logs here)

| node | gap to the cutoff (plain DNN, 2,500 iterations) | best gain in 10 rounds | largest violation at the plain DNN solution |
|---|---|---|---|
| depth 2 `0:13,6:14` (`c204.txt`) | 75.5 | +0.22 | 0.010 |
| depth 2 `0:2,23:6` (`c63.txt`) | 28.6 | +0.41 | 0.020 |
| root (`croot.txt`) | 174.6 | none | 0.004 |

nug12 root: +0.13 against a gap of 10. On the random n = 8 instances the DNN
bound is exact on 9 of 10; on the tenth the cuts close 19% of a gap of 39.

At a converged DNN solution the triangle inequalities are violated by at
most 0.01-0.02, and only a few hundred are violated: with the gangster
constraints and nonnegativity, the DNN relaxation nearly implies them. A
better multiplier method (bundle) would not change this by much. Triangle
inequalities are not the lever for fathoming higher in the tree; a stronger
relaxation family (level-2 RLT) would be the next thing to try.

Note: runs with only 300 ADMM iterations per round (an earlier attempt) look
like gains of 6-10 units, but those come from the extra ADMM iterations, not
the cuts: compare against the converged plain bound.
