# Experiment: face-restricted child bounds for lookahead

`face_bound.patch` (against 0efdd2f) tries idea 5 of the "bigger ideas" report:
bound each lookahead child with the parent's convex QP restricted to the
face x_ij = 1, instead of building the child's own QP. Not merged.

## Why it is a valid bound

For X = E_ij + Y (Y zero in row i and column j), the parent's QP
f(X) = tr(A X B X' - S X X' - X T X') + C.X + const
equals, as a polynomial in Y, the child's reduced problem (A_c, B_c, C_c as
reduceQAP/reduceC build them) with:
- S and T replaced by the parent's S and T with row/column i and j deleted
  (the cross terms of X X' and X'X vanish), and
- the constant lowered by S_ii + T_jj.

The face lies in the parent's affine hull, so the restricted QP is still
convex, and it equals the objective at permutations. No eigendecomposition
or S,T construction is needed per child.

## Results

Run with `QAP_FACE=2`, which computes both bounds for every lookahead child
and keeps the current one, so the search is unchanged:

| instance | children fixed out by bound, current | face | face - current, other children |
|---|---|---|---|
| nug15 | 10,959 | 7,577 | -10.8 |
| nug16a | 41,347 | 30,958 | -8.2 |
| had16 | 35,293 | 31,182 | -4.0 |
| chr12a | 1,711 | 1,183 | (much lower) |
| tai15a | 281,944 | 151,264 | (much lower) |

The face bound is also no cheaper: 1.69 s against 1.46 s of CPU on nug15.
The parent's S,T are a poor fit for the child, and the setup it saves is
small.

## Where lookahead time goes (gprof, nug16a, default.par)

- lookahead (select_ind_lookahead): 83% of the run
- Frank-Wolfe on the children (fw): almost all of that
- dlap (the LAP inside FW): 50% of the whole run
- init_ST, eigendecompositions and gradients together: under 10%

So a child bound only gets cheaper with fewer FW iterations. Cutting the
lookahead child budget ("Frank-Wolfe iterations, for row selection") from
50/25 to 25/12 or 12/6 keeps nodes about the same and saves 20-38% of the
time on nug16a, nug17, had16 and chr12a. On the suites (bench/LOG.csv; default, 25/12 and 12/6 run side by side):

| lookahead child FW iterations | standard nodes / time | large nodes / time |
|---|---|---|
| 50/25 (default.par) | 1 / 1 | 1 / 1 |
| 25/12 | 0.971 / 0.756 | 1.033 / 0.916 (every instance faster) |
| 12/6  | 0.958 / 0.674 | 1.137 / 0.925 (chr15a 1.07, tai15a 1.05 slower) |
