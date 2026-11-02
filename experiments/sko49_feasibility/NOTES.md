# Is sko49 reachable with the certified DNN bound? A numbers-based feasibility note (Nov 2 2026)

Companion to `../wil50_feasibility/NOTES.md`, written on a 4-vCPU cloud Xeon (OpenBLAS
`scipy-openblas32`, 1 thread per process, 4 processes at a time) in about 2.5 hours of
compute. Every bound below is a certified value printed by `dnntest` or `dnncuts` unless
it is marked *estimate*. Scripts and logs: `root/`, `depth3/`, `depth4/`, `cuts/`, `sym/`;
`summarize.py` builds the tables from the logs.

## 0. The instance

- `qaplib/sko49.dat` (already in the repo): n = 49, **the first matrix is the rectilinear
  distance on a 7 x 7 grid** (checked entry by entry: `d((r,c),(r',c')) = |r-r'| + |c-c'|`
  with row-major cells), the second is a symmetric random flow matrix (entries 0-10, 68%
  nonzero, zero diagonal). Same layout as wil50 (grid first), so the regional cuts take
  `-grid 7:7 -side 0`.
- **Best known 23,386**, confirmed from QAPLIB's `sko49.sln` (fetched from the Lehigh
  mirror; the permutation evaluates to 23,386 with the first matrix as the distance).
  The file is now in `QAPSOLN/sko49.sln`. Never proved optimal.
- `bnbdrv` prints `objective values are even integers: fathoming tolerance 1.95`, so the
  **cutoff is 23,384.05**. Root QPB bound 21,003.7 (gap 2,382, 10.2%).

## 5. Tree structure and symmetry (section numbers follow the task's measurement order)

Branching fixes one facility per level, so a depth-d node has 49 - d children (fewer where
the QPB reduced costs eliminate a location; `dnntest` runs on the plain node). The 7 x 7
grid distance is invariant under the dihedral group of order 8 (rotations and reflections),
against order 4 for wil50's 5 x 10 grid. Counting ordered d-tuples of distinct cells up to
symmetry (Burnside; `sym/make_sym.py` prints this):

| depth | plain nodes | inequivalent nodes | ratio | wil50 (order 4, estimate) |
|---|---|---|---|---|
| 1 | 49 | 10 | 4.9 | 15 |
| 2 | 2,352 | 315 | 7.5 | 635 |
| 3 | 110,544 | 13,923 | 7.9 | 29,580 |
| 4 | 5,085,024 | 636,048 | 8.0 | 1.4e6 |
| 5 | 2.29e8 | 2.86e7 | 8.0 | 6.4e7 |

So sko49's tree is **2.1 x narrower than wil50's at every level from depth 3**: the
symmetry (8 against 4) is worth 2 x, and n = 49 against 50 another 6%.

No `SYM/sko49.sym` existed. The format (`latest/branch.C`, `read_symmetry`) is J1 = the
locations to branch on at the root, then J2/J3 pairs: when every assigned location lies in
J2 (the fixed cells of a reflection), the children need only cover J3 (a fundamental domain
of that reflection); `SYM/nug30.sym` does exactly this for the 5 x 6 grid (J1 = the 3 x 3
corner block, J2 = the middle row, J3 = the top three rows). `sym/make_sym.py` writes the
7 x 7 version in a few lines (`sym/sko49.sym`: J1 = the 10 cells with row <= col <= 3; four
J2/J3 pairs, one per reflection, |J2| = 7, |J3| = 28). One caveat: **the solver's symmetry
acts on the columns (locations, the second matrix) and sko49's grid is its first matrix**,
so the file applies to the transposed instance `sym/sko49t.dat` (B then A; the optimal
value is unchanged, q = p^-1), on which the regional cuts would take `-side 1`. The
mechanism uses only the first matching J2 set at a node and gives nothing once a location
off every axis is assigned, so it realises the order-8 reduction at depth 1 (10 children)
and most of it at depth 2; below that the tree is the plain one under each inequivalent
depth-2 node, which is where the factor 8 of the table comes from.
