# wil50 symmetry file (Nov 3, cloud)

Item 1 of `../roadmap/README.md`. The estimates in `../NOTES.md` section 5 and
the roadmap already assume the order-4 symmetry of wil50's 5 x 10 grid
(3.98x fewer nodes from depth 3 down), but no symmetry file existed and every
run so far has been on the plain tree. This directory builds the file and
checks it.

## What the file is

`make_sym.py` writes `latest/SYM/wil50t.sym` in the solver's J1/J2/J3 format
(`latest/branch.C`, `read_symmetry`; `nug30.sym` is the 5 x 6 version):

- J1 = rows 0-2 x columns 0-4 (15 cells): one cell per orbit of the group
  {identity, row reflection, column reflection, both}, the only children the
  root needs.
- one J2/J3 pair: J2 = row 2 (the cells fixed by the row reflection), J3 =
  rows 0-2. When every assigned cell lies in the middle row the remaining
  column reflection halves the children.

Burnside gives (P(50,d) + P(10,d))/4 orbits: 15, 635, 29,580, 1.38e6, 6.36e7
at depths 1-5 against 50, 2,450, 117,600, 5.53e6, 2.54e8 plain. The mechanism
realises these exactly at depths 1 and 2 (15; 10 x 49 + 5 x 29 = 635).

**The file is for the transposed instance.** The solver's symmetry acts on
the locations, the second matrix, and wil50's grid is its first matrix.
`make_sym.py --transpose qaplib/wil50.dat` writes `qaplib/wil50t.dat` (B then
A), which has the same optimal value (sum A_ij B_p(i)p(j) = sum B_kl A_q(k)q(l)
with q = p^-1). Consequences for every wil50t run:

- a `-F i:j` / `-x i:j` spec on wil50.dat becomes `j:i` on wil50t.dat;
- the regional cuts take `QAP_DNN_CUT_SIDE=1` in the search and `-side 1` in
  `dnncuts` (the grid is now the second matrix);
- at nodes whose assigned cells all lie in row 2 the branching is restricted
  to rows (`rows_only`), which the marginal-line rules already handle.

The solver finds the file by name: `./SYM/wil50t.sym` relative to `latest/`,
or `-s path`.

## Checks

| check | command (from `latest/`) | expected | result | log |
|---|---|---|---|---|
| depth-1 children with the file | `bnbdrv -p PARAM/nug_deep.par -d 1 -u 48816 ../qaplib/wil50t.dat` | 1 + 15 = 16 nodes | **16** (67 s) | `logs/d1_sym.txt` |
| depth-1 children without it | same with `-s /nonexistent.sym` | 1 + 50 = 51 | **51** (294 s) | `logs/d1_plain.txt` |
| depth-2 nodes with the file | same with `-d 2` | 1 + 15 + 635 = 651 | **651** (15 at depth 1, 635 at depth 2; 55 min) | `logs/d2_sym.txt` |
| node 10 uncut, 400 iterations | `dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -delta 50 -nocap -r 1 -i 300 -x 39:27,30:48,25:49 qaplib/wil50.dat` against `-side 1 -x 27:39,48:30,49:25 qaplib/wil50t.dat` | equal to rounding | 48,612.890 / 48,612.981 | `logs/n10_*.log` |
| node 10 after one regional-cut round | same runs | equal to rounding | **48,650.809 / 48,650.785**; 9,051 / 9,052 violated triangles found, same CPU (553 / 554 s) | `logs/n10_*.log` |

The node counts are exactly the Burnside orbit counts at depths 1 and 2, so
the file realises the full order-4 reduction. The bound checks say the transposition and the side switch are exact to the
PRSM's floating-point path: the same cut population is found on the same
grid, and the certified values agree to 0.02-0.09 at 400 and 700 iterations.

## What it buys and what it does not

It buys the 3.98x already in the estimates: 635 depth-2 subtrees instead of
2,450, and the same ratio at every level below. It changes nothing inside a
subtree such as node 10's (`-F 39:27,30:48,25:49`: cells in rows 3, 3 and 2,
so no J2 set applies), which is why the node-10 A/Bs did not need it and are
unaffected. The full run, the sampling plan, and any depth-2 or depth-3
sample must use `wil50t.dat` with the file; the plain `wil50.dat` tree is
3.98x larger for the same answer.
