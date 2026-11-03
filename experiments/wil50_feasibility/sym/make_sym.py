# Symmetry file for wil50's 5 x 10 grid (group of order 4: identity, row
# reflection r -> 4-r, column reflection c -> 9-c, and their product), in the
# solver's SYM format (latest/branch.C, read_symmetry; nug30.sym is the 5 x 6
# version). J1 = one cell per orbit of the group at the root: rows 0-2 x
# columns 0-4 (15 cells). One J2/J3 pair: when every assigned cell lies in the
# middle row (J2 = row 2, the fixed cells of the row reflection), the children
# need only cover a fundamental domain of the column reflection within the
# remaining symmetry, J3 = rows 0-2 (30 cells). Burnside: (P(50,d) + P(10,d))/4
# orbits at depth d: 15, 635, 29,580, ...; the J1/J2/J3 mechanism realises
# exactly these at depths 1 and 2 (15; 10*49 + 5*29 = 635).
#
# The solver's symmetry acts on the locations (the second matrix) and wil50's
# grid is its FIRST matrix, so the file is for the transposed instance
# qaplib/wil50t.dat (B then A), which has the same optimal value (q = p^-1):
# a `-F i:j` spec on wil50.dat becomes `-F j:i` on wil50t.dat and the cut
# phase takes QAP_DNN_CUT_SIDE=1 (dnncuts -side 1).
#
#   python3 make_sym.py > ../../../latest/SYM/wil50t.sym
#   python3 make_sym.py --transpose ../../../qaplib/wil50.dat > ../../../qaplib/wil50t.dat
import sys
R, C = 5, 10
def cell(r, c): return r * C + c + 1          # 1-based, row-major

if len(sys.argv) > 1 and sys.argv[1] == "--transpose":
    t = open(sys.argv[2]).read().split(); n = int(t[0]); v = t[1:]
    A, B = v[:n * n], v[n * n:2 * n * n]
    print(n); print()
    for M in (B, A):
        for i in range(n): print(" ".join(M[i * n:(i + 1) * n]))
        print()
    sys.exit(0)

J1 = sorted(cell(r, c) for r in range(R) for c in range(C) if r <= 2 and c <= 4)
J2 = sorted(cell(2, c) for c in range(C))                    # the middle row
J3 = sorted(cell(r, c) for r in range(R) for c in range(C) if r <= 2)
out = [str(len(J1)), " ".join(map(str, J1)), "1",
       f"{len(J2)} {len(J3)}", " ".join(map(str, J2)), " ".join(map(str, J3))]
print("\n".join(out))
if "--count" in sys.argv:
    from math import perm
    for d in range(1, 6):
        print(f"# depth {d}: plain {perm(50, d)}, orbits {(perm(50, d) + perm(10, d)) // 4}", file=sys.stderr)
