# Symmetry file for the 7x7 grid (dihedral group of order 8), in the solver's SYM format
# (latest/branch.C read_symmetry): J1 = the locations to branch on at the root (one per orbit
# of the group: the 10 cells with row <= col <= 3), then J2/J3 pairs: when every assigned
# location lies in J2 (the fixed cells of a reflection), the children need only cover J3 (a
# fundamental domain of that reflection). Cells are 1-based, row-major, 7 per row.
# The solver's symmetry acts on the columns (locations, the second matrix), and sko49's grid
# is its first matrix, so this file is for the transposed instance (sko49t.dat: B, A), which
# has the same optimal value (sum A_ij B_p(i)p(j) = sum B_kl A_q(k)q(l), q = p^-1).
import sys
R = C = 7
def cell(r, c): return r * C + c + 1
J1 = sorted(cell(r, c) for r in range(R) for c in range(C) if r <= c <= 3)
pairs = []
# diagonal reflection (r,c)->(c,r): fixed r==c, domain r<=c
pairs.append((sorted(cell(i, i) for i in range(7)), sorted(cell(r, c) for r in range(R) for c in range(C) if r <= c)))
# anti-diagonal reflection (r,c)->(6-c,6-r): fixed r+c==6, domain r+c<=6
pairs.append((sorted(cell(i, 6 - i) for i in range(7)), sorted(cell(r, c) for r in range(R) for c in range(C) if r + c <= 6)))
# horizontal mirror (r,c)->(6-r,c): fixed r==3, domain r<=3
pairs.append((sorted(cell(3, c) for c in range(7)), sorted(cell(r, c) for r in range(R) for c in range(C) if r <= 3)))
# vertical mirror (r,c)->(r,6-c): fixed c==3, domain c<=3
pairs.append((sorted(cell(r, 3) for r in range(7)), sorted(cell(r, c) for r in range(R) for c in range(C) if c <= 3)))
out = [str(len(J1)), " ".join(map(str, J1)), str(len(pairs))]
for J2, J3 in pairs:
    out += [f"{len(J2)} {len(J3)}", " ".join(map(str, J2)), " ".join(map(str, J3))]
print("\n".join(out))
