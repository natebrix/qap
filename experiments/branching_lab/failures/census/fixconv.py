"""Which matrix the trace's rows (facilities) index, and what a fix row:col means.
Checks: (1) the best known objective under sum A[i][j] B[p i][p j] vs the transpose role;
(2) the trace's root path fixes against the .sln permutation (fix f:l with p[f] == l)."""
import sys
from common import *
for inst, run in [('tai60b','tai60b-261117'),('tai50b','tai50b-261116'),('tai35b','tai35b_base'),('tai40b','tai40b_base'),('tho40','tho40ab-32'),('nug30','nug30-trace-a')]:
    n,A,B = read_dat(inst); obj,p = read_sln(inst)
    q=[0]*n
    for i in range(n): q[p[i]]=i
    o1=sum(A[i][j]*B[p[i]][p[j]] for i in range(n) for j in range(n))
    o2=sum(B[i][j]*A[p[i]][p[j]] for i in range(n) for j in range(n))
    symA = all(A[i][j]==A[j][i] for i in range(n) for j in range(n)); symB=all(B[i][j]==B[j][i] for i in range(n) for j in range(n))
    # count trace nodes consistent with p (fix f:l, p[f]==l) and with q
    nodes,_=load(run)
    cp=sum(1 for x in nodes.values() if x.fixes and all(p[f]==l for f,l in x.fixes))
    cq=sum(1 for x in nodes.values() if x.fixes and all(q[f]==l for f,l in x.fixes))
    maxd_p=max([x.depth for x in nodes.values() if x.fixes and all(p[f]==l for f,l in x.fixes)]+[0])
    print(f'{inst}: sln obj {obj}; sum A_ij B_p(i)p(j) = {o1}; with roles swapped {o2}; A sym {symA} B sym {symB}; '
          f'trace nodes on path with f:l=p[f]: {cp} (deepest {maxd_p}); with the inverse: {cq}')
