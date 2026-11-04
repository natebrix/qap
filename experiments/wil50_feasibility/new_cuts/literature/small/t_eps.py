import sys, time
sys.path.insert(0, '/home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small')
from relax_small import *
n, A, B, C = read_qaplib(sys.argv[1]); L = lifted_cost(n, A, B, C)
for eps in map(float, sys.argv[2:]):
    val, st, t, yv = solve(n, L, ["dnn"], None, solver="SCS", eps=eps)
    print(eps, val, st, round(t, 1), flush=True)
