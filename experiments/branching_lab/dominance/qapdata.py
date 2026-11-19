"""QAPLIB loading and cost. cost(p) = sum_{u,v} A[u,v] * B[p[u], p[v]], A the file's first matrix.
A trace fix 'r:c' means p[r] = c (checked in check_fixes.py against QAPSOLN)."""
import numpy as np, gzip
QAPLIB = '/home/user/wt-cognac/qaplib'
SOLN = '/home/user/wt-cognac/QAPSOLN'

def read_dat(name):
    path = name if '/' in name else f'{QAPLIB}/{name}.dat'
    tok = open(path).read().split()
    n = int(tok[0]); v = np.array(tok[1:1 + 2 * n * n], dtype=np.int64)
    return v[:n * n].reshape(n, n), v[n * n:].reshape(n, n)

def read_sln(name):
    tok = open(f'{SOLN}/{name}.sln').read().split()
    n, val = int(tok[0]), int(tok[1])
    return val, np.array([int(t) - 1 for t in tok[2:2 + n]])

def cost(A, B, p):
    p = np.asarray(p)
    return int((A * B[np.ix_(p, p)]).sum())

def read_trace(path):
    """yield (depth, fixes list of (r,c), outcome, line, cpu)"""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for ln in f:
            t = ln.split()
            if not t or t[0] != 'node':
                continue
            fx = [] if t[2] == '-' else [tuple(map(int, s.split(':'))) for s in t[2].split(',')]
            yield int(t[1]), fx, t[7], t[9], float(t[11]), ln
