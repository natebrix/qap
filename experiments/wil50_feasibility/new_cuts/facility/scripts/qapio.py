import numpy as np
def read_qaplib(path):
    toks = open(path).read().split()
    n = int(toks[0]); v = np.array(toks[1:1 + 2 * n * n], dtype=float)
    A = v[:n * n].reshape(n, n); B = v[n * n:].reshape(n, n)
    return n, A, B
def parse_fixes(spec):
    if spec in ('-', ''): return []
    return [tuple(int(t) for t in s.split(':')) for s in spec.split(',')]
