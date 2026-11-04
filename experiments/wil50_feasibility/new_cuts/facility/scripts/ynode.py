"""Load a dumped DNN iterate of a wil50 node and the node's reduced data.
Conventions (latest/dnn.C, assign.C reduceQAP): the first index of the QAP is
the row/column of A (wil50: the 5 x 10 grid, 'locations' here) and the second
the row/column of B (the flows, 'facilities' here); a fix a:b puts facility b on
location a; local indices keep the original order; lifted index of x_(a,j)
(location a, facility j, local) is 1 + a*m + j."""
import numpy as np
from qapio import read_qaplib, parse_fixes
WIL = '/home/user/qap/qaplib/wil50.dat'

class Node:
    def __init__(self, spec, path=WIL):
        n, A, B = read_qaplib(path)
        A = (A + A.T) / 2; B = (B + B.T) / 2
        fx = parse_fixes(spec)
        self.n, self.fixes = n, fx
        fl = {a for a, b in fx}; ff = {b for a, b in fx}
        self.loc = [a for a in range(n) if a not in fl]     # free locations (original)
        self.fac = [b for b in range(n) if b not in ff]     # free facilities (original)
        self.m = m = len(self.loc)
        self.A, self.B = A[np.ix_(self.loc, self.loc)], B[np.ix_(self.fac, self.fac)]
        self.Afull, self.Bfull = A, B
        # linear term (assign.C): C_aj = sum_k fixed A_(a,k) B_(j,p(k)) + A_(k,a) B_(p(k),j)
        C = np.zeros((m, m))
        for k, pk in fx:
            C += np.outer(A[self.loc, k], B[self.fac, pk]) + np.outer(A[k, self.loc], B[pk, self.fac])
        self.C = C
        self.shift = sum(A[i, k] * B[pi, pk] for i, pi in fx for k, pk in fx)
        self.rc = np.array([[divmod(a, 10)[0], divmod(a, 10)[1]] for a in self.loc])  # grid coordinates

    def load(self, path):
        m = self.m; n1 = m * m + 1
        Y = np.fromfile(path, dtype=np.float64)
        assert Y.size == n1 * n1, (Y.size, n1)
        Y = Y.reshape(n1, n1)
        self.Y = Y
        self.x = Y[0, 1:].reshape(m, m)                   # [a, j]
        self.T = Y[1:, 1:].reshape(m, m, m, m)            # [a, j, b, l]
        return self

    def delta(self):
        """expected distance between facilities j, l: sum_ab A_ab Y_(aj),(bl)"""
        return np.einsum('ab,ajbl->jl', self.A, self.T, optimize=True)

    def objective(self):
        D = self.delta()
        return float((self.B * D).sum() + (self.C * self.x).sum() + self.shift)
