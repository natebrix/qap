"""lsq.py: Python side of ls.c (built on import into ls.so beside this file), the instances, best known
solutions and the node traces.

Conventions, checked in check_conventions() against QAPSOLN and the traces: the instance file holds A then B;
a permutation p is 0-based with p[i] = the location of facility i (QAPLIB's .sln, minus 1); the cost is
sum_{i,k} A[i,k] B[p[i],p[k]]; a fix 'r:v' in a node trace means p[r] = v (row r of the assignment matrix =
facility r = an index of A, column v = an index of B); line 'rN' branches on facility N (children N:v for every
free v), 'cN' on location N (children r:N for every free r)."""
import ctypes, gzip, os, subprocess
import numpy as np

H = os.path.dirname(os.path.abspath(__file__))
Q = H.rsplit('/', 3)[0]                      # /home/user/qap
WT = '/home/user/wt-cognac'
_so = H + '/ls.so'
if not os.path.exists(_so) or os.path.getmtime(_so) < os.path.getmtime(H + '/ls.c'):
    subprocess.check_call(['gcc', '-O3', '-march=native', '-shared', '-fPIC', '-o', _so, H + '/ls.c'])
_L = ctypes.CDLL(_so)
_dp = ctypes.POINTER(ctypes.c_double); _ip = ctypes.POINTER(ctypes.c_int); _lp = ctypes.POINTER(ctypes.c_long)
_L.qap_cost.restype = ctypes.c_double
_L.descent.restype = ctypes.c_double
_L.multistart.restype = ctypes.c_double
_L.qap_cost.argtypes = [ctypes.c_int, _dp, _dp, _ip]
_L.descent.argtypes = [ctypes.c_int, _dp, _dp, _ip, _ip, _lp]
_L.multistart.argtypes = [ctypes.c_int, _dp, _dp, _ip, _ip, ctypes.c_int, ctypes.c_uint, _lp]


def read_qaplib(path):
    t = open(path).read().split()
    n = int(t[0])
    v = np.array(t[1:1 + 2 * n * n], dtype=float)
    return n, v[:n * n].reshape(n, n).copy(), v[n * n:].reshape(n, n).copy()


class Inst:
    def __init__(self, name, dat=None, sln=None, perm=None, best=None):
        self.name = name
        self.n, self.A, self.B = read_qaplib(dat or f'{Q}/qaplib/{name}.dat')
        self.A = np.ascontiguousarray(self.A); self.B = np.ascontiguousarray(self.B)
        self._a = self.A.ctypes.data_as(_dp); self._b = self.B.ctypes.data_as(_dp)
        if perm is None:
            t = open(sln or f'{Q}/QAPSOLN/{name}.sln').read().split()
            best, perm = int(t[1]), [int(x) - 1 for x in t[2:2 + self.n]]
        self.p0 = np.array(perm, dtype=np.int32)
        self.best = best
        assert abs(self.cost(self.p0) - best) < 0.5, (name, self.cost(self.p0), best)
        self.moves = ctypes.c_long(0)

    def cost(self, p):
        p = np.ascontiguousarray(p, dtype=np.int32)
        return _L.qap_cost(self.n, self._a, self._b, p.ctypes.data_as(_ip))

    def adjust(self, p, fix):
        """p with the fixes forced, each by one exchange"""
        p = np.array(p, dtype=np.int32)
        for r, v in fix.items():
            if p[r] != v:
                j = int(np.where(p == v)[0][0]); p[j], p[r] = p[r], p[j]
        return p

    def free_mask(self, fix):
        f = np.ones(self.n, dtype=np.int32)
        for r in fix: f[r] = 0
        return f

    def descend(self, p, fix):
        """2-exchange descent from p (already satisfying fix); returns (cost, p)"""
        p = np.ascontiguousarray(p, dtype=np.int32).copy(); fr = self.free_mask(fix)
        c = _L.descent(self.n, self._a, self._b, p.ctypes.data_as(_ip), fr.ctypes.data_as(_ip), ctypes.byref(self.moves))
        return c, p

    def multistart(self, p, fix, nstart, seed=1):
        p = np.ascontiguousarray(p, dtype=np.int32).copy(); fr = self.free_mask(fix)
        c = _L.multistart(self.n, self._a, self._b, p.ctypes.data_as(_ip), fr.ctypes.data_as(_ip), nstart, seed,
                          ctypes.byref(self.moves))
        return c, p


def parse_fixes(s):
    if s == '-': return []
    return [tuple(map(int, x.split(':'))) for x in s.split(',')]


def read_trace(path):
    """the node trace as a list of dicts; a node traced twice (a resumed unit) keeps its last line"""
    op = gzip.open if path.endswith('.gz') else open
    nodes = {}
    for ln in op(path, 'rt'):
        f = ln.split()
        if len(f) < 12 or f[0] != 'node': continue
        fv = lambda s: None if s == '-' else float(s)
        nodes[f[2]] = dict(depth=int(f[1]), fixes=f[2], fx=parse_fixes(f[2]), qpb=fv(f[3]), dnn=fv(f[4]), cut=fv(f[5]),
                           cutoff=float(f[6]), out=f[7], rule=f[8], line=f[9],
                           nch=None if f[10] == '-' else int(f[10]), cpu=float(f[11]))
    return list(nodes.values())


def lb(nd):
    """the node's final bound: the cut phase's if it ran, else the first DNN solve's, else QPB"""
    for k in ('cut', 'dnn', 'qpb'):
        if nd[k] is not None: return nd[k]
    return None


def mirror(I):
    """tai60b's mirror: p0 with facilities 43 and 44 exchanged (43:45, 44:58; cost best + 2,786). Until Nov 19 18:00
    excess.py, truth.py and chain.py built it wrongly (the facilities AT locations 43 and 44 exchanged); see README."""
    m = I.p0.copy(); m[43], m[44] = I.p0[44], I.p0[43]
    assert abs(I.cost(m) - I.best - 2786) < 0.5
    return m


def model_cost(M, x):
    """model.py's predicted subtree of a child, relative: P(branch | x) max(x, xmin)^c"""
    import numpy as np
    x = np.asarray(x, dtype=float)
    p = 1 / (1 + np.exp(-(M['p0'] + M['p1'] * np.log(np.clip(x, 0.01, None)))))
    return p * np.clip(x, M['xmin'], None) ** M['c']
