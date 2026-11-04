"""Sparse cut files for fdrv -pool (.cb): int32 ncuts; per cut int64 key,
double h, int32 nnz, nnz x (int32 r, int32 c, double g): sum g Y_rc <= h on
the lifted matrix (each unordered pair (r, c) once, Y symmetric). Also: the
dense -Ladd term mu * sum G (each off-diagonal entry split in halves) and a
validity check on permutation matrices."""
import struct
import numpy as np

def write_cb(path, cuts):
    """cuts: list of (key, h, r, c, g) with r, c int arrays, g float array"""
    with open(path, 'wb') as f:
        f.write(struct.pack('<i', len(cuts)))
        for key, h, r, c, g in cuts:
            r = np.asarray(r, np.int32); c = np.asarray(c, np.int32); g = np.asarray(g, np.float64)
            lo = np.minimum(r, c); hi = np.maximum(r, c)
            f.write(struct.pack('<qdi', int(key), float(h), len(g)))
            rec = np.zeros(len(g), dtype=[('r', '<i4'), ('c', '<i4'), ('g', '<f8')])
            rec['r'] = lo; rec['c'] = hi; rec['g'] = g
            f.write(rec.tobytes())

def read_cb(path):
    data = open(path, 'rb').read(); off = 4
    nc, = struct.unpack_from('<i', data, 0); out = []
    for _ in range(nc):
        key, h, nnz = struct.unpack_from('<qdi', data, off); off += 20
        rec = np.frombuffer(data, dtype=[('r', '<i4'), ('c', '<i4'), ('g', '<f8')], count=nnz, offset=off); off += 16 * nnz
        out.append((key, h, rec['r'].astype(np.int64), rec['c'].astype(np.int64), rec['g'].copy()))
    return out

def value(cut, Y):
    key, h, r, c, g = cut
    return float((g * Y[r, c]).sum() - h)

def dense_term(cuts, mus, n1):
    L = np.zeros((n1, n1)); const = 0.0
    for (key, h, r, c, g), mu in zip(cuts, mus):
        d = r == c
        np.add.at(L, (r[d], c[d]), mu * g[d])
        np.add.at(L, (r[~d], c[~d]), mu * g[~d] / 2); np.add.at(L, (c[~d], r[~d]), mu * g[~d] / 2)
        const += mu * h
    return L, const

def perm_check(cuts, m, nperm=2000, seed=0, fixed_loc_free=None):
    """max over random permutations (location a -> facility pi[a], local) of
    sum g Y_rc - h with Y = xx' (x includes the constant 1): must be <= 0"""
    rng = np.random.default_rng(seed); worst = -np.inf
    for _ in range(nperm):
        pi = rng.permutation(m)
        xt = np.zeros(m * m + 1); xt[0] = 1; xt[1 + np.arange(m) * m + pi] = 1
        for key, h, r, c, g in cuts:
            worst = max(worst, float((g * xt[r] * xt[c]).sum() - h))
    return worst
