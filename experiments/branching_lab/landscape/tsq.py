"""tsq.py: ts.c's robust taboo search for lsq.Inst objects."""
import ctypes, os, numpy as np, lsq
_L = ctypes.CDLL(lsq.H + '/ts.so')
_dp = ctypes.POINTER(ctypes.c_double); _ip = ctypes.POINTER(ctypes.c_int)
_L.rots.restype = ctypes.c_double
_L.rots.argtypes = [ctypes.c_int, _dp, _dp, _ip, _ip, ctypes.c_int, ctypes.c_uint]
def rots(I, p, fix, iters, seed=1):
    p = np.ascontiguousarray(p, dtype=np.int32).copy(); fr = I.free_mask(fix)
    c = _L.rots(I.n, I._a, I._b, p.ctypes.data_as(_ip), fr.ctypes.data_as(_ip), iters, seed)
    return c, p
