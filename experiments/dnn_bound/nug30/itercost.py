"""Profile one OWX-ADMM iteration at m = 28 (lifted 785, V 730 columns): eigh vs matmuls vs clip.
Also times the dual-bound evaluation.  Run alone on one core."""
import sys, time
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn')
from qaputil import *
from dnn import build_lifted, proj_psd, dual_bound
n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug30.dat')
for m in (29, 28, 27):
    p = [-1] * n
    for i in range(n - m):
        p[i] = i
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    Lh, free, Vh = build_lifted(AFF, BLL, Cp)
    N1, K = Vh.shape
    rng = np.random.default_rng(0)
    Y = rng.random((N1, N1)); Y = (Y + Y.T) / 2
    Z = rng.random((N1, N1)); Z = (Z + Z.T) / 2
    reps = 5
    t = time.process_time()
    for _ in range(reps):
        W = Vh.T @ (Y + Z / 3) @ Vh
    t_in = (time.process_time() - t) / reps
    t = time.process_time()
    for _ in range(reps):
        w, U = np.linalg.eigh(W)
    t_eigh = (time.process_time() - t) / reps
    t = time.process_time()
    for _ in range(reps):
        R = (U * np.maximum(w, 0)) @ U.T
    t_rec = (time.process_time() - t) / reps
    t = time.process_time()
    for _ in range(reps):
        VRV = Vh @ R @ Vh.T
    t_out = (time.process_time() - t) / reps
    t = time.process_time()
    for _ in range(reps):
        Yn = np.clip(VRV - (Lh + Z) / 3, 0, 1); Yn[~free] = 0; Yn = (Yn + Yn.T) / 2; Z2 = Z + 1.618 * 3 * (Yn - VRV); d = np.linalg.norm(Yn - Y)
    t_clip = (time.process_time() - t) / reps
    t = time.process_time()
    for _ in range(reps):
        lb, _ = dual_bound(Lh, free, Vh, Z)
    t_lb = (time.process_time() - t) / reps
    # how many eigenvalues of W are positive (the needed part of the spectrum)
    npos = int((w > 0).sum())
    print('m=%d lifted=%d K=%d  V^T(.)V %.3fs  eigh(K) %.3fs  U diag U^T %.3fs  V(.)V^T %.3fs  clip/Z/norm %.3fs  | iter total %.3fs | dual_bound %.3fs | pos eigs of random W %d' % (
        m, N1, K, t_in, t_eigh, t_rec, t_out, t_clip, t_in + t_eigh + t_rec + t_out + t_clip, t_lb, npos), flush=True)
