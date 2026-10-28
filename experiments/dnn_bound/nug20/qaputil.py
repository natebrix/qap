"""Instance / state parsing and subproblem reduction for the QAP bound test."""
import numpy as np


def read_qaplib(fn):
    with open(fn) as f:
        vals = f.read().split()
    n = int(vals[0])
    a = np.array(vals[1:1 + n * n], dtype=float).reshape(n, n)
    b = np.array(vals[1 + n * n:1 + 2 * n * n], dtype=float).reshape(n, n)
    rest = vals[1 + 2 * n * n:]
    c = np.zeros((n, n))
    if len(rest) >= n * n:
        c = np.array(rest[:n * n], dtype=float).reshape(n, n)
    return n, a, b, c


def qap_obj(a, b, c, p):
    """p[i] = location of facility i."""
    p = np.asarray(p)
    return float(np.sum(a * b[np.ix_(p, p)]) + c[np.arange(len(p)), p].sum())


def read_state_nodes(fn):
    """Return list of dicts for each node line: depth, bound, pred, nfix, p (list), hist, warm."""
    nodes = []
    n = None
    with open(fn) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == 'n':
                n = int(t[1])
            elif t[0] == 'node':
                depth = int(t[1]); bound = float(t[2]); pred = float(t[3]); nfix = int(t[4])
                p = [int(x) for x in t[5:5 + n]]
                hist = [int(x) for x in t[5 + n:5 + 2 * n]]
                warm = int(t[5 + 2 * n])
                nodes.append(dict(depth=depth, bound=bound, pred=pred, nfix=nfix, p=p, hist=hist, warm=warm))
    return n, nodes


def parent_key(p, hist, depth):
    """Key identifying the depth-(depth-1) parent: the fixed pairs except the last fixed one.
    history[i] = order in which facility i was fixed (guess); fall back to using hist to drop the
    most recently fixed facility."""
    # hist[d] = facility branched on at depth d; the node at depth `depth` has
    # facilities hist[0..depth-1] fixed; its parent has hist[0..depth-2] fixed.
    fixed = [(i, p[i]) for i in range(len(p)) if p[i] >= 0]
    assert len(fixed) == depth, (fixed, depth)
    assert sorted(hist[:depth]) == sorted(i for i, _ in fixed)
    last = hist[depth - 1]
    return tuple(sorted((i, l) for i, l in fixed if i != last))


def reduce_subproblem(a, b, c, p):
    """Fixed assignment p (p[i]=loc or -1). Returns F, L, A_FF, B_LL, Cp (|F|x|L|), const such that
    full objective for any completion = tr(A_FF Y B_LL Y') + <Cp, Y> + const."""
    n = len(p)
    F = [i for i in range(n) if p[i] < 0]
    fixed = [i for i in range(n) if p[i] >= 0]
    L = [l for l in range(n) if l not in set(p[i] for i in fixed)]
    AFF = a[np.ix_(F, F)]
    BLL = b[np.ix_(L, L)]
    Cp = np.zeros((len(F), len(L)))
    for ii, i in enumerate(F):
        for ll, l in enumerate(L):
            s = c[i, l]
            for f in fixed:
                s += a[i, f] * b[l, p[f]] + a[f, i] * b[p[f], l]
            Cp[ii, ll] = s
    const = 0.0
    for f in fixed:
        const += c[f, p[f]]
        for g in fixed:
            const += a[f, g] * b[p[f], p[g]]
    return F, L, AFF, BLL, Cp, const


def verify_reduction(a, b, c, p, rng, trials=5):
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    for _ in range(trials):
        perm = rng.permutation(len(F))
        full = list(p)
        for ii, i in enumerate(F):
            full[i] = L[perm[ii]]
        v_full = qap_obj(a, b, c, full)
        v_red = qap_obj(AFF, BLL, Cp, perm) + const
        assert abs(v_full - v_red) < 1e-6, (v_full, v_red)
    return True
