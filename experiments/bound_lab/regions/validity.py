#!/usr/bin/env python3
"""Validity of the cut bound with arbitrary regions (bound_lab/regions).

  python3 -I validity.py [-n NODES] [-seed S] [-free M]

For tai12b, tai15b and small synthetic tai-b-like instances (clustered points in the plane,
heavy-tailed asymmetric flows), draw random nodes with M free elements (default 8), find
each node's optimum by enumeration, and run dnncuts at the node with the cutoff far above it
(-c 2 x optimum), with random regions and with the structured families, on both sides, in the
CP form and the Polyak rounds form. Every certified value (FINAL best_lb, and every "cert"
of the CP trace) must be <= the node's optimum. Prints one line per run and a summary.
"""
import itertools, os, re, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DC = '/home/user/wt-regions/experiments/dnn_cuts/dnncuts'
Q = '/home/user/qap/qaplib'
sys.path.insert(0, HERE)
import gen_regions as G


def node_opt(A, B, fixes):
    n = len(A)
    p = -np.ones(n, int)
    for i, j in fixes:
        p[i] = j
    fi = [i for i in range(n) if p[i] < 0]
    fj = [j for j in range(n) if j not in set(p[p >= 0].tolist())]
    perms = np.array(list(itertools.permutations(fj)), dtype=np.int16)
    P = np.empty((len(perms), n), dtype=np.int16)
    P[:] = p
    P[:, fi] = perms
    obj = np.zeros(len(perms))
    for i in range(n):
        for k in range(n):
            if A[i, k] != 0:
                obj += A[i, k] * B[P[:, i], P[:, k]]
    return obj.min()


def synth(rng, n):
    k = rng.integers(2, 4)
    centers = rng.uniform(0, 1000, size=(k, 2))
    lab = rng.integers(0, k, size=n)
    X = centers[lab] + rng.normal(0, 30, size=(n, 2))
    D = np.rint(np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1)))
    F = np.rint(np.exp(rng.uniform(0, 8, size=(n, n))) * (rng.random((n, n)) < 0.6))
    np.fill_diagonal(F, 0)
    return D, F


def write_dat(path, A, B):
    n = len(A)
    with open(path, 'w') as f:
        f.write(f'{n}\n\n')
        for M in (A, B):
            for r in M:
                f.write(' '.join(str(int(x)) for x in r) + '\n')
            f.write('\n')


def run(dat, fixes, cutoff, regfile, side, form, delta):
    spec = ','.join(f'{i}:{j}' for i, j in fixes) or '-'
    cmd = [DC, '-x', spec, '-c', str(cutoff), '-regions', regfile, '-side', str(side), '-T', '1']
    if form == 'cp':
        cmd += ['-I', '200', '-k', '1000', '-cpsep', '50:1000', '-cpiters', '800']
    else:
        cmd += ['-I', '200', '-i', '200', '-r', '5', '-k', '1000', '-delta', str(delta)]
    out = subprocess.run(cmd + [dat], capture_output=True, text=True).stdout
    m = re.search(r'FINAL best_lb (\S+) first_lb (\S+)', out)
    if not m:
        return None, None, None
    certs = [float(x) for x in re.findall(r'^cp it +\d+ cert (\S+)', out, re.M)]
    certs += [float(x) for x in re.findall(r'^round +\d+ +lb (\S+)', out, re.M)]
    return float(m.group(1)), float(m.group(2)), max(certs + [float(m.group(1))])


def main():
    a = sys.argv[1:]
    nodes = int(a[a.index('-n') + 1]) if '-n' in a else 4
    seed = int(a[a.index('-seed') + 1]) if '-seed' in a else 1
    mfree = int(a[a.index('-free') + 1]) if '-free' in a else 8
    rng = np.random.default_rng(seed)
    tmp = tempfile.mkdtemp(prefix='regval_')
    insts = []
    for name in ('tai12b', 'tai15b'):
        n, A, B = G.read_qaplib(f'{Q}/{name}.dat')
        insts.append((name, f'{Q}/{name}.dat', A, B))
    for s in range(3):
        n = mfree + int(rng.integers(0, 3))
        D, F = synth(rng, n)
        path = f'{tmp}/syn{s}.dat'
        write_dat(path, D, F)
        insts.append((f'syn{s}(n={n})', path, D, F))
    runs = bad = gains = 0
    worst = -np.inf
    for name, dat, A, B in insts:
        n = len(A)
        for t in range(nodes):
            d = n - mfree
            fi = rng.choice(n, d, replace=False)
            fj = rng.choice(n, d, replace=False)
            fixes = sorted(zip(fi.tolist(), fj.tolist()))
            opt = node_opt(A, B, fixes)
            for fam, side in (('random:12', 0), ('random:12', 1), ('line:2', 0), ('hier', 0), ('ball:3', 0),
                              ('fhier', 1), ('fspec:1', 1), ('line:1,fhub:6', 0)):
                if fam == 'line:1,fhub:6':   # deliberately mismatched index sets: still valid (any subset is)
                    pass
                regs = []
                for spec in fam.split(','):
                    regs += G.family(spec, n, A, B, rng)
                regs = G.canon(regs, n)
                rf = f'{tmp}/r.txt'
                with open(rf, 'w') as f:
                    for r in regs:
                        f.write(' '.join(map(str, r)) + '\n')
                for form in ('cp', 'rounds'):
                    best, first, mx = run(dat, fixes, 2 * opt + 1000, rf, side, form, max(1, int(opt) // 100))
                    if best is None:
                        print(f'{name} node {fixes} {fam} side {side} {form}: NO RESULT')
                        bad += 1
                        continue
                    runs += 1
                    ok = mx <= opt + 1e-6 * max(1.0, abs(opt))
                    bad += not ok
                    gains += best > first + 1e-6 * abs(opt)
                    worst = max(worst, (mx - opt) / max(1.0, abs(opt)))
                    print(f'{name} m={mfree} node {t} {fam:14s} side {side} {form:6s} opt {opt:.0f} first {first:.1f} '
                          f'best {best:.1f} max cert {mx:.1f} gap left {opt - mx:.1f} {"ok" if ok else "VIOLATION"}', flush=True)
    print(f'SUMMARY runs {runs} violations {bad} runs with a gain over the uncut 200-iteration bound {gains} '
          f'largest (cert - opt)/opt {worst:.3e}')


if __name__ == '__main__':
    main()
