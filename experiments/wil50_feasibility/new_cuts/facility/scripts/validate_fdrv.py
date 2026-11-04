#!/usr/bin/env python3
"""Validity check of fdrv's extra-region path (regions on either side, each
with its own side; dnncut_f.C) and of the phase/cont split: random n = 8
instances (tests/brute_force.py families), root and one-fix nodes, random
location and facility regions; every certified value printed (each round,
the continuation) must be <= the node optimum by enumeration. Cutoff = the
optimum, so a value above it would also print as a fathom.
  python3 validate_fdrv.py [-N 8] [--seed 3]"""
import argparse, importlib.util, os, random, re, subprocess, tempfile, itertools
ROOT = '/home/user/qap'
FDRV = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'fdrv')
spec = importlib.util.spec_from_file_location('brute_force', os.path.join(ROOT, 'tests', 'brute_force.py'))
bf = importlib.util.module_from_spec(spec); spec.loader.exec_module(bf)

def node_opt(A, B, C, fixes):
    n = len(A); best = None
    fx = dict(fixes)
    for p in itertools.permutations(range(n)):
        if any(p[i] != j for i, j in fx.items()): continue
        v = sum(A[i][k] * B[p[i]][p[k]] for i in range(n) for k in range(n)) + (sum(C[i][p[i]] for i in range(n)) if C else 0)
        best = v if best is None or v < best else best
    return best

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('-N', type=int, default=8); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('-n', type=int, default=8)
    a = ap.parse_args(); rng = random.Random(a.seed); tmp = tempfile.mkdtemp()
    checks = bad = 0; gains = []
    for t in range(a.N):
        fam = bf.FAMILIES[t % len(bf.FAMILIES)]
        A, B, C = bf.make_instance(fam, a.n, rng)
        dat = os.path.join(tmp, 'i%d.dat' % t); bf.write_dat(dat, A, B, C)
        reg = os.path.join(tmp, 'r%d.txt' % t)
        with open(reg, 'w') as f:
            for side in (0, 1):
                for k in range(4):
                    S = rng.sample(range(a.n), rng.randint(2, a.n - 2))
                    f.write('%d s%d_%d %s\n' % (side, side, k, ' '.join(map(str, S))))
        for fx in ('-', '%d:%d' % (rng.randrange(a.n), rng.randrange(a.n))):
            fixes = [] if fx == '-' else [tuple(map(int, fx.split(':')))]
            opt = node_opt(A, B, C, fixes)
            pre = os.path.join(tmp, 'p%d_%s' % (t, fx.replace(':', '_')))
            for nogrid in (['-nogrid'],):
                base = [FDRV, '-x', fx, '-c', str(opt), '-e', '0', '-g', '1.618', '-G', '0']
                out = subprocess.run(base + ['-regions', reg] + nogrid + ['-k', '30', '-I', '600', '-i', '300', '-rounds', '5', '-delta', '2', '-nocap', '-freeze', '400', '-phase', pre, dat], capture_output=True, text=True).stdout
                lbs = [float(x) for x in re.findall(r'round\s+\d+\s+lb (\S+)', out)]
                first = lbs[0] if lbs else None
                if os.path.exists(pre + '.meta'):
                    out2 = subprocess.run(base + ['-cont', pre, '-iters', '400', dat], capture_output=True, text=True).stdout
                    lbs += [float(x) for x in re.findall(r'RESULT lb=(\S+)', out2)]
                for lb in lbs:
                    checks += 1
                    if lb > opt + 1e-6: bad += 1; print('VIOLATION', dat, fx, lb, opt)
                if lbs: gains.append((max(lbs) - first, opt - first))
                print('%-8s node %-5s opt %6d  first %10.3f  best %10.3f  (%d values)' % (fam, fx, opt, first or float('nan'), max(lbs) if lbs else float('nan'), len(lbs)))
    print('%d certified values checked, %d above the optimum' % (checks, bad))
    g = [x / y for x, y in gains if y > 1e-6]
    print('nodes with a gap: %d, mean fraction closed %.1f%%' % (len(g), 100 * sum(g) / max(1, len(g))))

if __name__ == '__main__':
    main()
