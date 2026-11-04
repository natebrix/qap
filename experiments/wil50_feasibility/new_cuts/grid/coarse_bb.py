#!/usr/bin/env python3
"""One level of branch-and-bound on a node's COARSE problem (the K-part of the
objective alone, K = a set of thresholds; inst/wil50_mid5only.dat for
K = mid5), to measure how fast a certified lower bound LB_K on the coarse
optimum rises. Coarse children of a facility f: one child per super-cell s
with a free cell, represented by fixing f to ONE free cell of s (all free
cells of a super-cell are interchangeable for the K-part: the rows of D_K
are identical within a super-cell, so the child's optimum and its DNN bound
do not depend on the representative). LB_K = min over children of their
certified bounds is a valid lower bound of OPT_K(node).

  python3 coarse_bb.py SPEC FACILITY [--iters 600] [--inst inst/wil50_mid5only.dat]
Writes runs/coarse_<tag>_f<F>_s<k>.log per child (resumable: RESULT line = done).
"""
import argparse, os, re, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl

HERE = os.path.dirname(os.path.abspath(__file__))
ROWB = lambda r: 0 if r < 2 else (1 if r < 3 else 2)            # mid5 row bands: rows 0-1, 2, 3-4
COLB = lambda c: 0 if c < 4 else (1 if c < 5 else (2 if c < 6 else 3))   # cols 0-3, 4, 5, 6-9


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('fac', type=int)
    ap.add_argument('--iters', type=int, default=600); ap.add_argument('--cutoff', type=float, default=1e12); ap.add_argument('--inst', default=HERE + '/inst/wil50_mid5only.dat')
    a = ap.parse_args()
    fixes = gl.parse_spec(a.spec)
    used = {x for x, _ in fixes}
    assert a.fac not in {y for _, y in fixes}
    reps = {}
    for c in range(50):
        if c in used:
            continue
        s = (ROWB(c // 10), COLB(c % 10))
        reps.setdefault(s, c)
    tag = a.spec.replace(':', '_').replace(',', '-')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', DNN_AMX='1', DNN_CHOLQR_PASSES='1')
    res = {}
    for s, c in sorted(reps.items()):
        log = f'{HERE}/runs/coarse_{tag}_f{a.fac}_s{s[0]}{s[1]}.log'
        if not (os.path.exists(log) and 'RESULT' in open(log).read()):
            cmd = ['taskset', '-c', '2', gl.ROOT + '/latest/dnntest', '-x', f'{a.spec},{c}:{a.fac}', '-c', str(a.cutoff), '-T', '1',
                   '-d', '3', '-L', '2', '-r', '40', '-i', str(a.iters), '-H', '0', a.inst]
            with open(log + '.tmp', 'w') as f:
                f.write(' '.join(cmd) + '\n'); f.flush()
                subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
            os.replace(log + '.tmp', log)
        m = re.search(r'RESULT .* lb=(\S+)', open(log).read())
        res[s] = float(m.group(1))
        print(f'super-cell {s} (rep cell {c}): lb {res[s]:.2f}', flush=True)
    lo = min(res.values())
    print(f'RESULT LB_K after one coarse level on facility {a.fac}: {lo:.2f} (children {len(res)})')


if __name__ == '__main__':
    main()
