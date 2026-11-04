#!/usr/bin/env python3
"""Cold vs warm starts on the same 16 depth-4 children (stronger_bounds/depth4/nodes.txt).

cold: stronger_bounds/depth4/uncut_<k>.log (dnntest -v, cold, 1,500 iterations, no give-up)
warm: warm_y/runs/child_<k>.log (dnntest -K, parent given up at 300, then the child from
      Z only (z), Z + Y restricted (a), Z + Y conditioned (b); give-up after 300)

Prints, per child, the first check whose certified lb exceeds the cutoff (fathom iteration)
and the lb at iterations 50/100/200/300 in each run.
usage: cold_vs_warm.py <wil50_feasibility dir>
"""
import re, sys, os

base = sys.argv[1]
CUT = 48814.05


def blocks(path):
    """split a dnntest -v log into per-solve lists of (it, lb)."""
    out, cur = [], []
    for line in open(path):
        m = re.match(r'it\s+(\d+) lb (\S+)', line)
        if m:
            it = int(m.group(1))
            if cur and it <= cur[-1][0]:
                out.append(cur); cur = []
            cur.append((it, float(m.group(2))))
        elif line.startswith('RESULT') and cur:
            out.append(cur); cur = []
    if cur:
        out.append(cur)
    return out


def first_cross(b):
    for it, lb in b:
        if lb > CUT:
            return it
    return None


def at(b, k):
    v = [lb for it, lb in b if it <= k]
    return v[-1] if v else float('nan')


print(f"{'child':>5} {'cold':>6} {'z':>6} {'a':>6} {'b':>6} | lb@50 cold/z/b | lb@100 cold/z/b | lb@200 cold/z/b")
tot = {'cold': 0, 'z': 0, 'b': 0}
n = 0
for k in range(1, 17):
    cold = blocks(os.path.join(base, f'stronger_bounds/depth4/uncut_{k}.log'))
    warm = blocks(os.path.join(base, f'warm_y/runs/child_{k}.log'))
    c = cold[0]
    parent, z, a, b = warm[0], warm[1], warm[2], warm[3]
    fc, fz, fa, fb = first_cross(c), first_cross(z), first_cross(a), first_cross(b)
    def s(x):
        return f"{x:6d}" if x is not None else f"{'-':>6}"
    print(f"{k:5d} {s(fc)} {s(fz)} {s(fa)} {s(fb)} | {at(c,50):8.0f} {at(z,50):8.0f} {at(b,50):8.0f} | "
          f"{at(c,100):8.0f} {at(z,100):8.0f} {at(b,100):8.0f} | {at(c,200):8.0f} {at(z,200):8.0f} {at(b,200):8.0f}")
    if fc and fz and fb:
        tot['cold'] += fc; tot['z'] += fz; tot['b'] += fb; n += 1
print(f"\nfathomed in all of cold, z, b: {n} children; iterations cold {tot['cold']}, Z only {tot['z']}, Z+Y cond {tot['b']}")
