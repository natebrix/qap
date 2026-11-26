#!/usr/bin/env python3
"""Which size model fits the estimate's shares, and what each says of r45 against c0 (see predict.py)."""
import math, importlib.util, io, contextlib, os
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('p', os.path.join(here, 'predict.py')); m = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
L, est = m.L, m.est
def fit(w):
    pts = [(math.log(w(L[n][l])), math.log(sh)) for n, (l, sh) in est.items()]
    c = sum(y - x for x, y in pts) / len(pts)
    return c, math.sqrt(sum((y - x - c) ** 2 for x, y in pts) / len(pts))
def show(name, w):
    c, rms = fit(w)
    r = lambda n, a, b: w(L[n][a]) / w(L[n][b])
    print(f"{name:10s} rms x{math.exp(rms):.2f} | 21:10 r45/c0 {r('21-10','r45','c0'):.2f} | 21:1 r45/c0 {r('21-1','r45','c0'):.2f} | 21:0 best " +
          ' '.join(f"{l}:{w(x)/w(L['21-0']['r45']):.2f}" for l, x in sorted(L['21-0'].items(), key=lambda t: w(t[1]))[:3]) +
          " | 21:10 best " + ' '.join(f"{l}:{w(x)/w(L['21-10']['c0']):.2f}" for l, x in sorted(L['21-10'].items(), key=lambda t: w(t[1]))[:3]))
for k in (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1):
    show(f"exp k={k}", lambda x, k=k: sum(math.exp(k * s) for s in x if s > 0))
for P in (1, 2, 3, 4, 6, 8, 12, 16):
    show(f"power P={P}", lambda x, P=P: sum(s ** P for s in x))
