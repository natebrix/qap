"""tai35b's twin: the certified bound after the phase (runs/hi_mu*, cutoff raised so nothing stops at the real one)
against H(mu) = min over the completions found (ceiling_tai35b.out) of cost + mu Delta, an upper bound on the true
minimum of the term the bound bounds. H - bound is then an upper bound on the relaxation's gap at that mu.
usage: gaplaw.py > gaplaw_tai35b.out"""
import re, os
H = os.path.dirname(os.path.abspath(__file__)); CUT = 283315445.05
h = {float(m.group(1)): float(m.group(2).replace(',', '')) for m in re.finditer(r'mu\s+(\d+): min over S of cost \+ mu Delta - cutoff = ([+-][\d,]+)', open(f'{H}/ceiling_tai35b.out').read())}
print('| mu | bound after the phase - cutoff | H(mu) - cutoff | H - bound (>= the relaxation gap) | bound gain over mu=0 | H gain over mu=0 | realised |')
print('|---|---|---|---|---|---|---|')
b0 = None
for mu in [0, 8, 16, 32, 64]:
    t = open(f'{H}/runs/hi_mu{mu}/nodes.txt').read().split(); b = float(t[5]) - CUT
    if b0 is None: b0, h0 = b, h[mu]
    print(f'| {mu} | {b:+,.0f} | {h[mu]:+,.0f} | {h[mu] - b:,.0f} | {b - b0:+,.0f} | {h[mu] - h0:+,.0f} | {(b - b0) / (h[mu] - h0) if mu else float("nan"):.2f} |')
