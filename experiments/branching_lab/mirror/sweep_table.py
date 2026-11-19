"""The mu sweep at tai35b's twin (runs/tw_*): per run the trace's QPB, first DNN and cut-phase values (of cost + mu Delta),
the iterations of each (DNN trace), Delta at the first solve's and the cut phase's final marginals (hook's stderr),
and the gain over mu = 0 against the node's shortfall. usage: sweep_table.py [PREFIX=tw]"""
import sys, glob, re, os
H = os.path.dirname(os.path.abspath(__file__)); pre = sys.argv[1] if len(sys.argv) > 1 else 'tw'
CUT = 283315445.05
rows = []
for d in glob.glob(f'{H}/runs/{pre}_*'):
    if not os.path.exists(f'{d}/nodes.txt'): continue
    nl = [l.split() for l in open(f'{d}/nodes.txt') if l.startswith('node')]
    t = nl[0]; qpb, dnn, cut = t[3], t[4], t[5]
    dl = open(f'{d}/dnn.txt').read().splitlines()[0]
    it1 = re.search(r'first lb \S+ primal (\S+) status (\S+) iters (\d+) cpu (\S+)', dl)
    itc = re.search(r'cut 1 best \S+ rounds \d+ iters (\d+) cpu (\S+) status (\S+)', dl)
    mu = 0.0; dx1 = dxc = '-'
    for l in open(f'{d}/err.txt'):
        if l.startswith('mirror: '):
            mu = float(re.search(r' mu (\S+)', l).group(1))
            m = re.search(r'Delta\(x first\) (\S+) Delta\(x cut\) (\S+) rowsum_err (\S+) scale (\S+)', l); dx1, dxc = m.group(1), m.group(2); sc = m.group(4)
    rows.append((mu, os.path.basename(d), float(qpb), float(dnn), float(cut) if cut != '-' else float('nan'),
                 it1.group(3), it1.group(2), itc.group(1) if itc else '-', itc.group(3) if itc else '-', t[7], dx1, dxc, it1.group(1)))
rows.sort(key=lambda r: (r[0], r[1] != f'{pre}_orig'))
base = [r for r in rows if r[0] == 0][0]
sf = CUT - base[4]
print(f'cutoff {CUT}; the mu = 0 node: first DNN {base[3]:.0f}, after the phase {base[4]:.0f}, shortfall {sf:.0f}; a quarter of it {sf/4:.0f}')
print('| run | mu | QPB | first DNN | its | after phase | phase its | phase status | fate | gain over mu=0 (phase) | gain / shortfall | Delta(x) first | Delta(x) phase end |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    g = r[4] - base[4]
    print(f'| {r[1]} | {r[0]:g} | {r[2]:,.0f} | {r[3]:,.0f} | {r[5]} | {r[4]:,.0f} | {r[7]} | {r[8]} | {r[9]} | {g:+,.0f} | {g/sf:+.3f} | {r[10]} | {r[11]} |')
