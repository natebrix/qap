#!/usr/bin/env python3
"""summ.py: the T2 table. For each arm directory t2/PROBLEM_ARM: nodes (the node trace, qapctl export), the
CPU (results.txt), the rule's CPU (ls.log: the landscape signal and the lookahead; la.log for QAP_LA arms),
the CPU net of the rule, firings and gate passes by depth, and the lines chosen.
  summ.py [PROBLEM]   writes out/t2_PROBLEM.md"""
import sys, os, re, gzip, glob, collections
H = os.path.dirname(os.path.abspath(__file__))
P = sys.argv[1] if len(sys.argv) > 1 else 'tai35b'
rows = []
for d in sorted(glob.glob(f'{H}/t2/{P}_*')):
    arm = os.path.basename(d)[len(P) + 1:]
    res = open(f'{d}/results.txt').read() if os.path.exists(f'{d}/results.txt') else ''
    cpu = re.search(r'CPU ([0-9.]+) s', res)
    best = re.search(r'best objective ([0-9]+)', res)
    nodes, bydepth, opn = 0, collections.Counter(), 0
    nf = f'{d}/nodes.txt.gz'
    if os.path.exists(nf):
        for ln in gzip.open(nf, 'rt'):
            t = ln.split()
            if t and t[0] == 'node':
                nodes += 1
                if t[7] == 'B': bydepth[int(t[1])] += 1
    ls_cpu = la_cpu = 0.0
    fire, gate = collections.Counter(), collections.Counter()
    chosen = []
    if os.path.exists(f'{d}/ls.log'):
        for ln in open(f'{d}/ls.log'):
            dep = int(re.search(r' depth (\d+)', ln).group(1))
            fire[dep] += 1
            g = re.search(r'gate (\d)', ln)
            if g and g.group(1) == '1': gate[dep] += 1
            m = re.search(r'cpu_ls ([0-9.]+)', ln); ls_cpu += float(m.group(1)) if m else 0
            m = re.search(r'cpu_la ([0-9.]+)', ln); la_cpu += float(m.group(1)) if m else 0
            m = re.search(r'chose (\S+)', ln)
            if m: chosen.append(f'{dep}:{m.group(1)}')
    if os.path.exists(f'{d}/la.log'):
        for ln in open(f'{d}/la.log'):
            dep = int(re.search(r' depth (\d+)', ln).group(1))
            fire[dep] += 1; gate[dep] += 1
            m = re.search(r'([0-9.]+) CPU s;', ln); la_cpu += float(m.group(1)) if m else 0
            m = re.search(r'chose (\S+)', ln); chosen.append(f'{dep}:{m.group(1)}')
    c = float(cpu.group(1)) if cpu else float('nan')
    rule = ls_cpu + la_cpu
    rows.append((arm, nodes, c, ls_cpu, la_cpu, c - rule, best.group(1) if best else '-',
                 ' '.join(f'{k}:{gate[k]}/{fire[k]}' for k in sorted(fire)), ' '.join(chosen),
                 ' '.join(f'{k}:{bydepth[k]}' for k in sorted(bydepth) if k <= 8)))
out = [f'# T2 {P} (summ.py)', '',
       '| arm | nodes | CPU s | rule: LS CPU s | rule: lookahead CPU s | CPU net of the rule | best | gate passes / firings by depth | lines chosen (depth:line) | branched nodes by depth 0-8 |',
       '|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    out.append('| ' + ' | '.join([r[0], f'{r[1]:,}', f'{r[2]:,.0f}', f'{r[3]:,.0f}', f'{r[4]:,.0f}', f'{r[5]:,.0f}', r[6], r[7], r[8], r[9]]) + ' |')
os.makedirs(f'{H}/out', exist_ok=True)
open(f'{H}/out/t2_{P}.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
