#!/usr/bin/env python3
"""t3rank.py: at every node where line_signals' (c) arm (QAP_LA=4:30:0:2) fired on tai35b and tai40b,
run the C signal (t0, 30 kicks, no refinement) and set its line ranking beside the DNN lookahead's scores.
Slack s = cutoff - the node's best bound (qpb, dnn, cut) in that arm's node trace.
  t3rank.py   writes t3out/rank_*.txt and out/t3_rank.md"""
import os, re, gzip, subprocess
H = os.path.dirname(os.path.abspath(__file__))
LS = '/home/user/wt-cognac/experiments/line_signals/runs'
W = '/home/user/wt-cognac'
T0 = '/tmp/claude-0/t0'
fires = [l for l in open(f'{LS}/la.log') if ' 4:30:0:2 ' in l]
# the first six are tai35b's, the next six tai40b's (in the log's order; tai40b's root chose c27)
probs = ['tai35b'] * 6 + ['tai40b'] * 6
trace = {}
def nodes(p):
    if p not in trace:
        d = {}
        for ln in gzip.open(f'{LS}/{p}_la4p/nodes.txt.gz', 'rt'):
            t = ln.split()
            if t and t[0] == 'node':
                key = frozenset(t[2].split(',')) if t[2] != '-' else frozenset()
                v = [float(x) for x in t[3:6] if x != '-']
                if not v: continue
                b = max(v)
                d[key] = (b, float(t[6]), t[7], t[9])
        trace[p] = d
    return trace[p]
out = ['# T3: the signal\'s ranking against the DNN lookahead\'s at (c)\'s firings (t3rank.py)', '',
       'Signal: t0, 30 kicks per child, no refinement, seed 0. DNN: line_signals la.log, 30 pooled iterations per child,',
       'score = the sum of shortfalls (lower is better). Lines are those the lookahead judged (4 rows, 4 columns).', '',
       '| instance | depth | fixes | s | ex_node (gate) | signal CPU s | DNN choice | DNN order (score, M) | signal order of the same lines (survivors:off) | signal rank of DNN choice / DNN rank of signal top | signal top line overall (all lines) |',
       '|---|---|---|---|---|---|---|---|---|---|---|']
for i, (ln, p) in enumerate(zip(fires, probs)):
    m = re.search(r'depth (\d+) fixes ?(.*?) ?\| chose (\S+)', ln)
    dep, fx, chose = int(m.group(1)), m.group(2).split(), m.group(3)
    sc = re.findall(r' ([rc]\d+):([0-9.]+):([0-9.]+):(\d+)', ln.split(';', 1)[1])
    dnn = sorted(((float(s), name) for name, _, s, _ in sc))
    b, cut, outc, line = nodes(p)[frozenset(fx)]
    s = cut - b
    lines = ','.join(n for _, n in dnn)
    f = f'{H}/t3out/rank_{p}_{i}.txt'
    if not os.path.exists(f):
        r = subprocess.run([T0, f'{W}/qaplib/{p}.dat', f'{W}/QAPSOLN/{p}.sln', ','.join(fx), f'{s:.4f}', '30', '0', '2.5', '0', lines],
                           capture_output=True, text=True)
        open(f, 'w').write(r.stdout + r.stderr)
    txt = open(f).read()
    exn = float(re.search(r'ex_node (\S+)', txt).group(1))
    cpu = re.search(r'cpu gate ([0-9.]+) table ([0-9.]+)', txt)
    allL = re.search(r'by cnt then off: (.*)', txt).group(1).split()
    rank = {e.split(':')[0]: (k, e) for k, e in enumerate(allL)}
    sig = sorted((rank[n][0], n) for _, n in dnn)
    dnames = [n for _, n in dnn]
    out.append(f'| {p} | {dep} | {" ".join(fx) or "-"} | {s:,.0f} | {exn:,.0f} ({"pass" if exn <= s else "fail"}) | {float(cpu.group(1)) + float(cpu.group(2)):.1f} | {chose} | '
               + ' '.join(f'{n} {x / 1e6:.0f}' for x, n in dnn) + ' | '
               + ' '.join(f'{n} {rank[n][1].split(":")[2]}:{float(rank[n][1].split(":")[3]):.0f}' for _, n in sig) + ' | '
               + f'{[n for _, n in sig].index(chose) + 1} / {dnames.index(sig[0][1]) + 1} | {allL[0].split(":")[0]} ({allL[0].split(":")[2]}) |')
os.makedirs(f'{H}/out', exist_ok=True)
open(f'{H}/out/t3_rank.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
