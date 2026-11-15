#!/usr/bin/env python3
"""tab.py OPT PLAIN_LOG LOG...: per log, regions, violations at the plain 300-iteration iterate, the certified
bound at chosen CP iterations, and the share of (OPT - plain bound at the same total iterations) closed."""
import re, sys
opt = float(sys.argv[1])
def trace(path, key):
    out = {}
    for line in open(path):
        m = re.match(r'cp it +(\d+) cert (\S+)', line) if key == 'cp' else re.match(r'it +(\d+) lb (\S+)', line)
        if m: out[int(m.group(1))] = float(m.group(2))
    return out
plain = trace(sys.argv[2], 'it')
its = [100, 300, 500, 1000, 2000]
print('| regions file | side | regions | violated at round 0 (max) | ' + ' | '.join(f'CP it {t}' for t in its) + ' | final | gap closed vs plain at equal iterations |')
print('|---|---|---|---|' + '---|' * (len(its) + 2))
print('| plain (no cuts) | | | | ' + ' | '.join(f'{plain[t]:,.0f}' if t in plain else '' for t in its) + f' | {plain[max(plain)]:,.0f} (it {max(plain)}) | |')
for path in sys.argv[3:]:
    txt = open(path).read()
    tr = trace(path, 'cp')
    m = re.search(r'-regions regions/\S+?\.(\S+)\.txt -side (\d)', txt)
    nr = re.search(r'(\d+) regions\n', txt)
    v = re.search(r'violated (\d+) \(max (\S+)\)', txt)
    best = {}; b = -1e300
    for t in sorted(tr):
        b = max(b, tr[t]); best[t] = b
    last = max(tr) if tr else 0
    share = (best[last] - plain[last]) / (opt - plain[last]) if last in plain else float('nan')
    print(f'| {m.group(1)} | {m.group(2)} | {nr.group(1)} | {v.group(1)} ({v.group(2)}) | ' +
          ' | '.join(f'{best[t]:,.0f}' if t in best else '' for t in its) + f' | {best[last]:,.0f} (it {last}) | {100 * share:.1f}% |')
