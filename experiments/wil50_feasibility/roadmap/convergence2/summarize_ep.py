#!/usr/bin/env python3
"""Tests E, P and S of run_tests.sh / run_tests2.sh.

E: node 10's equip phase (E_n10_phase.log) and its depth-4 children solved
   once with that pool frozen (E_c<f>.log), against the same children in the
   node-10 subtree run through the coordinator (coord_search/trace.txt: cold,
   no pool; first solve, then a cut phase when within the zone).
P: parents run 1,000 iterations (P_<parent>_J1000.log) against the parent at
   its give-up (300) in warm_y/runs; children start from the conditioned Y (b).
S: parent 13 at 300 with its state rounded to bf16 (S_13_J300_bf16.log)
   against warm_y/runs (unrounded).
usage: summarize_ep.py <repo root>
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trace_breakdown import parse as parse_trace

root = sys.argv[1]
W = os.path.join(root, 'experiments/wil50_feasibility')
RUNS = os.path.join(W, 'roadmap/convergence2/runs')
CUT = 48814.05


def checks(lines):
    return [(int(m.group(1)), float(m.group(2))) for m in
            (re.match(r'it\s+(\d+) lb (\S+)', l) for l in lines) if m]


def at(tr, k):
    v = [lb for it, lb in tr if it <= k]
    return v[-1] if v else float('nan')


# ---------------- E ----------------
print('== E: equip phase at node 10')
p = os.path.join(RUNS, 'E_n10_phase.log')
if os.path.exists(p):
    for l in open(p):
        if l.startswith(('round', 'quiet', 'stalled', 'FINAL', 'primal gate')):
            print('  ' + l.rstrip()[:170])
    pool = os.path.join(RUNS, 'E_n10.pool')
    if os.path.exists(pool):
        print(f'  pool file {os.path.getsize(pool)} bytes')
trace = {r['fixes']: r for r in parse_trace(os.path.join(W, 'coord_search/trace.txt'))}
tot_new = tot_old = 0.0
print(f"  {'child':>6} | {'search (cold, no pool)':48} | {'with node-10 pool (cold)':40}")
for f in [48, 49, 38, 29, 28, 19, 18, 9, 8, 0, 46, 37, 47]:
    lp = os.path.join(RUNS, f'E_c{f}.log')
    key = None
    for k in trace:
        parts = k.split()
        if len(parts) == 4 and f'{f}:45' in parts:
            key = k
    if key is None or not os.path.exists(lp):
        continue
    t = trace[key]
    lines = open(lp).read().splitlines()
    ps = [l for l in lines if l.startswith('pool solve')]
    if not ps:
        print(f'  {f:>4}:45 incomplete')
        continue
    m = re.search(r'lb (\S+)\s+primal (\S+)\s+status (\S+)\s+iters (\d+)\s+cpu (\S+)\s+kept (\d+)', ps[0])
    lb, st, it, cpu, kept = float(m.group(1)), m.group(3), int(m.group(4)), float(m.group(5)), int(m.group(6))
    old_it = t['it'] + t.get('cit', 0)
    old_cpu = t['cpu'] + t.get('ccpu', 0)
    old = f"{t['status'][:4]} {t['it']:4d} it {CUT - t['lb']:6.1f} short"
    if t['cut']:
        old += f" + phase {t['cit']:4d} it {t['cstatus'][:4]}"
    print(f"  {f:>4}:45 | {old:48} | {st[:4]} {it:4d} it {CUT - lb:6.1f} short, kept {kept}")
    tot_old += old_it
    tot_new += it
print(f'  iterations, search (first solve + phase) {tot_old:.0f}; first solve with the pool {tot_new:.0f}'
      ' (children still short need a phase on top)')

# ---------------- P and S ----------------
fam = {'13': [5, 6, 7, 8], '7': [13, 14, 15, 16], '10': [1, 2, 3, 4], '2': [9, 10, 11, 12]}


def child_blocks(path):
    """parent checks, then per child: (checks, result line) in order."""
    lines = open(path).read().splitlines()
    out, cur, parent, res = [], [], None, []
    for l in lines:
        m = re.match(r'it\s+(\d+) lb (\S+)', l)
        if m:
            cur.append((int(m.group(1)), float(m.group(2))))
        elif l.startswith('RESULT'):
            res.append(l)
            if parent is None:
                parent = cur
            else:
                out.append(cur)
            cur = []
    return parent, out, res


def res_it(l):
    return int(re.search(r'iters=(\d+)', l).group(1)), re.search(r'status=(\S+)', l).group(1)


for tag, fn, base_pass in [('P', 'P_{}_J1000.log', 3), ('S', 'S_{}_J300_bf16.log', 3)]:
    print(f'\n== {tag}: ' + ('parent 1,000 iterations vs parent 300 (warm_y)' if tag == 'P'
                            else 'parent state rounded to bf16 vs exact (warm_y)'))
    for par, kids in fam.items():
        lp = os.path.join(RUNS, fn.format(par))
        if not os.path.exists(lp):
            continue
        parent, mine, res = child_blocks(lp)
        if not res:
            continue
        print(f'  parent {par}: lb at 300 {at(parent, 300):.4f}, at end {parent[-1][1]:.4f} ({parent[-1][0]} it)')
        for k, (tr, rl) in enumerate(zip(mine, res[1:])):
            c = kids[k]
            wl = os.path.join(W, f'warm_y/runs/child_{c}.log')
            wparent, wkids, wres = child_blocks(wl)
            wtr, wr = wkids[2], wres[3]   # passes z, a, b: b is the third
            if k == 0:
                print(f'    (warm_y parent lb at 300: {at(wparent, 300):.4f})')
            i1, s1 = res_it(wr)
            i2, s2 = res_it(rl)
            print(f'    child {c:2d}: parent@300 {s1[:4]} {i1:4d} it  lb@50/100/200 {at(wtr,50):8.1f} {at(wtr,100):8.1f} {at(wtr,200):8.1f}'
                  f'  |  {tag} {s2[:4]} {i2:4d} it  {at(tr,50):8.1f} {at(tr,100):8.1f} {at(tr,200):8.1f}')
