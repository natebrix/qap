#!/usr/bin/env python3
"""subtree_summary.py DIR [DIR ...]: one complete subtree's size and cost from its
trace.txt (subtree.sh): DNN nodes, fathomed, iterations (first solves + cut phases),
DNN CPU, the root's first-solve and phase shortfall, nodes per depth."""
import collections, re, sys, os
def summary(d):
    rows = []
    for ln in open(os.path.join(d, "trace.txt")):
        m = re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb ([\d.\-]+) primal [\d.\-]+ status (\w+) iters (\d+) cpu ([\d.]+)', ln)
        if not m: continue
        ph = re.search(r'\| cut 1 best ([\d.\-]+) rounds \d+ iters (\d+) cpu ([\d.]+) status (\w+)', ln)
        cut = float(re.search(r'\| cutoff ([\d.\-]+)', ln).group(1))
        rows.append(dict(d=int(m.group(1)), nf=len(m.group(2).split()), lb=float(m.group(3)), st=m.group(4), it=int(m.group(5)), cpu=float(m.group(6)),
                         ph=(float(ph.group(1)), int(ph.group(2)), float(ph.group(3)), ph.group(4)) if ph else None, fath=' fathomed 1' in ln, cut=cut))
    if not rows: return f"{d}: no trace lines"
    root = min(rows, key=lambda r: r['nf'])
    fi = sum(r['it'] for r in rows); pi = sum(r['ph'][1] for r in rows if r['ph'])
    cpu = sum(r['cpu'] + (r['ph'][2] if r['ph'] else 0) for r in rows)
    byd = collections.Counter(r['d'] for r in rows)
    out = [f"{os.path.basename(d.rstrip('/'))}: {len(rows):,} DNN nodes, {sum(r['fath'] for r in rows):,} fathomed, "
           f"{sum(1 for r in rows if r['ph']):,} cut phases ({sum(1 for r in rows if r['ph'] and r['ph'][3] == 'fathomed')} fathomed), "
           f"{fi + pi:,} DNN iterations ({fi:,} + {pi:,}), DNN CPU {cpu:,.0f} s ({cpu / len(rows):.1f} s per node)",
           f"  root (depth {root['d']}): first solve {root['cut'] - root['lb']:.1f} short after {root['it']} iterations"
           + (f", phase to {root['cut'] - root['ph'][0]:.1f} short in {root['ph'][1]} ({root['ph'][3]})" if root['ph'] else ", no phase"),
           "  nodes by depth: " + ", ".join(f"{k}: {v:,}" for k, v in sorted(byd.items()))]
    return "\n".join(out)
for d in sys.argv[1:]:
    print(summary(d))
