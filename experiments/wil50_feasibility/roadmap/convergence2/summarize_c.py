#!/usr/bin/env python3
"""Test C: the Y12q arm (warm Y, 12 rounds of 100, quiet 500) against the
budget study's B / Y12 / Y12s arms (same_node_continuation/budget_study).

Iterations are compared, not CPU (different host and load). Round 0 (the
400 uncut iterations) must agree with the study's to the printed decimal,
which checks that the arithmetic is the same.
usage: summarize_c.py <repo root>
"""
import json, os, re, sys

root = sys.argv[1]
bs = os.path.join(root, 'experiments/same_node_continuation/budget_study/results/runs')
mine = os.path.join(root, 'experiments/wil50_feasibility/roadmap/convergence2/runs')
CUT = 48814.05


def parse(path):
    """rounds: list of (label, lb, primal, iters); final bound; total iterations; first crossing."""
    rounds, total, first = [], 0, None
    final = None
    for line in open(path):
        m = re.match(r'round\s+(\d+)\s+lb (\S+) .*?primal (\S+)\s+status (\S+)\s+iters (\d+)', line)
        q = re.match(r'quiet round\s+lb (\S+) .*?primal (\S+)\s+status (\S+)\s+iters (\d+)', line)
        if m:
            r, lb, pr, st, it = int(m.group(1)), float(m.group(2)), float(m.group(3)), m.group(4), int(m.group(5))
            total += it
            rounds.append((f'r{r}', lb, pr, it))
            if lb > CUT and first is None:
                first = total
        elif q:
            lb, pr, st, it = float(q.group(1)), float(q.group(2)), q.group(3), int(q.group(4))
            total += it
            rounds.append(('quiet', lb, pr, it))
            if lb > CUT and first is None:
                first = total
        f = re.match(r'FINAL best_lb (\S+)', line)
        if f:
            final = float(f.group(1))
    return rounds, final, total, first


summ = json.load(open(os.path.join(bs, 'summary.json')))
ref = {(r['node'], r['arm']): r for r in summ['rows']}
nodes = ['n14', 'n6', 'n9', 'n13', 'd2', 'n10']
print(f"{'node':5} {'arm':5} {'round0':>11} {'final':>11} {'vs B':>7} {'fathom@':>8} {'iters':>6}  per-round lb (lag = primal - lb)")
for n in nodes:
    p = os.path.join(mine, f'C_{n}_Y12q.log')
    if not os.path.exists(p):
        continue
    rounds, final, total, first = parse(p)
    if final is None:
        print(f'{n}: incomplete ({len(rounds)} rounds so far)')
        continue
    for arm in ['B', 'Y12', 'Y12s']:
        r = ref.get((n, arm))
        if r:
            br, _, _, _ = parse(os.path.join(bs, f'{n}_{arm}.log'))
            print(f"{n:5} {arm:5} {br[0][1]:11.4f} {r['final_bound']:11.4f} {r['final_bound'] - ref[(n, 'B')]['final_bound']:7.2f} "
                  f"{str(r['first_crossing_iteration'] or '-'):>8} {r['iterations']:6d}")
    rb = ref[(n, 'B')]['final_bound']
    lags = ' '.join(f"{lb:.0f}({pr - lb:+.0f})" for lab, lb, pr, it in rounds[1:])
    print(f"{n:5} {'Y12q':5} {rounds[0][1]:11.4f} {final:11.4f} {final - rb:7.2f} {str(first or '-'):>8} {total:6d}  {lags}")
    print()
