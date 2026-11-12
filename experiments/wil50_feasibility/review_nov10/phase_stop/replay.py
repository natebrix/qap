#!/usr/bin/env python3
"""Replay stop rules on CP cut phases from a trace with histories
(QAP_DNN_TRACE_HIST=1: "... | hist it:best it:best ..."; review_nov10 item 6).

    python3 replay.py TRACE [RULE ...]

A rule is one of
    stall:W:F      the production rule: from iteration W, stop when
                   cutoff - best > F * (best(it) - best(it - W))   (production: 600:3)
    hopeless:T:K:W the first solve's give-up rule from iteration T: with d1, d2 the
                   gains over the last two windows of W iterations and r = d1/d2
                   (at most 0.95), stop when best + K d1 r / (1 - r) < cutoff
                   (the first solve: T = dnn.giveup, K = 5, W = 50)
    floor:T        stop at iteration T whatever the history (a reference)
    need:T:M:W     from iteration T: stop when the gain per iteration over the last W
                   iterations, held constant, would need more than M further
                   iterations to reach the cutoff (M as a multiple of W: M=3 means
                   three more windows)
    gap:T=G,T=G    a gap schedule: from iteration T, stop while the gap to the cutoff is
                   still above G (the solver's cut.giveup, written T:G,T:G there)
Rules combine with "+" (stall:600:3+hopeless:150:5:50): the first to fire stops.

For every rule: the phases it stops before their natural end, the iterations
saved (of the phases' total), the fathoms it would lose (a fathoming phase
stopped before the iteration at which its best crossed the cutoff), and the
bound given up at the stop against the phase's final best. What it cannot see:
a stopped phase hands its children a weaker pool and branch bound, which the
A/B of runs_gu (the floor leaked into the phases) priced at 930 extra nodes
against 1.072x in DNN CPU.
"""
import re
import sys


def parse(path):
    phases = []
    for line in open(path):
        if '| hist' not in line:
            continue
        seg = line.rstrip('\n').split('|')
        head, first, cut, tail, hist = seg[0], seg[1], seg[2], seg[3], seg[4]
        m = re.match(r'pid \d+ depth (\d+) fixes (.*)', head.strip())
        depth, fixes = int(m.group(1)), m.group(2).strip()
        c = cut.split()
        kv = {c[i]: c[i + 1] for i in range(2, len(c) - 1, 2)}
        t = tail.split()
        tv = {t[i]: t[i + 1] for i in range(0, len(t) - 1, 2)}
        f = first.split()
        fv = {f[i]: f[i + 1] for i in range(1, len(f) - 1, 2)}
        h = [(int(a), float(b)) for a, b in (p.split(':') for p in hist.split()[1:])]
        phases.append(dict(depth=depth, fixes=fixes, first_lb=float(fv['lb']), cutoff=float(tv['cutoff']),
                           best=float(kv['best']), iters=int(kv['iters']), status=kv['status'],
                           cpu=float(kv['cpu']), hist=h))
    return phases


def best_at(h, first_lb, it):
    """the running best at iteration it (first_lb before the first check)"""
    b = first_lb
    for i, v in h:
        if i > it:
            break
        b = v
    return b


def rule_stop(rule, ph):
    """the iteration at which the rule fires on this phase's history, or None"""
    h, first_lb, cutoff = ph['hist'], ph['first_lb'], ph['cutoff']
    parts = rule.split('+')
    for idx, (it, b) in enumerate(h):
        if b > cutoff:
            return None   # fathomed first
        for r in parts:
            f = r.split(':')
            if f[0] == 'floor':
                if it >= int(f[1]):
                    return it
            elif f[0] == 'stall':
                W, F = int(f[1]), float(f[2])
                if it >= W:
                    b0 = best_at(h, first_lb, it - W)
                    if cutoff - b > F * max(0.0, b - b0):
                        return it
            elif f[0] == 'hopeless':
                T, K, W = int(f[1]), float(f[2]), int(f[3])
                if it >= T and it >= 2 * W:
                    d1 = b - best_at(h, first_lb, it - W)
                    d2 = best_at(h, first_lb, it - W) - best_at(h, first_lb, it - 2 * W)
                    r = min(d1 / d2, 0.95) if d2 > 0 else (0.95 if d1 > 0 else 0.0)
                    if b + K * d1 * r / (1 - r) < cutoff:
                        return it
            elif f[0] == 'gap':
                for tg in r[4:].split(','):
                    T, G = tg.split('=')
                    if it >= int(T) and cutoff - b > float(G):
                        return it
            elif f[0] == 'need':
                T, M, W = int(f[1]), float(f[2]), int(f[3])
                if it >= T and it >= W:
                    g = (b - best_at(h, first_lb, it - W)) / W
                    if g <= 0 or (cutoff - b) / g > M * W:
                        return it
            else:
                sys.exit(f'unknown rule {r}')
    return None


def fathom_it(ph):
    for it, b in ph['hist']:
        if b > ph['cutoff']:
            return it
    return ph['iters'] if ph['status'] == 'fathomed' else None


def evaluate(rule, phases, verbose=False):
    tot = sum(p['iters'] for p in phases)
    tot_cpu = sum(p['cpu'] for p in phases)
    stopped = saved = saved_cpu = 0
    lost = []
    give = []
    rows = []
    for p in phases:
        s = rule_stop(rule, p)
        if s is None or s >= p['iters']:
            continue
        fi = fathom_it(p)
        is_lost = fi is not None and s < fi
        stopped += 1
        saved += p['iters'] - s
        saved_cpu += p['cpu'] * (p['iters'] - s) / p['iters']
        g = p['best'] - best_at(p['hist'], p['first_lb'], s)
        give.append(g)
        if is_lost:
            lost.append((p, s, fi))
        rows.append((p, s, g, is_lost))
    print(f"{rule:34s} stops {stopped:3d} of {len(phases)}  saves {saved:6d} it ({100*saved/tot:4.1f}%), "
          f"{saved_cpu:7.0f} s ({100*saved_cpu/tot_cpu:4.1f}%)  fathoms lost {len(lost):2d}  "
          f"bound given up: mean {sum(give)/len(give) if give else 0:5.1f} max {max(give) if give else 0:5.1f}")
    if verbose:
        for p, s, g, is_lost in sorted(rows, key=lambda r: -r[0]['iters']):
            print(f"    depth {p['depth']} {p['fixes'][:40]:40s} {p['status']:8s} {p['iters']:5d} it -> stop {s:5d}"
                  f"  gives up {g:6.1f}  {'FATHOM LOST' if is_lost else ''}")
    for p, s, fi in lost:
        print(f"    lost: depth {p['depth']} {p['fixes']} fathomed at {fi}, rule at {s}, "
              f"best there {best_at(p['hist'], p['first_lb'], s):.1f} of cutoff {p['cutoff']:.1f}")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    phases = parse(sys.argv[1])
    verbose = '-v' in sys.argv
    rules = [a for a in sys.argv[2:] if a != '-v']
    fath = [p for p in phases if p['status'] == 'fathomed']
    other = [p for p in phases if p['status'] != 'fathomed']
    print(f"{len(phases)} phases, {sum(p['iters'] for p in phases)} iterations, {sum(p['cpu'] for p in phases):.0f} s: "
          f"{len(fath)} fathom ({sum(p['iters'] for p in fath)} it, {sum(p['cpu'] for p in fath):.0f} s), "
          f"{len(other)} do not ({sum(p['iters'] for p in other)} it, {sum(p['cpu'] for p in other):.0f} s)")
    if fath:
        fi = sorted(fathom_it(p) for p in fath)
        print(f"fathoming phases cross the cutoff at iteration {fi[0]}..{fi[-1]} (median {fi[len(fi)//2]})")
    if not rules:
        rules = ['stall:600:3', 'floor:150', 'hopeless:150:5:50', 'hopeless:200:5:50', 'hopeless:300:5:50',
                 'stall:300:3', 'stall:400:3', 'need:200:3:100', 'need:300:3:100',
                 'gap:200=30', 'gap:300=20', 'gap:200=30,300=20,400=15,600=10',
                 'gap:200=28,300=18,400=14,600=9', 'gap:200=35,300=25,400=20,600=12']
    for r in rules:
        evaluate(r, phases, verbose)


if __name__ == '__main__':
    main()
