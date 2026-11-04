#!/usr/bin/env python3
"""CPU and iteration breakdown of a QAP_DNN_TRACE file by depth, pool, outcome.

usage: trace_breakdown.py trace.txt
"""
import re, sys, collections


def parse(path):
    rows = []
    for line in open(path):
        d = {}
        m = re.search(r'depth (\d+)', line)
        d['depth'] = int(m.group(1))
        m = re.search(r'first lb (\S+) primal (\S+) status (\S+) iters (\d+) cpu (\S+) '
                      r'warm (\d) pool (\d) kept (\d+)', line)
        d['lb'] = float(m.group(1)); d['primal'] = float(m.group(2)); d['status'] = m.group(3)
        d['it'] = int(m.group(4)); d['cpu'] = float(m.group(5))
        d['warm'] = int(m.group(6)); d['pool'] = int(m.group(7)); d['kept'] = int(m.group(8))
        m = re.search(r'\| cut (\d)(.*?)\|', line)
        d['cut'] = int(m.group(1))
        if d['cut']:
            mm = re.search(r'best (\S+) rounds (\d+) iters (\d+) cpu (\S+) status (\S+)', m.group(2))
            d['cbest'] = float(mm.group(1)); d['crounds'] = int(mm.group(2))
            d['cit'] = int(mm.group(3)); d['ccpu'] = float(mm.group(4)); d['cstatus'] = mm.group(5)
        m = re.search(r'cutoff (\S+) fathomed (\d)', line)
        d['cutoff'] = float(m.group(1)); d['fath'] = int(m.group(2))
        d['fixes'] = line.split('fixes')[1].split('|')[0].strip()
        rows.append(d)
    return rows


if __name__ == '__main__':
    rows = parse(sys.argv[1])
    tot_first = sum(r['cpu'] for r in rows)
    tot_cut = sum(r.get('ccpu', 0) for r in rows)
    print(f"{len(rows)} nodes, first-solve CPU {tot_first:.0f} s, cut CPU {tot_cut:.0f} s")
    g = collections.defaultdict(list)
    for r in rows:
        g[(r['depth'], r['pool'], r['warm'], r['status'])].append(r)
    print("depth pool warm status      n   it_mean it_min it_max  cpu_sum  ms/it  mean_short_at_end")
    for k in sorted(g):
        v = g[k]
        its = [r['it'] for r in v]
        cpu = sum(r['cpu'] for r in v)
        gap = sum(r['cutoff'] - r['lb'] for r in v) / len(v)
        print(f"{k[0]:5d} {k[1]:4d} {k[2]:4d} {k[3]:10s} {len(v):4d} {sum(its)/len(its):8.0f} "
              f"{min(its):6d} {max(its):6d} {cpu:8.0f} {1000*cpu/sum(its):6.0f} {gap:10.1f}")
    print("\ncut phases:")
    for r in rows:
        if r['cut']:
            print(f"depth {r['depth']} pool {r['pool']} warm {r['warm']} first {r['lb']:.1f} "
                  f"({r['cutoff']-r['lb']:.1f} short, it {r['it']}, {r['status']}) -> best {r['cbest']:.1f} "
                  f"rounds {r['crounds']} iters {r['cit']} cpu {r['ccpu']:.0f} {r['cstatus']}")
