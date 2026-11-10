import sys, re, collections
def parse(path):
    rows = []
    for ln in open(path):
        m = re.search(r'depth (\d+) fixes ([\d: ]+)\| first lb ([\d.\-]+) primal \S+ status (\w+) iters (\d+) cpu ([\d.]+) warm (\d) pool (\d)', ln)
        if not m: continue
        d = dict(depth=int(m.group(1)), fixes=m.group(2).strip(), lb=float(m.group(3)), st=m.group(4),
                 it=int(m.group(5)), cpu=float(m.group(6)), warm=int(m.group(7)), pool=int(m.group(8)))
        c = re.search(r'\| cut 1 best ([\d.\-]+) rounds (\d+) iters (\d+) cpu ([\d.]+) status (\w+)', ln)
        d['cut'] = None
        if c:
            d['cut'] = dict(best=float(c.group(1)), it=int(c.group(3)), cpu=float(c.group(4)), st=c.group(5))
        d['fath'] = int(re.search(r'fathomed (\d) pool_out', ln).group(1))
        rows.append(d)
    return rows
if __name__ == '__main__':
    for p in sys.argv[1:]:
        r = parse(p)
        n = len(r); f = sum(x['fath'] for x in r)
        fi = sum(x['it'] for x in r); fc = sum(x['cpu'] for x in r)
        ci = sum(x['cut']['it'] for x in r if x['cut']); cc = sum(x['cut']['cpu'] for x in r if x['cut'])
        np_ = sum(1 for x in r if x['cut']); npf = sum(1 for x in r if x['cut'] and x['cut']['st'] == 'fathomed')
        print(f"{p.split('/')[-3] if 'slot0' not in p else p.split('/')[-6]}: nodes {n} fathomed {f} | first it {fi} cpu {fc:.0f} | phases {np_} ({npf} fath) it {ci} cpu {cc:.0f} | total it {fi+ci} cpu {fc+cc:.0f}")
        by = collections.Counter()
        for x in r:
            by[(x['depth'], x['st'])] += 1
        print('   ', dict(sorted(by.items())))
        hop = [x for x in r if x['st'] != 'fathomed']
        print('    not fathomed in first solve:', len(hop), 'iters', sum(x['it'] for x in hop), 'cpu %.0f' % sum(x['cpu'] for x in hop))
        fa = [x for x in r if x['st'] == 'fathomed']
        its = sorted(x['it'] for x in fa)
        print('    fathomed first solves:', len(fa), 'iters', sum(its), 'cpu %.0f' % sum(x['cpu'] for x in fa), 'median', its[len(its)//2] if its else None, 'hist', dict(sorted(collections.Counter(min(i, 60) // 5 * 5 for i in its).items())))
