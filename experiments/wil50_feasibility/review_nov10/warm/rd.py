import re, sys, glob, os
W = os.path.dirname(os.path.abspath(__file__))
def parse(fn):
    nodes = []; cur = None
    for ln in open(fn):
        if re.match(r'node (\d+) fixed=', ln) or ln.startswith('CHILD-PASS'):
            cur = {'lb': {}, 'res': None}; nodes.append(cur); continue
        m = re.match(r'\s*it\s+(\d+) lb ([\d.\-]+)', ln)
        if m and cur is not None:
            cur['lb'][int(m.group(1))] = float(m.group(2)); continue
        m = re.match(r'RESULT node=(\d+) m=\d+ lb=([\d.\-]+) .*status=(\w+) iters=(\d+) cpu=([\d.]+)', ln)
        if m and cur is not None:
            cur['res'] = (float(m.group(2)), m.group(3), int(m.group(4)), float(m.group(5)))
    return nodes
its = [1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30]
for fam in sys.argv[1:]:
    arms = sorted(os.path.basename(f)[:-len('_%s.log' % fam)] for f in glob.glob(os.path.join(W, '*_%s.log' % fam)))
    order = ['base'] + [a for a in arms if a != 'base']
    data = {a: parse(os.path.join(W, '%s_%s.log' % (a, fam))) for a in order}
    par = data['base'][0]['res'][0]
    print('family', fam, 'parent final', '%.2f' % par, ' parents:', ' '.join('%s=%s' % (a, data[a][0]['res'][:3]) for a in order if data[a] and data[a][0]['res']))
    nk = max(len(d) for d in data.values())
    for k in range(1, nk):
        print(' child', k)
        print('  %-6s' % 'arm' + ''.join('%8d' % i for i in its) + '   cross status')
        for a in order:
            d = data[a]
            if k >= len(d): continue
            lb = d[k]['lb']; r = d[k]['res']
            row = ''.join(('%8.1f' % (lb[i] - par)) if i in lb else '       .' for i in its)
            print('  %-6s' % a + row + '   ' + (('%d %s' % (r[2], r[1])) if r else 'running'))
