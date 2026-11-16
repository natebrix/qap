#!/usr/bin/env python3
"""validity.py V_*.log: nodes on the path to the best known solution (optimum 240,516 if it is optimal): the largest
certified value over every check of the phase must not pass 240,516."""
import sys, re
worst=-1e30; n=0
print('| run | depth | m | first solve | phase its | certified at 100 / 300 / 1000 | max certified, all checks | 240,516 - max | checks | |')
print('|---|---|---|---|---|---|---|---|---|---|')
for f in sys.argv[1:]:
    first=None; vals=[]; spec=''; m=None
    for l in open(f):
        if l.startswith('# /'): spec=re.search(r'-x (\S+)',l)[1]
        x=re.match(r'first solve\s+lb ([\d.]+)',l) or re.match(r'first solve loaded.*lb ([\d.]+)',l)
        if x: first=float(x[1])
        x=re.search(r' m=(\d+) ',l)
        if x and m is None: m=int(x[1])
        x=re.match(r'cp it\s+(\d+) cert ([\d.]+)',l)
        if x: vals.append((int(x[1]),float(x[2])))
    if not vals: continue
    mx=max(v for _,v in vals); worst=max(worst,mx); n+=len(vals)
    at=lambda k: ('%.2f'%[v for i,v in vals if i==k][0]) if any(i==k for i,_ in vals) else '-'
    print('| %s | %d | %d | %.1f | %d | %s / %s / %s | %.4f | %.4f | %d | %s |'%(f.split('/')[-1][:-4],spec.count(':'),m,first,vals[-1][0],at(100),at(300),at(1000),mx,240516-mx,len(vals),'PASS' if mx<=240516 else 'FAIL'))
print('\nworst: max certified - 240516 = %.6f over %d certified values'%(worst-240516,n))
