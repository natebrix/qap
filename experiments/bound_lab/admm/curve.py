#!/usr/bin/env python3
# usage: curve.py its(comma) log... : certified running best (lb) at the given iterations, plus pres/dY at the last
import sys,re
its=[int(x) for x in sys.argv[1].split(',')]
for f in sys.argv[2:]:
    d={}
    for l in open(f):
        m=re.match(r'it\s+(\d+) lb ([-\d.e+]+) .*pres ([\d.e+-]+) dY ([\d.e+-]+) npos (\d+)',l)
        if m: d[int(m.group(1))]=(float(m.group(2)),float(m.group(3)),float(m.group(4)),int(m.group(5)))
    res=[l for l in open(f) if l.startswith('RESULT')]
    st=re.search(r'status=(\w+) iters=(\d+)',res[-1]) if res else None
    row=[]
    for i in its:
        k=[j for j in d if j<=i]
        row.append('%.0f'%d[max(k)][0] if k and (i in d or max(d)>=i or True) and (i<=max(d)) else '-')
    last=d[max(d)]
    print(f.split('/')[-1].replace('.log',''), *row, 'pres %.2e dY %.2e npos %d'%(last[1],last[2],last[3]), st.group(1)+'@'+st.group(2) if st else '', sep=' | ')
