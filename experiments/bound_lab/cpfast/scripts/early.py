#!/usr/bin/env python3
"""early.py DNN.txt ...: the first phase of each trace, shortfall to the cutoff of its running best at marks."""
import sys,re
for f in sys.argv[1:]:
    l=open(f).read()
    h=re.search(r'\| hist (.*)',l)
    if not h: print(f,'no hist'); continue
    d=dict((int(a),float(b)) for a,b in (x.split(':') for x in h[1].split()))
    print('%-34s'%f.split('/')[-2],' '.join('%d:%.0f'%(k,240514.05-d[k]) for k in [20,40,50,60,80,100,120,150,200,250,300,400,500,600,800] if k in d))
