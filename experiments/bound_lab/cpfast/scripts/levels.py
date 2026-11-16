#!/usr/bin/env python3
"""levels.py [-g GAIN] [-c CUTOFF] BASE.log VARIANT.log ...: iterations of the phase to reach 25/50/75/90% of
the reference gain (default: the base log's final gain over the first solve) and the cutoff (240,514.05),
each variant's factor against the base, and CPU per iteration. Markdown table."""
import sys, re
sys.path.insert(0, __file__.rsplit('/',1)[0])
from curve import read
args=sys.argv[1:]; gain=None; cutoff=240514.05
while args and args[0] in('-g','-c'):
    if args[0]=='-g': gain=float(args[1])
    else: cutoff=float(args[1])
    args=args[2:]
R=[read(f) for f in args]; names=[f.split('/')[-1].replace('.log','') for f in args]
b=R[0]; first=b['first']
if gain is None: gain=b['lb'][-1]-first
def reach(r,L):
    x=[i for i,v in zip(r['it'],r['lb']) if v>=L]
    return x[0] if x else None
levels=[('25%',first+.25*gain),('50%',first+.5*gain),('75%',first+.75*gain),('90%',first+.9*gain),('100%',first+gain),('cutoff',cutoff)]
print('reference gain %.1f over the first solve %.1f; cutoff %.2f'%(gain,first,cutoff))
print('| run | its | end bound | '+' | '.join('%s (%.0f)'%(n,L) for n,L in levels)+' | s/it |')
print('|---|---|---|'+'---|'*len(levels)+'---|')
base=[reach(b,L) for _,L in levels]
for n,r in zip(names,R):
    cells=[]
    for (ln,L),b0 in zip(levels,base):
        x=reach(r,L)
        if x is None: cells.append('> %d'%r['it'][-1])
        elif r is b or b0 is None: cells.append('%d'%x)
        else: cells.append('%d (%.1fx)'%(x,b0/x))
    print('| %s | %d | %.1f | %s | %.3f |'%(n,r['it'][-1],r['lb'][-1],' | '.join(cells),r['cpu']/r['it'][-1]))
