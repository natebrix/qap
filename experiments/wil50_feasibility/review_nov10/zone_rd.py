import re, statistics as st, itertools, sys
import os
S=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','sampling','results','')
cut=48814.05
A={}
for ln in open(S+'stageA_trace.txt'):
    m=re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb ([\d.]+) .*? status (\w+) iters (\d+) cpu ([\d.]+).*?\| cut (\d)(.*?)\| cutoff', ln)
    if not m: continue
    d=int(m.group(1)); fx=frozenset(m.group(2).split()); first=float(m.group(3))
    best=first; ph=None
    if m.group(7)=='1':
        mm=re.search(r'best ([\d.]+) rounds (\d+) iters (\d+) cpu ([\d.]+) status (\w+)', m.group(8))
        best=max(first,float(mm.group(1))); ph=(float(mm.group(1)),int(mm.group(3)),float(mm.group(4)),mm.group(5))
    A[fx]=dict(d=d,first=first,best=best,ph=ph,st=m.group(4),it=int(m.group(5)),cpu=float(m.group(6)),fx=m.group(2))
rows=[]
for ln in open(S+'stageB/summary.md'):
    if not ln.startswith('| ') or ln.startswith('| node') or ln.startswith('|---'): continue
    c=[x.strip() for x in ln.strip().strip('|').split('|')]
    if len(c)<11: continue
    try:
        spec=c[0]; inh=float(c[2]); flb=float(c[3]); status=c[5]; it=int(c[6]); cpu=float(c[7]); pool=int(c[9])
    except: continue
    fx=spec.split(',')
    par=None
    for sub in itertools.combinations(fx,2):
        a=A.get(frozenset(sub))
        if a and a['d']==2 and abs(a['best']-inh)<0.06: par=a
    rows.append(dict(spec=spec,inh=inh,flb=flb,status=status,it=it,cpu=cpu,pool=pool,par=par,c=c))
M=[r for r in rows if r['par']]
def pshort(r): return cut-r['par']['first']
def q(v,p):
    v=sorted(v); return v[int(p*(len(v)-1))] if v else float('nan')
print('parent first-solve shortfall window | n children | pooled | fathomed in first solve | mean iters | mean cpu | hopeless: n, child short q25/med/q75/max | parents')
for a,b in ((180,200),(200,215),(215,230),(230,244.08),(244.08,258),(258,272),(272,286),(286,300),(300,325)):
    g=[r for r in M if a<=pshort(r)<b]
    if not g: continue
    h=[cut-r['flb'] for r in g if r['status']!='fathomed']
    npar=len(set(r['par']['fx'] for r in g))
    print(' %6.1f-%6.1f | %3d | %3d | %2d (%.2f) | %5.0f | %5.1f | %2d: %s | %d'%(a,b,len(g),sum(r['pool'] for r in g),sum(r['status']=='fathomed' for r in g),sum(r['status']=='fathomed' for r in g)/len(g),sum(r['it'] for r in g)/len(g),sum(r['cpu'] for r in g)/len(g),len(h),'%.0f/%.0f/%.0f/%.0f'%(q(h,.25),q(h,.5),q(h,.75),max(h)) if h else '-',npar))
# detail: hopeless children by parent window: shortfall classes
print()
for a,b in ((200,244.08),(244.08,300)):
    g=[r for r in M if a<=pshort(r)<b]
    cl=[0,0,0,0,0]
    for r in g:
        s=cut-r['flb']
        if r['status']=='fathomed': cl[0]+=1
        elif s<50: cl[1]+=1
        elif s<130: cl[2]+=1
        elif s<244.08: cl[3]+=1
        else: cl[4]+=1
    print('parents %.0f-%.0f: n=%d  fathomed %d | 0-50 %d | 50-130 %d | 130-244 %d | 244+ %d'%(a,b,len(g),*cl))
print(rows[0]['c'])
