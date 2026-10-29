import re, os, itertools
D=os.path.dirname(os.path.abspath(__file__))+'/'
ch=[l.split() for l in open(D+'children.txt')]
out=open(D+'console_out.txt').read()
qpb=[float(x) for x in re.findall(r"^(?:bnb> )?QPBp\s*=\s*(\S+)", out, re.M)]
assert len(qpb)==len(ch), (len(qpb), len(ch))
dnn=[None]*len(ch)
for q in range(4):
    rs=[l for l in open(D+f'dnn_out_{q}.txt') if l.startswith('RESULT')]
    for t,l in enumerate(rs):
        d=dict(kv.split('=',1) for kv in l.split()[1:])
        dnn[q+4*t]=(float(d['lb']), d['fathomed']=='1', float(d['cpu']))
meta=[l.split() for l in open(D+'meta.txt')]
# threshold for "QPB predicts a fathom": from this data, best single threshold over all children
allq=list(zip(qpb,[d[1] for d in dnn]))
best=max(((sum((q>=t)==f for q,f in allq)/len(allq)),t) for t in sorted(set(qpb)))
T=best[1]
print(f"children: {len(ch)}; DNN fathoms {sum(d[1] for d in dnn)}; child-level QPB threshold {T:.1f} (accuracy {best[0]:.3f})")
def spearman(a,b):
    def rk(v):
        o=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v)
        for k,i in enumerate(o): r[i]=k
        return r
    ra,rb=rk(a),rk(b); n=len(a); ma=sum(ra)/n; mb=sum(rb)/n
    c=sum((x-ma)*(y-mb) for x,y in zip(ra,rb)); va=sum((x-ma)**2 for x in ra); vb=sum((y-mb)**2 for y in rb)
    return c/(va*vb)**.5 if va and vb else float('nan')
tot={'sum':[0,0],'min':[0,0],'count':[0,0]}
for pi,m in enumerate(meta):
    lines=m[5:]
    print(f"\nparent {m[1]} (fixes {m[3]}, DNN lb {float(m[2]):.1f}, most eccentric free location {m[4]})")
    print(f"  {'line':6s} {'misses':>6s} {'DNN cpu':>8s}  {'sum QPB':>9s} {'min QPB':>8s} {'#>=T':>5s}")
    stats=[]
    for li,name in enumerate(lines):
        idx=[k for k,c in enumerate(ch) if int(c[0])==pi and int(c[1])==li]
        miss=sum(1 for k in idx if not dnn[k][1]); cpu=sum(dnn[k][2] for k in idx)
        qs=[qpb[k] for k in idx]
        st=(name,miss,cpu,sum(qs)/len(qs),min(qs),sum(1 for x in qs if x>=T),len(idx))
        stats.append(st)
        print(f"  {name:6s} {miss:6d} {cpu:8.0f}  {st[3]:9.1f} {st[4]:8.1f} {st[5]:5d}/{st[6]}")
    misses=[s[1] for s in stats]; bestm=min(misses)
    for key,col in (('sum',3),('min',4),('count',5)):
        sc=[s[col] for s in stats]; choice=max(range(len(stats)),key=lambda i:sc[i])
        regret=misses[choice]-bestm
        tot[key][0]+=regret; tot[key][1]+=misses[choice]
        print(f"  score {key:5s}: Spearman(score, -misses) {spearman(sc,[-x for x in misses]):+.2f}; picks {stats[choice][0]} with {misses[choice]} misses (best line {bestm}, mean {sum(misses)/len(misses):.1f})")
print("\nover all parents: misses of the chosen line minus the best line's ->", {k:v[0] for k,v in tot.items()})
