import re, sys, glob, os
P=os.path.dirname(os.path.abspath(__file__))+'/'
nodes=[l.split() for l in open(P+'nodes.txt')]
out=open(P+'console_out.txt').read()
def seq(tag): return [float(x) for x in re.findall(r"^(?:bnb> )?"+tag+r"\s*=\s*(\S+)", out, re.M)]
qpb, qpbp, pb, glb = seq("QPB "), seq("QPBp"), seq("PB "), seq("GLB")
assert len(qpb)==len(qpbp)==len(pb)==len(glb)==225, (len(qpb),len(qpbp),len(pb),len(glb))
dnn=[None]*225
for q in range(4):
    rs=[l for l in open(P+f'dnn_out_{q}.txt') if l.startswith('RESULT')]
    for t,l in enumerate(rs):
        d=dict(kv.split('=',1) for kv in l.split()[1:])
        dnn[q+4*t]=(float(d['lb']), d['fathomed']=='1', int(d['iters']))
y=[1 if d[1] else 0 for d in dnn]
preds={'GLB':glb,'PB':pb,'QPB':qpb,'QPB (parametric)':qpbp,'lookahead estimate':[float(n[2]) for n in nodes],'inherited bound':[float(n[1]) for n in nodes]}
def auc(x,y):
    pos=[a for a,b in zip(x,y) if b]; neg=[a for a,b in zip(x,y) if not b]
    s=0
    for a in pos:
        for b in neg: s+= 1 if a>b else (0.5 if a==b else 0)
    return s/(len(pos)*len(neg))
def best_acc(x,y):
    best=(0,None)
    for t in sorted(set(x)):
        acc=sum((a>=t)==bool(b) for a,b in zip(x,y))/len(y)
        best=max(best,(acc,t))
    return best
def spearman(a,b):
    def rank(v):
        o=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v)
        for k,i in enumerate(o): r[i]=k
        return r
    ra,rb=rank(a),rank(b); n=len(a)
    ma=sum(ra)/n; mb=sum(rb)/n
    cov=sum((x-ma)*(z-mb) for x,z in zip(ra,rb)); va=sum((x-ma)**2 for x in ra); vb=sum((z-mb)**2 for z in rb)
    return cov/(va*vb)**0.5
lb=[d[0] for d in dnn]
print(f"{sum(y)} of 225 depth-2 nodes fathomed by the DNN bound within 1000 iterations")
print(f"{'predictor':22s} {'AUC':>6s} {'best acc':>9s} {'threshold':>10s} {'Spearman(DNN lb)':>17s}  range")
for k,v in preds.items():
    a,t=best_acc(v,y)
    print(f"{k:22s} {auc(v,y):6.3f} {a:9.3f} {t:10.1f} {spearman(v,lb):17.3f}  {min(v):.0f}..{max(v):.0f}")
with open(P+'table.csv','w') as f:
    f.write("k,fixes,inherited,lookahead,QPB,QPBp,PB,GLB,dnn_lb,fathomed,iters\n")
    for k in range(225):
        f.write(f"{k},{nodes[k][3].replace(',',';')},{nodes[k][1]},{nodes[k][2]},{qpb[k]},{qpbp[k]},{pb[k]},{glb[k]},{lb[k]},{y[k]},{dnn[k][2]}\n")
