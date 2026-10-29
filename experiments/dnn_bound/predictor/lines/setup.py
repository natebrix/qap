import csv, random, os
D=os.path.dirname(os.path.abspath(__file__))+'/'
rows=list(csv.DictReader(open(D+'../pred/table.csv')))
miss=sorted([r for r in rows if r['fathomed']=='0'], key=lambda r: -float(r['dnn_lb']))
pick=[miss[0], miss[len(miss)//4], miss[len(miss)//2], miss[3*len(miss)//4], miss[-1]]
n=30
def ecc(j, freecols):
    r,c=divmod(j,6)
    return max(abs(r-(k//6))+abs(c-(k%6)) for k in freecols if k!=j)
rng=random.Random(11)
children=[]  # (parent idx, line id, fixes string)
meta=[]
for pi,p in enumerate(pick):
    fx=[tuple(map(int,t.split(':'))) for t in p['fixes'].split(';')]
    fr=[i for i in range(n) if i not in [a for a,b in fx]]
    fc=[j for j in range(n) if j not in [b for a,b in fx]]
    # nug30 locations are a 5 x 6 grid, location j at row j//6, col j%6 (cf. the 0-based .dat order)
    ecol=max(fc, key=lambda j:(ecc(j,fc), -j))
    cols=rng.sample([j for j in fc if j!=ecol],4)+[ecol]
    rws=rng.sample(fr,4)
    lines=[('row',i) for i in rws]+[('col',j) for j in cols]
    meta.append((pi,p['k'],p['dnn_lb'],p['fixes'],lines,ecol))
    for li,(kind,x) in enumerate(lines):
        for y in (fc if kind=='row' else fr):
            a,b=(x,y) if kind=='row' else (y,x)
            f=fx+[(a,b)]
            children.append((pi,li,",".join(f"{i}:{j}" for i,j in sorted(f))))
open(D+'meta.txt','w').write("\n".join(f"{pi} {k} {lb} {fx} {ecol} "+" ".join(f"{a}{b}" for a,b in lines) for pi,k,lb,fx,lines,ecol in meta)+"\n")
open(D+'children.txt','w').write("\n".join(f"{a} {b} {c}" for a,b,c in children)+"\n")
cmds=[]
for _,_,c in children:
    cmds.append("c")
    for t in c.split(','):
        i,j=t.split(':'); cmds+=["f",i,j]
    cmds.append("b")
cmds.append("q")
open(D+'console_cmds.txt','w').write("\n".join(cmds)+"\n")
for q in range(4):
    with open(D+f'dnn_nodes_{q}.txt','w') as f:
        for k,(_,_,c) in enumerate(children):
            if k%4==q: f.write(c+" 6122.05\n")
print(len(children),"children;", [(m[1],round(float(m[2]),1)) for m in meta])
