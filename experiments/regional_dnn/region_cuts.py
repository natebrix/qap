"""Region Boolean cuts on qap's native X[grid_location, flow_facility].

Separation is heuristic evidence; only region_probe's certified lb is a bound.
The saved Y is an approximate iterate, not necessarily DNN feasible.
"""
import argparse, itertools, json
from pathlib import Path
import numpy as np

def node_indices(spec,n=50):
    fixes=[] if spec=='-' else [tuple(map(int,t.split(':'))) for t in spec.split(',')]
    return ([i for i in range(n) if i not in {a for a,b in fixes}],
            [j for j in range(n) if j not in {b for a,b in fixes}])

def grid_regions(rows, mode='threshold'):
    rr=np.array(rows)//10;cc=np.array(rows)%10
    regions=[(f'row<{r}',rr<r) for r in range(1,5)]
    regions += [(f'col<{c}',cc<c) for c in range(1,10)]
    if mode=='rectangles':
        regions += [(f'box_{r0}_{r1}_{c0}_{c1}',(rr>=r0)&(rr<r1)&(cc>=c0)&(cc<c1))
                    for r0 in range(5) for r1 in range(r0+1,6)
                    for c0 in range(10) for c1 in range(c0+1,11)]
    result=[];seen=set()
    for name,mask in regions:
        key=tuple(np.flatnonzero(mask));comp=tuple(np.flatnonzero(~mask))
        if not key or not comp or key in seen or comp in seen:continue
        seen.add(key);result.append((name,mask))
    return result

def aggregate(Y,mask,m):
    # diag(Y) is used for x: identical to row zero on the exact DNN face.
    p=np.diag(Y)[1:].reshape(m,m)[mask].sum(axis=0)
    W=Y[1:,1:].reshape(m,m,m,m)[mask,:,:,:][:,:,mask,:].sum(axis=(0,2))
    return p,W

def triangle_separate(Y,rows,mode='threshold',keep=200):
    m=len(rows);ijk=np.array(list(itertools.combinations(range(m),3)))
    a,b,c=ijk.T;found=[];summary=[]
    for name,mask in grid_regions(rows,mode):
        p,W=aggregate(Y,mask,m)
        vals=np.array([W[a,b]+W[a,c]-W[b,c]-p[a],
                       W[a,b]+W[b,c]-W[a,c]-p[b],
                       W[a,c]+W[b,c]-W[a,b]-p[c],
                       p[a]+p[b]+p[c]-W[a,b]-W[a,c]-W[b,c]-1])
        summary.append(dict(region=name,max_violation=float(vals.max()),
                            violated=int((vals>1e-5).sum()),capacity=int(mask.sum())))
        ids=np.argwhere(vals>1e-5)
        if len(ids)>keep:
            v=vals[ids[:,0],ids[:,1]];ids=ids[np.argpartition(v,-keep)[-keep:]]
        for typ,t in ids:
            found.append(dict(region=name,locations=np.flatnonzero(mask).tolist(),
                facilities=ijk[t].tolist(),type=int(typ),violation=float(vals[typ,t])))
    found.sort(key=lambda c:-c['violation'])
    return found[:keep],summary

def cut_value(Y,cut,m):
    mask=np.zeros(m,bool);mask[cut['locations']]=True;p,W=aggregate(Y,mask,m)
    a,b,c=cut['facilities'];t=cut['type']
    if t==0:return W[a,b]+W[a,c]-W[b,c]-p[a]
    if t==1:return W[a,b]+W[b,c]-W[a,c]-p[b]
    if t==2:return W[a,c]+W[b,c]-W[a,b]-p[c]
    return p[a]+p[b]+p[c]-W[a,b]-W[a,c]-W[b,c]-1

def cut_matrix(cut,m):
    N=m*m+1;G=np.zeros((N,N));loc=np.array(cut['locations'],dtype=int);a,b,c=cut['facilities'];t=cut['type']
    ids=lambda f:1+loc*m+f
    def pair(i,j,v):
        I,J=ids(i),ids(j);G[np.ix_(I,J)]+=v/2;G[np.ix_(J,I)]+=v/2
    def lin(i,v):
        I=ids(i);G[I,I]+=v
    if t<3:
        apex=[a,b,c][t];other=[f for f in [a,b,c] if f!=apex]
        pair(apex,other[0],1);pair(apex,other[1],1);pair(other[0],other[1],-1);lin(apex,-1)
        h=0
    else:
        for i in [a,b,c]:lin(i,1)
        pair(a,b,-1);pair(a,c,-1);pair(b,c,-1);h=1
    return G,h

def mixed_separate(Y,rows,keep=200):
    """Allow each of the three Boolean observables its own grid region."""
    from scipy.sparse import csr_matrix
    m=len(rows);regions=grid_regions(rows);obs=[];ii=[];jj=[]
    for name,mask in regions:
        loc=np.flatnonzero(mask).tolist()
        for f in range(m):
            t=len(obs);obs.append(dict(region=name,locations=loc,facility=f))
            ii.extend([t]*len(loc));jj.extend([a*m+f for a in loc])
    T=csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(obs),m*m))
    p=T@np.diag(Y)[1:];W=(T@(T@Y[1:,1:]).T).T
    found=[];count=0;mx=-np.inf;K=len(obs)
    for a in range(K-2):
        b,c=np.triu_indices(K-a-1,1);b=b+a+1;c=c+a+1
        vals=np.array([W[a,b]+W[a,c]-W[b,c]-p[a],
                       W[a,b]+W[b,c]-W[a,c]-p[b],
                       W[a,c]+W[b,c]-W[a,b]-p[c],
                       p[a]+p[b]+p[c]-W[a,b]-W[a,c]-W[b,c]-1])
        mx=max(mx,float(vals.max()));count+=int((vals>1e-5).sum())
        ids=np.argwhere(vals>1e-5)
        if len(ids)>keep:
            v=vals[ids[:,0],ids[:,1]];ids=ids[np.argpartition(v,-keep)[-keep:]]
        for typ,t in ids:
            found.append(dict(observables=[obs[a],obs[b[t]],obs[c[t]]],
                              type=int(typ),violation=float(vals[typ,t])))
        if len(found)>keep*3:
            found.sort(key=lambda c:-c['violation']);found=found[:keep]
    found.sort(key=lambda c:-c['violation'])
    return found[:keep],dict(observables=K,inequalities=4*K*(K-1)*(K-2)//6,
                           max_violation=mx,violated=count)

def mixed_matrix(cut,m):
    N=m*m+1;G=np.zeros((N,N));ids=[1+np.array(o['locations'],dtype=int)*m+o['facility'] for o in cut['observables']]
    def pair(a,b,v):
        # np.add.at handles possible overlap of the two observable supports.
        I,J=np.meshgrid(ids[a],ids[b],indexing='ij')
        np.add.at(G,(I,J),v/2);np.add.at(G,(J,I),v/2)
    def lin(a,v):np.add.at(G,(ids[a],ids[a]),v)
    t=cut['type']
    if t<3:
        b,c=[i for i in range(3) if i!=t];pair(t,b,1);pair(t,c,1);pair(b,c,-1);lin(t,-1);h=0
    else:
        for i in range(3):lin(i,1)
        pair(0,1,-1);pair(0,2,-1);pair(1,2,-1);h=1
    return G,h

def feasibility(Y,m):
    X=np.diag(Y)[1:].reshape(m,m);T=Y[1:,:].reshape(m,m,-1)
    return dict(symmetry=float(abs(Y-Y.T).max()),minimum_entry=float(Y.min()),
      diagonal_link=float(abs(Y[0]-np.diag(Y)).max()),
      row_moment=float(abs(T.sum(axis=1)-Y[0]).max()),
      col_moment=float(abs(T.sum(axis=0)-Y[0]).max()),
      row_assignment=float(abs(X.sum(axis=1)-1).max()),
      col_assignment=float(abs(X.sum(axis=0)-1).max()))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('prefix');ap.add_argument('--fixes',default='-')
    ap.add_argument('--mode',choices=['threshold','rectangles','mixed'],default='threshold')
    ap.add_argument('--keep',type=int,default=200);args=ap.parse_args()
    rows,cols=node_indices(args.fixes);m=len(rows)
    Y=np.fromfile(args.prefix+'.Y.bin').reshape(m*m+1,m*m+1)
    if args.mode=='mixed':
        cuts,stats=mixed_separate(Y,rows,args.keep)
        result=dict(fixes=args.fixes,rows=rows,cols=cols,feasibility=feasibility(Y,m),cuts=cuts,stats=stats)
        Path(args.prefix+'.mixed.json').write_text(json.dumps(result,indent=2))
        print(json.dumps(dict(feasibility=result['feasibility'],stats=stats,top=cuts[:3]),indent=2))
        raise SystemExit
    cuts,stats=triangle_separate(Y,rows,args.mode,args.keep)
    result=dict(fixes=args.fixes,rows=rows,cols=cols,feasibility=feasibility(Y,m),cuts=cuts,regions=stats)
    Path(args.prefix+'.'+args.mode+'.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(feasibility=result['feasibility'],regions=len(stats),
       max_violation=max(t['max_violation'] for t in stats),
       violated=sum(t['violated'] for t in stats),top=cuts[:3]),indent=2))
