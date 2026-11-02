"""Optional selective RLT2 separator for three facilities, using SciPy HiGHS.

Cut potentials are rounded to multiples of 1/1024 and checked on every
injective placement, so validity does not rely on the LP solver's tolerances.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from region_cuts import node_indices

ap=argparse.ArgumentParser();ap.add_argument('prefix');ap.add_argument('--fixes',default='-')
ap.add_argument('--facilities',required=True,help='three LOCAL free B-matrix indices')
ap.add_argument('--seconds',type=float,default=120)
args=ap.parse_args();rows,cols=node_indices(args.fixes);m=len(rows);N=m*m+1
fac=list(map(int,args.facilities.split(',')));assert len(set(fac))==3
Y=np.fromfile(args.prefix+'.Y.bin').reshape(N,N)
ids=[1+np.arange(m)*m+f for f in fac]
P=[Y[np.ix_(ids[0],ids[1])],Y[np.ix_(ids[0],ids[2])],Y[np.ix_(ids[1],ids[2])]]
a,b,c=np.indices((m,m,m));valid=(a!=b)&(a!=c)&(b!=c);a=a[valid];b=b[valid];c=c[valid]
K=len(a);rr=np.tile(np.arange(K),3)
jj=np.r_[a*m+b,m*m+a*m+c,2*m*m+b*m+c]
A=coo_matrix((-np.ones(3*K),(rr,jj)),shape=(K,3*m*m)).tocsr()
t=time.monotonic();lp=linprog(np.concatenate([p.ravel() for p in P]),A_ub=A,
    b_ub=np.zeros(K),bounds=(-1,1),method='highs',options={'time_limit':args.seconds})
if lp.x is None:
    print(json.dumps(dict(success=False,status=lp.message,seconds=time.monotonic()-t)));raise SystemExit(1)
v=np.rint(lp.x.reshape(3,m,m)*1024).astype(np.int64)
lo=int((v[0,a,b]+v[1,a,c]+v[2,b,c]).min())
v[0]+=max(0,-lo)
assert (v[0,a,b]+v[1,a,c]+v[2,b,c]).min()>=0
pot=v/1024.;G=np.zeros_like(Y)
for mat,(f,g) in zip(pot,[(0,1),(0,2),(1,2)]):
    G[np.ix_(ids[f],ids[g])]-=mat/2
    G[np.ix_(ids[g],ids[f])]-=mat.T/2
out=args.prefix+'.triplet_'+'_'.join(map(str,fac));G.tofile(out+'.G.bin')
np.savez_compressed(out+'.npz',integer_potentials=v,denominator=1024,facilities=fac)
result=dict(success=bool(lp.success),status=lp.message,facilities_local=fac,
    facilities_original=[cols[f] for f in fac],lp_objective=float(lp.fun),
    rounded_cut_violation=float((G*Y).sum()),integer_constraint_min=int((v[0,a,b]+v[1,a,c]+v[2,b,c]).min()),
    placement_checks=K,seconds=time.monotonic()-t)
Path(out+'.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
