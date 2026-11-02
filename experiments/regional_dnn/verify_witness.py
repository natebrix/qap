"""Exact rational witness: DNN + all atomic triangles != regional triangles."""
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path
import numpy as np
import sympy as sp
from region_cuts import cut_matrix, cut_value, mixed_matrix

n=6;N=37
Q=sp.eye(n)
for i in range(n):
    for j in range(n):
        if i!=j and i//3==j//3:Q[i,j]=-sp.Rational(1,2)
W=(Q+sp.ones(n))/4
Y=sp.zeros(N);Y[0,0]=1
for a in range(n):
 for i in range(n):
  u=1+a*n+i;Y[0,u]=Y[u,0]=sp.Rational(1,6)
  for b in range(n):
   for j in range(n):
    v=1+b*n+j
    if a==b and i==j:Y[u,v]=sp.Rational(1,6)
    elif a==b or i==j:Y[u,v]=0
    elif (a<3)==(b<3):Y[u,v]=W[i,j]/6
    else:Y[u,v]=(sp.Rational(1,2)-W[i,j])/9

assert all(v>=0 for v in Y)
assert Y.diagonal()==Y[0,:]
for k in range(N):
 for i in range(n):
  assert sum(Y[1+a*n+i,k] for a in range(n))==Y[0,k]
  assert sum(Y[1+i*n+a,k] for a in range(n))==Y[0,k]
eigs=Y.eigenvals();assert all(v>=0 for v in eigs)
max_atomic=-sp.oo
for a,b,c in combinations(range(1,N),3):
    vals=[Y[a,b]+Y[a,c]-Y[b,c]-Y[a,a],
          Y[a,b]+Y[b,c]-Y[a,c]-Y[b,b],
          Y[a,c]+Y[b,c]-Y[a,b]-Y[c,c],
          Y[a,a]+Y[b,b]+Y[c,c]-Y[a,b]-Y[a,c]-Y[b,c]-1]
    max_atomic=max(max_atomic,*vals)
assert max_atomic<=0
test_cut=dict(locations=[0,1,2],facilities=[0,1,2],type=3)
Yn=np.array(Y,dtype=float);G,h=cut_matrix(test_cut,n)
v=cut_value(Yn,test_cut,n)
assert abs(v-1/8)<1e-12 and abs((G*Yn).sum()-h-v)<1e-12
maximum=-1.;checks=0
for perm in permutations(range(n)):
    x=np.zeros((n,n));x[np.arange(n),perm]=1;y=np.r_[1,x.ravel()];YY=np.outer(y,y)
    for typ in range(4):
        cut=dict(test_cut,type=typ);GG,hh=cut_matrix(cut,n)
        violation=cut_value(YY,cut,n)
        assert abs((GG*YY).sum()-hh-violation)<1e-12
        assert violation<=1e-12;maximum=max(maximum,violation);checks+=1
# Each Boolean triangle slack g is 0 or 1, hence g^2=g modulo z_i^2=z_i.
z=sp.symbols('z0:3');g=1-sum(z)+sum(z[i]*z[j] for i,j in combinations(range(3),2))
poly=sp.Poly(sp.expand(g*g-g),*z);reduced=0
for powers,coeff in poly.terms():
    reduced+=coeff*sp.prod(z[i] for i,e in enumerate(powers) if e)
assert sp.expand(reduced)==0
rng=np.random.default_rng(310);mixed_checks=0
for trial in range(25):
    obs=[dict(facility=int(rng.integers(n)),locations=np.flatnonzero(rng.integers(0,2,size=n)).tolist()) for _ in range(3)]
    for typ in range(4):
        GG,hh=mixed_matrix(dict(observables=obs,type=typ),n)
        for _ in range(20):
            perm=rng.permutation(n);x=np.zeros((n,n));x[np.arange(n),perm]=1;y=np.r_[1,x.ravel()]
            a,b,c=[x[o['locations'],o['facility']].sum() for o in obs]
            v=[a*b+a*c-b*c-a,a*b+b*c-a*c-b,a*c+b*c-a*b-c,a+b+c-a*b-a*c-b*c-1][typ]
            assert v<=0 and abs(y@GG@y-hh-v)<1e-10;mixed_checks+=1
result=dict(exact_eigenvalues={str(k):v for k,v in eigs.items()},
   max_atomic_triangle_violation=str(max_atomic),regional_triangle_violation='1/8',
   permutation_cut_checks=checks,max_permutation_violation=maximum,moment_order2_SOS_identity=True,
   mixed_support_coefficient_checks=mixed_checks)
Path(__file__).with_name('witness_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
