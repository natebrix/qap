"""Cluster structure of the tai-b instances, from the matrices.
The trace's facilities (rows, the first index of a fix) index the first matrix A (symmetric; fixconv.py),
locations (columns) the second matrix B (asymmetric). Facilities: single-linkage clusters of A, cut at the
largest relative gap among the minimum spanning tree's edges. Locations: B has no geometry; we report the
same MST statistic on a profile distance between locations (1 - correlation of their B rows+columns) to
show whether any grouping exists."""
import sys, json, math
from common import *

def mst_edges(D, n):
    inT=[False]*n; best=[math.inf]*n; par=[-1]*n; best[0]=0; E=[]
    for _ in range(n):
        u=min((best[i],i) for i in range(n) if not inT[i])[1]; inT[u]=True
        if par[u]>=0: E.append((D[u][par[u]],u,par[u]))
        for v in range(n):
            if not inT[v] and D[u][v]<best[v]: best[v]=D[u][v]; par[v]=u
    return sorted(E)

def components(n, E):
    par=list(range(n))
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for _,u,v in E: par[f(u)]=f(v)
    g={}
    for i in range(n): g.setdefault(f(i),[]).append(i)
    return sorted(g.values(), key=lambda c:c[0])

def cluster(D, n):
    E=mst_edges(D,n)
    w=[e[0] for e in E]
    # largest ratio between consecutive sorted MST edge lengths (ignore zeros)
    best=(0,None)
    for k in range(1,len(w)):
        if w[k-1]>0:
            r=w[k]/w[k-1]
            if r>best[0]: best=(r,k)
    r,k=best
    keep=E[:k]
    return components(n,keep), r, (w[k-1], w[k])

def corrdist(B,n):
    prof=[B[i]+[B[j][i] for j in range(n)] for i in range(n)]
    def z(v):
        m=sum(v)/len(v); s=math.sqrt(sum((x-m)**2 for x in v)) or 1; return [(x-m)/s for x in v]
    Z=[z(p) for p in prof]
    return [[1-sum(a*b for a,b in zip(Z[i],Z[j])) for j in range(n)] for i in range(n)]

if __name__=='__main__':
    res={}
    for inst in ['tai35b','tai40b','tai50b','tai60b']:
        n,A,B=read_dat(inst)
        cl,r,(lo,hi)=cluster(A,n)
        if inst=='tai60b':
            # the largest gap (572/186) leaves a 36-facility cluster that row 0 of A shows is three groups
            # (A[0][14..25] ~1,390, [26..47] ~1,100, [48..51] ~1,530); cut the five MST edges >= 162 (next: 106)
            E=mst_edges(A,n); cl=components(n,E[:n-6]); lo,hi=E[n-7][0],E[n-6][0]; r=hi/lo
        within=[A[i][j] for c in cl for i in c for j in c if i<j]
        between=[A[i][j] for a in range(len(cl)) for b in range(a+1,len(cl)) for i in cl[a] for j in cl[b]]
        print(f'{inst}: facilities (A) {len(cl)} clusters, sizes {[len(c) for c in cl]}; cut between MST edges {lo} and {hi} (ratio {r:.1f}); '
              f'within-cluster A max {max(within)}, between min {min(between)}')
        for c in cl: print('   ', c)
        Dl=corrdist(B,n); cl2,r2,(lo2,hi2)=cluster(Dl,n)
        print(f'   locations (B profile 1-corr): largest MST gap ratio {r2:.2f} ({lo2:.3f} -> {hi2:.3f}); split sizes {sorted(len(c) for c in cl2)[::-1][:8]}')
        res[inst]={'fac_clusters':cl}
    json.dump(res,open('out/clusters.json','w'))
