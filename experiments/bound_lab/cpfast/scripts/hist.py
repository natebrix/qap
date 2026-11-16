#!/usr/bin/env python3
"""hist.py DNN.txt ... : the phases of bnbdrv DNN traces (QAP_DNN_TRACE, QAP_DNN_TRACE_HIST=1): one line per
phase: depth, fixes (last two), first-solve iterations and shortfall, pool, phase iterations and outcome, and the
phase's running best (as shortfall to the cutoff) at marks. -n: only node lines whose fixes end with the given token."""
import sys, re
marks=[50,100,200,300,500,750,1000,1500,2000,3000]
def phases(f):
    out=[]
    for l in open(f):
        m=re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb (\S+) primal \S+ status (\S+) iters (\d+) cpu (\S+) warm \d pool (\d+) kept (\d+).*?(?:\| cut 1 best (\S+) rounds \d+ iters (\d+) cpu (\S+) status (\S+) cp \d active \d+ stalled (\d+) gaveup \d+ )?\| cutoff (\S+) fathomed (\d)',l)
        if not m: continue
        d=dict(depth=int(m[1]),fixes=m[2],first=float(m[3]),fstat=m[4],fit=int(m[5]),fcpu=float(m[6]),pool=int(m[8]),cutoff=float(m[14]),fath=int(m[15]))
        if m[9]:
            d.update(cut=float(m[9]),cit=int(m[10]),ccpu=float(m[11]),cstat=m[12],stalled=int(m[13]))
            h=re.search(r'\| hist (.*)',l)
            d['hist']=[(int(a),float(b)) for a,b in (x.split(':') for x in h[1].split())] if h else []
        out.append(d)
    return out
if __name__=='__main__':
    for f in sys.argv[1:]:
        P=phases(f)
        nit=sum(p['fit']+p.get('cit',0) for p in P); cpu=sum(p['fcpu']+p.get('ccpu',0) for p in P)
        print('%s: %d DNN nodes, %d first-solve + %d phase iterations = %d, %.0f CPU s, %d phases (%d fathomed)'%(f,len(P),sum(p['fit'] for p in P),sum(p.get('cit',0) for p in P),nit,cpu,sum('cit' in p for p in P),sum(p.get('cstat')=='fathomed' for p in P)))
        for p in P:
            if 'cit' not in p: continue
            row=[]
            for mk in marks:
                v=[b for i,b in p['hist'] if i<=mk]
                row.append('%7.0f'%(p['cutoff']-v[-1]) if v and p['cit']>=mk else '      -')
            print('  d%d %-22s first %4d it %7.1f short pool %4d | phase %4d it %-9s%s end %7.1f short | short at %s: %s'%(p['depth'],' '.join(p['fixes'].split()[-2:]),p['fit'],p['cutoff']-p['first'],p['pool'],p['cit'],p['cstat'],' stalled' if p['stalled'] else '',max(0,p['cutoff']-p['cut']),'/'.join(map(str,marks)),' '.join(row)))
