#!/usr/bin/env python3
"""curve.py LOG...: certified bound of the phase vs iteration (running best), one column per log.
   -t L1,L2,..: also iterations to reach the given bound levels."""
import re, sys
def read(f):
    it=[]; lb=[]; best=-1e30; first=None; cpu=None; extra={}
    for l in open(f, errors='replace'):
        m=re.match(r'first solve.*lb ([\d.]+)',l) or re.match(r'first solve loaded.*lb ([\d.]+)',l)
        if m: first=float(m[1])
        m=re.match(r'it\s+(\d+) lb ([\d.]+) .* cpu ([\d.]+)',l)
        if m:
            best=max(best,float(m[2])); it.append(int(m[1])); lb.append(best); cpu=float(m[3])
        m=re.match(r'cp it\s+(\d+) .*cuts (\d+) active (\d+) sum_nu (\S+) viol (\d+) \(max (\S+)/',l)
        if m: extra[int(m[1])]=(int(m[2]),int(m[3]),float(m[4]),int(m[5]),float(m[6]))
    return dict(it=it,lb=lb,first=first,cpu=cpu,extra=extra)
if __name__=='__main__':
    args=sys.argv[1:]; levels=None; marks=[50,100,200,300,500,750,1000,1500,2000,2500,3000,3400]
    if args and args[0]=='-t': levels=[float(x) for x in args[1].split(',')]; args=args[2:]
    R=[read(f) for f in args]
    names=[f.split('/')[-1].replace('.log','') for f in args]
    print('it'.rjust(6),' '.join(n.rjust(14) for n in names))
    print('first'.rjust(6),' '.join(('%.1f'%r['first'] if r['first'] else '-').rjust(14) for r in R))
    for mk in marks:
        row=[]
        for r in R:
            v=[b for i,b in zip(r['it'],r['lb']) if i<=mk]
            row.append(('%.1f'%v[-1]).rjust(14) if v and r['it'][-1]>=mk else '-'.rjust(14))
        print(str(mk).rjust(6),' '.join(row))
    print('s/it'.rjust(6),' '.join(('%.3f'%(r['cpu']/r['it'][-1]) if r['it'] else '-').rjust(14) for r in R))
    if levels:
        for L in levels:
            row=[]
            for r in R:
                x=[i for i,b in zip(r['it'],r['lb']) if b>=L]
                row.append((str(x[0]) if x else '>%d'%(r['it'][-1] if r['it'] else 0)).rjust(14))
            print(('%.0f'%L).rjust(8),' '.join(row))
