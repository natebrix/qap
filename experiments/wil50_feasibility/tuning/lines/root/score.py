import glob,re,math,sys,collections
S=sys.argv[1]; cutoff=48814.05
res=collections.defaultdict(list)
for p in glob.glob(S+'/children/f*_c*.log'):
    m=re.search(r'f(\d+)_c(\d+)\.log',p); f,c=int(m.group(1)),int(m.group(2))
    for l in open(p):
        if l.startswith('RESULT'):
            d=dict(kv.split('=') for kv in l.split()[1:]); res[f].append((c,float(d['lb']),d['status'],int(d['iters']),float(d['cpu'])))
for f in sorted(res, key=lambda f: sum(47**(max(0,cutoff-lb)/130) for _,lb,_,_,_ in res[f])):
    ch=sorted(res[f]); fath=sum(1 for _,lb,st,_,_ in ch if lb>=cutoff or st=='fathomed')
    gaps=sorted(max(0,cutoff-lb) for _,lb,_,_,_ in ch)
    cost=sum(47**(g/130) for g in gaps); cpu=sum(x[4] for x in ch)
    print(f"facility {f:2d}: {len(ch):2d} children, fathomed {fath:2d}, cost {cost:9.1f}, median gap {gaps[len(gaps)//2]:6.1f}, max gap {gaps[-1]:6.1f}, cpu {cpu:6.0f}")
    print("   gaps:", ' '.join(f"{g:.0f}" for g in gaps))
