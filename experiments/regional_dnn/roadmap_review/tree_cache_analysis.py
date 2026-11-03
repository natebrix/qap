import json,pathlib,math,sys,numpy as np
from scipy.stats import spearmanr
D=pathlib.Path(sys.argv[1])/'experiments/wil50_feasibility/tuning/lines/work'
parents={1:{34:40,28:22,14:29},3:{27:44,5:16,17:9},10:{39:27,30:48,25:49}}
allrows=[]
for pid,fix in parents.items():
 rows=json.loads((D/f'report_{pid}.json').read_text());qmap=json.loads((D/f'qpb_{pid}.json').read_text()); pdata=json.loads((D/f'parent_{pid}.json').read_text()); pLB=pdata['lb']
 print('\nPARENT',pid)
 for row in rows:
  kind,x=row['line'].split();x=int(x); ch=[]
  for y in range(50):
   i,j=(x,y) if kind=='row' else (y,x)
   if i in fix or j in fix.values():continue
   f=dict(fix);f[i]=j;s=','.join(f'{a}:{b}' for a,b in sorted(f.items()))
   r=json.loads((D/f"child_{s.replace(':','-').replace(',','_')}.json").read_text());ch.append((max(0,48814.05-r['lb']),r['fathomed'],qmap[s],s,r['cpu']))
  row.update(pid=pid,ch=ch)
  for g in [100,130,160]:
   for cut in [0,80,130]:row[f'C{g}_{cut}']=sum(0 if f or gap<=cut else 47**((gap-cut)/g) for gap,f,q,s,t in ch)
  row['qtail']=sum(math.exp(-(q-47000)/100) for gap,f,q,s,t in ch)
  row['qmin']=min(c[2] for c in ch)
  row['qbottom5']=np.mean(sorted(c[2] for c in ch)[:5])
  row['maxgap']=max(c[0] for c in ch)
  row['costtop5share']=sum(sorted((47**(c[0]/130) if not c[1] else 0 for c in ch),reverse=True)[:5])/row['C130_0']
  print(row['line'],'fath',sum(c[1] for c in ch),'maxgap',round(row['maxgap'],1),'top5share',round(row['costtop5share'],2),'C130cuts',[round(row[f'C130_{cut}']) for cut in [0,80,130]],'qbottom5',round(row['qbottom5'],1))
 for key in ['mean child QPB','qmin','qbottom5','qtail']:
  chosen=(min if key=='qtail' else max)(rows,key=lambda r:r[key]);print('CHOOSE',key,chosen['line'],'regret',round(chosen['C130_0']/min(r['C130_0'] for r in rows),2),'rho',round(spearmanr([r[key] for r in rows],[r['C130_0'] for r in rows]).statistic,2))
 for g in [100,130,160]:
  for cut in [0,80,130]:
   key=f'C{g}_{cut}';best=min(rows,key=lambda r:r[key]);print('BEST',g,cut,best['line'],round(best[key],1))
 # Partition using QPB-predicted easy sets: how many actually fathom, rank disagreements.
 best=min(rows,key=lambda r:r['C130_0']);ch=sorted(best['ch'],key=lambda c:-c[2])
 print('BESTLINE_QPB_EASY_PREFIX',[(k,sum(c[1] for c in ch[:k]),round(max(c[0] for c in ch[:k]),1)) for k in [8,16,24,32]])
 print('HARDEST_BESTLINE',[(round(c[0],1),c[3]) for c in sorted(best['ch'],reverse=True)[:6]])
 allrows+=rows
# Leave-one-parent-out child QPB -> gap OLS, with parent fixed effect by centering gap,QPB within parent.
for test in [1,3,10]:
 train=[r for r in allrows if r['pid']!=test]; testrows=[r for r in allrows if r['pid']==test]
 pairs=[]
 for pid in parents:
  rr=[r for r in train if r['pid']==pid]
  if not rr:continue
  a=np.array([(c[2],c[0]) for r in rr for c in r['ch']]); a-=a.mean(axis=0);pairs.extend(a)
 a=np.array(pairs);b=(a[:,0]@a[:,1])/(a[:,0]@a[:,0]);
 scored=[]
 for r in testrows:scored.append((sum(math.exp(b*(c[2]-47000)*math.log(47)/130) for c in r['ch']),r))
 chosen=min(scored,key=lambda z:z[0])[1]
 print('LOPO_QEXP',test,'slope',b,'pick',chosen['line'],'regret',chosen['C130_0']/min(r['C130_0'] for r in testrows))
