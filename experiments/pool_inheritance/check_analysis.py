"""Small prelaunch checks for the queue, log parser and declared decision gates."""
from pathlib import Path
from run import plan, parse, CUTOFF
from summarize import gate, evaluate, first_target, prefix_cpu_interval

p=plan()
assert p['sampled_facilities']==[5,10,11]
assert len(p['jobs'])==len({(j['node'],j['arm']) for j in p['jobs']})==19
assert sum(j['dependency'] is not None for j in p['jobs'])==12
assert gate([41,42,43,44,45,46],48760)['decision']=='PASS'
assert gate([30,31,32,33,34,35],48760)['decision']=='KILL'
assert gate([41,42,43,44,45,46],48759)['decision']=='KILL'
assert gate([40]*6,48760)['decision']=='ROUNDING_AMBIGUOUS'
try:evaluate({'jobs':[]},p)
except ValueError:pass
else:raise AssertionError('Incomplete matrix accepted')
job=dict(rounds=[dict(round='round  0',lb=48780,status='iter_limit',cumulative_iterations=400,dnn_cpu=4),
                 dict(round='round  1',lb=48795,status='iter_limit',cumulative_iterations=700,dnn_cpu=6)],
         process_cpu_seconds=12,
         stage_observations=[dict(observed_process_cpu_seconds=5),dict(observed_process_cpu_seconds=9)])
assert first_target(job,48789.9)['iterations']==700
assert first_target(job,CUTOFF) is None and first_target(job,None) is None
interval=prefix_cpu_interval(job,0)
assert interval['lower_seconds']<4 and interval['upper_seconds']>6
base=Path(__file__).parent/'results/preflight'
if (base/'smoke-load.log').exists():
    q=parse(base/'smoke-load.log')
    assert len(q['rounds'])==1 and q['rounds'][0]['round']=='quiet round' and q['iterations']==30
    assert q['pool']['kept']==69 and q['pool']['dropped']==11
    old=parse(base/'smoke-old.log');new=parse(base/'smoke-new.log');saved=parse(base/'smoke-save.log')
    assert old['final']==new['final']==saved['final']
print('PASS queue=19 unique runs, sample/dependencies, both KILL conditions, PASS, rounding ambiguity, incomplete-matrix rejection, target-stage accounting, driver smoke parsing')
