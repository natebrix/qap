"""Replay the pair and one-fixed exchange-dominance tests over every node of a solver trace.
usage: replay.py <label> <instance> <trace> <out.tsv.gz>
Per node: key, depth, outcome, cpu, best pair min Delta (and the pair), best one-fixed value (and (k,i)), test seconds."""
import sys, gzip, time, numpy as np
from qapdata import read_dat, read_trace
from dom import Inst, onefixed_fast
label, name, trace, out = sys.argv[1:5]
inst = Inst(*read_dat(name))
with gzip.open(out, 'wt') as f:
    f.write('key\tdepth\toutcome\tcpu\tpairmin\tpair\tonefix\tonefix_ki\tt_pair\tt_one\n')
    t0 = time.time()
    for idx, (d, fx, oc, line, cpu, ln) in enumerate(read_trace(trace)):
        F = dict(fx); key = ln.split()[2]
        s = time.perf_counter(); pm, pr = inst.pair_test(F); tp = time.perf_counter() - s
        s = time.perf_counter(); of, ki = onefixed_fast(inst, F); to = time.perf_counter() - s
        f.write(f"{key}\t{d}\t{oc}\t{cpu}\t{'' if pm is None else pm}\t{'' if pr is None else '%d,%d' % (int(pr[0]), int(pr[1]))}\t"
                f"{'' if of is None else int(round(of))}\t{'' if ki is None else '%d,%d' % ki}\t{tp:.6f}\t{to:.6f}\n")
        if idx % 20000 == 0: print(label, idx, f'{time.time()-t0:.0f}s', flush=True)
print(label, 'done', f'{time.time()-t0:.0f}s')
