"""Assemble the nug30 break-even analysis from the result files (DNN jsonl + solver estimator logs)."""
import sys, json, re, os, glob
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn')
from qaputil import *
D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn30/'
FATH = 6124 - 1.95


def load(fn):
    return [json.loads(l) for l in open(D + fn)] if os.path.exists(D + fn) else []


def own_qpb(state):
    _, kids = read_state_nodes(D + state)
    o = {}
    for k in kids:
        o[parent_key(k['p'], k['hist'], k['depth'])] = k['bound']
    return o


def parse_est(fn):
    """Final estimator block: tmean, nmean, std errs, per-level time/nodes."""
    if not os.path.exists(fn):
        return None
    s = open(fn).read()
    m = re.search(r'tmean = ([\d.e+-]+)\s+nmean = ([\d.e+-]+)', s)
    if not m:
        return None
    r = dict(tmean=float(m.group(1)), nmean=float(m.group(2)))
    e = re.search(r'time stddev = ([\d.e+-]+), err = ([\d.e+-]+)', s); r['terr'] = float(e.group(2))
    e = re.search(r'node stddev = ([\d.e+-]+), err = ([\d.e+-]+)', s); r['nerr'] = float(e.group(2))
    # last per-level table
    blocks = s.split('level    nodes')
    tab = blocks[-1].split('Strategy breakdown')[0]
    lev = {}
    for line in tab.splitlines():
        t = line.split()
        if len(t) >= 5 and t[0].isdigit():
            lev[int(t[0])] = (float(t[1]), float(t[-1]))
    r['levels'] = lev
    return r


# ---------- depth-2 frontier and DNN outcome vs QPB ----------
q2 = own_qpb('n30_d2.state')
qv = np.array(sorted(q2.values()))
d2 = load('results_depth2.jsonl'); d2low = load('results_depth2_low.jsonl')
allrows = {}
for r in d2 + d2low:
    allrows[tuple(tuple(t) for t in r['key'])] = r
miss = [r for r in allrows.values() if r['lb'] <= FATH]
fath = [r for r in allrows.values() if r['lb'] > FATH]
hi_miss = max(r['qpb'] for r in miss); lo_fath = min(r['qpb'] for r in fath)
T = (hi_miss + lo_fath) / 2
print('depth-2: random sample %d (fathomed %d), lowest-QPB sample %d (fathomed %d); distinct nodes %d' % (
    len(d2), sum(r['lb'] > FATH for r in d2), len(d2low), sum(r['lb'] > FATH for r in d2low), len(allrows)))
print('highest-QPB miss %.1f, lowest-QPB fathom %.1f -> threshold T = %.1f' % (hi_miss, lo_fath, T))
n_surv = int((qv < T).sum()); n_f = len(qv) - n_surv
print('frontier: %d of 225 depth-2 nodes have QPB < T (predicted DNN survivors), %d predicted fathomed (%.1f%%)' % (n_surv, n_f, 100 * n_f / 225))
cpu_f = np.mean([r['cpu'] for r in fath]); cpu_m = np.mean([r['cpu'] for r in miss])
print('mean DNN cpu: fathomed %.0f s (n=%d), missed (run to convergence) %.0f s (n=%d)' % (cpu_f, len(fath), cpu_m, len(miss)))
dnn2_total = n_f * cpu_f + n_surv * cpu_m
print('DNN pass over all 225 depth-2 nodes: %.0f s = %.1f CPU-h' % (dnn2_total, dnn2_total / 3600))
gap = [r['obj'] - r['qpb'] for r in miss]
print('DNN value - QPB at missed depth-2 nodes: %.0f..%.0f (mean %.0f)' % (min(gap), max(gap), np.mean(gap)))

# ---------- depth 3 ----------
d3 = load('results_depth3.jsonl'); k3 = load('results_depth3_kids.jsonl'); kl3 = load('results_depth3_kidslow.jsonl')
for name, rs in [('random', d3), ('children of random-sample misses', k3), ('children of lowest-QPB misses', kl3)]:
    if rs:
        print('depth-3 %s: %d nodes, fathomed %d, cpu mean %.0f s, qpb %.0f..%.0f' % (
            name, len(rs), sum(r['lb'] > FATH for r in rs), np.mean([r['cpu'] for r in rs]),
            min(r['qpb'] for r in rs), max(r['qpb'] for r in rs)))
all3 = d3 + k3 + kl3
cpu3 = np.mean([r['cpu'] for r in all3])
q3 = own_qpb('n30_d3.state')
# depth-3 QPB-open children of the predicted depth-2 survivors
surv_keys = [k for k, v in q2.items() if v < T]
kids_of_surv = [v for k, v in q3.items() if any(set(k) >= set(sk) for sk in surv_keys)]
print('depth-3 QPB-open children of the %d predicted survivors: %d (of 6101); their QPB pct 0/25/50/75/100: %s' % (
    n_surv, len(kids_of_surv), np.percentile(kids_of_surv, [0, 25, 50, 75, 100]).round(0)))
m3 = [r for r in all3 if r['lb'] <= FATH]; f3 = [r for r in all3 if r['lb'] > FATH]
if m3:
    T3 = (max(r['qpb'] for r in m3) + min(r['qpb'] for r in f3)) / 2
else:
    T3 = min(r['qpb'] for r in f3)
    print('no depth-3 miss observed; lowest fathomed depth-3 QPB = %.1f' % T3)
n3_surv = int((np.array(kids_of_surv) < T3).sum())
print('depth-3 children of survivors with QPB < %.1f (possible depth-3 survivors, upper-ish bound): %d' % (T3, n3_surv))
dnn3_total = len(kids_of_surv) * cpu3
print('DNN pass over those depth-3 children: %d x %.0f s = %.0f s = %.1f CPU-h' % (len(kids_of_surv), cpu3, dnn3_total, dnn3_total / 3600))

# ---------- QPB cost estimates ----------
for fn in ['est/full_nug30.log', 'est/full_nug30_b.log']:
    e = parse_est(D + fn)
    if e:
        print('%s: est. nodes %.3g (+-%.2g), est. CPU %.3g s = %.0f CPU-days (+-%.0f); time at depth<=2: %.0f s; depth 3-4: %.0f s' % (
            fn, e['nmean'], e['nerr'], e['tmean'], e['tmean'] / 86400, e['terr'] / 86400,
            sum(e['levels'][l][1] for l in e['levels'] if l <= 2), sum(e['levels'][l][1] for l in e['levels'] if l in (3, 4))))
        print('   per-level nodes/time:', {l: ('%.3g' % e['levels'][l][0], '%.3g s' % e['levels'][l][1]) for l in sorted(e['levels'])})
full = parse_est(D + 'est/full_nug30.log')
sub_avg = (full['tmean'] - sum(full['levels'][l][1] for l in full['levels'] if l <= 2)) / 225
print('average QPB subtree below a depth-2 node (full estimate / 225): %.0f s = %.1f CPU-h' % (sub_avg, sub_avg / 3600))
print('average QPB subtree below a depth-3 node: %.0f s' % ((full['tmean'] - sum(full['levels'][l][1] for l in full['levels'] if l <= 3)) / 6101))
recs2 = load('results_depth2.jsonl')
for fn in sorted(glob.glob(D + 'est/est_d2_*.log')):
    i = int(re.search(r'est_d2_(\d+)', fn).group(1)); e = parse_est(fn)
    r = recs2[i]
    if e:
        print('reduced depth-2 node idx %d key %s qpb %.1f dnn %s: est. subtree nodes %.3g (+-%.2g), CPU %.3g s = %.1f CPU-h (+-%.1f)' % (
            i, r['key'], r['qpb'], 'MISS %.0f' % r['obj'] if r['lb'] <= FATH else 'fathomed', e['nmean'], e['nerr'], e['tmean'], e['tmean'] / 3600, e['terr'] / 3600))
    else:
        print('reduced depth-2 node idx %d: no final estimate in %s' % (i, fn))
