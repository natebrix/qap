#!/usr/bin/env python3
"""down.py [DSTART] [DMAX]: down the tree under the strong c47 arm (README "Down the tree under c47").

Level D bounds every open depth-D node from the states of level D-1 (level 10 starts from runs/d9_*/out.state),
with production's settings (infra/phase1/recipes/tai60b.par) and production's own line rule, setting aside depth
D+1 (bnbdrv -d D). Each input state is split into one chunk per parent (the nodes sharing one DNN warm start,
splitstate-style), so a process holds the DNN states of one parent's branching children only; 4 chunks run at a
time. A chunk that finished (runs/down/dD/done/K) is not rerun, so the script resumes after a container restart.
After each level: out/down_levels.tsv gets a row (with production's under c54 at that depth), the README table
is rewritten from it, the directory is committed and pushed, and the kill rule is checked: stop if the level's
branching count reaches half of production's at that depth. Stops when a level has no open nodes (the subtree
closed) or at the wall-time limit (WALL hours since the first level started, default 10)."""
import os, sys, glob, subprocess, time, collections, re
H = os.path.dirname(os.path.abspath(__file__))
BUILD = '/home/user/wt-linetest'
PAR = f'{BUILD}/infra/phase1/recipes/tai60b.par'
DAT = f'{BUILD}/qaplib/tai60b.dat'
CUT = 608215054.05
NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39,'
JOBS = int(os.environ.get('JOBS', 4)); WALL = float(os.environ.get('WALL', 10)) * 3600
ENV = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
PROD = {}
for l in open(H + '/out/production_by_depth.txt'):
    d, n, b, c = l.split(); PROD[int(d)] = (int(n), int(b), float(c))

def split_by_warm(src, outdir, tag):
    """one chunk per warm start (a parent's children), plus one for nodes without a warm start"""
    head, warms, nodes, tail = [], [], [], []
    with open(src) as fh:
        for ln in fh:
            if ln.startswith('node '): nodes.append(ln)
            elif ln.startswith('warm '): warms.append(ln)
            elif ln.startswith(('warm_starts ', 'nodes ')): pass
            elif nodes or ln == 'end\n': tail.append(ln)
            else: head.append(ln)
    if not nodes: return []
    assert tail == ['end\n'], (src, tail[:3])
    groups = collections.defaultdict(list)
    for ln in nodes:
        body, w = ln.rsplit(None, 1); groups[int(w)].append(body)
    out = []
    for w, bodies in sorted(groups.items()):
        p = f'{outdir}/{tag}_w{w}.state'
        with open(p, 'w') as o:
            o.writelines(head)
            if w >= 0: o.write('warm_starts 1\n'); o.write(warms[w]); nw = 0
            else: o.write('warm_starts 0\n'); nw = -1
            o.write(f'nodes {len(bodies)}\n')
            for b in bodies: o.write(f'{b} {nw}\n')
            o.write('end\n')
        out.append(p)
    return out

def trace(paths):
    d = {}
    for p in paths:
        if not os.path.exists(p): continue
        for ln in open(p):
            f = ln.split()
            if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d

def lb(f):
    for i in (5, 4, 3):
        if f[i] != '-': return float(f[i])

def level_row(D):
    T = trace(glob.glob(f'{H}/runs/down/d{D}/nodes/*.txt'))
    L = [f for f in T.values() if int(f[1]) == D]
    c = collections.Counter(f[7] for f in L); cpu = sum(float(f[11]) for f in L)
    sf = [CUT - lb(f) for f in L if f[7] == 'B']
    return dict(depth=D, nodes=len(L), D=c['D'], C=c['C'], B=c['B'], other=len(L) - c['D'] - c['C'] - c['B'],
                cpu=cpu, sfsum=sum(sf), sfmax=max(sf) if sf else 0)

def write_table(rows):
    with open(H + '/out/down_levels.tsv', 'w') as w:
        w.write('depth\tnodes\tD\tC\tB\tother\tcpu\tsf_sum\tsf_max\tprod_nodes\tprod_B\tprod_cpu\n')
        for r in rows:
            pn, pb, pc = PROD.get(r['depth'], (0, 0, 0))
            w.write(f"{r['depth']}\t{r['nodes']}\t{r['D']}\t{r['C']}\t{r['B']}\t{r['other']}\t{r['cpu']:.0f}\t{r['sfsum']:.0f}\t{r['sfmax']:.0f}\t{pn}\t{pb}\t{pc:.0f}\n")
    lines = ['| depth | nodes bounded | D / C / B (other) | CPU s [m, this box] | branching shortfall: sum / largest | production under c54: nodes / B / CPU s (cloud boxes) |',
             '|---|---|---|---|---|---|',
             '| 9 (from above) | 260 | 247 / 4 / 9 | 9,002 | 1,152k / 254k | 520 / 56 / 80,387 |']
    for r in rows:
        pn, pb, pc = PROD.get(r['depth'], (0, 0, 0))
        oth = f" ({r['other']})" if r['other'] else ''
        lines.append(f"| {r['depth']} | {r['nodes']:,} | {r['D']:,} / {r['C']:,} / {r['B']:,}{oth} | {r['cpu']:,.0f} | {r['sfsum']/1e3:,.0f}k / {r['sfmax']/1e3:,.0f}k | {pn:,} / {pb:,} / {pc:,.0f} |")
    tab = '\n'.join(lines) + '\n'
    R = open(H + '/README.md').read()
    a, b = '<!-- down-table -->\n', '<!-- /down-table -->\n'
    if a in R:
        R = R[:R.index(a) + len(a)] + tab + R[R.index(b):]
        open(H + '/README.md', 'w').write(R)

def commit(msg):
    subprocess.run(f'cd {H}/../../.. && git add experiments/branching_lab/linetest && git commit -q -m "{msg}\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_014vvMYRPMz8BP8Xc3WigdYm" && '
                   f'for i in 1 2 3 4; do git push -q origin claude/cool-babbage-ilnogj-linetest && break; sleep $((i*5)); done', shell=True)

def log(s):
    with open(H + '/runs/down/down.log', 'a') as w: w.write(time.strftime('%H:%M:%S ', time.gmtime()) + s + '\n')
    print(s, flush=True)

def main():
    D0 = int(sys.argv[1]) if len(sys.argv) > 1 else 10; DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    os.makedirs(H + '/runs/down', exist_ok=True)
    t0f = H + '/runs/down/t0'
    if not os.path.exists(t0f): open(t0f, 'w').write(str(time.time()))
    t0 = float(open(t0f).read())
    rows = []
    if os.path.exists(H + '/out/down_levels.tsv'):
        for l in open(H + '/out/down_levels.tsv').read().split('\n')[1:]:
            if l:
                t = l.split('\t')
                if int(t[0]) < D0:
                    rows.append(dict(depth=int(t[0]), nodes=int(t[1]), D=int(t[2]), C=int(t[3]), B=int(t[4]), other=int(t[5]),
                                     cpu=float(t[6]), sfsum=float(t[7]), sfmax=float(t[8])))
    for D in range(D0, DMAX + 1):
        L = f'{H}/runs/down/d{D}'
        for s in ('in', 'out', 'nodes', 'done', 'log'): os.makedirs(f'{L}/{s}', exist_ok=True)
        srcs = sorted(glob.glob(f'{H}/runs/d9_*/out.state')) if D == 10 else sorted(glob.glob(f'{H}/runs/down/d{D-1}/out/*.state'))
        if not os.path.exists(f'{L}/split.ok'):
            for p in glob.glob(f'{L}/in/*.state'): os.remove(p)
            n = 0
            for s in srcs:
                tag = os.path.basename(os.path.dirname(s)) if D == 10 else os.path.basename(s)[:-6]
                n += len(split_by_warm(s, f'{L}/in', tag))
            open(f'{L}/split.ok', 'w').write(str(n))
        chunks = sorted(glob.glob(f'{L}/in/*.state'))
        if not chunks:
            log(f'depth {D}: no open nodes; the subtree has closed'); return 'closed', D
        todo = [c for c in chunks if not os.path.exists(f'{L}/done/' + os.path.basename(c))]
        log(f'depth {D}: {len(chunks)} chunks ({len(todo)} to run), {sum(int(open(c).read().split("nodes ",1)[1].split()[0]) for c in chunks)} nodes')
        running = {}
        while todo or running:
            while todo and len(running) < JOBS:
                c = todo.pop(0); k = os.path.basename(c)[:-6]
                e = dict(ENV, QAP_NODE_TRACE=f'{L}/nodes/{k}.txt')
                if os.path.exists(f'{L}/nodes/{k}.txt'): os.remove(f'{L}/nodes/{k}.txt')
                p = subprocess.Popen(['./bnbdrv', '-p', PAR, '--allow-settings-change', '-r', c, '-d', str(D), '-o', f'{L}/out/{k}.state', DAT],
                                     cwd=BUILD + '/latest', env=e, stdout=open(f'{L}/log/{k}.log', 'w'), stderr=open(f'{L}/log/{k}.err', 'w'))
                running[p] = c
            time.sleep(5)
            for p in list(running):
                if p.poll() is not None:
                    c = running.pop(p); k = os.path.basename(c)
                    if p.returncode == 0: open(f'{L}/done/{k}', 'w').write('0')
                    else: log(f'depth {D}: chunk {k} exit {p.returncode}'); return 'error', D
            if time.time() - t0 > WALL:
                for p in running: p.terminate()
                log(f'depth {D}: wall-time limit reached'); return 'wall', D
        r = level_row(D); rows = [x for x in rows if x['depth'] != D] + [r]; rows.sort(key=lambda x: x['depth'])
        write_table(rows)
        pn, pb, pc = PROD.get(D, (0, 0, 0))
        log(f"depth {D}: {r['nodes']} nodes, D {r['D']} C {r['C']} B {r['B']} other {r['other']}, CPU {r['cpu']:.0f} s; production {pn} / {pb} B / {pc:.0f} s")
        # free the disk: this level's inputs are no longer needed (the previous level's outputs are its inputs)
        for p in glob.glob(f'{L}/in/*.state'): os.remove(p)
        if D > 10:
            for p in glob.glob(f'{H}/runs/down/d{D-1}/out/*.state'): os.remove(p)
        commit(f"branching_lab/linetest: down the tree under c47, depth {D}: {r['nodes']} nodes, {r['B']} branching (production {pb})")
        if pb > 0 and r['B'] >= pb / 2:
            log(f'depth {D}: kill rule: {r["B"]} branching >= half of production\'s {pb}'); return 'kill', D
    return 'dmax', DMAX

if __name__ == '__main__':
    print(main())
