#!/usr/bin/env python3
"""rebuild.py NAME SRC.state FIXES: rebuild a node's state by bounding its ancestors one at a time.
FIXES is the full fix list of the last node to bound (depth 16 for the depth-17 test; any depth for the scan). Starting from SRC (a state holding the depth-10 ancestor set aside),
for D = first depth in SRC .. 16: keep the one depth-D node whose fixes agree with FIXES (and its warm start),
bound it with production's settings (bnbdrv -d D, its children set aside), and go on from that output. The depth-16
output (runs/rebuild/NAME/d16.out.state) holds that node's depth-17 children set aside with its DNN state.
Each bound is compared with the down.py trace (runs/down/dD/nodes): they should be the same to the digits printed."""
import os, sys, subprocess, glob
H = os.path.dirname(os.path.abspath(__file__))
BUILD = '/home/user/wt-linetest'
name, src, fixes = sys.argv[1], sys.argv[2], sys.argv[3]
want = dict(tuple(map(int, f.split(':'))) for f in fixes.split(','))
R = f'{H}/runs/rebuild/{name}'; os.makedirs(R, exist_ok=True)
ref = {}
for p in glob.glob(H + '/runs/down/d*/nodes/*.txt') + glob.glob(H + '/runs/d9_*/nodes.txt'):
    for ln in open(p):
        f = ln.split()
        if len(f) >= 12 and f[0] == 'node': ref[f[2]] = f
def cut(srcf, dst, depth):
    head, warms, nodes, tail = [], [], [], []
    for ln in open(srcf):
        if ln.startswith('node '): nodes.append(ln)
        elif ln.startswith('warm '): warms.append(ln)
        elif ln.startswith(('warm_starts ', 'nodes ')): pass
        elif nodes or ln == 'end\n': tail.append(ln)
        else: head.append(ln)
    keep = []
    for ln in nodes:
        t = ln.split()
        if int(t[1]) != depth: continue
        p = [int(x) for x in t[5:65]]
        fx = {r: v for r, v in enumerate(p) if v >= 0}
        if all(want.get(r) == v for r, v in fx.items()): keep.append(ln)
    assert len(keep) == 1, (srcf, depth, len(keep))
    body, w = keep[0].rsplit(None, 1); w = int(w)
    with open(dst, 'w') as o:
        o.writelines(head)
        if w >= 0: o.write('warm_starts 1\n'); o.write(warms[w]); w = 0
        else: o.write('warm_starts 0\n')
        o.write(f'nodes 1\n{body} {w}\n'); o.writelines(tail)
first = min(int(l.split()[1]) for l in open(src) if l.startswith('node '))
cur = src
LAST = len(want)            # FIXES names the node of depth LAST; its ancestors from SRC's depth to LAST are bounded
for D in range(first, LAST + 1):
    inp = f'{R}/d{D}.in.state'; out = f'{R}/d{D}.out.state'
    if not os.path.exists(out):
        cut(cur, inp, D)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', QAP_NODE_TRACE=f'{R}/nodes.txt')
        subprocess.run(['./bnbdrv', '-p', f'{BUILD}/infra/phase1/recipes/tai60b.par', '--allow-settings-change', '-r', inp, '-d', str(D),
                        '-o', out, f'{BUILD}/qaplib/tai60b.dat'], cwd=BUILD + '/latest', env=env, check=True,
                       stdout=open(f'{R}/d{D}.log', 'w'), stderr=open(f'{R}/d{D}.err', 'w'))
        os.remove(inp)
    mine = [l.split() for l in open(f'{R}/nodes.txt') if l.split()[1] == str(D)][-1]
    r = ref.get(mine[2])
    print(f'depth {D} {mine[2].split(",")[-1]}: {mine[7]} dnn {mine[4]} cut {mine[5]} line {mine[9]} | run: {r[7] if r else "?"} dnn {r[4] if r else "?"} cut {r[5] if r else "?"} line {r[9] if r else "?"}', flush=True)
    cur = out
print('done', cur)
