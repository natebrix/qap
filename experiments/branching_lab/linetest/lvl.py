#!/usr/bin/env python3
"""lvl.py NAME DEPTH PAR SRC.state [SRC.state ...]: bound every open depth-DEPTH node of the given states under the
params file PAR, setting aside depth DEPTH+1 (bnbdrv -d DEPTH). The states are split as down.py splits them (a
parent's children in pieces of at most 13, each with a copy of the parent's warm start) and run JOBS (default 4) at
a time. Traces in runs/lvl/NAME/nodes/*.txt, set-aside states in runs/lvl/NAME/out/. Finished pieces are kept on a
rerun. Used for the 2x2 at depth 9 (README)."""
import os, sys, glob, subprocess, time
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import down
name, D, par, srcs = sys.argv[1], int(sys.argv[2]), os.path.abspath(sys.argv[3]), sys.argv[4:]
L = f'{H}/runs/lvl/{name}'
for s in ('in', 'out', 'nodes', 'done', 'log'): os.makedirs(f'{L}/{s}', exist_ok=True)
if not os.path.exists(f'{L}/split.ok'):
    for i, s in enumerate(srcs): down.split_by_warm(s, f'{L}/in', f's{i}')
    open(f'{L}/split.ok', 'w').write('1')
todo = [c for c in sorted(glob.glob(f'{L}/in/*.state')) if not os.path.exists(f'{L}/done/' + os.path.basename(c))]
print(name, len(todo), 'pieces to run', flush=True)
running = {}
while todo or running:
    while todo and len(running) < down.JOBS:
        c = todo.pop(0); k = os.path.basename(c)[:-6]
        e = dict(down.ENV, QAP_NODE_TRACE=f'{L}/nodes/{k}.txt')
        if os.path.exists(f'{L}/nodes/{k}.txt'): os.remove(f'{L}/nodes/{k}.txt')
        p = subprocess.Popen(['./bnbdrv', '-p', par, '--allow-settings-change', '-r', c, '-d', str(D), '-o', f'{L}/out/{k}.state', down.DAT],
                             cwd=down.BUILD + '/latest', env=e, stdout=open(f'{L}/log/{k}.log', 'w'), stderr=open(f'{L}/log/{k}.err', 'w'))
        running[p] = c
    time.sleep(5)
    for p in list(running):
        if p.poll() is not None:
            c = running.pop(p)
            if p.returncode == 0: open(f'{L}/done/' + os.path.basename(c), 'w').write('0')
            else: print('piece', c, 'exit', p.returncode, flush=True); sys.exit(1)
print(name, 'done', flush=True)
