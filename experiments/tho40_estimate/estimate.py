"""Tree-size estimate from QAP_DNN_TRACE files.

Each trace line is one DNN-bounded node: depth, fixes, first-solve lb, cut-phase
best lb, cutoff, fathomed flag. The tree is rebuilt from the fixes (a node's
parent is the node with one fix fewer), subtree sizes are summed bottom-up, and
log10(subtree size) is fitted against the relative gap (cutoff - lb)/cutoff,
per depth, on nodes whose subtree is complete. The fit from one instance is
then applied to another instance's traced nodes whose subtrees are incomplete.

usage: estimate.py fit <trace> [n]            # fit and print the model
       estimate.py predict <fit-trace> <n_fit> <target-trace> <n_target> <root_children>
       estimate.py both <traceA> <nA> <traceB> <nB>   # cross-validate A->B and B->A
"""
import sys, re, math, collections
import numpy as np

def read_trace(path):
    nodes = {}
    for l in open(path):
        m = re.match(r'pid \d+ depth (\d+) fixes((?: \d+:\d+)*) \| first lb ([-\d.]+) .*? \| cut (\d)(?: best ([-\d.]+))?.*\| cutoff ([-\d.]+) fathomed (\d)', l)
        if not m: continue
        depth = int(m.group(1)); fixes = frozenset(m.group(2).split())
        lb = float(m.group(3))
        if m.group(4) == '1' and m.group(5): lb = max(lb, float(m.group(5)))
        cutoff = float(m.group(6)); fath = m.group(7) == '1'
        if fixes in nodes:   # a node bounded twice (resume); keep the better
            o = nodes[fixes]; lb = max(lb, o['lb']); fath = fath or o['fath']
        nodes[fixes] = dict(depth=depth, lb=lb, cutoff=cutoff, fath=fath, children=[])
    # link parents
    for f, v in nodes.items():
        if v['depth'] <= 1: v['parent'] = None; continue
        v['parent'] = None
        for x in f:
            p = f - {x}
            if p in nodes: v['parent'] = p; break
    for f, v in nodes.items():
        if v['parent'] is not None: nodes[v['parent']]['children'].append(f)
    return nodes

def nchildren(n, depth, root_children):
    return root_children if depth == 0 else n - depth

def complete_sizes(nodes, n, root_children, finished=False):
    """Bottom-up subtree size; a node is 'complete' if fathomed or all its
    n-depth children are present and complete. With finished=True (a trace of
    a search that ran to the end) a missing child was fathomed by the QPB bound
    alone and is not traced: it counts as a leaf of size 1."""
    order = sorted(nodes, key=lambda f: -nodes[f]['depth'])
    for f in order:
        v = nodes[f]
        if v['fath']:
            v['size'] = 1; v['complete'] = True; continue
        kids = v['children']; nk = nchildren(n, v['depth'], root_children)
        missing = nk - len(kids)
        v['size'] = 1 + sum(nodes[k]['size'] for k in kids) + (missing if finished else 0)
        v['complete'] = finished or (missing == 0 and all(nodes[k]['complete'] for k in kids))
    return nodes

def fit(nodes, label=''):
    """Per-depth linear fit of log10(size) on relative gap, complete non-fathomed nodes.
    Also the fathomed fraction per depth as a function of gap is implicit (size 1)."""
    byd = collections.defaultdict(list)
    for v in nodes.values():
        if v['complete']:
            g = max(0.0, v['cutoff'] - v['lb']) / v['cutoff'] * 100
            byd[v['depth']].append((g, v['size'], v['fath']))
    model = {}
    print(f"-- fit {label}: depth, nodes (fathomed), gap% range, log10 size = a + b*gap%")
    for d in sorted(byd):
        pts = byd[d]; nf = sum(1 for _, _, f in pts if f)
        nonf = [(g, s) for g, s, f in pts if not f and s > 1]
        if len(nonf) >= 3:
            x = np.array([g for g, _ in nonf]); y = np.log10([s for _, s in nonf])
            b, a = np.polyfit(x, y, 1); r = np.corrcoef(x, y)[0, 1]
            model[d] = (a, b)
            print(f"   d={d}: {len(pts):6d} ({nf:6d} fathomed) gap {min(g for g,_,_ in pts):.2f}-{max(g for g,_,_ in pts):.2f}%  a={a:.3f} b={b:.3f} r={r:.2f}  (non-fathomed {len(nonf)})")
        else:
            print(f"   d={d}: {len(pts):6d} ({nf:6d} fathomed), too few open nodes to fit")
    # pooled fit with depth as a covariate, for depths missing from the model
    xs, ds, ys = [], [], []
    for d, pts in byd.items():
        for g, s, f in pts:
            if not f and s > 1: xs.append(g); ds.append(d); ys.append(math.log10(s))
    A = np.c_[np.ones(len(xs)), xs, ds]; coef, *_ = np.linalg.lstsq(A, np.array(ys), rcond=None)
    print(f"   pooled: log10 size = {coef[0]:.3f} + {coef[1]:.3f}*gap% + {coef[2]:.3f}*depth  ({len(xs)} nodes)")
    return model, coef

def model_size(model, coef, depth, gap_pct):
    if depth in model:
        a, b = model[depth]; return 10 ** (a + b * gap_pct)
    return 10 ** (coef[0] + coef[1] * gap_pct + coef[2] * depth)

def predict(nodes, n, root_children, model, coef, label=''):
    """Estimated whole-tree size: actual where complete, model where not.
    Missing children of an incomplete node get the parent's model size shared out."""
    order = sorted(nodes, key=lambda f: -nodes[f]['depth'])
    for f in order:
        v = nodes[f]
        if v['complete']: v['est'] = v['size']; continue
        g = max(0.0, v['cutoff'] - v['lb']) / v['cutoff'] * 100
        m = max(1.0, model_size(model, coef, v['depth'], g))
        kids = v['children']; nk = nchildren(n, v['depth'], root_children)
        missing = nk - len(kids)
        per_child = max(0.0, (m - 1) / nk)
        v['est'] = 1 + sum(nodes[k]['est'] for k in kids) + missing * per_child
        v['model_only'] = m
    roots = [f for f, v in nodes.items() if v['parent'] is None]
    opn = [(v['depth'], max(0.0, v['cutoff'] - v['lb']) / v['cutoff'] * 100) for v in nodes.values() if not v['complete']]
    if opn:
        byd = collections.defaultdict(list)
        for d, g in opn: byd[d].append(g)
        print("   open (model-estimated) nodes by depth: " + '; '.join(f"d={d}: {len(g)} nodes, gap {min(g):.2f}-{max(g):.2f}%" for d, g in sorted(byd.items())))
    # the QPB-only root is not traced: 1 + the depth-1 roots
    missing_roots = root_children - len(roots)
    mean_root = sum(nodes[f]['est'] for f in roots) / max(1, len(roots))
    tot = 1 + sum(nodes[f]['est'] for f in roots) + missing_roots * mean_root
    if missing_roots: print(f"   {missing_roots} of {root_children} depth-1 nodes not yet bounded: each counted at the mean of the traced ones, {mean_root:.3e}")
    traced = len(nodes)
    model_only = 1 + sum(max(1.0, model_size(model, coef, nodes[f]['depth'], max(0.0, nodes[f]['cutoff'] - nodes[f]['lb']) / nodes[f]['cutoff'] * 100)) for f in roots)
    print(f"-- predict {label}: traced {traced}, depth-1 roots {len(roots)}, estimate {tot:.3e} nodes (model alone from depth 1: {model_only:.3e})")
    # per-depth summary of what the estimate rests on
    byd = collections.defaultdict(lambda: [0, 0, 0])
    for v in nodes.values():
        byd[v['depth']][0] += 1; byd[v['depth']][1] += v['fath']; byd[v['depth']][2] += v['complete']
    for d in sorted(byd): print(f"   depth {d}: traced {byd[d][0]}, fathomed {byd[d][1]}, complete subtrees {byd[d][2]}")
    return tot

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'fit':
        tr, n, rc = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        nodes = complete_sizes(read_trace(tr), n, rc, finished=True); fit(nodes, tr)
        print("total traced", len(nodes), "root-complete:", all(nodes[f]['complete'] for f in nodes if nodes[f]['parent'] is None))
    elif cmd == 'predict':
        ftr, fn, frc, ttr, tn, trc = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6]), int(sys.argv[7])
        fit_finished = ftr != ttr          # the same unfinished trace as fit and target: complete subtrees only
        fnodes = complete_sizes(read_trace(ftr), fn, frc, finished=fit_finished); model, coef = fit(fnodes, ftr)
        finished_target = len(sys.argv) > 8   # truncation test on a finished run
        tnodes = complete_sizes(read_trace(ttr), tn, trc, finished=finished_target)
        # hide the target's completeness below a depth to test the model: --from-depth D
        if len(sys.argv) > 8:
            D = int(sys.argv[8])
            for v in tnodes.values():
                if v['depth'] >= D and not v['fath']: v['complete'] = False; v['children'] = []
            for v in tnodes.values():
                if v['depth'] > D: v['parent'] = None
            tnodes = {f: v for f, v in tnodes.items() if v['depth'] <= D}
            tnodes = complete_sizes(tnodes, tn, trc, finished=False)   # sizes and completeness from the truncated trace only
            print(f"   (target truncated to depth {D}: only nodes at depths <= {D} known, open ones estimated)")
        predict(tnodes, tn, trc, model, coef, ttr)
