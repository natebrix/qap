"""Shared loaders for the failure census (branching_lab/failures/census).

Trace format (latest/bnb.C node_trace):
  node <depth> <fixes> <qpb> <dnn> <cut> <cutoff> <outcome> <level> <line> <kids> <cpu>
The 'level' field is the par file's parameter LEVEL (pick_param), not the rule number.
RUNS maps each run's levels to the rule number from that run's own par file (see README,
"Which rule chose each line").
"""
import gzip, os, re
from collections import defaultdict

WT = '/home/user/wt-cognac'
LT = '/home/user/qap/experiments/branching_lab/linetest'
LAP = '/home/user/qap/experiments/branching_lab/laptop'

# level -> rule, and how it was established
RUNS = {
    'tai60b-261117': dict(inst='tai60b', files=[f'{WT}/infra/phase1/runs/tai60b-261117/nodes.txt.gz'],
                          rules={0: 9, 1: 9, 2: 9, 3: 2}, complete=True,
                          how='recipe.sh: commit 316361fc; git show 316361fc:infra/phase1/recipes/tai60b.par: levels 0-2 rule 9, bottom block rule 2'),
    'tai50b-261116': dict(inst='tai50b', files=[f'{WT}/infra/phase1/runs/tai50b-261116/nodes.txt.gz'],
                          rules={0: 9, 1: 9, 2: 9, 3: 2}, complete=True,
                          how='recipe.sh: commit 749b36b2; git show 749b36b2:infra/phase1/recipes/tai50b.par: levels 0-2 rule 9, bottom rule 2'),
    'tai50b-r9': dict(inst='tai50b', files=[f'{WT}/infra/phase1/runs/tai50b-r9/nodes.txt.gz'],
                      rules={0: 9, 1: 9, 2: 9, 3: 9}, complete=True,
                      how='recipe.sh: commit e406a047, PARAMS tai50b_r9.par: rule 9 in all four blocks'),
    'tai35b_base': dict(inst='tai35b', files=[f'{WT}/experiments/line_signals/runs/tai35b_base/nodes.txt.gz'],
                        rules={0: 9, 1: 9, 2: 9, 3: 2}, complete=True,
                        how='run.log: params experiments/line_signals/tai35b.par (= tai/tai35b.par), commit 8df95737: levels 0-2 rule 9, bottom rule 2'),
    'tai40b_base': dict(inst='tai40b', files=[f'{WT}/experiments/line_signals/runs/tai40b_base/nodes.txt.gz'],
                        rules={0: 9, 1: 9, 2: 9, 3: 2}, complete=True,
                        how='run.log: params experiments/line_signals/tai40b.par, commit 8df95737: levels 0-2 rule 9, bottom rule 2'),
    'tho40ab-32': dict(inst='tho40', files=[f'{WT}/infra/phase1/runs/tho40ab-32/nodes.txt.gz'],
                       rules={0: 9, 1: 9, 2: 9, 3: 2}, complete=True,
                       how='launches.tsv: commit f27dd94, recipe tho40ab-ship32 -> tho40-ship32.par at f27dd94: levels 0-2 rule 9, bottom rule 2'),
    'nug30-trace-a': dict(inst='nug30', files=[f'{WT}/infra/phase1/runs/nug30-trace-a/nodes.txt.gz'],
                          rules={0: 9, 1: 9, 2: 9, 3: 2}, complete=True,
                          how='launches.tsv: commit c79368d, recipe nug30trace -> nug30.par at c79368d: levels 0-2 rule 9, bottom rule 2'),
    # the c47 descent: depth 7 forced to c47, depths 8-17 lines by rule 9 (tai60b.par at 0e6c2ecf)
    'tai60b-c47': dict(inst='tai60b', files=None, rules={0: 9, 1: 9, 2: 9, 3: 9}, complete=False,
                       how='linetest README "Correction": tai60b.par at 0e6c2ecf (bottom block rule 9); depth 7 forced c47; depth 18 partly bounded'),
    # the landscape-signal run: signal where ls.log says "chose", rule 9 elsewhere; stopped part-way
    'tai50b-ls': dict(inst='tai50b', files=[f'{LAP}/runs/tai50b-ls/nodes.txt.gz'], rules={0: 9, 1: 9, 2: 9, 3: 9},
                      complete=False,
                      how='laptop/runs/tai50b-ls/README.md: recipes/tai50b.par at 0e6c2ecf (rule 9 everywhere); QAP_LS=a:2:0:10: the signal chose where ls.log says "chose"'),
}


def c47_files():
    fs = [f'{LT}/runs/sc47/nodes.txt']
    for d in ['d9_29_47', 'd9_45_47_31_47', 'd9_46_47', 'd9_47_47']:
        fs.append(f'{LT}/runs/{d}/nodes.txt')
    for D in range(10, 19):
        nd = f'{LT}/runs/down/d{D}/nodes'
        if os.path.isdir(nd):
            fs += [os.path.join(nd, f) for f in sorted(os.listdir(nd))]
    return fs


RUNS['tai60b-c47']['files'] = c47_files()


def opener(f):
    return gzip.open(f, 'rt') if f.endswith('.gz') else open(f)


def fval(s):
    return None if s == '-' else float(s)


class Node:
    __slots__ = ('depth', 'fixes', 'key', 'qpb', 'dnn', 'cut', 'cutoff', 'out', 'level', 'line', 'kids', 'cpu',
                 'kidsl', 'parent', 'sub_n', 'sub_cpu', 'b', 'dn', 'dc', 'ex', '__dict__')

    def __init__(self, t):
        self.depth = int(t[1])
        self.fixes = () if t[2] == '-' else tuple(tuple(int(x) for x in p.split(':')) for p in t[2].split(','))
        self.key = t[2]
        self.qpb, self.dnn, self.cut, self.cutoff = fval(t[3]), fval(t[4]), fval(t[5]), fval(t[6])
        self.out = t[7]
        self.level = None if t[8] == '-' else int(t[8])
        self.line = None if t[9] == '-' else t[9]
        self.kids = None if t[10] == '-' else int(t[10])
        self.cpu = float(t[11])
        self.kidsl = []
        self.parent = None
        self.sub_n = 0
        self.sub_cpu = 0.0

    @property
    def bound(self):
        v = [x for x in (self.qpb, self.dnn, self.cut) if x is not None]
        return max(v) if v else None

    @property
    def short(self):
        """cutoff minus the node's best bound (positive: not fathomed by that margin)."""
        b = self.bound
        return None if (b is None or self.cutoff is None) else self.cutoff - b


def load(run):
    """Return (nodes dict key->Node, root) with children linked and subtree totals (descendants only)."""
    nodes = {}
    for f in RUNS[run]['files']:
        with opener(f) as fh:
            for line in fh:
                t = line.split()
                if not t or t[0] != 'node' or len(t) < 12:
                    continue
                n = Node(t)
                if n.key in nodes:
                    # duplicate (a node bounded twice, e.g. a rerun chunk): keep the later line
                    pass
                nodes[n.key] = n
    for n in nodes.values():
        if n.fixes:
            pk = ','.join(f'{a}:{b}' for a, b in n.fixes[:-1]) or '-'
            p = nodes.get(pk)
            if p is not None:
                n.parent = p
                p.kidsl.append(n)
    # subtree totals, deepest first
    for n in sorted(nodes.values(), key=lambda x: -x.depth):
        if n.parent is not None:
            n.parent.sub_n += 1 + n.sub_n
            n.parent.sub_cpu += n.cpu + n.sub_cpu
    tops = [n for n in nodes.values() if n.parent is None]
    return nodes, tops


def read_dat(inst):
    toks = open(f'{WT}/qaplib/{inst}.dat').read().split()
    n = int(toks[0])
    v = [int(x) for x in toks[1:1 + 2 * n * n]]
    A = [v[i * n:(i + 1) * n] for i in range(n)]
    B = [v[n * n + i * n:n * n + (i + 1) * n] for i in range(n)]
    return n, A, B


def read_sln(inst):
    toks = open(f'{WT}/QAPSOLN/{inst}.sln').read().split()
    n, obj = int(toks[0]), int(toks[1])
    p = [int(x) - 1 for x in toks[2:2 + n]]
    return obj, p


def rule_of(run, n):
    if n.level is None:
        return None
    return RUNS[run]['rules'][n.level]
