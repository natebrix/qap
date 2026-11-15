"""A node trace (QAP_NODE_TRACE, format in latest/bnb.C) as a tree, for the
analyses of docs/tho40b_results_plan.md.

    node <depth> <fixes> <qpb> <dnn> <cut> <cutoff> <outcome> <rule> <line> <kids> <cpu>

fixes are facility:location pairs in branching order, so a node's parent is
the node whose fixes are its own less the last pair. The root is the node
whose parent is not traced: the root of a whole run, the -F node of a subtree
run (its frontier pass bounds nothing). Outcomes: Q, D, C fathomed by QPB, the
DNN bound, the cut phase; P discarded when popped; E enumerated; B branched,
with kids children.

load() returns a Tree: one entry per node in trace order, with the subtree
totals (nodes, CPU, deepest depth) summed bottom-up over the traced nodes.
A complete trace (qapctl export, qapctl nodes, a finished run) has every
child of every branched node; a prefix of one (a run part-way) has not, and
Tree.complete says which subtrees are whole.
"""
import gzip
import math

FATHOMED = set("QDCP")


def _num(tok):
    return math.nan if tok == "-" else float(tok)


def open_trace(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


class Tree:
    def __init__(self, lines):
        self.depth, self.fixes, self.qpb, self.dnn, self.cut, self.cutoff = [], [], [], [], [], []
        self.outcome, self.kids, self.cpu, self.line = [], [], [], []
        for ln in lines:
            f = ln.split()
            if len(f) < 12 or f[0] != "node":
                continue
            self.depth.append(int(f[1]))
            self.fixes.append(() if f[2] == "-" else tuple(f[2].split(",")))
            self.qpb.append(_num(f[3]))
            self.dnn.append(_num(f[4]))
            self.cut.append(_num(f[5]))
            self.cutoff.append(_num(f[6]))
            self.outcome.append(f[7])
            self.line.append(f[9])
            self.kids.append(0 if f[10] == "-" else int(f[10]))
            self.cpu.append(float(f[11]))
        n = len(self.depth)
        index = {fx: i for i, fx in enumerate(self.fixes)}
        self.parent = [index.get(fx[:-1], -1) if fx else -1 for fx in self.fixes]
        self.children = [[] for _ in range(n)]
        for i, p in enumerate(self.parent):
            if p >= 0:
                self.children[p].append(i)
        self.roots = [i for i, p in enumerate(self.parent) if p < 0]
        # subtree totals, deepest first
        self.size = [1] * n
        self.sub_cpu = list(self.cpu)
        self.max_depth = list(self.depth)
        self.complete = [False] * n
        for i in sorted(range(n), key=lambda i: -self.depth[i]):
            ch = self.children[i]
            self.complete[i] = self.outcome[i] != "B" or (len(ch) == self.kids[i] and all(self.complete[c] for c in ch))
            p = self.parent[i]
            if p >= 0:
                self.size[p] += self.size[i]
                self.sub_cpu[p] += self.sub_cpu[i]
                self.max_depth[p] = max(self.max_depth[p], self.max_depth[i])

    def __len__(self):
        return len(self.depth)

    def lb(self, i):
        """the best bound computed at the node"""
        vals = [v for v in (self.qpb[i], self.dnn[i], self.cut[i]) if not math.isnan(v)]
        return max(vals) if vals else math.nan

    def shortfall(self, i):
        """cutoff - bound: how far the node's best bound fell short of fathoming it"""
        return self.cutoff[i] - self.lb(i)

    def gap_pct(self, i):
        """the shortfall as a percentage of the cutoff (estimate.py's relative gap)"""
        return max(0.0, self.shortfall(i)) / self.cutoff[i] * 100

    def fathomed(self, i):
        return self.outcome[i] in FATHOMED

    def ancestors(self, i):
        p = self.parent[i]
        while p >= 0:
            yield p
            p = self.parent[p]

    def fixes_flag(self, i):
        """the node as a -F spec (bnbdrv -F i:j,...)"""
        return ",".join(self.fixes[i])


def load(path, limit=None):
    """the tree of a node-trace file (.gz or plain); limit: only its first lines"""
    with open_trace(path) as f:
        lines = []
        for ln in f:
            lines.append(ln)
            if limit is not None and len(lines) >= limit:
                break
    return Tree(lines)
