"""The live tree of a --trace run, for brixius-web's tree viewer: every node
the run has merged, rebuilt from the node traces (QAP_NODE_TRACE, format in
latest/bnb.C) and published under public/runs/<run>/tree/ with the status.

The nodes come in a fixed order: the frontier pass's trace, then each merged
unit's (its recorded attempt's nodes.txt, as qapctl nodes reads them) by merge
sequence. A node's id is its place in that order, so ids never change while
the run goes on, and a restarted coordinator, which reads them all again,
gives every node the same id. A node's parent is the node whose fixes are its
own less the last pair (the trace lists them in branching order).

The tree is cut into shards the way brixius-web's viewer reads them
(lib/qap/tree.ts: Shard, SHARD_COLUMNS, rowToNode): a shard is a node at an
even depth, its children and its grandchildren, and the whole subtree of each
grandchild with at most INLINE nodes; a larger grandchild is marked as a
boundary, and its children are in its own shard, named by its id. So a node's
children are always in one shard, and a new node changes only the shards of
its even-depth ancestors (their subtree totals), which are the only ones
published again. A subtree only grows, so a grandchild once a boundary stays
one, and a rebuild cuts the tree as the run did. On nug17 (41,443 nodes) that
is 285 shards, 103 rows at the median, 1,667 at most. meta.json (TreeMeta in
lib/qap/tree.ts, with layout, generatedAt and final) is written last.
"""
import json
import math
import time
from array import array
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

from .common import run_key, unit_key

FATHOMED = frozenset(b"QDCP")
OUTCOMES = b"QDCPEB"
SHARD_COLUMNS = ["id", "parent", "depth", "fix", "outcome", "lb", "qpb", "dnn", "cut", "cutoff", "iterations", "cpu",
                 "status", "rule", "line", "kids", "size", "fathomed", "subCpu", "maxDepth", "shard"]
CACHE = "max-age=60"
INLINE = 50     # a grandchild's subtree up to this many nodes is in its grandparent's shard
NAN = float("nan")


def _value(tok):
    if tok == b"-":
        return NAN
    try:
        return float(tok)
    except ValueError:
        return NAN


def _round(v, places):
    return None if math.isnan(v) else float(f"{v:.{places}f}")


class LiveTree:
    """the merged nodes, added a trace at a time (ingest), with each subtree's
    totals kept up to date and the shards a change touches remembered"""

    def __init__(self, prefix=()):
        # a subtree run (-F): every node's fixes start with the -F node's, which is the root
        self.prefix = tuple(prefix)
        self.par = array("i", [-1])
        self.dep = array("h", [len(self.prefix)])
        self.fix = array("i", [-1])             # facility << 10 | location; -1 at the root
        self.rec = array("b", [0])              # the node's own line read (only the root starts without)
        self.out = array("B", [ord("B")])
        self.qpb, self.dnn, self.cut, self.cutoff = (array("d", [NAN]) for _ in range(4))
        self.rule = array("h", [-1])
        self.line = array("h", [0])             # r<k>: k + 1, c<k>: -(k + 1), none: 0
        self.kids = array("h", [-1])
        self.cpu = array("d", [0.0])
        self.size = array("i", [1])
        self.fat = array("i", [0])
        self.sub = array("d", [0.0])
        self.mx = array("h", [len(self.prefix)])
        self.first = array("i", [-1])           # first child, next sibling
        self.nxt = array("i", [-1])
        self.child = {}                         # parent << 20 | facility << 10 | location -> id
        self.waiting = {}                       # a parent's fixes text -> lines that came before it
        self.duplicates = 0
        self.outcomes = {}
        self.dirty = {0}

    def __len__(self):
        return len(self.par)

    def ingest(self, data):
        """add a trace's node lines (bytes); other lines are skipped"""
        for line in data.split(b"\n"):
            self._line(line)

    def _line(self, line):
        f = line.split()
        if len(f) < 12 or f[0] != b"node" or len(f[7]) != 1 or f[7][0] not in OUTCOMES:
            return
        fixes = [] if f[2] == b"-" else [tuple(int(x) for x in p.split(b":")) for p in f[2].split(b",")]
        if self.prefix:
            if tuple(fixes[:len(self.prefix)]) != self.prefix:
                return                      # not in the subtree (not written by a -F run)
            fixes = fixes[len(self.prefix):]
        if not fixes:
            if self.rec[0]:
                self.duplicates += 1
            else:
                self._fill(0, f)
                self.dirty.add(0)
            return
        pid = 0
        for a, b in fixes[:-1]:
            pid = self.child.get(pid << 20 | a << 10 | b)
            if pid is None:
                self.waiting.setdefault(f[2].rsplit(b",", 1)[0], []).append(line)
                return
        a, b = fixes[-1]
        k = pid << 20 | a << 10 | b
        if k in self.child:
            self.duplicates += 1
            return
        i = len(self.par)
        self.child[k] = i
        self.par.append(pid)
        self.dep.append(int(f[1]))
        self.fix.append(a << 10 | b)
        for col, v in ((self.rec, 0), (self.out, 0), (self.qpb, NAN), (self.dnn, NAN), (self.cut, NAN),
                       (self.cutoff, NAN), (self.rule, -1), (self.line, 0), (self.kids, -1), (self.cpu, 0.0),
                       (self.size, 1), (self.fat, 0), (self.sub, 0.0), (self.first, -1)):
            col.append(v)
        self.mx.append(int(f[1]))
        self.nxt.append(self.first[pid])
        self.first[pid] = i
        self._fill(i, f)
        x, d = pid, self.dep[i]
        while x != -1:
            self.size[x] += 1
            self.fat[x] += self.fat[i]
            self.sub[x] += self.cpu[i]
            if self.mx[x] < d:
                self.mx[x] = d
            if self.dep[x] % 2 == 0:
                self.dirty.add(x)
            x = self.par[x]
        later = self.waiting.pop(f[2], None)
        for w in later or ():
            self._line(w)

    def _fill(self, i, f):
        o = f[7][0]
        self.rec[i] = 1
        self.out[i] = o
        self.qpb[i], self.dnn[i], self.cut[i], self.cutoff[i] = (_value(t) for t in f[3:7])
        if o == ord("B"):
            self.rule[i] = -1 if f[8] == b"-" else int(f[8])
            ln = f[9].decode()
            self.line[i] = (int(ln[1:]) + 1) * (1 if ln[0] == "r" else -1) if ln[:1] in ("r", "c") else 0
            self.kids[i] = -1 if f[10] == b"-" else int(f[10])
        cpu = _value(f[11])
        self.cpu[i] = 0.0 if math.isnan(cpu) else cpu
        # += : the root's line may come after its children's (a frontier of depth 0)
        self.fat[i] += 1 if o in FATHOMED else 0
        self.sub[i] += self.cpu[i]
        self.outcomes[chr(o)] = self.outcomes.get(chr(o), 0) + 1

    # ---------------------------------------------------------------- reading
    def children(self, i):
        out, c = [], self.first[i]
        while c != -1:
            out.append(c)
            c = self.nxt[c]
        # by location, then facility (lib/qap/tree.ts assemble: siblings share a row or a column)
        out.sort(key=lambda c: (self.fix[c] & 1023, self.fix[c] >> 10))
        return out

    def fix_text(self, i):
        return "" if self.fix[i] < 0 else f"{self.fix[i] >> 10}:{self.fix[i] & 1023}"

    def row(self, i, boundary):
        qpb, dnn, cut = self.qpb[i], self.dnn[i], self.cut[i]
        bounds = [v for v in (qpb, dnn, cut) if not math.isnan(v)]
        ln = self.line[i]
        return [i, self.par[i], self.dep[i], self.fix_text(i), chr(self.out[i]),
                _round(max(bounds), 2) if bounds else None, _round(qpb, 2), _round(dnn, 2), _round(cut, 2),
                _round(self.cutoff[i], 2), 0, _round(self.cpu[i], 4), None,
                None if self.rule[i] < 0 else self.rule[i],
                None if ln == 0 else (f"r{ln - 1}" if ln > 0 else f"c{-ln - 1}"),
                None if self.kids[i] < 0 else self.kids[i],
                self.size[i], self.fat[i], _round(self.sub[i], 2), self.mx[i], 1 if boundary else 0]

    def is_root(self, i):
        """a shard of its own: the root, or an even depth with a subtree past INLINE"""
        return i == 0 or (self.dep[i] % 2 == 0 and self.first[i] != -1 and self.size[i] > INLINE)

    def shard(self, r):
        """node r, its children, its grandchildren and the small ones' subtrees"""
        rows = [self.row(r, False)]
        kids = self.children(r)
        rows += [self.row(c, False) for c in kids]
        for c in kids:
            for g in self.children(c):
                rows.append(self.row(g, self.is_root(g)))
                if not self.is_root(g):
                    stack = self.children(g)[::-1]
                    while stack:
                        x = stack.pop()
                        rows.append(self.row(x, False))
                        stack += self.children(x)[::-1]
        return {"root": r, "columns": SHARD_COLUMNS, "rows": rows}

    def shard_roots(self):
        return [i for i in range(len(self.par)) if self.is_root(i)]

    def take_dirty(self):
        d, self.dirty = self.dirty, set()
        return sorted(i for i in d if self.is_root(i))

    def outcome_counts(self):
        out = dict(self.outcomes)
        if not self.rec[0]:
            out["B"] = out.get("B", 0) + 1     # the root, before its own line (as lib/qap/tree.ts)
        return out

    def cutoff_value(self):
        for i in range(min(len(self.par), 1000)):
            if not math.isnan(self.cutoff[i]):
                return _round(self.cutoff[i], 2)
        return None


def grid_of(n, *matrices):
    """the rows x cols grid whose Manhattan distances one of the matrices is
    (location l at row l // cols, column l % cols), or None"""
    for m in matrices:
        if not m:
            continue
        for cols in range(2, n + 1):
            if all(m[i][j] == abs(i // cols - j // cols) + abs(i % cols - j % cols)
                   for i in range(n) for j in range(n)):
                return {"rows": -(-n // cols), "cols": cols}
    return None


def iso(t):
    return datetime.fromtimestamp(t, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class TreePublisher:
    """keeps a run's LiveTree up to date with the merged units and publishes
    the shards that changed, then meta.json"""

    def __init__(self, b, run_id, instance, n, grid, source, log=print, threads=16, prefix=()):
        self.b, self.run_id, self.log, self.threads = b, run_id, log, threads
        self.prefix = tuple(prefix)
        self.meta = {"run": run_id, "title": run_id, "instance": instance, "n": n, "grid": grid, "format": "nodes",
                     "source": source, "notes": None, "layout": "depth2",
                     # a subtree run: the -F node's fixes, the root's own (facility:location, 0-based)
                     "rootFixes": ",".join(f"{a}:{b}" for a, b in self.prefix) or None}
        self.tree = None
        self.seq = 0
        self.missing = 0

    def _fetch(self, u):
        return self.b.get_bytes(unit_key(self.run_id, u["uid"], f"a{u['attempt']}", "nodes.txt"))

    def update(self, pool):
        """read the frontier's trace (once) and the units merged since; False
        when the run has no trace"""
        if self.tree is None:
            # a frontier pass that bounded no node (a -F run at depth 0) leaves no trace
            frontier = self.b.get_bytes(run_key(self.run_id, "frontier", "nodes.txt")) or b""
            self.tree = LiveTree(self.prefix)
            self.tree.ingest(frontier)
        units = sorted(self.b.units_merged_after(self.run_id, self.seq), key=lambda u: u["merged_seq"])
        for u, data in zip(units, pool.map(self._fetch, units)):
            if data is None:
                self.missing += 1
            else:
                self.tree.ingest(data)
            self.seq = u["merged_seq"]
        return True

    def publish(self, final=False, now=None):
        """True once the tree is published (meta.json written)"""
        t0 = time.time()
        with ThreadPoolExecutor(self.threads) as pool:
            if not self.update(pool):
                return False
            t = self.tree
            roots = t.take_dirty()
            key = lambda name: f"public/runs/{self.run_id}/tree/{name}"
            put = lambda r: self.b.put_bytes(key(f"{r}.json"), json.dumps(t.shard(r), separators=(",", ":")).encode(),
                                             "application/json", CACHE)
            list(pool.map(put, roots))
        meta = dict(self.meta, cutoff=t.cutoff_value(), nodes=len(t), maxDepth=t.mx[0], shards=len(t.shard_roots()),
                    outcomes=t.outcome_counts(), orphans=sum(len(v) for v in t.waiting.values()),
                    duplicates=t.duplicates, missingUnits=self.missing, generatedAt=iso(now or time.time()),
                    final=final)
        self.b.put_bytes(key("meta.json"), json.dumps(meta, indent=1).encode(), "application/json", CACHE)
        self.log(f"tree: {len(t)} nodes, {len(roots)} of {meta['shards']} shards written in {time.time() - t0:.1f} s")
        return True


def write_dir(trace, out, run, instance=None, n=None, grid=None, source=None):
    """a whole trace's tree as the publisher would leave it, into a directory
    (meta.json and every shard): a fixture for brixius-web, or a look at a tree"""
    import os
    t = LiveTree()
    with open(trace, "rb") as f:
        t.ingest(f.read())
    os.makedirs(out, exist_ok=True)
    roots = t.shard_roots()
    for r in roots:
        with open(os.path.join(out, f"{r}.json"), "w") as f:
            json.dump(t.shard(r), f, separators=(",", ":"))
    meta = {"run": run, "title": run, "instance": instance, "n": n, "grid": grid, "format": "nodes",
            "source": source or os.path.basename(trace), "notes": None, "layout": "depth2",
            "cutoff": t.cutoff_value(), "nodes": len(t), "maxDepth": t.mx[0], "shards": len(roots),
            "outcomes": t.outcome_counts(), "orphans": sum(len(v) for v in t.waiting.values()),
            "duplicates": t.duplicates, "missingUnits": 0, "generatedAt": None, "final": True}
    with open(os.path.join(out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    return meta


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="a node trace's tree as coord/tree.py publishes it, into a directory")
    p.add_argument("trace")
    p.add_argument("out")
    p.add_argument("--run", required=True)
    p.add_argument("--instance")
    p.add_argument("--n", type=int)
    p.add_argument("--grid", help="ROWSxCOLS")
    a = p.parse_args()
    g = None
    if a.grid:
        rows, cols = map(int, a.grid.split("x"))
        g = {"rows": rows, "cols": cols}
    m = write_dir(a.trace, a.out, a.run, a.instance, a.n, g)
    print(f"{a.out}: {m['nodes']} nodes, {m['shards']} shards")
