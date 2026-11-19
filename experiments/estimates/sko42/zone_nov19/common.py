"""Shared by the scripts here: a whole-tree run's node trace (QAP_NODE_TRACE, latest/bnb.C: "node depth fixes
qpb dnn cut cutoff outcome rule line kids cpu"), each node once, with its subtree's CPU and nodes and whether
the subtree is closed (a branched node is closed when every child it made is traced and closed)."""
import collections
import gzip
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
RUNS = os.path.join(ROOT, "infra", "phase1", "runs")
ARMS = {"z10": os.path.join(RUNS, "sko42-zone-261119-z10", "nodes.txt.gz"),
        "z25": os.path.join(RUNS, "sko42-zone-261119-z25", "nodes.txt.gz")}


def parent(key):
    return key.rsplit(",", 1)[0] if "," in key else ""


def read_lines(path):
    return [l.split() for l in gzip.open(path, "rt") if l.startswith("node ")]


def build(lines):
    """{fixes: node}, {fixes: [children]}, subtree cpu, subtree size, closed; the root's fixes are ''"""
    nodes = {}
    for f in lines:
        k = "" if f[2] == "-" else f[2]
        if k in nodes:                 # a node traced twice (a retried unit) counts once
            continue
        v = lambda x: None if x == "-" else float(x)
        nodes[k] = dict(depth=int(f[1]), qpb=v(f[3]), dnn=v(f[4]), cut=v(f[5]), out=f[7],
                        kids=-1 if f[10] == "-" else int(f[10]), cpu=float(f[11]))
    ch = collections.defaultdict(list)
    for k in nodes:
        if k:
            ch[parent(k)].append(k)
    cpu, size, closed = {}, {}, {}
    for k in sorted(nodes, key=lambda k: -nodes[k]["depth"]):
        n, c = nodes[k], ch.get(k, [])
        cpu[k] = n["cpu"] + sum(cpu[x] for x in c)
        size[k] = 1 + sum(size[x] for x in c)
        closed[k] = n["out"] != "B" or (len(c) == n["kids"] and all(closed[x] for x in c))
    return dict(nodes=nodes, ch=ch, cpu=cpu, size=size, closed=closed)


def load(path):
    return build(read_lines(path))
