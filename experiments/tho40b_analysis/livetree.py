"""A running run's tree, from its published live tree (qap coord/tree.py), as a node
trace the other scripts here read: the analyses on what is already finished, while the
run goes on.

  crawl RUN OUT.jsonl.gz   every merged node of a public --trace run, through the site's
                           tree endpoint (brixius-web app/projects/qap/tree/[run]/data:
                           the bundle's shards, read without naming the bucket); about
                           20 min for tho40-b at 1.6e6 nodes, a few requests at a time
  trace IN.jsonl.gz OUT.txt.gz [--complete-at D] [--min-nodes N]
                           the rows as QAP_NODE_TRACE lines (fixes rebuilt from the
                           parents). --complete-at D keeps only the complete subtrees
                           whose roots are at depth D, each its own root: every child a
                           branched node made is merged, down to the leaves
  summary IN.jsonl.gz      by depth: nodes, branched, how many of those are complete,
                           the nodes and CPU in complete subtrees
  levels IN.jsonl.gz       by depth, over every merged node: fathoms by bound, children
                           kept per branched node of those the line allowed, the share
                           branched on a facility (row), CPU per node and the depth's
                           share of all CPU, a branched node's shortfall from the cutoff

A shard row: id, parent, depth, fix (the node's last facility:location), outcome, lb,
qpb, dnn, cut, cutoff, iterations, cpu, status, rule, line, kids, size, fathomed,
subCpu, maxDepth, shard (1: a boundary node, its children are in the shard named by
its id). Subtrees finish small ones first, so what is complete part-way leans to the
small and the easy: per-node findings carry over, totals do not.
"""
import argparse
import collections
import gzip
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

SITE = "https://nathanbrixius.com/projects/qap/tree/{run}/data/{id}.json"


def crawl(run, out, workers=4):
    def get(sid):
        err = None
        for a in range(5):
            try:
                return json.loads(urllib.request.urlopen(SITE.format(run=run, id=sid), timeout=120).read())
            except Exception as e:   # transient: retry with backoff
                err = e
                time.sleep(2 * (a + 1))
        raise err
    s0 = get(0)
    cols = s0["columns"]
    C = {c: i for i, c in enumerate(cols)}
    rows, seen = {}, {0}

    def add(shard, sid):
        new = []
        for r in shard["rows"]:
            nid = r[C["id"]]
            rows[nid] = r
            if r[C["shard"]] == 1 and nid != sid and nid not in seen:
                seen.add(nid)
                new.append(nid)
        return new
    todo, n = add(s0, 0), 1
    with ThreadPoolExecutor(workers) as ex:
        while todo:
            batch, todo = todo, []
            for sh, sid in zip(ex.map(get, batch), batch):
                todo += add(sh, sid)
                n += 1
            print(f"{n} shards, {len(rows):,} rows", flush=True)
    with gzip.open(out, "wt") as f:
        f.write(json.dumps({"run": run, "columns": cols, "fetched": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()),
                            "shards": n}) + "\n")
        for r in rows.values():
            f.write(json.dumps(r) + "\n")


def load(path):
    f = gzip.open(path, "rt")
    hdr = json.loads(f.readline())
    C = {c: i for i, c in enumerate(hdr["columns"])}
    rows = {}
    for line in f:
        r = json.loads(line)
        rows[r[C["id"]]] = r
    kids = collections.defaultdict(list)
    for r in rows.values():
        p = r[C["parent"]]
        if p is not None and p >= 0:
            kids[p].append(r[C["id"]])
    comp = {}
    for i in sorted(rows, key=lambda i: -rows[i][C["depth"]]):   # deepest first
        r = rows[i]
        comp[i] = r[C["outcome"]] != "B" or (len(kids[i]) == (r[C["kids"]] or 0) and all(comp[k] for k in kids[i]))
    return hdr, C, rows, kids, comp


def to_trace(path, out, complete_at=None, min_nodes=1):
    hdr, C, rows, kids, comp = load(path)
    fixes = {}

    def fx(i):
        if i not in fixes:
            r = rows[i]
            p = r[C["parent"]]
            fixes[i] = (fx(p) + [r[C["fix"]]]) if p is not None and p >= 0 else []
        return fixes[i]
    sys.setrecursionlimit(10000)
    if complete_at is None:
        keep = list(rows)
    else:
        keep = []
        for i, r in rows.items():
            if r[C["depth"]] == complete_at and comp[i] and r[C["outcome"]] == "B" and r[C["size"]] >= min_nodes:
                stack = [i]
                while stack:
                    j = stack.pop()
                    keep.append(j)
                    stack += kids[j]
    v = lambda x: "-" if x is None else f"{x:.4f}"
    n = 0
    with gzip.open(out, "wt") as f:
        for i in sorted(keep, key=lambda i: (rows[i][C["depth"]], i)):
            r = rows[i]
            b = r[C["outcome"]] == "B"
            fs = fx(i)
            f.write(f"node {r[C['depth']]} {','.join(fs) if fs else '-'} {v(r[C['qpb']])} {v(r[C['dnn']])} "
                    f"{v(r[C['cut']])} {v(r[C['cutoff']])} {r[C['outcome']]} {r[C['rule']] if b else '-'} "
                    f"{r[C['line']] if b else '-'} {r[C['kids']] if b else '-'} {r[C['cpu']]}\n")
            n += 1
    print(f"{out}: {n:,} nodes (of {len(rows):,} merged, fetched {hdr['fetched']})")


def summary(path):
    hdr, C, rows, kids, comp = load(path)
    by = collections.defaultdict(lambda: [0, 0, 0, 0, 0.0])
    for i, r in rows.items():
        b = by[r[C["depth"]]]
        b[0] += 1
        if r[C["outcome"]] == "B":
            b[1] += 1
            if comp[i]:
                b[2] += 1
                b[3] += r[C["size"]]
                b[4] += r[C["subCpu"]]
    print(f"{hdr.get('run', '')} merged nodes {len(rows):,}, fetched {hdr['fetched']}")
    print(f"{'depth':>5} {'nodes':>9} {'branched':>9} {'complete':>9} {'nodes in complete subtrees':>27} {'their CPU s':>13}")
    for d in sorted(by):
        b = by[d]
        print(f"{d:>5} {b[0]:>9,} {b[1]:>9,} {b[2]:>9,} {b[3]:>27,} {b[4]:>13,.0f}")


def levels(path):
    hdr, C, rows, kids, comp = load(path)
    by = collections.defaultdict(list)
    for r in rows.values():
        by[r[C["depth"]]].append(r)
    total_cpu = sum(r[C["cpu"]] for r in rows.values())
    print(f"{hdr.get('run', '')} merged nodes {len(rows):,}, CPU {total_cpu:,.0f} s, fetched {hdr['fetched']}")
    print(f"{'depth':>5} {'nodes':>9} {'D':>6} {'Q':>6} {'C+P':>5} {'E':>5} {'branched':>8} {'kids/free':>9} "
          f"{'rows':>5} {'CPU/node':>9} {'CPU %':>6} {'shortfall: median (10%-90%)':>30}")
    for d in sorted(by):
        rs = by[d]
        cnt = collections.Counter(r[C["outcome"]] for r in rs)
        br = [r for r in rs if r[C["outcome"]] == "B"]
        pct = lambda k: f"{100 * cnt[k] / len(rs):5.1f}%"
        free = max(1, 40 - d)   # tho40: n = 40
        kf = sum(r[C["kids"]] or 0 for r in br) / max(1, len(br)) / free
        rowshare = sum(1 for r in br if (r[C["line"]] or "").startswith("r")) / max(1, len(br))
        cpu = sum(r[C["cpu"]] for r in rs)
        sf = sorted(r[C["cutoff"]] - r[C["lb"]] for r in br if r[C["lb"]] is not None)
        q = (lambda f: sf[min(len(sf) - 1, int(f * len(sf)))]) if sf else (lambda f: float("nan"))
        print(f"{d:>5} {len(rs):>9,} {pct('D')} {pct('Q')} {100 * (cnt['C'] + cnt['P']) / len(rs):4.1f}% "
              f"{100 * cnt['E'] / len(rs):4.1f}% {len(br):>8,} {kf:>9.2f} {rowshare:>5.2f} {cpu / len(rs):>9.2f} "
              f"{100 * cpu / total_cpu:>5.1f}% {q(0.5):>12,.0f} ({q(0.1):,.0f}-{q(0.9):,.0f})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("crawl")
    c.add_argument("run")
    c.add_argument("out")
    t = sub.add_parser("trace")
    t.add_argument("inp")
    t.add_argument("out")
    t.add_argument("--complete-at", type=int)
    t.add_argument("--min-nodes", type=int, default=1)
    s = sub.add_parser("summary")
    s.add_argument("inp")
    lv = sub.add_parser("levels")
    lv.add_argument("inp")
    a = ap.parse_args()
    if a.cmd == "crawl":
        crawl(a.run, a.out)
    elif a.cmd == "trace":
        to_trace(a.inp, a.out, a.complete_at, a.min_nodes)
    elif a.cmd == "levels":
        levels(a.inp)
    else:
        summary(a.inp)


if __name__ == "__main__":
    main()
