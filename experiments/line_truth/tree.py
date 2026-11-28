#!/usr/bin/env python3
"""tree.py SHARDDIR: the published tree as one table (a row per node, the shards' columns), and whether each
node's subtree is closed: a fathomed or enumerated node is; a branched node is when every child it made (kids) is
in the tree and closed. load(dir) returns (cols, nodes by id, children by id, closed set, meta)."""
import json, os, glob, sys
FATH = set('QDCPE')
def load(d):
    meta = json.load(open(os.path.join(d, 'meta.json')))
    nodes = {}; cols = None
    for f in glob.glob(os.path.join(d, '*.json')):
        if f.endswith('meta.json'): continue
        s = json.load(open(f)); cols = s['columns']; i_id = cols.index('id')
        for r in s['rows']:
            nodes.setdefault(r[i_id], r)
    ix = {c: i for i, c in enumerate(cols)}
    kids = {}
    for r in nodes.values(): kids.setdefault(r[ix['parent']], []).append(r[ix['id']])
    closed = set()
    order = sorted(nodes.values(), key=lambda r: -r[ix['depth']])
    for r in order:
        i = r[ix['id']]
        if r[ix['outcome']] in FATH: closed.add(i)
        elif r[ix['outcome']] == 'B':
            ch = kids.get(i, [])
            if len(ch) == r[ix['kids']] and all(c in closed for c in ch): closed.add(i)
    return ix, nodes, kids, closed, meta
if __name__ == '__main__':
    ix, nodes, kids, closed, meta = load(sys.argv[1])
    import collections
    print(f"{meta['run']}: {len(nodes)} nodes read of {meta['nodes']} (generated {meta['generatedAt']})")
    byd = collections.defaultdict(lambda: [0, 0, 0, 0.0, 0.0])
    for r in nodes.values():
        b = byd[r[ix['depth']]]; b[0] += 1
        if r[ix['outcome']] == 'B':
            b[1] += 1; b[2] += r[ix['id']] in closed
        b[3] += r[ix['cpu']] or 0
    print("depth  nodes  branched  closed subtrees  own CPU h")
    for d in sorted(byd): print(f"{d:5d} {byd[d][0]:7d} {byd[d][1]:7d} {byd[d][2]:7d} {byd[d][3]/3600:10.1f}")
