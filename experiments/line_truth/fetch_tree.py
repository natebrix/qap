#!/usr/bin/env python3
"""fetch_tree.py RUN OUTDIR [WORKERS]: the published tree of a public --trace run (qap coord/tree.py), through the
site's proxy (brixius-web app/projects/qap/tree/[run]/data/[file]): meta.json, the root shard 0.json, and the
shard <id>.json of every boundary node (column `shard` = 1), breadth first. Resumable: a shard already in OUTDIR
is read from there. A few requests at a time."""
import json, os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
run, out = sys.argv[1], sys.argv[2]; workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
base = f"https://nathanbrixius.com/projects/qap/tree/{run}/data/"
os.makedirs(out, exist_ok=True)

def get(name, refresh=False):
    p = os.path.join(out, name)
    if not refresh and os.path.exists(p) and os.path.getsize(p) > 50:
        return json.load(open(p))
    for attempt in range(5):
        try:
            req = urllib.request.Request(base + name, headers={"User-Agent": "tree-fetch"})
            with urllib.request.urlopen(req, timeout=60) as r: body = r.read()
            d = json.loads(body)
            with open(p + ".tmp", "wb") as f: f.write(body)
            os.replace(p + ".tmp", p)
            return d
        except Exception as e:
            err = e; time.sleep(2 * (attempt + 1))
    print("failed", name, err, flush=True); return None

meta = get("meta.json", refresh=True)
print(f"{run}: {meta['nodes']} nodes, {meta['shards']} shards, generated {meta['generatedAt']}, final {meta['final']}", flush=True)
todo = ["0.json"]; done = 0; t0 = time.time()
with ThreadPoolExecutor(workers) as ex:
    while todo:
        res = list(ex.map(get, todo)); nxt = []
        for d in res:
            if not d: continue
            done += 1
            ci = {c: i for i, c in enumerate(d["columns"])}
            root = d["root"]
            for r in d["rows"]:
                if r[ci["shard"]] == 1 and r[ci["id"]] != root: nxt.append(f"{r[ci['id']]}.json")
        print(f"{done} shards, {len(nxt)} next, {time.time()-t0:.0f} s", flush=True)
        todo = nxt
print("done", done)
