import json, sys, urllib.request, os
from concurrent.futures import ThreadPoolExecutor
run, out = sys.argv[1], sys.argv[2]
base = f"https://nathanbrixius.com/projects/qap/tree/{run}/data/"
os.makedirs(out, exist_ok=True)
def get(name):
    p = os.path.join(out, name)
    if not os.path.exists(p):
        with urllib.request.urlopen(base + name, timeout=60) as r: open(p, "wb").write(r.read())
    return json.load(open(p))
meta = get("meta.json")
seen, todo = set(), ["0"]
with ThreadPoolExecutor(8) as ex:
    while todo:
        todo = [t for t in todo if t not in seen]; seen.update(todo)
        shards = list(ex.map(lambda s: get(f"{s}.json"), todo))
        todo = []
        for sh in shards:
            cols = sh["columns"] if isinstance(sh, dict) and "columns" in sh else None
            rows = sh["rows"] if cols else sh
            for r in rows:
                d = dict(zip(cols, r)) if cols else r
                if d.get("shard") not in (None, False, 0) and str(d["id"]) not in seen:
                    todo.append(str(d["id"]))
print(meta["nodes"], len(seen), "shards")
