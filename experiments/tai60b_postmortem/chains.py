"""tai60b-261117's finished tree, from the website's public tree shards: the two chains below the
root (44:45, the best known solution's path, and 44:58, its mirror), and at every chain depth the
subtrees of the chain node's other children (nodes, CPU, deepest node, branched children).
    python3 -I chains.py CACHE_DIR > chains.txt
Reads https://nathanbrixius.com/projects/qap/tree/tai60b-261117/data/<shard>.json (fetched once
into CACHE_DIR): only the shards along the chains."""
import json, os, sys, urllib.request
BASE = "https://nathanbrixius.com/projects/qap/tree/tai60b-261117/data/"
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
cache = {}
def shard(name):
    if name not in cache:
        p = os.path.join(out, f"{name}.json")
        if not os.path.exists(p):
            with urllib.request.urlopen(BASE + f"{name}.json", timeout=60) as r: open(p, "wb").write(r.read())
        s = json.load(open(p)); c = s["columns"]
        cache[name] = [dict(zip(c, r)) for r in s["rows"]]
    return cache[name]
def kids_of(node, rows):
    k = [r for r in rows if r["parent"] == node["id"]]
    if not k and node.get("shard"):          # 1: the node's subtree is in the shard named by its id
        rows = shard(str(node["id"])); k = [r for r in rows if r["parent"] == node["id"]]
    return k, rows
top = shard("0"); root = [r for r in top if r["depth"] == 0][0]
per = {}
for start in ("44:45", "44:58"):
    node = [r for r in top if r["depth"] == 1 and r["fix"] == start][0]; R = top
    while True:
        k, R = kids_of(node, R)
        if not k: break
        k.sort(key=lambda r: -(r["size"] or 0)); nxt, off = k[0], k[1:]
        e = per.setdefault(node["depth"], dict(n=0, cpu=0.0, br=0, deep=0, rule=set(), lines=[]))
        e["n"] += sum(r["size"] or 0 for r in off); e["cpu"] += sum(r["subCpu"] or 0 for r in off)
        e["br"] += sum(r["outcome"] == "B" for r in off); e["deep"] = max([e["deep"]] + [r["maxDepth"] or 0 for r in off])
        e["rule"].add(node["rule"]); e["lines"].append(node["line"])
        if nxt["outcome"] != "B": break
        node = nxt
N, C = root["size"], root["subCpu"]
print(f"tree: {N} nodes, {C/3600:.0f} CPU h; root line {root['line']}")
print("| chain depth | strategy (2: rule 9, 3: rule 2) | lines (path, mirror) | off-chain children branched | their nodes | % of nodes | CPU h | % of CPU | deepest |")
print("|---|---|---|---|---|---|---|---|---|")
for d in sorted(per):
    e = per[d]
    if e["cpu"] / C < 0.002 and d > 12: continue
    print(f"| {d} | {','.join(map(str, sorted(e['rule'])))} | {', '.join(map(str, e['lines']))} | {e['br']} | {e['n']:,} | {100*e['n']/N:.1f} | {e['cpu']/3600:,.0f} | {100*e['cpu']/C:.1f} | {e['deep']} |")
