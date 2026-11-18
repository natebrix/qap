"""Keep the nodes of a bnbdrv -o state file at depth <= DEPTH and only the warm starts they use (renumbered):
  python3 trimstate.py IN.state OUT.state DEPTH
run.sh trim uses it: an arm's set-aside grandchildren (-d 3) are never bounded, but each branched child's
DNN state stays in memory for them (~57 MB on disk, ~0.4 GB of RAM at n = 60)."""
import sys
src, dst, depth = sys.argv[1], sys.argv[2], int(sys.argv[3])
head, warms, nodes, tail = [], [], [], []
with open(src) as fh:
    for ln in fh:
        if ln.startswith("node "): nodes.append(ln)
        elif ln.startswith("warm "): warms.append(ln)
        elif ln.startswith(("warm_starts ", "nodes ")): pass
        elif nodes or ln == "end\n": tail.append(ln)
        else: head.append(ln)
assert tail == ["end\n"], tail[:3]
ids, kept = {}, []
for ln in nodes:
    if int(ln.split(None, 2)[1]) > depth: continue
    body, w = ln.rsplit(None, 1)
    w = int(w)
    if w >= 0: w = ids.setdefault(w, len(ids))
    kept.append(f"{body} {w}\n")
with open(dst, "w") as out:
    out.writelines(head)
    out.write(f"warm_starts {len(ids)}\n")
    for old in sorted(ids, key=ids.get): out.write(warms[old])
    out.write(f"nodes {len(kept)}\n")
    out.writelines(kept)
    out.write("end\n")
print(f"{len(kept)} of {len(nodes)} nodes kept, {len(ids)} of {len(warms)} warm starts", file=sys.stderr)
