"""Keep the depth-DEPTH nodes of a bnbdrv -o state file that hold any of the fixes FAC:LOC,... and only the warm
starts they use (renumbered; trimstate.py's logic with a filter on the permutation):
  python3 splitstate.py IN.state OUT.state DEPTH FAC:LOC[,FAC:LOC...]"""
import sys
src, dst, depth = sys.argv[1], sys.argv[2], int(sys.argv[3])
fx = [tuple(map(int, s.split(':'))) for s in sys.argv[4].split(',')]
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
    t = ln.split()
    if int(t[1]) != depth or not any(int(t[5 + f]) == l for f, l in fx): continue
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
