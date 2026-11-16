"""Split or filter a bnbdrv -o state file by its open nodes (after tho40_search/cutstate.py):
  python3 splitstate.py IN.state OUTPREFIX K            -> OUTPREFIX0.state .. OUTPREFIX{K-1}.state, node lines dealt round-robin
  python3 splitstate.py IN.state OUT.state -pick F:L,.. -> the nodes whose LAST fix is one of the listed facility:location pairs
  python3 splitstate.py IN.state OUT.state -depth D     -> the nodes at depth D only
  python3 splitstate.py -merge OUT.state IN1.state IN2.state ...   -> one file with every node line (headers from IN1)
Node line: "node DEPTH ... p[0..n-1] ..." as cutstate.py reads it (t[5+f] = location of facility f)."""
import sys
def read(src):
    head, nodes, tail = [], [], []
    with open(src) as fh:
        for ln in fh:
            if ln.startswith("nodes "): continue
            if ln.startswith("node "): nodes.append(ln)
            elif nodes: tail.append(ln)
            else: head.append(ln)
    return head, nodes, tail
def write(dst, head, nodes, tail):
    with open(dst, "w") as w:
        w.writelines(head); w.write("nodes %d\n" % len(nodes)); w.writelines(nodes); w.writelines(tail)
if sys.argv[1] == "-merge":
    head, nodes, tail = read(sys.argv[3])
    for s in sys.argv[4:]: nodes += read(s)[1]
    write(sys.argv[2], head, nodes, tail); print(len(nodes), "nodes")
elif sys.argv[3] == "-depth":
    head, nodes, tail = read(sys.argv[1])
    keep = [ln for ln in nodes if int(ln.split()[1]) == int(sys.argv[4])]
    write(sys.argv[2], head, keep, tail); print(len(keep), "nodes")
elif sys.argv[3] == "-pick":
    head, nodes, tail = read(sys.argv[1])
    want = [tuple(map(int, x.split(":"))) for x in sys.argv[4].split(",")]
    keep = [ln for ln in nodes if any(int(ln.split()[5 + f]) == l for f, l in want)]
    write(sys.argv[2], head, keep, tail); print(len(keep), "nodes")
else:
    head, nodes, tail = read(sys.argv[1]); k = int(sys.argv[3])
    for i in range(k): write("%s%d.state" % (sys.argv[2], i), head, nodes[i::k], tail)
    print(len(nodes), "nodes in", k, "files")
