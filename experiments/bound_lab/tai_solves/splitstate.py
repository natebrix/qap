"""Split or filter a bnbdrv -o state file by its open nodes (after tho40_search/cutstate.py):
  python3 splitstate.py IN.state OUTPREFIX K            -> OUTPREFIX0.state .. OUTPREFIX{K-1}.state, node lines dealt round-robin
  python3 splitstate.py IN.state OUT.state -pick F:L,.. -> the nodes with one of the listed facility:location pairs fixed (-drop: without any)
  python3 splitstate.py IN.state OUT.state -depth D     -> the nodes at depth D only
  python3 splitstate.py IN.state OUT.state -onpath x    -> the deepest node on the path of the best known solution (-offpath: every node off it); a file name in place of x: that permutation's path
  python3 splitstate.py -merge OUT.state IN1.state IN2.state ...   -> one file with every node line (headers from IN1)
Node line: "node DEPTH ... p[0..n-1] ..." as cutstate.py reads it (t[5+f] = location of facility f)."""
import sys, os
def read(src):
    head, nodes, tail = [], [], []
    with open(src) as fh:
        for ln in fh:
            if ln.startswith("nodes "): continue
            if ln.startswith("node "): nodes.append(ln)
            elif nodes or ln.strip() == "end": tail.append(ln)   # (a file with no open node still ends with "end")
            else: head.append(ln)
    return head, nodes, tail
def write(dst, head, nodes, tail):
    with open(dst, "w") as w:
        w.writelines(head); w.write("nodes %d\n" % len(nodes)); w.writelines(nodes); w.writelines(tail)
if __name__ == "__main__":
    if sys.argv[1] == "-merge":
        # node lines end with the index of their parent's "warm" line (-1: none), so the warm lines of every
        # input are kept, in order, and each input's indices are shifted by the number of warm lines before it
        H, W, NODES, tail = None, [], [], None
        for src in sys.argv[3:]:
            head, nodes, tl = read(src)
            w = [ln for ln in head if ln.startswith("warm ")]
            if H is None: H = [ln for ln in head if not ln.startswith("warm ") and not ln.startswith("warm_starts ")]; tail = tl
            for ln in nodes:
                t = ln.split(); k = int(t[-1])
                if k >= 0: t[-1] = str(k + len(W))
                NODES.append(" ".join(t) + "\n")
            W += w
        write(sys.argv[2], H + ["warm_starts %d\n" % len(W)] + W, NODES, tail); print(len(NODES), "nodes,", len(W), "warm lines")
    elif sys.argv[3] == "-depth":
        head, nodes, tail = read(sys.argv[1])
        keep = [ln for ln in nodes if int(ln.split()[1]) == int(sys.argv[4])]
        write(sys.argv[2], head, keep, tail); print(len(keep), "nodes")
    elif sys.argv[3] in ("-onpath", "-offpath"):
        # the deepest node all of whose fixes agree with QAPSOLN/tai50b.sln (0-based, facility -> location)
        BK = [5,23,14,12,21,4,29,40,13,8,16,19,25,18,46,11,15,37,17,30,47,20,33,26,36,48,22,10,3,31,2,7,0,32,6,42,24,28,44,27,39,35,43,45,1,34,49,38,9,41]
        if os.path.exists(sys.argv[4]): BK = [int(x) for x in open(sys.argv[4]).read().split()]   # another permutation's path
        head, nodes, tail = read(sys.argv[1])
        def ok(ln):
            p = [int(x) for x in ln.split()[5:55]]; return all(x < 0 or x == BK[i] for i, x in enumerate(p))
        keep = sorted([ln for ln in nodes if ok(ln)], key=lambda ln: -int(ln.split()[1]))[:1]
        if sys.argv[3] == "-offpath": keep = [ln for ln in nodes if not ok(ln)]
        write(sys.argv[2], head, keep, tail); print(len(keep), "nodes", "depth", keep[0].split()[1] if keep else "-")
    elif sys.argv[3] == "-drop":
        head, nodes, tail = read(sys.argv[1])
        want = [tuple(map(int, x.split(":"))) for x in sys.argv[4].split(",")]
        keep = [ln for ln in nodes if not any(int(ln.split()[5 + f]) == l for f, l in want)]
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
