"""The open nodes of a state file as text: depth, inherited bound, fixes (facility:location, 0-based) -- enough to
restart a node cold with bnbdrv -F when the state files are gone.   python3 openlist.py FILE.state"""
import sys
for ln in open(sys.argv[1]):
    if not ln.startswith("node "): continue
    t = ln.split(); p = [int(x) for x in t[5:55]]
    print(t[1], t[2], ",".join("%d:%d" % (i, x) for i, x in enumerate(p) if x >= 0))
