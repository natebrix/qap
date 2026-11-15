"""Cut a bnbdrv -o state file to one of its open nodes (as experiments/tho40b_analysis/path_check/dive.py --warm):
  python3 cutstate.py IN.state OUT.state DEPTH FACILITY:LOCATION
keeps the node line at DEPTH whose permutation has FACILITY at LOCATION, with its shipped DNN state."""
import sys
src, dst, depth, fl = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
f, l = map(int, fl.split(":"))
kept = 0
with open(src) as fh, open(dst, "w") as w:
    for ln in fh:
        if ln.startswith("nodes "): continue
        if ln.startswith("node "):
            t = ln.split()
            if int(t[1]) == depth and int(t[5 + f]) == l:
                w.write("nodes 1\n" + ln); kept += 1
            continue
        w.write(ln)
sys.exit(0 if kept == 1 else 3)
