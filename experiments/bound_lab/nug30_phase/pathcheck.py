"""Sanity on the path to the optimum.   python3 pathcheck.py PATHRUN [ARM...]
PATHRUN: a whole search started from incumbent 6126 (cutoff 6124.05), which must find 6124 itself. The permutation
it reports (final.state) is checked against the data, then:
  in PATHRUN: every traced node whose fixes agree with it must branch (or be enumerated), with every certified
              value (first solve, phase) at most 6124;
  in each ARM (incumbent 6124, cutoff 6122.05): a node on that path may fathom, but no certified value on it may pass 6124.
"""
import os, re, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(H, "../../../qaplib/nug30.dat")).read().split(); n = 30
v = np.array(t[1:1 + 2 * n * n], float).reshape(2, n, n)
run = sys.argv[1]
st = open(os.path.join(H, "runs", run, "final.state")).read()
inc = re.search(r"^incumbent (\S+) (\d+) (.*)$", st, re.M)
p = np.array(inc.group(3).split(), int)
obj = (v[0] * v[1][np.ix_(p, p)]).sum()
log = open(os.path.join(H, "runs", run, "log.txt")).read()
print(f"{run}: incumbent {inc.group(1)}, found flag {inc.group(2)}, its permutation costs {obj:.0f} on qaplib/nug30.dat; 'New objective value' lines: {re.findall(r'New objective value: (\S+)', log)}; open nodes at the end: {st.count(chr(10) + 'node ')}")
ok = float(inc.group(1)) == 6124.0 and obj == 6124.0
def path_nodes(arm):
    out = []
    for ln in open(os.path.join(H, "runs", arm, "nodes.txt")):
        f = ln.split()
        if len(f) < 12 or f[0] != "node": continue
        fx = [] if f[2] == "-" else [tuple(map(int, x.split(":"))) for x in f[2].split(",")]
        if all(p[a] == b for a, b in fx): out.append(f)
    return sorted(out, key=lambda f: int(f[1]))
for arm in sys.argv[1:]:
    strict = arm == run
    print(f"{arm}: nodes on the path ({'must branch' if strict else 'may fathom; values must not pass 6124'})")
    for f in path_nodes(arm):
        vals = [float(x) for x in (f[4], f[5]) if x != "-"]
        good = all(x <= 6124.0 + 1e-6 for x in vals) and (not strict or f[7] in "BE")
        ok &= good
        print("  depth %s %-28s first solve %-10s phase %-10s cutoff %s outcome %s  %s" % (f[1], f[2], f[4], f[5], f[6], f[7], "ok" if good else "FAIL"))
print("PASS" if ok else "FAIL")
