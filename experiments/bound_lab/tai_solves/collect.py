"""A bfs.sh search summed up:  python3 collect.py TAG INSTANCE OPTIMUM  (runs/TAG/d*_c*/{nodes.txt,dnn.txt,log.txt,out.state})
Per depth from the node traces (one line per node: outcome Q/D/C fathomed by QPB / the DNN first solve / the cut
phase, B branched, E enumerated, P discarded when popped) and the DNN traces (iterations, CPU); the path check
against logs/INSTANCE.perm.txt (a node all of whose fixes agree with the optimal permutation must not fathom and
no certified bound there may exceed the optimum); the incumbents found; the open nodes left."""
import sys, glob, os, re, collections
H = os.path.dirname(os.path.abspath(__file__)); tag, inst, opt = sys.argv[1], sys.argv[2], int(sys.argv[3])
P = [int(x) for x in open(f"{H}/logs/{inst}.perm.txt").read().split()]
dirs = sorted(glob.glob(f"{H}/runs/{tag}/d*_c*")) + ([f"{H}/runs/{sys.argv[4]}"] if len(sys.argv) > 4 else [])
lev = collections.defaultdict(lambda: collections.Counter()); onp = []; bad = []
for d in dirs:
    if os.path.exists(d + "/nodes.txt"):
        for ln in open(d + "/nodes.txt"):
            t = ln.split()
            if t[0] != "node": continue
            dep = int(t[1]); out = t[7]; L = lev[dep]; L[out] += 1; L["cpu"] += float(t[-1])
            fx = [] if t[2] == "-" else [tuple(map(int, x.split(":"))) for x in t[2].split(",")]
            if all(P[f] == l for f, l in fx):
                vals = [float(x) for x in t[3:6] if x != "-"]
                onp.append((dep, t[2], max(vals) if vals else None, out))
                if out in "QDCP" or (vals and max(vals[1:] or [0]) > opt): bad.append(ln.strip())
    if os.path.exists(d + "/dnn.txt"):
        for ln in open(d + "/dnn.txt"):
            m = re.search(r"depth (\d+) .*\| first lb \S+ primal \S+ status \S+ iters (\d+) cpu (\S+)", ln)
            if not m: continue
            L = lev[int(m[1])]; L["dnn"] += 1; L["it"] += int(m[2]); L["dcpu"] += float(m[3])
            c = re.search(r"\| cut 1 best \S+ rounds \d+ iters (\d+) cpu (\S+)", ln)
            if c: L["phases"] += 1; L["it"] += int(c[1]); L["dcpu"] += float(c[2])
found = sorted({ln.strip() for d in dirs for ln in open(d + "/log.txt") if ln.startswith("New objective value")})
openn = sum(1 for d in dirs if os.path.exists(d + "/out.state") for ln in open(d + "/out.state") if ln.startswith("node "))
last = max((int(os.path.basename(d).split("_")[0][1:]) for d in dirs if "/d" in d), default=0)
left = sum(1 for d in dirs if os.path.basename(d).startswith(f"d{last}_") or "/d" not in d and True for ln in (open(d + "/out.state") if os.path.exists(d + "/out.state") else []) if ln.startswith("node "))
unfinished = [os.path.basename(d) for d in dirs if "exit 0" not in open(d + "/cmd.txt").read()]
print(f"{inst} {tag}: optimum {opt:,}; incumbents found: {found or 'none'}")
print(f"open nodes left in the last level's state files: {left}; processes not ended with exit 0: {unfinished or 'none'}")
T = collections.Counter()
print("\n| depth | nodes | branched | fathomed: QPB / DNN first solve / cut phase | enumerated | discarded at pop | DNN bounds | phases | DNN iterations | CPU-s |\n|---|---|---|---|---|---|---|---|---|---|")
for dep in sorted(lev):
    L = lev[dep]; n = sum(L[k] for k in "BQDCE"); T.update(L); T["n"] += n
    print(f"| {dep} | {n} | {L['B']} | {L['Q']} / {L['D']} / {L['C']} | {L['E']} | {L['P']} | {L['dnn']} | {L['phases']} | {L['it']:,} | {L['cpu']:,.0f} |")
print(f"| all | {T['n']} | {T['B']} | {T['Q']} / {T['D']} / {T['C']} | {T['E']} | {T['P']} | {T['dnn']} | {T['phases']} | {T['it']:,} | {T['cpu']:,.0f} |")
print(f"\npath of the optimum: {len(onp)} nodes on it; violations (fathomed, or a DNN/cut bound above the optimum): {len(bad)}")
for dep, fx, v, out in sorted(onp)[:60]: print(f"  depth {dep} {out} bound {v:,.4f} ({opt - v:,.1f} under)  {fx}" if v is not None else f"  depth {dep} {out} {fx}")
for b in bad: print("  VIOLATION:", b)
