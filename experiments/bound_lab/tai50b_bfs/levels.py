"""After each level: the table by depth from the node and DNN traces (runs/dD_cK/), the path check, the stop
conditions, and README.md = README.head.md + the tables.   python3 levels.py   (exit 2: a stop condition)
Stop conditions: a node on the path of the best known solution fathomed or with a certified bound above
458,821,517; an incumbent below 458,821,517; a level's nodes or CPU more than 3x the high prediction (E3 at
gap scale 2.0, ../tai50b_size/logs/webcost.txt)."""
import sys, glob, os, re, collections
H = os.path.dirname(os.path.abspath(__file__)); BKV = 458821517
P = [int(x) for x in open(f"{H}/logs/tai50b.perm.txt").read().split()]
# predictions by depth, ../tai50b_size (made without wc): E2 (nodes, CPU-s) from logs/estimate.txt; E3 at gap scale 1.0
# and 2.0 (nodes, open, CPU-s) from logs/webcost.txt (depth 2 with its 588 children of the 12 depth-1 near misses)
E2 = {1: (50, 36321), 2: (654, 74635), 3: (141, 7510), 4: (135, 46558), 5: (3651, 136438), 6: (4043, 64626), 7: (129, 5288), 8: (130, 2057)}
E3 = {1: (50, 5, 29416), 2: (245 + 588, 14, 20700 + 17384), 3: (672, 11, 34732), 4: (517, 43, 70500), 5: (1978, 65, 100950), 6: (2925, 7, 221324), 7: (308, 5, 12015), 8: (215, 6, 3811)}
E3HI = {1: (50, 27121), 2: (980 + 588, 98444 + 17384), 3: (5568, 296717), 4: (5828, 655751), 5: (8142, 486533), 6: (21735, 1637615), 7: (1672, 60621), 8: (731, 13470)}
lev = collections.defaultdict(collections.Counter); onp = []; bad = []
for d in sorted(glob.glob(f"{H}/runs/d*_c*")):
    if os.path.exists(d + "/nodes.txt"):
        for ln in open(d + "/nodes.txt"):
            t = ln.split()
            if t[0] != "node": continue
            dep = int(t[1]); out = t[7]; L = lev[dep]; L[out] += 1; L["cpu"] += float(t[-1])
            fx = [] if t[2] == "-" else [tuple(map(int, x.split(":"))) for x in t[2].split(",")]
            if all(P[f] == l for f, l in fx):
                vals = [float(x) for x in t[4:6] if x != "-"]          # the DNN first solve's and the phase's certified values
                onp.append((dep, t[2], max(vals) if vals else None, out))
                if out in "QDCP" or (vals and max(vals) > BKV): bad.append(ln.strip())
    if os.path.exists(d + "/dnn.txt"):
        for ln in open(d + "/dnn.txt"):
            m = re.search(r"depth (\d+) .*\| first lb \S+ primal \S+ status \S+ iters (\d+) cpu (\S+)", ln)
            if not m: continue
            L = lev[int(m[1])]; L["dnn"] += 1; L["it"] += int(m[2])
            c = re.search(r"\| cut 1 best \S+ rounds \d+ iters (\d+) cpu (\S+)", ln)
            if c: L["phases"] += 1; L["it"] += int(c[1])
    dep = int(os.path.basename(d).split("_")[0][1:])
    if os.path.exists(d + "/err.txt"):
        e = open(d + "/err.txt").read(); lev[dep]["pretry"] += e.count("pool retry:"); lev[dep]["cretry"] += e.count("cold retry:")
found = sorted({ln.strip() for d in glob.glob(f"{H}/runs/d*_c*") for ln in open(d + "/log.txt") if ln.startswith("New objective value")})
better = [f for f in found if float(f.split(":")[1]) < BKV]
opens = {}
for f in glob.glob(f"{H}/open/open_*.txt"): opens[int(re.search(r"open_(\d+)", f)[1])] = sum(1 for _ in open(f))
path = {dep: (v, out) for dep, fx, v, out in onp}
out = ["| depth | nodes bounded | branched | fathomed: QPB / DNN first solve / cut phase | open nodes carried to the next depth | DNN bounds, phases | DNN iterations | CPU-s | retries: pool, cold | on-path node: short of 458,821,517 | E2 predicted: nodes, CPU-s | E3 (scale 1) predicted: nodes, branched, CPU-s |",
       "|---|---|---|---|---|---|---|---|---|---|---|---|"]
T = collections.Counter(); stop = []
for dep in sorted(lev):
    L = lev[dep]; n = sum(L[k] for k in "BQDCE"); T.update(L); T["n"] += n
    if dep not in opens: continue          # (a level still running is not shown)
    sh = f"{BKV - path[dep][0]:,.0f} ({path[dep][1]})" if dep in path and path[dep][0] is not None else "-"
    e2 = "%s, %s" % (f"{E2[dep][0]:,}", f"{E2[dep][1]:,}") if dep in E2 else ""
    e3 = "%s, %d, %s" % (f"{E3[dep][0]:,}", E3[dep][1], f"{E3[dep][2]:,}") if dep in E3 else ""
    out.append(f"| {dep} | {n:,} | {L['B']} | {L['Q']} / {L['D']} / {L['C']} | {opens.get(dep + 1, '?')} | {L['dnn']}, {L['phases']} | {L['it']:,} | {L['cpu']:,.0f} | {L['pretry']}, {L['cretry']} | {sh} | {e2} | {e3} |")
    if dep in E3HI and (n > 3 * E3HI[dep][0] or L["cpu"] > 3 * E3HI[dep][1]): stop.append(f"depth {dep}: {n} nodes, {L['cpu']:.0f} CPU-s against the high prediction {E3HI[dep]}")
done = [d for d in sorted(lev) if d in opens]
tn = sum(sum(lev[d][k] for k in "BQDCE") for d in done); tc = sum(lev[d]["cpu"] for d in done); ti = sum(lev[d]["it"] for d in done)
summ = [f"Levels complete: depths {done[0] if done else '-'} to {done[-1] if done else '-'}; {tn:,} nodes below the root, {ti:,} DNN iterations, {tc:,.0f} CPU-s = {tc/3600:.1f} core-hours (the root: 6,470 iterations, 4,341 CPU-s, run in `../tai50b_size`).",
        f"Open nodes now: {opens.get(done[-1] + 1, '?') if done else 50} at depth {done[-1] + 1 if done else 1} (`open/open_{done[-1] + 1 if done else 1}.txt`).",
        f"Incumbents found: {found or 'none yet (the search starts from 458,821,518 = best known + 1)'}.",
        f"Path of the best known solution: {len(onp)} nodes bounded on it, **{len(bad)} violations** (fathomed, or a certified bound above 458,821,517)."
        + (f" Largest certified bound on it: {max(v for _, _, v, _ in onp if v is not None):,.4f}." if any(v is not None for _, _, v, _ in onp) else "")]
if bad: stop += ["PATH VIOLATION: " + b for b in bad]
if better: stop += ["BETTER SOLUTION: " + b for b in better]
prog = open(f"{H}/progress.txt").read() if os.path.exists(f"{H}/progress.txt") else ""
txt = open(f"{H}/README.head.md").read().rstrip("\n") + "\n\n## Status\n\n" + "\n".join("- " + s for s in summ) + ("\n- **STOPPED:** " + "; ".join(stop) if stop else "") + \
      "\n\n## The levels [m]\n\n" + "\n".join(out) + "\n\nThe predictions are `../tai50b_size`'s, made without `wc=2.88` and from a 2,400-iteration phase limit under the same root; E3's depth 2 includes the 588 children of the 12 depth-1 near misses it prices separately.\n" + \
      "\n## progress.txt\n\n```\n" + prog + "```\n"
open(f"{H}/README.md", "w").write(txt); print("\n".join(summ + out + stop))
sys.exit(2 if stop else 0)
