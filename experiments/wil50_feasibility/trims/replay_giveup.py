# Replay the search's give-up rule (latest/dnn.C: extrapolation of the running
# best certified value from hopeless_after on, checked every 10 iterations) on
# the depth-4 children's uncut traces (stronger_bounds/depth4/uncut_*.log,
# dnntest -v -H 0 -i 1500), with and without the zone deferral
# (QAP_DNN_GU_ZONE: a give-up within zone*|cutoff| of the cutoff is deferred
# hopeless_delay iterations, once; it stands only if the rule fires again).
# For each node: when the rule fires, when the deferred give-up would stand or
# the node fathoms, and the CPU of each policy against the phase that ran.
# Usage: python3 replay_giveup.py H zone delay [delay...]
import re, sys, glob, os
H = int(sys.argv[1]); zone = float(sys.argv[2]); delays = [int(d) for d in sys.argv[3:]]
D = os.path.join(os.path.dirname(__file__), "..", "stronger_bounds", "depth4")
def trace(path):
    cut = None; pts = []
    for line in open(path, errors="replace"):
        m = re.search(r"cutoff=([0-9.]+)", line)
        if m: cut = float(m.group(1))
        m = re.match(r"it\s+(\d+)\s+lb\s+([0-9.]+).*?cpu\s+([0-9.]+)", line)
        if m: pts.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return cut, pts
def run(pts, cut, H, zone, delay):
    """returns (outcome, iteration, cpu): 'fathom', 'give-up' or 'limit'"""
    best = []; confirm_at = -1; w = 5
    for it, lb, cpu in pts:
        best.append(lb)
        if lb > cut: return "fathom", it, cpu
        if H > 0 and it >= H and it >= confirm_at:
            h = len(best); fire = False
            if h > 2 * w:
                d1 = best[-1] - best[-1 - w]; d2 = best[-1 - w] - best[-1 - 2 * w]
                r = min(d1 / d2, 0.95) if d2 > 0 else (0.95 if d1 > 0 else 0.0)
                if lb + 5 * d1 * r / (1 - r) < cut: fire = True
            if fire:
                if zone > 0 and confirm_at < 0 and cut - lb <= zone * max(1.0, abs(cut)):
                    confirm_at = it + delay
                else:
                    return "give-up", it, cpu
    it, lb, cpu = pts[-1]
    return "limit", it, cpu
phase = {}
for line in open(os.path.join(D, "nodes.txt")):
    pass
# phase cost from the cuts logs: total cpu of the rounds that ran
for f in glob.glob(os.path.join(D, "cuts_*.log")):
    k = int(re.search(r"cuts_(\d+)", f).group(1))
    txt = open(f, errors="replace").read()
    cpus = [float(x) for x in re.findall(r"^round\s+\d+.*?cpu\s+([0-9.]+)", txt, re.M)]
    q = re.search(r"^quiet round.*?cpu\s+([0-9.]+)", txt, re.M)
    if q: cpus.append(float(q.group(1)))
    fath = "fathomed" in txt
    phase[k] = (sum(cpus), fath, len(cpus))
hdr = f"{'node':>4} {'no zone':>22} | " + " | ".join(f"{'zone '+str(zone)+' delay '+str(d):>24}" for d in delays) + f" | {'phase from 400 (cuts_k.log)':>28}"
print(hdr)
tot = {None: 0.0}; tot.update({d: 0.0 for d in delays}); tot_phase = 0.0
for f in sorted(glob.glob(os.path.join(D, "uncut_*.log")), key=lambda p: int(re.search(r"uncut_(\d+)", p).group(1))):
    k = int(re.search(r"uncut_(\d+)", f).group(1)); cut, pts = trace(f)
    base = run(pts, cut, H, 0, 0)
    pcpu, pf, pr = phase.get(k, (0, False, 0))
    # policy cost: the uncut iterations to the decision; a give-up then pays the phase
    def cost(o): return o[2] + (pcpu if o[0] == "give-up" else 0)
    cells = [f"{base[0]:7s}@{base[1]:4d} {cost(base):7.0f}s"]
    tot[None] += cost(base)
    for d in delays:
        o = run(pts, cut, H, zone, d); cells.append(f"{o[0]:7s}@{o[1]:4d} {cost(o):7.0f}s")
        tot[d] += cost(o)
    tot_phase += pcpu
    print(f"{k:>4} {cells[0]:>22} | " + " | ".join(f"{c:>24}" for c in cells[1:]) + f" | {pcpu:7.0f}s {'fathom' if pf else 'miss  '} {pr:2d} rounds")
print(f"{'sum':>4} {tot[None]:>21.0f}s | " + " | ".join(f"{tot[d]:>23.0f}s" for d in delays) + f" | {tot_phase:7.0f}s")
print("(a give-up at iteration t costs its uncut CPU plus the logged phase; a fathom costs its uncut CPU; the phase column is the logged cuts run from 400)")
