import json, os, sys
# cmp_level.py LEVEL RUN_A RUN_B: the draws at one level of two probes of the same root, side by side.
# RUN_x is a probe's own directory in a probe.py run (it reads L<level>/draws.json and t*/result.json,
# which are not committed). Used for the matched depth-4 children of Cp396 and Zp396.
L = int(sys.argv[1]); dirs = sys.argv[2:]
def load(d):
    dr = json.load(open(os.path.join(d, f"L{L}", "draws.json")))
    out = []
    for j, u in enumerate(dr["draws"]):
        p = os.path.join(d, f"L{L}", f"t{j}", "result.json")
        if not os.path.exists(p):
            out.append((u, None)); continue
        doc = json.load(open(p)); tr = doc["trace"]
        first = tr[0]["cutoff"] - tr[0]["lb"] if tr else float("nan")
        best = min(t["short"] for t in tr) if tr else float("nan")
        ph = sum(1 for t in tr if t["cut"])
        out.append((u, (first, best, ph, int(not doc["children"]), doc["dnn_cpu"], len(doc["children"]))))
    return dr, out
R = [load(d) for d in dirs]
print("draws equal:", R[0][0]["draws"] == R[1][0]["draws"])
tot = [[0, 0, 0, 0] for _ in R]
for j in range(len(R[0][1])):
    cells = []
    for k, (dr, out) in enumerate(R):
        u, v = out[j]
        if v is None: cells.append(f"{str(u):>14} (running)"); continue
        f, b, ph, fa, cpu, ch = v
        tot[k][0] += fa; tot[k][1] += ph; tot[k][2] += cpu; tot[k][3] += ch
        cells.append(f"{str(u):>14} first {f:6.1f} best {b:6.1f} ph {ph} fath {fa} cpu {cpu:6.0f} kids {ch:2d}")
    print(" | ".join(cells))
for k, d in enumerate(dirs):
    print(os.path.basename(d), "fathomed %d phases %d cpu %.0f children %d" % tuple(tot[k]))
