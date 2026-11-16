"""Third study: the hand-over cost model refit, and the rule replayed on phase histories that include the
wc arm (the phase at penalty factor wc = 2.88: ../../cpfast/search/w18_wc2.88, w18_wc2.88_gap300, w38_wc2.88).

  python3 replay3.py            the hand-overs measured so far against three models, then the rules replayed
Model (per child, g = gap left in 0.1% of the cutoff):  quad3 = 20 + 2 g + 3 g^2 (the second study's key),
quad1.5 = 20 + 1.5 g^2, cubic = 20 + 0.22 g^3.
Rules: 'old' = cut.handover 20:2:3 (W 200, condition (i) only); 'new' = 20:0:0.22:1.5:100:3:1.25: from iteration
200, rate over the last 100; stop when (i) 1.5 x gap / rate > C(gap) and (ii) rate x dC/dgap < 1.25.
"""
import os, sys
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, H); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load
import hist as histmod
UNIT = 240.514
MODELS = {"quad3": lambda u: 20 + 2 * u + 3 * u * u, "quad1.5": lambda u: 20 + 1.5 * u * u, "cubic": lambda u: 20 + 0.22 * u ** 3}

def phases_from(d, name):
    nodes, dnn = load(d)
    lines = [l for l in open(os.path.join(d, "dnn.txt")) if " | first lb " in l]
    below = {}
    for x, e in zip(nodes, dnn):
        f = x["fixes"].split(",")
        for k in range(1, len(f)):
            v = below.setdefault(",".join(f[:k]), [0, 0]); v[0] += 1; v[1] += e["iters"] + e["cit"]
    out = []
    for x, e, ln in zip(nodes, dnn, lines):
        if not e["cut"]: continue
        h = [(int(a), e["cutoff"] - float(b)) for a, b in (t.split(":") for t in ln.split("| hist")[1].split())] if "| hist" in ln else []
        fath = e["cstatus"] == "fathomed"
        out.append(dict(run=name, depth=x["depth"], m=40 - x["depth"], pooled=e["kept"] > 0, short0=e["cutoff"] - e["lb"], hist=h, status=e["cstatus"],
                        its=e["cit"], left=0.0 if fath else e["cutoff"] - e["cbest"], after=None if fath else below.get(x["fixes"], [0, None])[1],
                        kids=below.get(x["fixes"], [0, 0])[0]))
    return out

def all_phases():
    P = []
    import glob
    for d in sorted(glob.glob(os.path.join(H, "runs", "*_*"))):
        n = os.path.basename(d)
        if n.startswith("path") or not os.path.exists(os.path.join(d, "dnn.txt")): continue
        P += phases_from(d, n)
    for n in ("w18_wc2.88", "w18_wc2.88_gap300", "w38_wc2.88"):
        P += phases_from(os.path.join(H, "..", "cpfast", "search", n), "cpfast/" + n)
    return P

def rule(p, n0, a, b, phi, W, pw, kappa):
    """the iteration at which cut.handover = n0:a:b:phi:W:pw:kappa stops the phase, its gap there; None if it does not"""
    h = p["hist"]
    for k, (it, g) in enumerate(h):
        if g < 0: return None
        if it >= p["its"] and p["status"] != "fathomed": return None
        if it < max(W, 200): continue
        g0 = p["short0"]
        for jt, gj in h:
            if jt > it - W: break
            g0 = gj
        rate = max(0.0, g0 - g) / W; u = g / UNIT
        C = p["m"] * (n0 + a * u + b * u ** pw); dC = p["m"] * (a + b * pw * u ** (pw - 1)) / UNIT
        if phi * g > rate * C and (kappa <= 0 or rate * dC < kappa): return it, g
    return None

if __name__ == "__main__":
    P = all_phases()
    if "--table" in sys.argv:
        print("| run | depth | gap left | children's iterations per child | quad3 | quad1.5 | cubic 20 + 0.22 g^3 |")
        print("|---|---|---|---|---|---|---|")
        for p in sorted(P, key=lambda p: p["left"]):
            if p["after"] and p["status"] != "fathomed" and p["kids"] >= p["m"] - 3:
                u = p["left"] / UNIT
                print("| `%s` | %d | %s | %.1f | %.0f | %.0f | %.0f |" % (p["run"], p["depth"], format(round(p["left"]), ","), p["after"] / p["m"], *[MODELS[k](u) for k in ("quad3", "quad1.5", "cubic")]))
        sys.exit(0)
    print("hand-overs whose children ran complete: gap left, children's iterations per child, and the three models")
    print("%-28s %5s %8s %3s %9s | %7s %7s %7s" % ("run", "depth", "gap left", "m", "per child", "quad3", "quad1.5", "cubic"))
    for p in sorted(P, key=lambda p: p["left"]):
        if p["after"] and p["status"] != "fathomed" and p["kids"] >= p["m"] - 3:
            u = p["left"] / UNIT
            print("%-28s %5d %8.1f %3d %9.1f | %7.1f %7.1f %7.1f" % (p["run"], p["depth"], p["left"], p["m"], p["after"] / p["m"], *[MODELS[k](u) for k in ("quad3", "quad1.5", "cubic")]))
    cub = MODELS["cubic"]
    print("\nthe rules replayed on every phase with a history; cost = phase iterations + (cubic model) children after a stop")
    for label, sel in (("phases at w (this directory's runs)", [p for p in P if p["hist"] and not p["run"].startswith("cpfast") and "hw" not in p["run"]]),
                       ("phases at wc = 2.88 (cpfast's z030 runs)", [p for p in P if p["hist"] and p["run"].startswith("cpfast/w18_wc2.88") or p["run"] == "cpfast/w38_wc2.88"]),
                       ("phases of the combined arms here (zallhw*)", [p for p in P if p["hist"] and "hw" in p["run"]])):
        if not sel: continue
        print("%s: %d phases, %d fathomed as recorded" % (label, len(sel), sum(p["status"] == "fathomed" for p in sel)))
        for name, args in (("recorded", None), ("old 20:2:3 (W 200)", (20, 2, 3, 1.5, 200, 2, 0)), ("cubic, (i) only, W 100", (20, 0, 0.22, 1.5, 100, 3, 0)),
                           ("new 20:0:0.22:1.5:100:3:1.25", (20, 0, 0.22, 1.5, 100, 3, 1.25)), ("new with W 200", (20, 0, 0.22, 1.5, 200, 3, 1.25))):
            tp = th = 0.0; lost = 0; det = []
            for p in sel:
                r = rule(p, *args) if args else None
                if r is None:
                    t, g = p["its"], p["left"]
                else:
                    t, g = r; lost += p["status"] == "fathomed"
                    det.append("%s d%d %.0f: stop %d at gap %.0f (recorded %s %d)" % (p["run"].split("/")[-1], p["depth"], p["short0"], t, g, p["status"], p["its"]))
                tp += t; th += p["m"] * cub(g / UNIT) if g > 0 else 0
            print("   %-30s phase it %6.0f + children %6.0f = %6.0f, fathoms lost %d" % (name, tp, th, tp + th, lost))
            if "-v" in sys.argv:
                for d in det: print("        " + d)
