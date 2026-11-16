"""The dives of a set, node by node, from the DNN traces:  python3 dives.py runs/SET [-nodes]
A dive = consecutive trace lines of one process with increasing depth. Prints one line per dive
(depth:outcome/iterations/final shortfall in M; D, C fathomed by the first solve / the phase, B branched),
then per depth: nodes seen, fathomed, mean iterations and CPU-s of fathomed and of branched nodes."""
import sys, glob, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N
def dives(setdir):
    out = []
    for f in sorted(glob.glob(setdir + "/s*/dnn.txt")):
        cur = []
        for d in N.load([f]):
            if cur and d["depth"] <= cur[-1]["depth"]: out.append((f, cur)); cur = []
            cur.append(d)
        if cur: out.append((f, cur))
    return out
if __name__ == "__main__":
    ds = dives(sys.argv[1]); by = collections.defaultdict(list)
    for f, c in ds:
        done = c[-1]["fath"] == 1
        print(os.path.basename(os.path.dirname(f)), "" if done else "(unfinished)", " -> ".join(
            "%d:%s%s/%d/%.2f" % (d["depth"], d["out"], "*" if d["onpath"] else "", d["fit"] + d["cit"], (d["cutoff"] - d["final"]) / 1e6) for d in c),
            "| CPU-s %.0f" % sum(d["fcpu"] + d["ccpu"] for d in c))
        for d in c: by[d["depth"]].append(d)
    print("\n| depth | nodes | fathomed: first solve / phase | branched | fathomed: mean iterations, CPU-s | branched: mean iterations, CPU-s, mean final shortfall |")
    print("|---|---|---|---|---|---|")
    for dep in sorted(by):
        L = by[dep]; F = [d for d in L if d["fath"]]; B = [d for d in L if not d["fath"]]
        m = lambda X, k: sum(k(d) for d in X) / len(X) if X else float("nan")
        print("| %d | %d | %d / %d | %d | %.0f, %.0f | %.0f, %.0f, %s |" % (dep, len(L), sum(d["out"] == "D" for d in L), sum(d["out"] == "C" for d in L), len(B),
              m(F, lambda d: d["fit"] + d["cit"]), m(F, lambda d: d["fcpu"] + d["ccpu"]), m(B, lambda d: d["fit"] + d["cit"]), m(B, lambda d: d["fcpu"] + d["ccpu"]),
              f"{m(B, lambda d: d['cutoff'] - d['final']):,.0f}"))
    if "-nodes" in sys.argv:
        for f, c in ds:
            for d in c: print(d["depth"], ",".join(d["fixes"]), d["out"], d["fit"], d["cit"], round(d["cutoff"] - d["lb"]), round(d["cutoff"] - d["final"]), round(d["fcpu"] + d["ccpu"]))
