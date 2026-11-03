#!/usr/bin/env python3
"""Summarize warm_y/runs: per child, the three starts' outcome (fathomed or
not), iterations, CPU, and the certified bound at iterations 50, 100, 200;
then totals over the fathoms all three starts share and over the misses.
Kill rule (roadmap/admm/splitting.md): a Y start is dead if it saves under
10% of the iterations to fathom AND gains under +10 at iteration 100."""
import glob, os, re, sys
D = os.path.dirname(os.path.abspath(__file__))
START = {"z": "Z only", "a": "Y restricted", "b": "Y conditioned"}
rows = {}
for f in sorted(glob.glob(os.path.join(D, "runs", "child_*.log")),
                key=lambda p: int(re.search(r"child_(\d+)", p).group(1))):
    k = int(re.search(r"child_(\d+)", f).group(1))
    cur = None; out = {}
    for line in open(f, errors="replace"):
        m = re.match(r"CHILD-PASS \d+ start=(\w)", line)
        if m: cur = m.group(1); out[cur] = {"trace": {}}; continue
        if cur is None: continue
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m: out[cur]["trace"][int(m.group(1))] = float(m.group(2)); continue
        if line.startswith("RESULT"):
            kv = dict(t.split("=", 1) for t in line.split()[1:] if "=" in t)
            out[cur].update(status=kv["status"], iters=int(kv["iters"]), cpu=float(kv["cpu"]))
    if all(s in out and "status" in out[s] for s in "zab"): rows[k] = out
def lb_at(r, it):
    t = r["trace"]
    return t.get(it, t[max(t)] if t else float("nan")) if r["status"] != "fathomed" or it <= r["iters"] else float("inf")
print(f"{len(rows)} children done")
print(f"{'child':>5} " + "  ".join(f"{START[s]:>30}" for s in "zab"))
for k, out in rows.items():
    cells = []
    for s in "zab":
        r = out[s]
        tag = "F" if r["status"] == "fathomed" else "-"
        cells.append(f"{tag} {r['iters']:5d} it {r['cpu']:6.0f}s  lb100 {lb_at(r,100):8.1f}")
    print(f"{k:>5} " + "  ".join(f"{c:>30}" for c in cells))
fath = [k for k, o in rows.items() if all(o[s]["status"] == "fathomed" for s in "zab")]
miss = [k for k, o in rows.items() if all(o[s]["status"] != "fathomed" for s in "zab")]
mixed = [k for k in rows if k not in fath and k not in miss]
print(f"\nfathomed by all three: {len(fath)}; missed by all three: {len(miss)}; outcome differs: {mixed}")
for s in "zab":
    it = sum(rows[k][s]["iters"] for k in fath); cpu = sum(rows[k][s]["cpu"] for k in fath)
    mcpu = sum(rows[k][s]["cpu"] for k in miss)
    def avg(at):
        v = [lb_at(rows[k][s], at) for k in miss]
        return sum(v) / len(v) if v else float("nan")
    print(f"{START[s]:>14}: fathoms {it:6d} iterations {cpu:7.0f} s;  misses {mcpu:7.0f} s, mean lb at 50/100/200 {avg(50):9.2f} {avg(100):9.2f} {avg(200):9.2f}")
