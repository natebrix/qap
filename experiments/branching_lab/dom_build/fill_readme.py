#!/usr/bin/env python3
"""Fill README.md's @@...@@ placeholders from the logs (run once, at the end)."""
import json, os, re, subprocess
D = os.path.dirname(os.path.abspath(__file__))
os.chdir(D)
r = open("README.md").read()


def rows(path):
    by = {}
    for l in open(path):
        x = json.loads(l); by.setdefault(x["instance"], {})[x["key"]] = x
    return by


o1 = rows("logs/valid_opt1.jsonl")
for k in ("lipa10a", "lipa10b"):
    o1.pop(k, None)
o1.update(rows("logs/valid_opt1_lipa.jsonl"))
po = rows("logs/valid_poor.jsonl")
closed = {k: d for k, d in o1.items() if 0 in d and 1 in d and d[0]["valid"] and d[1]["valid"]}
notclosed = sorted(k for k, d in o1.items() if k not in closed)
bad = [k for k, d in list(o1.items()) + list(po.items()) if 1 in d and d[1]["status"] == "ok" and not d[1]["valid"]]
assert not bad, bad
names = set(closed) | {k for k, d in po.items() if d[0]["valid"] and d[1]["valid"]}
fired = {k for k in names if (closed.get(k, {}).get(1, {}).get("dropped") or 0) + (po.get(k, {}).get(1, {}).get("dropped") or 0) > 0}
nruns = len(closed) + sum(1 for d in po.values() if d[1]["valid"])


def table(by, keep=None):
    out = ["| instance | n | optimum (source) | start | sym | nodes off | nodes on | on/off | dropped | final off / on | open | valid |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    tot = [0, 0, 0]
    for k, d in by.items():
        if 0 not in d or 1 not in d: continue
        a, b = d[0], d[1]
        if keep and not keep(k): continue
        ok = a["valid"] and b["valid"]
        src = a["opt_src"] + (f", {a['n_optima']} optima" if a.get("n_optima") else "")
        no, nn = a.get("nodes"), b.get("nodes")
        st = "yes" if ok else f"timeout ({a['status']} / {b['status']})"
        out.append(f"| {k} | {a['n']} | {a['opt']} ({src}) | {a.get('ub_used') or 'none':.0f} | {'yes' if a['sym_file'] else ''} | "
                   f"{no if no is not None else '-'} | {nn if nn is not None else '-'} | {(nn / no if no and nn else float('nan')):.3f} | "
                   f"{b.get('dropped')} | {a.get('final_incumbent')} / {b.get('final_incumbent')} | {a.get('open_nodes')} / {b.get('open_nodes')} | {st} |")
        if ok: tot[0] += no; tot[1] += nn; tot[2] += b.get("dropped") or 0
    return "\n".join(out), tot


car = lambda k: not k.startswith("car10")
t1, tot1 = table(o1, car)
_, tot1all = table(o1)
t1 += (f"\n\nThe 30 car10 instances (linear term C) are left out of the rows above: all valid, the test firing on {sum(1 for k, d in closed.items() if k.startswith('car10') and d[1]['dropped'])} of them. "
       f"Over the {len(closed)} closed instances: {tot1all[0]:,} nodes off, {tot1all[1]:,} on ({tot1all[1] / tot1all[0]:.3f}), "
       f"{tot1all[2]:,} dropped.")
sel = lambda k: k in ("had10", "had12", "had14", "nug12", "nug15", "rou12", "rou15", "tai10a", "tai10b", "tai12a", "tai12b",
                      "chr12b", "chr12c", "esc08b", "esc08d", "scr15", "lipa10a", "lipa10b", "law07")
t2, _ = table(po, sel)
r = r.replace("@@TABLE1@@", t1).replace("@@TABLE2@@", t2)
r = r.replace("@@NINST@@", str(len(names))).replace("@@NFIRED@@", str(len(fired))).replace("@@NRUNS@@", str(nruns))
r = r.replace("@@NOTCLOSED@@", (", ".join(notclosed) if notclosed else "none"))
r = r.replace("@@DESIGN@@", open("design.md.part").read().strip())
open("README.md", "w").write(r)
print("closed", len(closed), "instances", len(names), "fired", len(fired), "runs", nruns, "notclosed", notclosed)
print("left:", sorted(set(re.findall(r"@@\w+@@", r))))
