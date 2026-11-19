#!/usr/bin/env python3
"""The solver's pair test (domcheck, latest/dominance.C, on the symmetrised data
as bnbdrv holds it) against the Python reference (dominance/dom.py, on the
original, unsymmetrised data), on every node of a trace.

  python3 cmp_replay.py LABEL INSTANCE TRACE.gz OUT.txt

Checks for every node with >= 2 fixes: 4 x dom.py's best pair min Delta equals
domcheck's value exactly, and the dominated flags agree; prints a summary
(nodes, dominated, mismatches) and, for the highest dominated nodes, the
subtree sizes the replay predicts (as dominance/analyze.py: the trace's
nodes under each highest flagged node).
"""
import gzip, subprocess, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dominance"))
from qapdata import read_dat
from dom import Inst

label, name, trace, out = sys.argv[1:5]
dat = f"/home/user/wt-dombuild/qaplib/{name}.dat"
p = subprocess.Popen(f"zcat {trace} | {HERE}/bin/domcheck {dat} | grep -v -e '^objective values' -e '^symmetrized'", shell=True, stdout=subprocess.PIPE, text=True)
inst = Inst(*read_dat(name))
nodes = dom = mism = flagmism = 0
keys, flagged, cpu = [], {}, {}
op = gzip.open if trace.endswith(".gz") else open
tr = [ln.split() for ln in op(trace, "rt") if ln.startswith("node ")]
with open(out, "w") as f:
    for t, ln in zip(tr, p.stdout):
        fixes, v4, pair, d = ln.split()
        assert fixes == t[2]
        keys.append(fixes); cpu[fixes] = float(t[11])
        if v4 == "-":
            continue
        F = dict(tuple(map(int, s.split(":"))) for s in fixes.split(","))
        pm, pr = inst.pair_test(F)
        nodes += 1
        if 4 * pm != int(v4):
            mism += 1
            f.write(f"MISMATCH {fixes} py4={4*pm} c4={v4}\n")
        if (pm >= 1) != (d == "1") or (int(v4) > 0) != (d == "1"):
            flagmism += 1
            f.write(f"FLAGMISMATCH {fixes} py={pm} c4={v4} d={d}\n")
        if d == "1":
            dom += 1; flagged[fixes] = int(v4)
# the highest flagged nodes and the trace nodes under them
top = [k for k in flagged if not any(k.startswith(j + ",") for j in [",".join(k.split(",")[:i]) for i in range(1, len(k.split(",")))] if j in flagged)]
under = {k: 0 for k in top}
undercpu = {k: 0.0 for k in top}
topset = set(top)
for k in keys:
    parts = k.split(",") if k != "-" else []
    for i in range(1, len(parts) + 1):
        pre = ",".join(parts[:i])
        if pre in topset:
            under[pre] += 1; undercpu[pre] += cpu[k]; break
saved = sum(under.values()); savedcpu = sum(undercpu.values())
tot = len(keys); totcpu = sum(cpu.values())
summary = (f"{label}: {tot} trace nodes, {nodes} with >= 2 fixes compared, value mismatches {mism}, "
           f"flag mismatches {flagmism}; dominated {dom}, highest dominated {len(top)}, "
           f"trace nodes under them {saved} ({100*saved/tot:.2f}%), CPU {savedcpu/3600:.3f} h of {totcpu/3600:.3f} h "
           f"({100*savedcpu/max(totcpu,1e-9):.2f}%)")
with open(out, "a") as f:
    f.write(summary + "\n")
    for k in sorted(top, key=lambda k: -under[k])[:40]:
        f.write(f"top {len(k.split(','))} {k} 4min={flagged[k]} subtree={under[k]} cpu={undercpu[k]:.1f}\n")
print(summary)
