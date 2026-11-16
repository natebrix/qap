"""The table of nug30.sh: python3 nug30_summary.py OUT ARM [ARM ...]"""
import os
import re
import sys

out, arms = sys.argv[1], sys.argv[2:]
rows = {}
for a in arms:
    d = os.path.join(out, a)
    nodes = sum(1 for l in open(os.path.join(d, "nodes.txt")) if l.startswith("node ") and l.split()[7] != "P")
    its = 0
    for l in open(os.path.join(d, "dnn.txt")):
        m = re.search(r"\| first lb \S+ primal \S+ status \w+ iters (\d+)", l)
        if m:
            its += int(m.group(1))
        c = re.search(r"\| cut 1 best \S+ rounds \d+ iters (\d+)", l)
        if c:
            its += int(c.group(1))
    p = re.search(r"wall ([\d.]+) s, CPU ([\d.]+) s.*best objective (\d+)", open(os.path.join(d, "log.txt")).read())
    rows[a] = (nodes, its, float(p.group(2)), float(p.group(1)), p.group(3)) if p else (nodes, its, float("nan"), float("nan"), "no end")
b = rows.get("base")
print(f"{'arm':>6} {'nodes':>7} {'DNN iterations':>15} {'CPU s':>8} {'wall s':>7} {'best':>7} {'base wall / arm':>16}")
for a, (n, i, c, w, o) in rows.items():
    print(f"{a:>6} {n:>7,} {i:>15,} {c:>8,.0f} {w:>7,.0f} {o:>7} {b[3] / w if b else float('nan'):>16.2f}")
print("cloud box, 4 jobs: base 3,688 nodes, 110,926 iterations, 3,030 CPU s, 771 s wall; g5e 1,302, 52,041, 1,667, 424")
print("'best' must read 6124 and each log.txt must end with no open node")
