#!/usr/bin/env python3
"""Score the lines of the tai60b line study (README.md) by their children, one row per arm
(runs/<arm>/nodes.txt) and production's (production_r1_children.tsv):

  python3 experiments/tai60b_lines/score.py

For each line: the depth-2 node's bound (the same in every arm: one parent state), how many of its
children fathom (D: the first solve, C: the cut phase, Q: QPB), how many branch, the child that
holds QAPLIB's best known permutation (the path child: it cannot fathom, and its subtree is most of
the node's) with its shortfall (cutoff - lb), the other branching children's shortfalls, and the
children's CPU. A line is better when its path child is closer to the cutoff and fewer of its other
children branch, and those that do are closer.
"""
import glob
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
N = 60


def best_known():
    t = open(os.path.join(ROOT, "QAPSOLN", "tai60b.sln")).read().split()
    p = [int(x) - 1 for x in t[2:2 + N]]      # QAPLIB: 1-based locations, facility i at p[i]
    return p


def val(s):
    try:
        v = float(s)
        return None if math.isnan(v) else v
    except ValueError:
        return None


def read_arm(path):
    """the depth-2 node and its depth-3 children from a node trace (a node's last line: a resumed arm
    traces the nodes after its checkpoint again)"""
    parent, kids = None, {}
    for ln in open(path):
        f = ln.split()
        if len(f) < 12 or f[0] != "node":
            continue
        d = int(f[1])
        lbs = [v for v in map(val, f[3:6]) if v is not None]
        rec = {"fix": f[2].split(",")[-1], "fixes": f[2], "lb": max(lbs) if lbs else None, "cutoff": float(f[6]),
               "out": f[7], "line": f[9], "cpu": float(f[11])}
        if d == 2:
            parent = rec
        elif d == 3:
            kids[rec["fixes"]] = rec
    return parent, list(kids.values())


def read_production():
    parent, kids = None, []
    for ln in open(os.path.join(HERE, "production_r1_children.tsv")):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        rec = {"fix": f[0], "out": f[1], "lb": float(f[2]), "cutoff": float(f[2]) + float(f[3]), "cpu": float(f[4]),
               "size": int(f[5])}
        if f[0] == "parent":
            rec["line"] = f[7].split("=")[1]
            parent = rec
        else:
            kids.append(rec)
    return parent, kids


def on_path(fix, p):
    i, j = map(int, fix.split(":"))
    return p[i] == j


def row(name, parent, kids, p):
    br = [k for k in kids if k["out"] == "B"]
    path = [k for k in kids if on_path(k["fix"], p)]
    off = sorted((k["cutoff"] - k["lb"] for k in br if not on_path(k["fix"], p)), reverse=True)
    pc = path[0] if path else None
    fat = {o: sum(1 for k in kids if k["out"] == o) for o in "DCQ"}
    return (f"| {name} | {parent['line']} | {parent['cutoff'] - parent['lb']:,.0f} | {len(kids)} | "
            f"{fat['D']} / {fat['C']} / {fat['Q']} | {len(br)} | "
            f"{(pc['fix'] + ' ' + pc['out'] + ' ' + format(pc['cutoff'] - pc['lb'], ',.0f')) if pc else '-'} | "
            f"{len(off)}: {', '.join(format(v / 1e3, ',.0f') + 'k' for v in off)} | "
            f"{sum(k['cpu'] for k in kids):,.0f} |")


def main():
    p = best_known()
    print("| arm | line | depth-2 node short | children | fathomed D / C / Q | branched | path child, short | "
          "other branched: short | children CPU s |")
    print("|---|---|---|---|---|---|---|---|---|")
    pp, pk = read_production()
    print(row("production (warm chain from the root)", pp, pk, p))
    for path in sorted(glob.glob(os.path.join(HERE, "runs", "*", "nodes.txt"))):
        name = os.path.basename(os.path.dirname(path))
        if name == "parent" or "stopped" in open(os.path.join(os.path.dirname(path), "cmd.txt")).read():
            continue
        parent, kids = read_arm(path)
        if parent is None:
            continue
        done = "exit" in open(os.path.join(os.path.dirname(path), "cmd.txt")).read()
        print(row(name, parent, kids, p) + ("" if done else " (running)"))
    fidelity(pk)


def fidelity(pk):
    """the r1 arm against production's r1 children, child by child: the same fixes, bounded from a
    weaker parent state (production's 43:58 came warm down the chain from the root)"""
    path = os.path.join(HERE, "runs", "r1", "nodes.txt")
    if not os.path.exists(path):
        return
    prod = {k["fix"]: k for k in pk}
    _, kids = read_arm(path)
    print("\nr1 arm against production, child by child (short = cutoff - lb; negative: fathomed with room):\n")
    print("| child | arm | arm short | production | production short | arm - production | arm CPU s | production CPU s | production subtree |")
    print("|---|---|---|---|---|---|---|---|---|")
    same = 0
    for k in kids:
        q = prod[k["fix"]]
        same += (k["out"] == "B") == (q["out"] == "B")
        a, b = k["cutoff"] - k["lb"], q["cutoff"] - q["lb"]
        print(f"| {k['fix']} | {k['out']} | {a:,.0f} | {q['out']} | {b:,.0f} | {a - b:,.0f} | {k['cpu']:,.0f} | {q['cpu']:,.0f} | {q['size']:,} |")
    print(f"\n{same} of {len(kids)} children branch or fathom as in production.")


if __name__ == "__main__":
    main()
