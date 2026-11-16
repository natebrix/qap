"""A guided dive toward a known solution, from the incumbent a run was given.

tho40-b was given the best known value + 2 (240,518; tho40's objective values are
even), as the 2000 nug30 solve was given the optimum + 2, so it must find a solution
of 240,516 itself. It can only if no node on the way to one is fathomed: a bound
above a subtree's optimum anywhere on that path would lose the solution, and the
run would end without it. This follows that path one node at a time, with the run's
settings, and checks it.

One bnbdrv process per node, run to the node's own depth (-d), so the solver bounds
it with the recipe (the DNN bound at its depth, the branching rule of its level, its
own choice of line) and stops; its children go to the -o state file. The next node
is the child on the chosen line that agrees with the solution:

  --warm  (as the run): the parent's -o state, cut to the path child alone, is the
          next process's -r state, so each node starts from its parent's shipped DNN
          state (ship.bits) as a coordinator unit does
  default (cold): the next node is started by -F with the path so far, no parent state

The solution is QAPSOLN/<sln> (a permutation of the instance it names; --transposed
when the run's instance swaps A and B, as tho40t.dat does tho40.dat, so the solution
is inverted). The symmetry file keeps one image of it under the instance's grid
symmetries: the one whose location for the root's branching facility is in the
file's root set (J1). Only a grid of r x c locations, numbered row-major, is handled.

Passes when every node on the path branches (B), keeps the path child, and the last
one, enumerated (E), yields the solution's value. Writes OUTDIR/path.tsv and each
node's log, traces and state.

usage (from anywhere; tho40-b's settings are the defaults):
  python3 dive.py OUTDIR [--warm]
  python3 dive.py OUTDIR --instance qaplib/tho40t.dat --sln tho40.sln --transposed \
      --params ../infra/phase1/recipes/tho40.par --sym SYM/tho40t.sym --grid 5x8 --plus 2
"""
import argparse
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
LATEST = os.path.join(ROOT, "latest")


def read_sym_root(path):
    """the root set J1 of a symmetry file (0-based locations); see latest/branch.C read_symmetry"""
    t = open(path).read().split()
    k = int(t[0])
    return {int(x) - 1 for x in t[1:1 + k]}


def images(q, rows, cols):
    """the solution under the grid's four symmetries (row flip, column flip)"""
    out = []
    for fr in (0, 1):
        for fc in (0, 1):
            def g(l):
                r, c = divmod(l, cols)
                return (rows - 1 - r if fr else r) * cols + (cols - 1 - c if fc else c)
            out.append([g(l) for l in q])
    return out


def state_nodes(path, n):
    nodes, inc = [], None
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if line.startswith("node "):
                nodes.append((int(t[1]), float(t[3]), list(map(int, t[5:5 + n]))))
            elif line.startswith("incumbent "):
                inc = (float(t[1]), list(map(int, t[3:3 + n])))   # obj, ifs_found, perm
    return nodes, inc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--warm", action="store_true")
    ap.add_argument("--instance", default="qaplib/tho40t.dat", help="relative to the repository")
    ap.add_argument("--sln", default="tho40.sln", help="in QAPSOLN/")
    ap.add_argument("--transposed", action="store_true", default=None)
    ap.add_argument("--params", default="../infra/phase1/recipes/tho40.par", help="relative to latest/")
    ap.add_argument("--sym", default="SYM/tho40t.sym", help="relative to latest/")
    ap.add_argument("--grid", default="5x8")
    ap.add_argument("--plus", type=int, default=2, help="the incumbent given: the solution's value + this")
    ap.add_argument("--resume", action="store_true",
                    help="reuse the nodes OUT already holds (their trace, state and wall time), so an interrupted "
                         "dive goes on where it stopped; the same options as the first run")
    a = ap.parse_args()
    if a.transposed is None:
        a.transposed = a.instance == "qaplib/tho40t.dat"
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)

    s = open(os.path.join(ROOT, "QAPSOLN", a.sln)).read().split()
    n, best = int(s[0]), int(s[1])
    p = [int(x) - 1 for x in s[2:2 + n]]
    q = p
    if a.transposed:
        q = [0] * n
        for i, j in enumerate(p):
            q[j] = i
    rows, cols = map(int, a.grid.split("x"))
    J1 = read_sym_root(os.path.join(LATEST, a.sym))
    cands = images(q, rows, cols)
    qq = None   # the image the symmetry keeps, chosen once the root's line is known

    prefix = []
    summary = open(f"{out}/path.tsv", "w")
    summary.write("depth\tfix\toutcome\tqpb\tdnn\tlb\tslack\tline\tkids\tfree\tpath_child_kept\tchild_pred\tcpu\twall\n")
    verdict = 1
    while True:
        k = len(prefix)
        nt, dt, st, lg = (f"{out}/d{k:02d}.{x}" for x in ("nodes", "dnn", "state", "log"))
        wf = f"{out}/d{k:02d}.wall"   # written last: a node with it finished
        if a.resume and all(os.path.exists(f) for f in (nt, st, wf)):
            wall = float(open(wf).read())
        else:
            for f in (nt, dt, wf):
                if os.path.exists(f):
                    os.remove(f)
            cmd = ["./bnbdrv", "-p", a.params, "-s", a.sym, "-u", str(best + a.plus), "-d", str(k), "-o", st]
            if prefix and a.warm:
                cmd += ["-r", f"{out}/d{k:02d}.in.state"]
            elif prefix:
                cmd += ["-F", ",".join(f"{i}:{j}" for i, j in prefix)]
            cmd.append(os.path.join(ROOT, a.instance))
            t0 = time.time()
            with open(lg, "w") as log:
                subprocess.run(cmd, cwd=LATEST, env=dict(os.environ, QAP_NODE_TRACE=nt, QAP_DNN_TRACE=dt),
                               stdout=log, stderr=subprocess.STDOUT, check=False)
            wall = time.time() - t0
            with open(wf, "w") as w:
                w.write(f"{wall:.1f}\n")
        mine = [l.split() for l in open(nt) if l.startswith("node ") and int(l.split()[1]) == k] \
            if os.path.exists(nt) else []
        if len(mine) != 1:
            print(f"depth {k}: expected one node line, found {len(mine)}; see {lg}", flush=True)
            return 2
        _, _, fixes, qpb, dnn, cut, cutoff, outcome, rule, line, kids, cpu = mine[0][:12]
        num = lambda v: float("nan") if v == "-" else float(v)
        vals = [v for v in (num(qpb), num(dnn), num(cut)) if v == v]
        lb = max(vals) if vals else float("nan")   # an enumerated node (E) has no bound
        nodes, inc = state_nodes(st, n)
        kept = pred = ""
        nxt = None
        if outcome == "B":
            idx = int(line[1:])
            if qq is None:
                keep = [c for c in cands if line[0] == "r" and c[idx] in J1]
                if len(keep) != 1:
                    print(f"the root branched on {line}: {len(keep)} images of the solution in its root set; "
                          "a symmetry axis through it is not handled", flush=True)
                    return 2
                qq = keep[0]
            inv = [0] * n
            for f, l in enumerate(qq):
                inv[l] = f
            nxt = (idx, qq[idx]) if line[0] == "r" else (inv[idx], idx)
            match = [nd for nd in nodes if nd[0] == k + 1 and nd[2][nxt[0]] == nxt[1]]
            kept = "yes" if match else "NO"
            pred = f"{match[0][1]:.2f}" if match else ""
        row = [k, f"{prefix[-1][0]}:{prefix[-1][1]}" if prefix else "root", outcome, qpb, dnn,
               "-" if lb != lb else f"{lb:.4f}", "-" if lb != lb else f"{best - lb:.2f}", line, kids, n - k,
               kept or "-", pred or "-", cpu, f"{wall:.1f}"]
        summary.write("\t".join(map(str, row)) + "\n")
        summary.flush()
        print("\t".join(map(str, row)), flush=True)
        if outcome == "E":
            ok = inc is not None and inc[0] == best
            print(f"enumerated at depth {k}: incumbent {inc[0]:.0f}, "
                  f"{'the solution' if ok and inc[1] == qq else 'another permutation' if ok else 'NOT the solution value'}")
            print("PASS" if ok else "FAIL")
            verdict = 0 if ok else 1
            break
        if outcome != "B":
            print(f"FAIL: the path node at depth {k} was fathomed ({outcome}), bound {lb:.4f}, cutoff {cutoff}")
            break
        if kept != "yes":
            print(f"FAIL: the path child {nxt[0]}:{nxt[1]} was not kept at depth {k}")
            break
        prefix.append(nxt)
        if a.warm:   # the next process's state: this one's, with the path child as its only node
            with open(st) as fh, open(f"{out}/d{k + 1:02d}.in.state", "w") as w:
                for l in fh:
                    if l.startswith("nodes "):
                        continue
                    if l.startswith("node "):
                        t = l.split()
                        if int(t[1]) == k + 1 and int(t[5 + nxt[0]]) == nxt[1]:
                            w.write("nodes 1\n" + l)
                        continue
                    w.write(l)
    return verdict


if __name__ == "__main__":
    sys.exit(main())
