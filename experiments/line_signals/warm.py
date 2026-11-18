#!/usr/bin/env python3
"""The children's DNN bound warm from the node's state, capped: does a short warm solve rank the lines as
the cold first solve (a) does? One dnntest per node: the node solved (uncut, -J iterations), then every
child of the study's lines warm from its Z (-Q z), CAP iterations, no give-up, the certified bound printed
every 10 iterations (-v), so one run gives every cap up to CAP.

  python3 warm.py run tho30 5:2 [CAP]      -> dumps/warm_tho30_5-2.out
  python3 warm.py run tho40 26:3,31:7,37:0 [CAP]
"""
import glob, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.abspath(os.path.join(HERE, "..", ".."))
CFG = {"tho30": dict(par=f"{R}/experiments/tho30_solve/cutd2_cp.par", ub=149936, dat=f"{R}/qaplib/tho30.dat",
                     res=f"{R}/experiments/tho30_solve/line_truth/results"),
       "tho40": dict(par=f"{R}/infra/phase1/recipes/nug30.par", ub=240516, dat=f"{R}/qaplib/tho40t.dat",
                     res=f"{R}/experiments/tho40_line_truth/results")}


def study(inst, node):
    path = os.path.join(CFG[inst]["res"], node.replace(":", "-").replace(",", "_") + ".json")
    return json.load(open(path))


def run(inst, node, cap):
    c = CFG[inst]
    rec = study(inst, node)
    kids = []
    for L in rec["lines"].values():
        for k in L["children"]:
            if k not in kids:
                kids.append(k)
    out = os.path.join(HERE, "dumps", f"warm_{inst}_{node.replace(':', '-').replace(',', '-')}.out")
    cmd = [f"{R}/latest/dnntest", "-p", c["par"], "-x", node, "-u", str(c["ub"]), "-J", "1000", "-K", ",".join(kids),
           "-Q", "z", "-i", str(cap), "-H", "0", "-v", c["dat"]]
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with open(out, "w") as fh:
        fh.write("# " + " ".join(cmd) + "\n")
        fh.flush()
        subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=f"{R}/latest", check=True)


def parse(path):
    """{child 'i:j': [(iteration, lb, cpu), ...]} for the children (the parent's trajectory skipped)"""
    tr, cur = {}, None
    for ln in open(path):
        if ln.startswith("CHILD-SPEC "):
            cur = ln.split()[1]
            tr[cur] = []
        elif cur and ln.startswith("it "):
            f = ln.split()
            cpu = float(f[f.index("cpu") + 1])
            tr[cur].append((int(f[1]), float(f[3]), cpu))
        elif cur and ln.startswith("RESULT"):
            kv = dict(x.split("=", 1) for x in ln.split()[1:] if "=" in x)
            tr[cur].append((int(kv["iters"]), float(kv["lb"]), float(kv["cpu"])))
            cur = None
    return tr


def at_cap(traj, cap):
    """the certified bound and the CPU at iteration cap (a fathomed child stops earlier: its last values)"""
    best, cpu = -1e300, 0.0
    for it, lb, c in traj:
        if it > cap:
            break
        best, cpu = max(best, lb), c
    return best, cpu


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 200)


def analyze(caps=(10, 20, 30, 50, 100, 150, 200)):
    sys.path.insert(0, HERE)
    import signals
    for inst in ("tho30", "tho40"):
        for path in sorted(glob.glob(os.path.join(CFG[inst]["res"], "*.json"))):
            rec = json.load(open(path))
            node = rec["node"]
            out = os.path.join(HERE, "dumps", f"warm_{inst}_{node.replace(':', '-').replace(',', '-')}.out")
            if not os.path.exists(out):
                continue
            tr = parse(out)
            cutoff = rec["ub"] - (1.95 if rec["ub"] % 2 == 0 else 0.95)
            lines, truth, cold = [], [], []
            for line, L in rec["lines"].items():
                kids = list(L["children"])
                if not all(k in tr and tr[k] for k in kids):
                    continue
                lines.append(line)
                truth.append(1 + sum(L["children"][k]["truth"]["nodes"] for k in kids))
                cold.append(sum(max(0.0, cutoff - L["children"][k]["score"]["lb"]) for k in kids))
            if len(lines) < 3:
                print(f"{inst} {node}: {len(lines)} lines complete so far")
                continue
            neg = [-x for x in truth]
            row = [f"{inst} {node} ({len(lines)} lines): cold-300 {signals.spearman([-x for x in cold], neg):+.2f}"]
            best = truth.index(min(truth))
            for cap in caps:
                sh, cpu = [], []
                for line in lines:
                    kids = list(rec["lines"][line]["children"])
                    v = [at_cap(tr[k], cap) for k in kids]
                    sh.append(sum(max(0.0, cutoff - lb) for lb, _ in v))
                    cpu.append(sum(c for _, c in v))
                top = min(range(len(lines)), key=lambda i: sh[i])
                row.append(f"cap {cap}: {signals.spearman([-x for x in sh], neg):+.2f} top {lines[top]} "
                           f"{truth[top] / truth[best]:.2f}x, {sum(cpu) / len(lines):.0f} CPU s/line")
            print("\n  ".join(row))


if __name__ == "__main__" and sys.argv[1] == "analyze":
    analyze()
