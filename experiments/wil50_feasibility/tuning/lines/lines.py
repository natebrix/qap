#!/usr/bin/env python3
"""Which branching line at a node leaves the smallest subtree below it?

  python3 experiments/wil50_feasibility/tuning/lines/lines.py qaplib/wil50.dat \\
      --ub 48816 --parents experiments/wil50_feasibility/tuning/lines/parents.txt \\
      --jobs 8 --work /tmp/wil50_lines

For each parent node (one "id spec" per line of --parents):

1. Parent: the certified DNN bound with the primal marginals x_ij (dnntest -X,
   --parent-iters, no give-up), and QPB's dual matrix U (solver console, `u`).
2. Candidate lines (a row = branch on a facility, a column = on a location):
   - ecc:      the most eccentric free location (largest distance to another
               free location, as rule 6);
   - U:        the 2 rows and 2 columns with the largest U sums (rule 4's
               candidate set; its own choice among them is reported from the
               children's QPB bounds, as lookahead would make it);
   - xbar:     the row and the column whose marginals are most concentrated
               (largest sum of x_ij^2): most of their children are
               assignments the relaxation puts near zero;
   - random:   --random rows and as many columns, as a baseline.
3. Children: every child of every candidate line, bounded by the DNN bound
   with a fixed budget (--child-iters, no give-up, so all children are
   comparable) and by QPB (console `b`).
4. Score: a child short of the cutoff by g needs about g / --gain more levels,
   so it leaves a subtree of about --branch ** (g / --gain) nodes; a fathomed
   child leaves none. A line's cost is the sum over its children. Also
   reported: children fathomed, children within one level (g < gain), and the
   rank correlation of x_ij and of child QPB with the child's gap, pooled over
   all children (does x_ij, known before branching, predict a child's gap?).

Results are cached in --work (one JSON per child), so a rerun resumes. DNN
runs use --dnn-flags (default: the solver's mode with PRSM and degree 4).
"""
import argparse, concurrent.futures as cf, json, math, os, random, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))))
DNNTEST = os.path.join(ROOT, "latest", "dnntest")
BNBDRV = os.path.join(ROOT, "latest", "bnbdrv")
ENV = dict(os.environ, VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")


def read_qap(path):
    t = open(path).read().split()
    n = int(t[0]); v = [float(x) for x in t[1:1 + 2 * n * n]]
    return n, [v[i * n:(i + 1) * n] for i in range(n)], [v[n * n + i * n:n * n + (i + 1) * n] for i in range(n)]


def parse_fixes(spec):
    return {} if spec == "-" else {int(a): int(b) for a, b in (x.split(":") for x in spec.split(","))}


def fmt(fixed):
    return ",".join(f"{i}:{j}" for i, j in sorted(fixed.items())) or "-"


def dnn(args, spec, iters, extra=()):
    out = subprocess.run([DNNTEST] + args.dnn_flags.split() + ["-H", "0", "-i", str(iters),
                          "-c", repr(args.cutoff), "-x", spec, *extra, args.instance],
                         capture_output=True, text=True, env=ENV, cwd=ROOT).stdout
    res = [l for l in out.splitlines() if l.startswith("RESULT")]
    if not res:
        raise RuntimeError(f"dnntest failed on {spec}:\n{out[-2000:]}")
    r = dict(kv.split("=", 1) for kv in res[0].split()[1:])
    xbar = {(int(a), int(b)): float(c) for _, a, b, c in
            (l.split() for l in out.splitlines() if l.startswith("XBAR"))}
    return float(r["lb"]), r["fathomed"] == "1", float(r["cpu"]), xbar


def console(args, specs, want_u=False):
    """QPB (parametric) at each node, and U (local order) if want_u."""
    cmds = []
    for s in specs:
        cmds.append("c")
        for i, j in sorted(parse_fixes(s).items()):
            cmds += ["f", str(i), str(j)]
        cmds.append("u" if want_u else "b")
    cmds.append("q")
    out = subprocess.run([BNBDRV, "-p", "PARAM/nug_deep.par", "-u", repr(args.ub), "-c", args.instance],
                         input="\n".join(cmds) + "\n", capture_output=True, text=True,
                         cwd=os.path.join(ROOT, "latest"), env=ENV).stdout
    if want_u:
        blocks = out.split("U = ")[1:]
        us = []
        for blk in blocks:
            rows = []
            for line in blk.splitlines()[1:]:
                vals = line.split()
                if not vals or not re.match(r"^-?[0-9.]+$", vals[0]):
                    break
                rows.append([float(x) for x in vals])
            us.append(rows)
        q = [float(x) for x in re.findall(r"QPBp?\s*=\s*(\S+)", out)]
        return us, q
    return [float(x) for x in re.findall(r"^(?:bnb> )?QPBp\s*=\s*(\S+)", out, re.M)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("instance")
    p.add_argument("--ub", type=float, required=True)
    p.add_argument("--delta", type=float, default=1.95, help="fathoming tolerance (cutoff = ub - delta)")
    p.add_argument("--parents", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--parent-iters", type=int, default=1500)
    p.add_argument("--child-iters", type=int, default=300)
    p.add_argument("--random", type=int, default=2)
    p.add_argument("--gain", type=float, default=130.0, help="bound gained per level (wil50: ~130)")
    p.add_argument("--branch", type=float, default=47.0, help="children per node one level down")
    p.add_argument("--seed", type=int, default=20261031)
    p.add_argument("--dnn-flags", default="-T 1 -d 4 -g 1.0 -G 1.0")
    p.add_argument("--dry-run", action="store_true", help="list the candidate lines and stop")
    args = p.parse_args()
    args.instance = os.path.abspath(args.instance)
    args.cutoff = args.ub - args.delta
    os.makedirs(args.work, exist_ok=True)
    n, A, B = read_qap(args.instance)
    rng = random.Random(args.seed)
    parents = [l.split() for l in open(args.parents) if l.strip() and not l.startswith("#")]

    report = []
    pooled = []   # (xbar, child qpb, gap)
    for pid, spec in parents:
        fixed = parse_fixes(spec)
        fr = [i for i in range(n) if i not in fixed]
        fc = [j for j in range(n) if j not in fixed.values()]
        # 1. parent
        pf = os.path.join(args.work, f"parent_{pid}.json")
        if os.path.exists(pf):
            par = json.load(open(pf))
        else:
            lb, fath, cpu, xbar = dnn(args, spec, args.parent_iters, ["-X"])
            us, q = console(args, [spec], want_u=True)
            par = {"lb": lb, "cpu": cpu, "xbar": [[i, j, x] for (i, j), x in xbar.items()], "U": us[0]}
            json.dump(par, open(pf, "w"))
        xbar = {(i, j): x for i, j, x in par["xbar"]}
        U = par["U"]
        # 2. candidate lines
        lines = {}
        ecc = max(fc, key=lambda j: (max(B[j][k] for k in fc if k != j), -j))
        lines.setdefault(("col", ecc), []).append("ecc")
        rs = sorted(range(len(fr)), key=lambda a: -sum(U[a]))[:2]
        cs = sorted(range(len(fc)), key=lambda b: -sum(U[a][b] for a in range(len(fr))))[:2]
        for a in rs: lines.setdefault(("row", fr[a]), []).append("U")
        for b in cs: lines.setdefault(("col", fc[b]), []).append("U")
        xr = max(fr, key=lambda i: sum(xbar[(i, j)] ** 2 for j in fc))
        xc = max(fc, key=lambda j: sum(xbar[(i, j)] ** 2 for i in fr))
        lines.setdefault(("row", xr), []).append("xbar")
        lines.setdefault(("col", xc), []).append("xbar")
        for i in rng.sample([i for i in fr if ("row", i) not in lines], args.random):
            lines.setdefault(("row", i), []).append("random")
        for j in rng.sample([j for j in fc if ("col", j) not in lines], args.random):
            lines.setdefault(("col", j), []).append("random")
        print(f"parent {pid} {spec}: lb {par['lb']:.1f} ({args.cutoff - par['lb']:.1f} short); "
              f"lines {[(k, x, '+'.join(v)) for (k, x), v in lines.items()]}", flush=True)
        if args.dry_run:
            continue
        # 3. children
        todo = []
        for (kind, x) in lines:
            for y in (fc if kind == "row" else fr):
                i, j = (x, y) if kind == "row" else (y, x)
                child = dict(fixed); child[i] = j
                todo.append((kind, x, i, j, fmt(child)))
        def run(t):
            kind, x, i, j, cspec = t
            f = os.path.join(args.work, f"child_{cspec.replace(':', '-').replace(',', '_')}.json")
            if os.path.exists(f):
                return t, json.load(open(f))
            lb, fath, cpu, _ = dnn(args, cspec, args.child_iters)
            r = {"lb": lb, "fathomed": fath, "cpu": cpu}
            json.dump(r, open(f, "w"))
            return t, r
        res = {}
        with cf.ThreadPoolExecutor(args.jobs) as ex:
            for k, (t, r) in enumerate(ex.map(run, todo), 1):
                res[t] = r
                if k % 50 == 0 or k == len(todo):
                    print(f"  {k}/{len(todo)} children", flush=True)
        # child QPB, batched through the console (cached)
        qf = os.path.join(args.work, f"qpb_{pid}.json")
        qmap = json.load(open(qf)) if os.path.exists(qf) else {}
        missing = [t[4] for t in todo if t[4] not in qmap]
        chunks = [missing[k::args.jobs] for k in range(args.jobs)]
        with cf.ThreadPoolExecutor(args.jobs) as ex:
            for ch, qs in zip(chunks, ex.map(lambda c: console(args, c) if c else [], chunks)):
                qmap.update(dict(zip(ch, qs)))
        json.dump(qmap, open(qf, "w"))
        # 4. scores
        rows = []
        for (kind, x), sel in lines.items():
            ch = [(t, res[t]) for t in todo if t[0] == kind and t[1] == x]
            gaps = [0.0 if r["fathomed"] else max(0.0, args.cutoff - r["lb"]) for t, r in ch]
            cost = sum(0.0 if r["fathomed"] else args.branch ** (g / args.gain) for (t, r), g in zip(ch, gaps))
            qsum = sum(qmap[t[4]] for t, r in ch)
            for (t, r), g in zip(ch, gaps):
                pooled.append((xbar[(t[2], t[3])], qmap[t[4]], g))
            rows.append({"line": f"{kind} {x}", "selected by": "+".join(sel), "children": len(ch),
                         "fathomed": sum(r["fathomed"] for t, r in ch),
                         "within a level": sum(1 for g in gaps if g < args.gain),
                         "median gap": sorted(gaps)[len(gaps) // 2], "cost": cost,
                         "mean child QPB": qsum / len(ch), "dnn cpu": sum(r["cpu"] for t, r in ch)})
        # rule 4's own pick among its U candidates: the largest mean child QPB
        ucands = [r for r in rows if "U" in r["selected by"]]
        if ucands:
            max(ucands, key=lambda r: r["mean child QPB"])["selected by"] += "+lookahead"
        rows.sort(key=lambda r: r["cost"])
        report.append((pid, spec, par, rows))
        print(f"\nparent {pid} {spec}, DNN lb {par['lb']:.1f} ({args.cutoff - par['lb']:.1f} short)")
        print(f"  {'line':9s} {'selected by':22s} {'fathomed':>8s} {'<1 level':>8s} {'median gap':>10s} "
              f"{'cost':>12s} {'mean QPB':>9s} {'cpu':>7s}")
        for r in rows:
            print(f"  {r['line']:9s} {r['selected by']:22s} {r['fathomed']:4d}/{r['children']:<3d} "
                  f"{r['within a level']:8d} {r['median gap']:10.1f} {r['cost']:12.4g} "
                  f"{r['mean child QPB']:9.1f} {r['dnn cpu']:7.0f}", flush=True)
        json.dump(rows, open(os.path.join(args.work, f"report_{pid}.json"), "w"), indent=1)
    if pooled and not args.dry_run:
        def spearman(a, b):
            def rk(v):
                o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
                for k, i in enumerate(o): r[i] = k
                return r
            ra, rb = rk(a), rk(b); m = len(a); ma, mb = sum(ra) / m, sum(rb) / m
            c = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
            return c / (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
        xs, qs, gs = zip(*pooled)
        print(f"\npooled over {len(pooled)} children: Spearman(x_ij, gap) {spearman(xs, gs):+.3f}, "
              f"Spearman(child QPB, gap) {spearman(qs, gs):+.3f} (negative: a larger x or QPB, a smaller gap)")
        print("per parent: cost of each selector's line relative to the best line and to the random lines' mean")
        for pid, spec, par, rows in report:
            best = rows[0]["cost"]
            rnd = [r["cost"] for r in rows if r["selected by"] == "random"]
            parts = []
            ratio = lambda a, b: f"{a / b:.2f}x" if b > 0 else ("1.00x" if a == 0 else "inf")
            for selname in ("ecc", "lookahead", "xbar"):
                cs = [r["cost"] for r in rows if selname in r["selected by"].split("+")]
                if cs:
                    parts.append(f"{selname} {ratio(min(cs), best)} best" +
                                 (f", {ratio(min(cs), sum(rnd) / len(rnd))} random" if rnd else ""))
            print(f"  parent {pid}: " + "; ".join(parts))


if __name__ == "__main__":
    main()
