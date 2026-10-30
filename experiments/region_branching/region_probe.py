#!/usr/bin/env python3
"""Region-branching probe: is the DNN bound of "facility f in block R" close to
the best of the single-location children "f at location l", l in R?

  python3 experiments/region_branching/region_probe.py qaplib/nug30.dat \\
      --grid 5x6 --block 2x3 --fixes 0:13,6:14 --ub 6124 [--facility F] [--jobs 4]

The instance's distance matrix must be a rows x cols grid in row-major order
(location l at row l // cols, column l % cols); the script checks this.
For each block of the free locations (blocks of block_rows x block_cols grid
cells), it bounds the region node (the node, with f forbidden at every free
location outside the block: dnntest's elimination field) and every cell child
(the node plus f -> l), with the solver's DNN settings (-H 300). A block whose
region bound is close to the minimum over its cells can be fathomed with one
DNN solve instead of one per cell.
"""
import argparse, os, subprocess, tempfile, concurrent.futures as cf
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DNNTEST = os.path.join(ROOT, "latest", "dnntest")


def read(path):
    t = open(path).read().split()
    n = int(t[0]); v = list(map(float, t[1:1 + 2 * n * n]))
    A = [v[i * n:(i + 1) * n] for i in range(n)]
    B = [v[n * n + i * n:n * n + (i + 1) * n] for i in range(n)]
    return n, A, B


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("instance"); p.add_argument("--grid", required=True); p.add_argument("--block", required=True)
    p.add_argument("--fixes", default="-"); p.add_argument("--ub", type=float, required=True)
    p.add_argument("--facility", type=int); p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--iters", type=int, default=1500); p.add_argument("--hopeless", type=int, default=300)
    a = p.parse_args()
    n, A, B = read(a.instance)
    gr, gc = map(int, a.grid.split("x")); br, bc = map(int, a.block.split("x"))
    assert gr * gc == n
    # B may carry a unit cost per axis (e.g. kra); require only the grid's equality pattern
    if not all((B[i][j] == 0) == (i == j) for i in range(n) for j in range(n)):
        raise SystemExit("distance matrix does not look like a grid")
    fixed = {} if a.fixes == "-" else {int(x.split(":")[0]): int(x.split(":")[1]) for x in a.fixes.split(",")}
    free_f = [i for i in range(n) if i not in fixed]
    free_l = [l for l in range(n) if l not in fixed.values()]
    f = a.facility if a.facility is not None else max(free_f, key=lambda i: sum(A[i]))
    blocks = {}
    for l in free_l:
        blocks.setdefault(((l // gc) // br, (l % gc) // bc), []).append(l)
    jobs = []
    for key, cells in sorted(blocks.items()):
        el = ",".join(f"{f}:{l}" for l in free_l if l not in cells)
        jobs.append((key, "region", a.fixes, el if el else "-"))
        for l in cells:
            fx = (a.fixes + "," if a.fixes != "-" else "") + f"{f}:{l}"
            jobs.append((key, l, fx, "-"))
    def run(job):
        key, what, fx, el = job
        with tempfile.NamedTemporaryFile("w", suffix=".nodes", delete=False) as fh:
            fh.write(f"{fx} {a.ub - 1.95!r} {el}\n"); nf = fh.name
        out = subprocess.run([DNNTEST, "-f", "-T", "1", "-H", str(a.hopeless), "-i", str(a.iters), "-N", nf, a.instance],
                             capture_output=True, text=True).stdout
        os.unlink(nf)
        r = [dict(kv.split("=", 1) for kv in l.split()[1:]) for l in out.splitlines() if l.startswith("RESULT")][0]
        return job, float(r["lb"]), r["status"], float(r["cpu"])
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(run, jobs))
    cutoff = a.ub - 1.95
    print(f"{a.instance} node {a.fixes}, facility {f}, blocks {a.block} of the {a.grid} grid, cutoff {cutoff}")
    print(f"{'block':8s} {'cells':>5s} {'region lb':>10s} {'status':>9s} {'cpu':>6s} | {'min cell lb':>11s} {'cells fathomed':>15s} {'cells cpu':>9s}")
    tot_r = tot_c = 0
    for key in sorted(blocks):
        reg = [x for x in res if x[0][0] == key and x[0][1] == "region"][0]
        cells = [x for x in res if x[0][0] == key and x[0][1] != "region"]
        tot_r += reg[3]; tot_c += sum(c[3] for c in cells)
        print(f"{str(key):8s} {len(cells):5d} {reg[1]:10.1f} {reg[2]:>9s} {reg[3]:6.0f} | {min(c[1] for c in cells):11.1f} "
              f"{sum(1 for c in cells if c[1] > cutoff):8d}/{len(cells):<6d} {sum(c[3] for c in cells):9.0f}")
    print(f"CPU: regions {tot_r:.0f} s, cells {tot_c:.0f} s")


if __name__ == "__main__":
    main()
