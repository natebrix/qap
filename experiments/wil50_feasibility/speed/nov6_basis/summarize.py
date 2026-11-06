#!/usr/bin/env python3
"""Warm set (run_warmset.sh) with and without DNN_SHIP_BASIS: per child and in total.

  summarize.py [--dir warm] BASE NEW [NEW2 ...] [--table]
  summarize.py [--dir warm] --profile [--short] CONFIG [CONFIG ...]

For each configuration against BASE, over the families both finished: the
parents (must be identical: a cold call is never seeded); per child the
iterations, fathom, certified running best at iterations 50, 100 and 200 (the
"it" lines; the last value when the solve stopped earlier), CPU, the first
R-step's CPU (the "basis profile" line, -v -v) and the seed's outcome (RESULT
seed=); then the totals: children's CPU and its saving, the start-up saving
(sum over children of the first R-step's difference), the fathom sets and the
kill line (less than 3% CPU saved on the warm children, or a child fathoming
later or a fathom lost, kills it).

--profile: where the children's CPU goes, per configuration (the -v -v profile
lines summed over the children of the families every listed configuration
finished): the first R-step, the full refreshes, the filter steps (count, ms
per step, mean block size kb), the rest of the eigen part (the seed, the power
steps, the fp32 copy of W), forward + expand + Y/Z-step, the certificates and
the unprofiled remainder; then, over the children that reach iteration 20, the
mean CPU to iteration 10 without the first R-step and of iterations 11-20 (the
"it" lines' cumulative CPU). --short: the children that fathomed only."""
import glob
import os
import re
import sys

D = os.path.dirname(os.path.abspath(__file__))


def parse(f):
    """-> list of 5 solves: parent then 4 children"""
    solves, cur = [], None
    for line in open(f, errors="replace"):
        if line.startswith("node ") or line.startswith("CHILD-PASS"):
            cur = {"trace": {}, "first": None, "seedres": None, "eig": None, "refresh": None}
        if cur is None:
            continue
        m = re.match(r"it\s+(\d+) lb (\S+) \(cert", line)
        if m and "extra" not in line:
            cur["trace"][int(m.group(1))] = float(m.group(2))
        m = re.match(r"basis profile: first R-step ([0-9.]+) ms .*?res (\S+),", line)
        if m:
            cur["first"] = float(m.group(1)); cur["seedres"] = float(m.group(2))
        m = re.match(r"it\s+(\d+) lb \S+ \(cert.* cpu (\S+)", line)
        if m and "extra" not in line:
            cur.setdefault("cum", {})[int(m.group(1))] = float(m.group(2))
        m = re.match(r"profile \(s\): forward (\S+)\s+eig (\S+)\s+expand (\S+)\s+Y/Z-step (\S+)\s+certify (\S+)", line)
        if m:
            cur["eig"] = float(m.group(2))
            cur["prof"] = [float(x) for x in m.groups()]
        m = re.match(r"eig split \(s\): full ssyevd (\S+) \((\d+) refreshes.*filter\+RR (\S+) \((\d+) steps, mean kb (\S+),"
                     r".*form R (\S+)", line)
        if m:
            cur["refresh"] = (float(m.group(1)), int(m.group(2)))
            cur["filt"] = (float(m.group(3)), int(m.group(4)), float(m.group(5)), float(m.group(6)))
        if line.startswith("RESULT"):
            kv = dict(x.split("=", 1) for x in line.split()[1:] if "=" in x)
            cur.update(node=int(kv["node"]), iters=int(kv["iters"]), lb=float(kv["lb"]), cpu=float(kv["cpu"]),
                       status=kv["status"], fathomed=kv["status"] == "fathomed", seed=int(kv.get("seed", 0)))
            solves.append(cur)
            cur = None
    return solves if len(solves) == 5 else None


def lb_at(s, it):
    """the running best certified value at iteration it (the last check at or before it)"""
    ks = [k for k in s["trace"] if k <= it]
    return s["trace"][max(ks)] if ks else float("nan")


def profile(runs, configs, short):
    order = ["f9", "f6", "p10", "p13", "p2", "p7"]
    fams = [f for f in order if all(f in runs.get(c, {}) for c in configs)]
    print(f"== profile of the children{' that fathomed' if short else ''}, families {','.join(fams)} (s unless stated)")
    print("%-7s %4s %6s %8s %7s %11s %18s %6s %8s %8s %7s %6s %18s" % (
        "config", "kids", "iters", "CPU", "1st R", "refreshes", "filter steps (ms)", "kb", "eig rest", "fwd+exp+YZ",
        "certify", "other", "to it10 / it11-20"))
    for c in configs:
        kids = [k for f in fams for k in runs[c][f][1:] if not short or k["fathomed"]]
        if not kids or any("prof" not in k or "filt" not in k for k in kids):
            print(f"{c}: no profile lines (-v -v)"); continue
        cpu = sum(k["cpu"] for k in kids); it = sum(k["iters"] for k in kids)
        first = sum(k["first"] or 0 for k in kids) / 1e3
        rs = sum(k["refresh"][0] for k in kids); rn = sum(k["refresh"][1] for k in kids)
        fs = sum(k["filt"][0] for k in kids); fn = sum(k["filt"][1] for k in kids)
        kb = sum(k["filt"][1] * k["filt"][2] for k in kids) / max(1, fn); fr = sum(k["filt"][3] for k in kids)
        eig = sum(k["prof"][1] for k in kids)
        fey = sum(k["prof"][0] + k["prof"][2] + k["prof"][3] for k in kids); cert = sum(k["prof"][4] for k in kids)
        other = cpu - sum(sum(k["prof"]) for k in kids)
        e20 = [k for k in kids if 10 in k.get("cum", {}) and 20 in k.get("cum", {})]
        a = sum(k["cum"][10] - (k["first"] or 0) / 1e3 for k in e20) / max(1, len(e20))
        b = sum(k["cum"][20] - k["cum"][10] for k in e20) / max(1, len(e20))
        print("%-7s %4d %6d %8.2f %7.2f %4d %6.2f %5d %6.2f %5.1f %6.0f %8.2f %8.2f %7.2f %6.2f %8.3f / %.3f (%d)" % (
            c, len(kids), it, cpu, first, rn, rs, fn, fs, 1e3 * fs / max(1, fn), kb, eig - rs - fs - fr, fey, cert,
            other, a, b, len(e20)))
    print("(eig rest: the eigen part without the refreshes, the filter steps and form R: the seed, the power steps,"
          " W's fp32 copy)")


def main():
    args = sys.argv[1:]
    table = "--table" in args
    if table:
        args.remove("--table")
    prof = "--profile" in args
    if prof:
        args.remove("--profile")
    short = "--short" in args
    if short:
        args.remove("--short")
    d = os.path.join(D, "warm")
    if "--dir" in args:
        k = args.index("--dir"); d = args[k + 1]; del args[k:k + 2]
    runs = {}
    for f in glob.glob(os.path.join(d, "*.log")):
        name, fam = re.match(r"(.+)_([a-z0-9]+)\.log$", os.path.basename(f)).groups()
        s = parse(f)
        if s:
            runs.setdefault(name, {})[fam] = s
    if prof:
        profile(runs, args, short)
        return
    base = args[0]
    order = ["f9", "f6", "p10", "p13", "p2", "p7"]
    for new in args[1:]:
        fams = [f for f in order if f in runs.get(base, {}) and f in runs.get(new, {})]
        print(f"== {new} against {base}, families {','.join(fams)}")
        B, N = runs[base], runs[new]
        same_par = all(B[f][0]["iters"] == N[f][0]["iters"] and B[f][0]["lb"] == N[f][0]["lb"] for f in fams)
        print(f"parents identical (iterations and lb): {'yes' if same_par else 'NO'}; "
              f"parents' CPU {sum(B[f][0]['cpu'] for f in fams):.1f} / {sum(N[f][0]['cpu'] for f in fams):.1f} s")
        rows = []
        for f in fams:
            for k in range(1, 5):
                rows.append((f, k, B[f][k], N[f][k]))
        if table:
            print("%-4s %s  %-11s %-9s %-23s %-23s %-23s %-15s %-15s %s" % (
                "fam", "k", "iters b/n", "fathom", "lb@50 b/n", "lb@100 b/n", "lb@200 b/n", "CPU b/n s",
                "1st R-step ms", "seed (res)"))
            for f, k, b, n in rows:
                print("%-4s %d  %4d/%-6d %-9s %-23s %-23s %-23s %6.2f/%-8.2f %6.0f/%-8.0f %d (%s)" % (
                    f, k, b["iters"], n["iters"], ("F" if b["fathomed"] else "m") + "/" + ("F" if n["fathomed"] else "m"),
                    "%.2f/%.2f" % (lb_at(b, 50), lb_at(n, 50)), "%.2f/%.2f" % (lb_at(b, 100), lb_at(n, 100)),
                    "%.2f/%.2f" % (lb_at(b, 200), lb_at(n, 200)), b["cpu"], n["cpu"],
                    b["first"] or float("nan"), n["first"] or float("nan"), n["seed"],
                    "%.3g" % n["seedres"] if n["seedres"] is not None else "-"))
        cb = sum(b["cpu"] for *_, b, n in rows); cn = sum(n["cpu"] for *_, b, n in rows)
        ib = sum(b["iters"] for *_, b, n in rows); inn = sum(n["iters"] for *_, b, n in rows)
        fb = {(f, k) for f, k, b, n in rows if b["fathomed"]}; fn = {(f, k) for f, k, b, n in rows if n["fathomed"]}
        later = [(f, k, b["iters"], n["iters"]) for f, k, b, n in rows if b["fathomed"] and n["fathomed"] and n["iters"] > b["iters"]]
        earlier = [(f, k, b["iters"], n["iters"]) for f, k, b, n in rows if b["fathomed"] and n["fathomed"] and n["iters"] < b["iters"]]
        fr = [(b["first"], n["first"]) for *_, b, n in rows if b["first"] is not None and n["first"] is not None]
        seeded = sum(1 for *_, b, n in rows if n["seed"] == 1); rejected = sum(1 for *_, b, n in rows if n["seed"] == -1)
        fcb = sum(b["cpu"] for f, k, b, n in rows if (f, k) in fb & fn); fcn = sum(n["cpu"] for f, k, b, n in rows if (f, k) in fb & fn)
        print(f"children: {len(rows)}, seeded {seeded}, rejected {rejected}")
        print(f"iterations {ib} / {inn}; fathomed {len(fb)} / {len(fn)}; lost {sorted(fb - fn)}, gained {sorted(fn - fb)}; "
              f"fathom later {later}, earlier {earlier}")
        if fr:
            sb = sum(x for x, _ in fr); sn = sum(y for _, y in fr)
            print(f"first R-step: mean {sb / len(fr):.0f} / {sn / len(fr):.0f} ms, saved {(sb - sn) / 1e3:.2f} s in all "
                  f"({(sb - sn) / len(fr):.0f} ms per child)")
        print(f"children's CPU {cb:.2f} / {cn:.2f} s: saved {cb - cn:+.2f} s, {100 * (cb - cn) / cb:+.2f}% "
              f"(per iteration {1e3 * cb / ib:.1f} / {1e3 * cn / inn:.1f} ms); "
              f"on the {len(fb & fn)} common fathoms {fcb:.2f} / {fcn:.2f} s, {100 * (fcb - fcn) / fcb:+.2f}%")
        if fr and len(fr) == len(rows):
            # the CPU outside the first R-step is the same work in both configurations
            # (same iterations give the same R-steps, refreshes and certificates), so
            # its difference per iteration is the noise of this comparison; the first
            # R-step's difference (timed in the process) is the switch's saving
            sb = sum(x for x, _ in fr) / 1e3; sn = sum(y for _, y in fr) / 1e3
            rb, rn = cb - sb, cn - sn
            print(f"split: first R-steps {sb:.2f} / {sn:.2f} s (saved {100 * (sb - sn) / cb:.2f}% of the children's CPU); "
                  f"the rest {rb:.2f} / {rn:.2f} s, per iteration {1e3 * rb / ib:.1f} / {1e3 * rn / inn:.1f} ms "
                  f"({100 * (rb / ib - rn / inn) / (rb / ib):+.2f}%: noise and trajectory)")
        misses = [(f, k, b, n) for f, k, b, n in rows if not b["fathomed"]]
        for at in (50, 100, 200):
            dd = [lb_at(n, at) - lb_at(b, at) for f, k, b, n in rows]
            dd = [x for x in dd if x == x]
            dm = [lb_at(n, at) - lb_at(b, at) for f, k, b, n in misses]
            dm = [x for x in dm if x == x]
            if dd:
                print(f"lb@{at} new - base: all children median {sorted(dd)[len(dd) // 2]:+.3f} "
                      f"[{min(dd):+.3f}, {max(dd):+.3f}]" + (f"; base's misses [{min(dm):+.3f}, {max(dm):+.3f}]" if dm else ""))
        kill = (cb - cn) / cb < 0.03 or bool(fb - fn) or bool(later)
        print(f"kill line (< 3% CPU saved, a fathom lost or later): {'KILLED' if kill else 'passes'}")
        print()


if __name__ == "__main__":
    main()
