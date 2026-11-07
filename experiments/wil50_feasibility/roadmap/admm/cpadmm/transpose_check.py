#!/usr/bin/env python3
"""Transpose equivalence of the CP solve (README.md, "Side 1"): the same node run
on wil50.dat with the regions on side 0 and on wil50t.dat with them on side 1 and
the fixes transposed (dnncuts -cpsep, run_sep.sh's settings). The DNN bound is
invariant under the transposition, so the certified values must agree at every
check up to the floating-point path's drift. Prints:

- the give-up points (the uncut round 0, the plain solve's own drift);
- the certified value at matching checks: max |difference| over windows of t;
- the running best: max |difference| and the first t at each level;
- the separations and the pool;
- the CPU per iteration of the CP solve and of its P-steps;
- with --ref LOG: whether the side-0 log's lines equal the reference log's up to
  its last check (the side-0 path unchanged), CPU fields stripped.

  python3 transpose_check.py SIDE0.log SIDE1.log [--ref REF.log] [--levels 48713 48743 48753]
"""
import argparse
import re


def parse(path):
    txt = open(path).read()
    d = {"cert": {}, "lines": []}
    m = re.search(r"^round  0  lb (\S+) .*? violated (\d+)", txt, re.M)
    d["r0"], d["r0_viol"] = (float(m.group(1)), int(m.group(2))) if m else (None, None)
    for mm in re.finditer(r"^cp it +(\d+) cert (\S+) .*?cuts (\d+) active (\d+).*?ms pstep (\S+) \((\d+)\) hpass (\S+)", txt, re.M):
        t = int(mm.group(1))
        d["cert"][t] = float(mm.group(2))
        d.setdefault("cuts", {})[t] = int(mm.group(3))
        d.setdefault("active", {})[t] = int(mm.group(4))
        d["pstep_ms"], d["hpass_ms"] = float(mm.group(5)), float(mm.group(7))
    m = re.search(r"^cp solve  lb (\S+) .*?iters (\d+)  cpu (\S+)", txt, re.M)
    if m:
        d["lb"], d["iters"], d["cpu"] = float(m.group(1)), int(m.group(2)), float(m.group(3))
    m = re.search(r"^profile \(s\): (.*)$", txt, re.M)
    d["profile"] = m.group(1) if m else ""
    m = re.search(r"^profile \(s\): cp projection (\S+)", txt, re.M)
    d["proj"] = float(m.group(1)) if m else None
    m = re.search(r"^cp separation: (\d+) separations, (\d+) cuts added, (\d+) dropped, pool (\d+)", txt, re.M)
    d["sep"] = tuple(int(x) for x in m.groups()) if m else None
    # the trace lines that do not depend on timing: "it" and "cp it" with the CPU fields cut
    for line in txt.splitlines():
        if line.startswith("it ") or line.startswith("cp it "):
            line = re.sub(r" cpu \S+", "", line)
            line = re.sub(r" cpcpu \S+", "", line)
            line = re.sub(r" \| ms pstep .*?diag \S+ s", "", line)
            line = re.sub(r" sep \S+ s", "", line)
            d["lines"].append(line)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("side0"); ap.add_argument("side1")
    ap.add_argument("--ref")
    ap.add_argument("--levels", type=float, nargs="+", default=[48686.4, 48713, 48743, 48753])
    a = ap.parse_args()
    s0, s1 = parse(a.side0), parse(a.side1)
    print(f"give-up point (uncut, 400 iterations): side 0 {s0['r0']:.4f} ({s0['r0_viol']} violated), "
          f"side 1 {s1['r0']:.4f} ({s1['r0_viol']} violated), difference {s1['r0'] - s0['r0']:+.4f}")
    ts = sorted(set(s0["cert"]) & set(s1["cert"]))
    print(f"matching checks: {len(ts)} (t = {ts[0]}..{ts[-1]})")
    print("| t window | max |cert1 - cert0| | at t | mean |diff| | max |best1 - best0| |")
    print("|---|---|---|---|---|")
    b0 = b1 = -1e300
    best = {}
    for t in ts:
        b0 = max(b0, s0["cert"][t]); b1 = max(b1, s1["cert"][t])
        best[t] = (b0, b1)
    for lo, hi in ((1, 100), (101, 300), (301, 500), (501, 1000), (1001, 1500), (1501, 10 ** 9), (1, 10 ** 9)):
        w = [t for t in ts if lo <= t <= hi]
        if not w:
            continue
        diffs = [abs(s1["cert"][t] - s0["cert"][t]) for t in w]
        k = max(range(len(w)), key=lambda i: diffs[i])
        bd = max(abs(best[t][1] - best[t][0]) for t in w)
        hi_s = "end" if hi >= 10 ** 9 else str(hi)
        print(f"| {lo}-{hi_s} | {diffs[k]:.4f} | {w[k]} | {sum(diffs) / len(diffs):.4f} | {bd:.4f} |")
    print("first t at a level (running best):")
    for lv in a.levels:
        f0 = next((t for t in ts if best[t][0] >= lv), None)
        f1 = next((t for t in ts if best[t][1] >= lv), None)
        print(f"  {lv:,.1f}: side 0 t = {f0}, side 1 t = {f1}")
    for t in (100, 300, 500, 1000, 1070, 1500, 2000, 3000):
        if t in best:
            print(f"  best at t = {t}: side 0 {best[t][0]:.4f}, side 1 {best[t][1]:.4f}, difference {best[t][1] - best[t][0]:+.4f}"
                  f" (pool {s0['cuts'][t]} / {s1['cuts'][t]}, active {s0['active'][t]} / {s1['active'][t]})")
    print(f"separations (count, added, dropped, final pool): side 0 {s0['sep']}, side 1 {s1['sep']}")
    for nm, s in (("side 0", s0), ("side 1", s1)):
        if "cpu" in s:
            print(f"{nm}: CP solve {s['iters']} iterations, {s['cpu']:.1f} s, {1000 * s['cpu'] / s['iters']:.1f} ms per iteration;"
                  f" P-step {s['pstep_ms']:.2f} ms per call, half-space pass {s['hpass_ms']:.2f} ms;"
                  f" projection {s['proj']:.1f} s ({1000 * s['proj'] / s['iters']:.1f} ms per iteration)")
            print(f"        profile: {s['profile']}")
    if "cpu" in s0 and "cpu" in s1:
        d = 1000 * (s1["cpu"] / s1["iters"] - s0["cpu"] / s0["iters"])
        print(f"side 1 - side 0 per iteration: {d:+.1f} ms ({100 * d / (1000 * s0['cpu'] / s0['iters']):+.1f}%); "
              f"projection {1000 * (s1['proj'] / s1['iters'] - s0['proj'] / s0['iters']):+.1f} ms")
    if a.ref:
        r = parse(a.ref)
        n = min(len(r["lines"]), len(s0["lines"]))
        same = sum(x == y for x, y in zip(r["lines"][:n], s0["lines"][:n]))
        first = next((i for i in range(n) if r["lines"][i] != s0["lines"][i]), None)
        print(f"side 0 against {a.ref.split('/')[-1]}: {same} of the first {n} trace lines identical (CPU stripped)"
              + ("" if first is None else f"; first difference at line {first}:\n  ref: {r['lines'][first][:160]}\n  new: {s0['lines'][first][:160]}"))


if __name__ == "__main__":
    main()
