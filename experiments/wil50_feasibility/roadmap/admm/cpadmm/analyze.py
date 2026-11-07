#!/usr/bin/env python3
"""Tables for the CP-ADMM study (README.md) from the logs in logs/.

  python3 -I analyze.py [--logs logs] [--ref polyak_ref_13.log] [--node 13] [--levels L ...] [--] [ARM ...]

Iterations are counted from the give-up point (after round 0's 400 uncut
iterations), for both protocols. Polyak: the best certified value so far over
its rounds and quiet round (every check of every solve is certified), and the
value at the check itself. CP-ADMM: the running best of its single solve.
"""
import argparse
import os
import re

D = os.path.dirname(os.path.abspath(__file__))
IT = re.compile(r"^it +(\d+) lb (\S+) \(cert (\S+),.*cpu ([0-9.]+)")
CPIT = re.compile(r"^cp it +(\d+) cert (\S+) hsum (\S+) cpcpu (\S+) cuts (\d+) active (\d+) sum_nu (\S+) viol (\d+) "
                  r"\(max (\S+)/\|\|G\|\|, rnorm (\S+)\) dstep (\S+)(?: \| ms pstep (\S+) \((\d+)\) hpass (\S+) \((\d+)\) diag (\S+) s)?(?: \| pool (\d+) sep (\d+) \(last: violated (\d+) max (\S+) added (\d+); total added (\d+) dropped (\d+)\) sep (\S+) s)?(?:  diag@(\d+) \((\d+) sweeps\): \|Y-Ylong\| (\S+) "
                  r"\(\|Ylong\| (\S+)\) \|G'\(nu-nulong\)\| (\S+) long maxviol (\S+))?")


def fl(x):
    return float("-inf") if x in ("-inf", "-nan", "nan") else float(x)


def polyak(path):
    """[(t since give-up, cert, best so far, call index)], per-call iteration counts"""
    out, calls = [], []
    off, last, call = 0, 0, 0
    cpu_off, last_cpu = 0.0, 0.0
    best = float("-inf")
    rounds = []   # (name, iters, cpu) of every solve, from the round lines
    for line in open(path):
        rm = re.match(r"^(round +\d+|quiet round) +lb .* iters (\d+) +cpu ([0-9.]+)", line)
        if rm:
            rounds.append((rm.group(1), int(rm.group(2)), float(rm.group(3))))
        m = IT.match(line)
        if not m:
            continue
        it = int(m.group(1))
        if it <= last:   # a new solve
            calls.append(last); off += last; call += 1
            cpu_off += last_cpu
        last = it
        last_cpu = float(m.group(4))
        # a DNN_CORE_LAZYCERT check that skipped Rump's test prints an estimate: not certified
        c = float("-inf") if ", lazy" in line else fl(m.group(3))
        best = max(best, c)
        out.append((off + it - 400, c, best, call, cpu_off + float(m.group(4))))
    calls.append(last)
    polyak.rounds = rounds
    return out, calls


def cp(path):
    """the CP solve's checks: [(t, cert, running best, cpu)], cp lines, profile, round-0 lb"""
    its, cps, prof, r0, solve = [], [], {}, None, None
    seen_round0 = False
    quiet, quiet_started, q_off, qline = [], False, 0, None
    for line in open(path):
        if line.startswith("round  0"):
            seen_round0 = True
            r0 = float(line.split()[3])
            continue
        m = IT.match(line)
        if m and seen_round0:
            c = float("-inf") if ", lazy" in line else fl(m.group(3))
            it = int(m.group(1))
            if its and it <= its[-1][0] and not quiet_started:
                quiet_started = True   # the -cpquiet round after the CP solve
                q_off = its[-1][0]
            if quiet_started:
                quiet.append((q_off + it, c, float(m.group(4))))
            else:
                its.append((it, c, None, float(m.group(4))))
            continue
        if line.startswith("cp quiet round"):
            qline = line.strip()
            continue
        m = CPIT.match(line)
        if m:
            g = m.groups()
            cps.append(dict(it=int(g[0]), cert=fl(g[1]), hsum=float(g[2]), cpcpu=float(g[3]), active=int(g[5]),
                            sum_nu=float(g[6]), nviol=int(g[7]), maxviol=float(g[8]), rnorm=float(g[9]),
                            dstep=float(g[10]),
                            ms_p=float(g[11]) if g[11] else None, n_p=int(g[12]) if g[12] else None,
                            ms_h=float(g[13]) if g[13] else None, n_h=int(g[14]) if g[14] else None,
                            t_diag=float(g[15]) if g[15] else None,
                            pool=int(g[16]) if g[16] else None, nsep=int(g[17]) if g[17] else None,
                            sep_viol=int(g[18]) if g[18] else None, sep_max=float(g[19]) if g[19] else None,
                            sep_added=int(g[21]) if g[21] else None, sep_dropped=int(g[22]) if g[22] else None,
                            t_sep=float(g[23]) if g[23] else None,
                            diag_it=int(g[24]) if g[24] else None,
                            diag_dy=float(g[26]) if g[26] else None, diag_y=float(g[27]) if g[27] else None,
                            diag_dnu=float(g[28]) if g[28] else None, diag_viol=float(g[29]) if g[29] else None))
            continue
        if line.startswith("profile (s): forward"):
            v = re.findall(r"([a-zA-Z/-]+(?:-step)?) ([0-9.]+)", line)
            for k, x in v: prof[k] = float(x)
        elif line.startswith("profile (s): cp projection"):
            prof["cp"] = float(line.split()[4])
        elif line.startswith("cp solve"):
            solve = line.strip()
    best, out = float("-inf"), []
    for it, c, _, cpu in its:
        best = max(best, c)
        out.append((it, c, best, cpu))
    if len(out) > 20:
        # CPU per iteration of the CP solve from the cumulative CPU of its checks, iterations 100..end
        i0 = next(k for k, r in enumerate(out) if r[0] >= 100)
        prof["s_per_it"] = (out[-1][3] - out[i0][3]) / (out[-1][0] - out[i0][0])
        prof["s_per_it_range"] = "%d-%d" % (out[i0][0], out[-1][0])
    qbest, qout = best, []
    for t, c, cpu in quiet:
        qbest = max(qbest, c)
        qout.append((t, c, qbest, cpu))
    if qline:
        solve = (solve or "") + "\n  " + qline
    return out + qout, cps, prof, r0, solve


def at(curve, t, key=2):
    """the value at the last check <= t"""
    v = None
    for row in curve:
        if row[0] <= t:
            v = row[key]
    return v


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--logs", default=os.path.join(D, "logs"))
    p.add_argument("--ref", default="polyak_ref_13.log")
    p.add_argument("--target", type=float, default=48763.0)
    p.add_argument("--levels", type=float, nargs="+", default=[48713.0, 48743.0, 48753.0, 48755.39, 48763.0])
    p.add_argument("--node", default="13", help="the node id in the log names (13, d2)")
    p.add_argument("arms", nargs="*", default=["s1p1", "s1p5", "s3p3", "plain"])
    a = p.parse_args()
    ref = os.path.join(a.logs, a.ref)
    pol, calls = polyak(ref) if os.path.exists(ref) else ([], [])
    arms = []
    for name in a.arms:
        f = os.path.join(a.logs, "cp_%s_%s.log" % (name, a.node))
        if not os.path.exists(f):
            f += ".tmp"
        if os.path.exists(f):
            arms.append((name, f) + cp(f))
    ts = [10, 50, 100, 200, 300, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3600, 4000, 5100]
    hdr = ["iterations after the give-up point", "Polyak best so far", "Polyak at the check"]
    hdr += ["CP %s" % n for n, *_ in arms]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for t in ts:
        row = [str(t)]
        pb = at(pol, t); pc = at(pol, t, 1)
        row.append("%.1f" % pb if pb is not None else "-")
        row.append("%.1f" % pc if pc is not None else "-")
        any_cp = False
        for name, f, its, cps, prof, r0, solve in arms:
            curve = [(it, c, b) for it, c, b, _ in its]
            v = at(curve, t) if curve and curve[-1][0] >= t else None
            any_cp = any_cp or v is not None
            row.append("%.1f" % v if v is not None else "-")
        if pb is None and not any_cp:
            continue
        print("| " + " | ".join(row) + " |")
    if pol:
        print("\nPolyak calls (iterations each):", calls, " final best %.2f" % pol[-1][2])
        rr = [r for r in polyak.rounds if not r[0].startswith("round  0") and r[0] != "round 0"]
        cut = [r for r in rr if r[0].startswith("round")]
        if cut:
            print("  Polyak cut rounds: %d iterations, %.1f CPU s, %.4f s per iteration (setup included)"
                  % (sum(r[1] for r in cut), sum(r[2] for r in cut), sum(r[2] for r in cut) / sum(r[1] for r in cut)))
        q = [r for r in rr if r[0].startswith("quiet")]
        if q:
            print("  Polyak quiet round: %d iterations, %.1f CPU s, %.4f s per iteration" % (q[0][1], q[0][2], q[0][2] / q[0][1]))
        cpu0 = next((r[4] for r in pol if r[0] >= 0), 0.0)
        for lvl in a.levels:
            first = next((r for r in pol if r[2] >= lvl and r[0] > 0), None)
            print("  Polyak first reaches %.1f at t = %s (CPU since the give-up point %s s)"
                  % (lvl, first[0] if first else None, "%.0f" % (first[4] - cpu0) if first else "-"))
    for name, f, its, cps, prof, r0, solve in arms:
        print("\n== %s (%s)" % (name, os.path.basename(f)))
        print("  round 0 lb %s; %s" % (r0, solve))
        if solve:
            sm = re.search(r"iters (\d+)  cpu ([0-9.]+)", solve)
            if sm:
                print("  solve: %.4f s per iteration (setup included; quiet round excluded)" % (float(sm.group(2)) / int(sm.group(1))))
        for lvl in a.levels:
            first = next((r for r in its if r[2] >= lvl), None)
            print("  first reaches %.1f at t = %s (CPU %s s)" % (lvl, first[0] if first else None,
                                                             "%.0f" % first[3] if first else "-"))
        if prof:
            print("  profile:", prof)
        if cps and cps[-1]["ms_p"] is not None:
            c = cps[-1]
            print("  projection at the end: %.2f ms per P-step (%d), %.2f ms per half-space pass (%d), diagnostics %.1f s, hook total %.1f s"
                  % (c["ms_p"], c["n_p"], c["ms_h"], c["n_h"], c["t_diag"], c["cpcpu"]))
            nit = c["it"]
            print("  per iteration without the diagnostics: %.2f ms" % (1e3 * (c["cpcpu"] - c["t_diag"]) / nit))
            if c["pool"] is not None:
                print("  separation: %d separations, pool %d at the end, %d added, %d dropped, %.1f s (%.2f ms per iteration)"
                      % (c["nsep"], c["pool"], c["sep_added"], c["sep_dropped"], c["t_sep"], 1e3 * c["t_sep"] / nit))
                rows = [x for x in cps if x["pool"] is not None and x["it"] in (300, 600, 900, 1200, 1500, 2100, 3000, 4500)]
                if rows:
                    print("  | it | pool | violated > 1e-4 at the last separation | its max violation | cert |")
                    print("  |---|---|---|---|---|")
                    for x in rows:
                        print("  | %d | %d | %d | %.4f | %.2f |" % (x["it"], x["pool"], x["sep_viol"], x["sep_max"], x["cert"]))
        if cps:
            print("  | it | cert | sum h'mu | active | violated > 1e-4 | max viol/||G|| | ||viol/G|| | dstep | diag: |Y-Ylong|/|Ylong| | diag |G'dnu| | diag long maxviol |")
            print("  |---|---|---|---|---|---|---|---|---|---|---|")
            for c in cps:
                if c["it"] in (10, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 4000) or c["it"] == cps[-1]["it"]:
                    dg = c["diag_it"] == c["it"]
                    print("  | %d | %.2f | %.1f | %d | %d | %.2e | %.2e | %.2e | %s | %s | %s |" % (
                        c["it"], c["cert"], c["hsum"], c["active"], c["nviol"], c["maxviol"], c["rnorm"], c["dstep"],
                        "%.2e" % (c["diag_dy"] / c["diag_y"]) if dg and c["diag_y"] else "-",
                        "%.2e" % c["diag_dnu"] if dg else "-", "%.2e" % c["diag_viol"] if dg else "-"))


if __name__ == "__main__":
    main()
