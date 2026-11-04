#!/usr/bin/env python3
"""Summarize DNN_CERTX logs (run_e1.sh / run_e2.sh).

Per dnn_bound call: iterations, status, the first check at which each variant
crosses the cutoff, and at each evaluated check the gain of every variant over
the standard certificate and the fraction of the (primal - standard) lag it
closes. For E2 also the per-round end values against the quiet round.
Usage: summarize.py log [log ...]
"""
import re, sys

PAT = re.compile(r"certx it\s+(\d+) std (\S+) exact (\S+) hg(.*?) ls (\S+) \(t (\S+)\) avg(\d+) (\S+) "
                 r"primal (\S+) cutoff (\S+) cpu (\S+)/(\S+)/(\S+)/(\S+)")


def f(x):
    return float(x) if x not in ("-inf", "inf", "nan") else float("-inf")


def parse(path):
    calls, cur = [], None
    for line in open(path):
        if line.startswith("certx call"):
            cur = {"checks": [], "summary": None, "head": line.strip()}
            calls.append(cur)
        elif line.startswith("certx it") and cur is not None:
            m = PAT.search(line)
            if not m:
                continue
            hg = [(int(a), float(b)) for a, b in re.findall(r" (\d+):(-?[\d.]+)", m.group(4))]
            mr = re.search(r" hgr(.*?) cpu (\S+)", line)
            mp = re.search(r" proj (\S+) cpup (\S+)", line)
            hgr = [(int(a), float(b)) for a, b in re.findall(r" (\d+):(-?[\d.]+)", mr.group(1))] if mr else []
            cur["checks"].append(dict(it=int(m.group(1)), std=f(m.group(2)), exact=f(m.group(3)), hg=hg, hgr=hgr,
                                      hgr_cpu=float(mr.group(2)) if mr else 0.0,
                                      proj=f(mp.group(1)) if mp else float("-inf"), proj_cpu=float(mp.group(2)) if mp else 0.0,
                                      ls=f(m.group(5)), t=float(m.group(6)), avg=f(m.group(8)),
                                      primal=f(m.group(9)), cutoff=f(m.group(10)),
                                      cpu=[float(m.group(i)) for i in range(11, 15)]))
        elif line.startswith("certx-summary") and cur is not None:
            cur["summary"] = line.strip()
    return calls


def main():
    for path in sys.argv[1:]:
        print(f"== {path}")
        for k, c in enumerate(parse(path)):
            s = c["summary"] or ""
            ms = re.search(r"status (\S+) iters (\d+).*first std (-?\d+) exact (-?\d+) hg (-?\d+) ls (-?\d+) avg (-?\d+)"
                           r"(?: hgr (-?\d+))?.*cpu_certx (\S+)", s)
            if ms:
                print(f" call {k}: {ms.group(1)} iters {ms.group(2)} first-cross std {ms.group(3)} exact {ms.group(4)}"
                      f" hg {ms.group(5)} hgr {ms.group(8)} ls {ms.group(6)} avg {ms.group(7)}"
                      f" (instrumentation cpu {ms.group(9)} s)")
            for ch in c["checks"]:
                hgl = ch["hg"][-1][1] if ch["hg"] else float("-inf")
                hg2 = ch["hg"][1][1] if len(ch["hg"]) > 1 else float("-inf")
                hg4 = ch["hg"][2][1] if len(ch["hg"]) > 2 else float("-inf")
                hgr = max(v for _, v in ch["hgr"]) if ch["hgr"] else float("-inf")
                lag = ch["primal"] - ch["std"]
                def frac(v):
                    return (v - ch["std"]) / lag if lag > 1e-9 and v > -1e300 else float("nan")
                print(f"   it {ch['it']:5d} short {ch['cutoff'] - ch['std']:8.2f} lag(primal-std) {lag:7.2f} | "
                      f"exact {ch['exact'] - ch['std']:+6.2f} hg2 {hg2 - ch['std']:+6.2f} hg4 {hg4 - ch['std']:+6.2f} "
                      f"hgP {hgl - ch['std']:+6.2f} ({100 * frac(hgl):3.0f}% of lag) hgr {hgr - ch['std']:+6.2f} "
                      f"ls {ch['ls'] - ch['std']:+7.2f} (t {ch['t']:.2f}) avg {ch['avg'] - ch['std']:+7.2f} "
                      f"| cpu ex {ch['cpu'][0]:.2f} hg {ch['cpu'][1]:.2f} hgr {ch['hgr_cpu']:.2f} "
                      f"ls {ch['cpu'][2]:.2f} avg {ch['cpu'][3]:.2f}")


if __name__ == "__main__":
    main()
