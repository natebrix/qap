#!/usr/bin/env python3
"""E1 table: per node, the iteration at which each certificate variant first
crosses the cutoff (same iterates), the iterations saved against the standard
certificate, the gains at fixed iterations, and the CPU per evaluation.
Usage: e1_table.py log [log ...]
"""
import re, statistics, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from summarize import parse  # noqa: E402


def main():
    rows = []
    for path in sys.argv[1:]:
        for c in parse(path):
            ch = c["checks"]
            if not ch:
                continue
            cut = ch[0]["cutoff"]
            def first(key):
                for x in ch:
                    v = key(x)
                    if v > cut:
                        return x["it"]
                return None
            fs = {"std": first(lambda x: x["std"]), "exact": first(lambda x: x["exact"]),
                  "hg": first(lambda x: max([v for _, v in x["hg"]] or [-1e300])),
                  "hgr": first(lambda x: max([v for _, v in x["hgr"]] or [-1e300])),
                  "avg": first(lambda x: x["avg"]), "ls": first(lambda x: x["ls"]), "proj": first(lambda x: x["proj"])}
            # gains at the last evaluated check before the standard crossing
            pre = [x for x in ch if x["std"] <= cut]
            last = pre[-1] if pre else ch[-1]
            g = lambda v: v - last["std"]
            hgv = max([v for _, v in last["hg"]] or [float("-inf")])
            hgrv = max([v for _, v in last["hgr"]] or [float("-inf")])
            cpu = {k: statistics.median([x["cpu"][i] for x in ch]) for k, i in (("ex", 0), ("hg", 1), ("ls", 2), ("avg", 3))}
            cpu["hgr"] = statistics.median([x["hgr_cpu"] for x in ch])
            cpu["proj"] = statistics.median([x["proj_cpu"] for x in ch])
            s = c["summary"] or ""
            ms = re.search(r"status (\S+) iters (\d+)", s)
            rows.append((path.rsplit("/", 1)[-1], ms.group(1) if ms else "?", int(ms.group(2)) if ms else -1, fs, last,
                         g(last["exact"]), g(hgv), g(hgrv), g(last["avg"]), g(last["ls"]), g(last["proj"]), cpu))
    print("| log | status | iters | first above cutoff: std / exact / hg / hgr / avg / ls / proj | at it (short, primal-std) | gain exact / hg / hgr / avg / ls / proj | cpu s per eval ex / hg / hgr / avg / ls / proj |")
    print("|---|---|---|---|---|---|---|")
    for name, st, it, fs, last, ge, gh, ghr, ga, gl, gp, cpu in rows:
        f = lambda k: str(fs[k]) if fs[k] is not None else "-"
        print(f"| {name} | {st} | {it} | {f('std')} / {f('exact')} / {f('hg')} / {f('hgr')} / {f('avg')} / {f('ls')} / {f('proj')} | "
              f"{last['it']} ({last['cutoff'] - last['std']:.2f}, {last['primal'] - last['std']:.2f}) | "
              f"{ge:+.2f} / {gh:+.2f} / {ghr:+.2f} / {ga:+.2f} / {gl:+.2f} / {gp:+.2f} | "
              f"{cpu['ex']:.2f} / {cpu['hg']:.2f} / {cpu['hgr']:.2f} / {cpu['avg']:.2f} / {cpu['ls']:.2f} / {cpu['proj']:.2f} |")


main()
