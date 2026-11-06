#!/usr/bin/env python3
"""The switch's saving on the production node mix [e]: a pipeline_ab node trace
(QAP_DNN_TRACE lines of the bestcore arm on the node-10 subtree) and the saving
per warm solve measured on the warm set, both in units of one ADMM iteration of
this box (the saving is a fixed cost per warm first solve; the rest of a solve is
unchanged), so that the machine's speed cancels:

  est_production.py SAVED_MS ITER_MS [trace.txt]

SAVED_MS: the mean start-up saving per warm child (summarize.py's first R-step line);
ITER_MS: the warm children's mean CPU per iteration on the same box (base config).
Default trace: ../../pipeline_ab/runs_wfix/bestcore/trace.txt (laptop, Nov 6)."""
import os
import re
import sys

D = os.path.dirname(os.path.abspath(__file__))
saved, per_it = float(sys.argv[1]), float(sys.argv[2])
path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(D, "..", "..", "pipeline_ab", "runs_wfix", "bestcore", "trace.txt")
rows = [l for l in open(path) if " fixes " in l]


def num(pat, l, conv=float, default=0):
    m = re.search(pat, l)
    return conv(m.group(1)) if m else default


first = [num(r"iters (\d+) cpu", l, int) for l in rows]
cut = [num(r"rounds \d+ iters (\d+)", l, int) for l in rows]
warm = [l for l in rows if " warm 1 " in l]
w_it = sorted(num(r"iters (\d+) cpu", l, int) for l in warm)
eq = saved / per_it   # iterations' worth saved per warm first solve
tot, tot_first = sum(first) + sum(cut), sum(first)
print(f"trace {os.path.relpath(path, D)}: {len(rows)} DNN nodes, {len(warm)} warm first solves "
      f"(median {w_it[len(w_it) // 2]} iterations, mean {sum(w_it) / len(w_it):.1f}); "
      f"{tot:,} DNN iterations ({tot_first:,} first solves, {sum(cut):,} cut rounds)")
print(f"saving {saved:.0f} ms per warm solve at {per_it:.0f} ms per iteration = {eq:.2f} iterations' worth")
print(f"[e] {len(warm) * eq:.0f} iterations' worth: {100 * len(warm) * eq / tot:.2f}% of the DNN CPU, "
      f"{100 * len(warm) * eq / tot_first:.2f}% of the first solves'")
short = [x for x in w_it if x <= 40]
print(f"[e] on a warm child of 20 iterations {100 * eq / 20:.1f}%, of 34 (the mean) {100 * eq / 34:.1f}%; "
      f"{len(short)} of {len(w_it)} warm first solves take 40 iterations or fewer")
