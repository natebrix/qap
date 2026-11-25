#!/usr/bin/env python3
"""Figure for dnn_plus.tex: the certified DNN+ value against cut-phase iterations,
as the gap to the best known value, with the plain DNN bound (4,000 iterations) as a line.
Reads runs_long/P/plus.trace and runs_serial/results.txt; writes fig_history.pdf.
  python3 docs/dnn_plus/fig_history.py"""
import os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
BEST = {"tho40": 240516, "tai50b": 458821517}
plain = {}
for l in open(f"{H}/runs_serial/results.txt"):
    f = l.split()
    if f[1] == "plain": plain[f[0]] = float(f[3])
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.6))
for ax, p in zip(axs, BEST):
    h = [(int(a), float(b)) for a, b in re.findall(r"(?<= )(\d+):([\d.]+)(?= |$)", open(f"{H}/runs_long/{p}/plus.trace").read().split("hist", 1)[1])]
    gap = lambda v: 100 * (BEST[p] - v) / BEST[p]
    ax.plot([i for i, _ in h], [gap(v) for _, v in h], color="black", lw=1.2)
    ax.axhline(gap(plain[p]), color="0.45", lw=0.8, ls="--")
    ax.text(20000, gap(plain[p]), "DNN ", ha="right", va="bottom", color="0.3")
    ax.text(20000, gap(h[-1][1]), "DNN$^+$ ", ha="right", va="bottom")
    ax.set_title(p, fontsize=9); ax.set_xlabel("cutting-plane iterations"); ax.set_ylim(0, None)
    ax.set_xticks(range(0, 20001, 5000)); ax.grid(axis="y", color="0.9", lw=0.5)
axs[0].set_ylabel("gap to best known value (%)")
fig.tight_layout(); fig.savefig(f"{H}/fig_history.pdf")
