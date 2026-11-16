"""Where a run's iterations are: first solves of the nodes that fathom there, first solves of the nodes that give
up (the 'entry' cost before a phase), and the phases.   python3 firstsplit.py RUN [RUN ...]   (cpfast/NAME for
../../cpfast/search/NAME)"""
import os, sys
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, H)
from summarize import load
print("| run | nodes | first solves that fathom: nodes, iterations (share) | first solves that give up: nodes, iterations (share), per node | phases: iterations (share) | total |")
print("|---|---|---|---|---|---|")
for r in sys.argv[1:]:
    d = os.path.join(H, "..", "cpfast", "search", r[7:]) if r.startswith("cpfast/") else os.path.join(H, "runs", r)
    nodes, dnn = load(d)
    fd = [e["iters"] for x, e in zip(nodes, dnn) if x["out"] == "D"]
    fg = [e["iters"] for x, e in zip(nodes, dnn) if x["out"] != "D"]
    ph = sum(e["cit"] for e in dnn); tot = sum(fd) + sum(fg) + ph
    print("| `%s` | %d | %d, %s (%.0f%%) | %d, %s (%.0f%%), %.0f | %s (%.0f%%) | %s |" % (r, len(nodes), len(fd), format(sum(fd), ","), 100 * sum(fd) / tot, len(fg), format(sum(fg), ","),
          100 * sum(fg) / tot, sum(fg) / max(len(fg), 1), format(ph, ","), 100 * ph / tot, format(tot, ",")))
