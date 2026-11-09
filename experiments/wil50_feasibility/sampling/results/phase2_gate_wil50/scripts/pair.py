"""two replays of the same task against each other: trace (pid, cpu aside), log counters, nodes, warm lines"""
import os, sys
G = "/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/gate"
sys.path.insert(0, G)
from compare import State, log_counts, norm_trace  # noqa: E402

a_side, b_side = sys.argv[1], sys.argv[2]
for tid in sys.argv[3:]:
    A, B = os.path.join(G, "runs", a_side, tid), os.path.join(G, "runs", b_side, tid)
    sa, sb = State(os.path.join(A, "out.state")), State(os.path.join(B, "out.state"))
    tr = norm_trace(os.path.join(A, "trace.txt")) == norm_trace(os.path.join(B, "trace.txt"))
    lg = log_counts(os.path.join(A, "solver.log")) == log_counts(os.path.join(B, "solver.log"))
    nodes = [t for t, _ in sa.nodes] == [t for t, _ in sb.nodes]
    warm = [w.text if w else None for _, w in sa.nodes] == [w.text if w else None for _, w in sb.nodes]
    hdr = (sa.root_bound, sa.obj, sa.perm) == (sb.root_bound, sb.obj, sb.perm)
    print(f"{a_side} vs {b_side} {tid}: trace {tr}, log counters {lg}, header {hdr} (root {sa.root_bound} / {sb.root_bound}), "
          f"{len(sa.nodes)}/{len(sb.nodes)} node lines identical {nodes}, warm lines identical {warm}")
