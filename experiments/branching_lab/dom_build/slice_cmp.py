#!/usr/bin/env python3
"""Node-for-node comparison of a key-off and a key-on trace of the same search,
and of both against the replay's prediction on the production trace.
  python3 slice_cmp.py OFF.trace.gz ON.trace.gz [INSTANCE PRODUCTION.trace.gz ROOTFIXES]

Prints: nodes (trace lines) and CPU per arm; the X nodes of the on arm, whether
each is a node of the off arm and how many off nodes lie under it; whether
on = off less the subtrees under the X nodes (exactly, as sets of fixes, with
the same outcome for every common node); and, with a production trace, that
subtree's nodes, the replay's predicted drop (domcheck), and the overlap."""
import gzip, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))


def load(p):
    op = gzip.open if p.endswith(".gz") else open
    rows = [l.split() for l in op(p, "rt") if l.startswith("node ")]
    return {r[2]: r for r in rows}, rows


def under(keys, tops):
    tops = set(tops); out = set()
    for k in keys:
        parts = k.split(",")
        if any(",".join(parts[:i]) in tops for i in range(1, len(parts) + 1)):
            out.add(k)
    return out


off, offr = load(sys.argv[1]); on, onr = load(sys.argv[2])
cpu = lambda rows: sum(float(r[11]) for r in rows)
X = [k for k, r in on.items() if r[7] == "X"]
print(f"off: {len(offr)} lines, {cpu(offr)/3600:.3f} CPU h; on: {len(onr)} lines ({len(X)} X), {cpu(onr)/3600:.3f} CPU h; "
      f"saved {len(offr) - len(onr) + len(X)} bounded nodes, {100*(1 - cpu(onr)/cpu(offr)):.1f}% CPU")
notin = [k for k in X if k not in off]
gone = under(off.keys(), X)
expect = set(off) - gone | set(X)
print(f"X nodes in the off tree: {len(X) - len(notin)} of {len(X)}; off nodes at or under them: {len(gone)}")
same = set(on) == expect
diffout = [k for k in set(on) & set(off) if on[k][7] != off[k][7] and on[k][7] != "X"]
print(f"on == off less the X subtrees (+ the X nodes themselves): {same}; common nodes with another outcome: {len(diffout)}")
if not same:
    print("  only in on:", sorted(set(on) - expect)[:5], " only in expected:", sorted(expect - set(on))[:5])
if len(sys.argv) > 3:
    name, prod, root = sys.argv[3:6]
    out = subprocess.run(f"zcat {prod} | {HERE}/bin/domcheck /home/user/wt-dombuild/qaplib/{name}.dat | grep -v -e '^objective' -e '^symm'",
                         shell=True, capture_output=True, text=True).stdout.split("\n")
    prow = [l.split() for l in gzip.open(prod, "rt") if l.startswith("node ")]
    sub = {r[2]: r for r in prow if r[2] == root or r[2].startswith(root + ",")}
    pdom = {t[2] for t, l in zip(prow, out) if t[2] in sub and len(l.split()) == 4 and l.split()[3] == "1"}
    ptops = [k for k in pdom if not any(",".join(k.split(",")[:i]) in pdom for i in range(1, len(k.split(","))))]
    pgone = under(sub.keys(), ptops)
    print(f"production subtree: {len(sub)} nodes, {cpu(sub.values())/3600:.3f} CPU h; replay: {len(ptops)} highest dominated, "
          f"{len(pgone)} nodes under them, {cpu([sub[k] for k in pgone])/3600:.3f} CPU h")
    print(f"  production nodes also in the off run: {len(set(sub) & set(off))} of {len(sub)} (off run {len(off)}); "
          f"replay's highest dominated also X in the on run: {len(set(ptops) & set(X))} of {len(ptops)}; "
          f"X in the on run not predicted: {len(set(X) - set(ptops))}")
