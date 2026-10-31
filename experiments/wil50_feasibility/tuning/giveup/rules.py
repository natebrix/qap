#!/usr/bin/env python3
"""Evaluate variants of the hopeless rule on -v traces: fire at the first check
>= H where best + F * d1 r / (1 - r) < cutoff (windows of w checks). For the
fathoming traces it must never fire before the fathom; for the misses, the
earlier it fires the better.

  python3 rules.py cutoff fathom_traces... -- miss_traces...
"""
import re, sys
cut = float(sys.argv[1]); args = sys.argv[2:]; i = args.index("--")
fath_files, miss_files = args[:i], args[i + 1:]
def load(f):
    pts = [(int(m.group(1)), float(m.group(2))) for m in re.finditer(r"^it\s+(\d+) lb (\S+)", open(f).read(), re.M)]
    out = []; b = -float("inf")
    for it, lb in pts:
        b = max(b, lb); out.append((it, b))
    return out
def fire(best, H, F, w=5):
    for h in range(2 * w + 1, len(best) + 1):
        it, bb = best[h - 1]
        if bb > cut: return None          # fathomed first
        if it < H: continue
        d1 = bb - best[h - 1 - w][1]; d2 = best[h - 1 - w][1] - best[h - 1 - 2 * w][1]
        r = min(d1 / d2, 0.95) if d2 > 0 else (0.95 if d1 > 0 else 0.0)
        if bb + F * d1 * r / (1 - r) < cut: return it
    return "never"
F_ = [load(f) for f in fath_files]; M_ = [load(f) for f in miss_files]
print(f"{'H':>4s} {'F':>4s} {'false give-ups':>14s}  give-up iteration on the misses")
for H in (100, 150, 200, 300):
    for F in (5, 10, 20, 40):
        false = sum(1 for b in F_ if fire(b, H, F) not in (None, "never"))
        miss = [fire(b, H, F) for b in M_]
        print(f"{H:4d} {F:4d} {false:9d}/{len(F_):<4d}  {miss}")
