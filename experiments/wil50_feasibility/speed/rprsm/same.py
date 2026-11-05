#!/usr/bin/env python3
"""Whether two dnntest logs agree line by line once the timings (cpu,
ms_per_iter, profile lines) are removed: the default-path identity check.
   usage: same.py a.log b.log [a.log b.log ...]"""
import re, sys


def norm(path):
    out = []
    for line in open(path):
        if line.startswith(("profile", "eig split", "step split", "form R")):
            continue
        out.append(re.sub(r"(cpu[= ]|ms_per_iter=)[0-9.]+", "", line))
    return out


ok = True
for a, b in zip(sys.argv[1::2], sys.argv[2::2]):
    same = norm(a) == norm(b)
    ok &= same
    print("%s %s %s (%d lines)" % (a, "==" if same else "DIFFERS", b, len(norm(a))))
sys.exit(0 if ok else 1)
