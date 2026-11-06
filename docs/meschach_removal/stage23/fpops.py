#!/usr/bin/env python3
"""For a function that changed (fndiff.py), list its floating-point
arithmetic instructions and its calls, in address order, in the old and new
object, and say whether the two lists are equal. Register allocation may
differ, so only the mnemonic is compared for arithmetic.
   python3 fpops.py OLD.o NEW.o FUNCTION_SUBSTRING"""
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from fndiff import functions  # noqa: E402

FP = re.compile(r"^(v?(add|sub|mul|div|sqrt|min|max|fmadd\w*|fmsub\w*|fnmadd\w*|fnmsub\w*)[sp]d"
                r"|v?cvt\w+|v?u?comi[sp]d|v?xorpd|v?andpd|v?andnpd)\b")


def summary(lines):
    ops, calls = [], []
    for i, line in enumerate(lines):
        m = FP.match(line)
        if m:
            ops.append(m.group(1))
        if line.startswith("R_X86_64_PLT32"):
            calls.append(line.split("\t", 1)[1].rsplit("-0x", 1)[0])
    return ops, calls


def main():
    a, b = functions(sys.argv[1]), functions(sys.argv[2])
    key = sys.argv[3]
    for name in sorted(set(a) & set(b)):
        if key not in name:
            continue
        oa, ca = summary(a[name])
        ob, cb = summary(b[name])
        print(name)
        print(f"  floating-point ops: {len(oa)} old, {len(ob)} new, "
              f"{'same sequence' if oa == ob else 'DIFFERENT'}")
        print(f"  old calls: {' '.join(ca)}")
        print(f"  new calls: {' '.join(cb)}")


if __name__ == "__main__":
    main()
