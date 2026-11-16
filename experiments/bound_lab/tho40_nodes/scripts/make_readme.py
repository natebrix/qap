#!/usr/bin/env python3
"""Assemble README.md from README.tmpl.md: a line "{{script.py args}}" is replaced by the
output of scripts/script.py (the tables are always those of the logs as they stand)."""
import os, re, subprocess, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
out = []
for line in open(os.path.join(D, "README.tmpl.md")):
    m = re.match(r"^\{\{(\S+)(.*)\}\}\s*$", line)
    if m:
        r = subprocess.run([sys.executable, os.path.join(D, "scripts", m.group(1))] + m.group(2).split(), capture_output=True, text=True)
        out.append(r.stdout.rstrip("\n") + "\n")
        if r.returncode != 0: out.append("(script failed: %s)\n" % r.stderr.strip()[-300:])
    else:
        out.append(line)
open(os.path.join(D, "README.md"), "w").write("".join(out))
