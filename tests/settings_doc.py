#!/usr/bin/env python3
"""RUNNING.md's settings table is the solver's (bnbdrv --settings-table): fail
when they differ, so the documentation cannot drift from latest/settings.C.

    python3 tests/settings_doc.py            # check (make test runs it)
    python3 tests/settings_doc.py --update   # rewrite the table in RUNNING.md

The table sits between the markers <!-- settings-table-begin --> and
<!-- settings-table-end --> in RUNNING.md."""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "latest", "bnbdrv")
DOC = os.path.join(ROOT, "RUNNING.md")
BEGIN, END = "<!-- settings-table-begin -->", "<!-- settings-table-end -->"


def main():
    table = subprocess.run([SOLVER, "--settings-table"], capture_output=True, text=True, check=True).stdout
    text = open(DOC).read()
    i, j = text.find(BEGIN), text.find(END)
    if i < 0 or j < 0 or j < i:
        sys.exit(f"{DOC}: the markers {BEGIN} / {END} are missing")
    head, tail = text[:i + len(BEGIN)], text[j:]
    new = head + "\n" + table + tail
    if new == text:
        print("settings table: RUNNING.md agrees with bnbdrv --settings-table")
        return 0
    if "--update" in sys.argv[1:]:
        open(DOC, "w").write(new)
        print("settings table: RUNNING.md updated")
        return 0
    old_lines = text[i + len(BEGIN):j].strip().splitlines()
    new_lines = table.strip().splitlines()
    for a, b in zip(old_lines, new_lines):
        if a != b:
            print(f"  doc:    {a}\n  solver: {b}")
            break
    else:
        print(f"  {len(old_lines)} lines in the doc, {len(new_lines)} from the solver")
    sys.exit("settings table: RUNNING.md differs from bnbdrv --settings-table; run tests/settings_doc.py --update")


if __name__ == "__main__":
    sys.exit(main())
