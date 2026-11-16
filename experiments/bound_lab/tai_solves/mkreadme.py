"""README.md from README.tmpl.md: @FILE:logs/x@ inserts a file, @TREE:tag@ the per-depth table of logs/solve_TAG.txt,
@PATH:tag@ the optimum's path as one line (depth: shortfall of the certified bound).   python3 mkreadme.py"""
import re, os
H = os.path.dirname(os.path.abspath(__file__)); rd = lambda f: open(os.path.join(H, f)).read().rstrip("\n")
def tree(tag): return "\n".join(l for l in rd(f"logs/solve_{tag}.txt").split("\n") if l.startswith("|"))
def path(tag):
    out = []
    for l in rd(f"logs/solve_{tag}.txt").split("\n"):
        m = re.match(r"  depth (\d+) (\S) bound \S+ \((\S+) under\)", l)
        if m: out.append(f"{m[1]}: {float(m[3].replace(',', '')):,.0f}" if float(m[3].replace(',', '')) >= 1 else f"{m[1]}: {m[3]}")
    return ", ".join(out)
t = rd("README.tmpl.md")
t = re.sub(r"@FILE:([^@]+)@", lambda m: rd(m[1]), t); t = re.sub(r"@TREE:([^@]+)@", lambda m: tree(m[1]), t); t = re.sub(r"@PATH:([^@]+)@", lambda m: path(m[1]), t)
open(os.path.join(H, "README.md"), "w").write(t + "\n"); print(len(t.splitlines()), "lines")
