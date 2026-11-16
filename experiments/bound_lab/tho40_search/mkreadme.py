"""Assemble README.md from README.tmpl.md: lines '@@ CMD' are replaced by the output of 'python3 CMD' run here."""
import os, subprocess
H = os.path.dirname(os.path.abspath(__file__))
out = []
for ln in open(os.path.join(H, "README.tmpl.md")):
    if ln.startswith("@@ "):
        out.append(subprocess.run("python3 " + ln[3:].strip(), shell=True, cwd=H, capture_output=True, text=True).stdout)
    else:
        out.append(ln)
open(os.path.join(H, "README.md"), "w").write("".join(out))
