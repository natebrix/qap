# Replay the primal gate on every logged cut phase: after round R, stop when the
# best (largest) Lagrangian primal value seen so far is more than M below the
# cutoff. A phase that the gate would stop but that later fathomed (any round or
# the quiet round) is a FALSE STOP, the kill criterion of roadmap/tree_size.md
# item 5. Logs: dnncuts output (round lines with "primal", a header with
# "cutoff", a FINAL line). Usage: python3 replay_pgate.py R M log...
import re, sys
R = int(sys.argv[1]); M = float(sys.argv[2])
false_stops = true_stops = passed_fathom = passed_miss = 0
rows = []
for path in sys.argv[3:]:
    txt = open(path, errors="replace").read()
    mc = re.search(r"cutoff\s+([0-9.]+)", txt)
    if not mc: continue
    cutoff = float(mc.group(1))
    rounds = re.findall(r"^round\s+(\d+)\s+lb\s+([0-9.]+).*?primal\s+([0-9.]+)\s+status\s+(\w+)", txt, re.M)
    rounds = [(int(a), float(b), float(c), d) for a, b, c, d in rounds if int(a) >= 1]
    if len(rounds) < R: continue
    quiet = re.search(r"^quiet round\s+lb\s+([0-9.]+).*?status\s+(\w+)", txt, re.M)
    final = re.search(r"FINAL best_lb\s+([0-9.]+)", txt)
    best_lb = float(final.group(1)) if final else max(r[1] for r in rounds)
    fathomed = best_lb > cutoff or any(r[3] == "fathomed" for r in rounds) or (quiet and quiet.group(2) == "fathomed")
    best_primal_R = max(r[2] for r in rounds[:R])
    fires = best_primal_R < cutoff - M
    nrounds = len(rounds) + (1 if quiet else 0)
    tag = ("FALSE STOP" if fires and fathomed else "stop" if fires else "pass/fathom" if fathomed else "pass/miss")
    if fires and fathomed: false_stops += 1
    elif fires: true_stops += 1
    elif fathomed: passed_fathom += 1
    else: passed_miss += 1
    rows.append((tag, path, cutoff - best_primal_R, cutoff - best_lb, nrounds))
for tag, path, dp, dl, nr in sorted(rows):
    print(f"{tag:12s} {path:75s} primal short by {dp:8.2f} after round {R}; final lb short by {dl:8.2f}; {nr} rounds")
print(f"\nR={R} M={M}: {false_stops} false stops, {true_stops} true stops (misses caught), "
      f"{passed_fathom} fathoms passed, {passed_miss} misses passed")
