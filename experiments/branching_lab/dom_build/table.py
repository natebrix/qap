#!/usr/bin/env python3
"""The validity table from run_valid.py's JSON lines (README).
  python3 table.py logs/valid_opt1.jsonl [logs/valid_poor.jsonl ...]"""
import json, sys
for path in sys.argv[1:]:
    rows = [json.loads(l) for l in open(path)]
    by = {}
    for r in rows:
        by.setdefault(r["instance"], {})[r["key"]] = r
    mode = rows[0]["mode"] if rows else "?"
    print(f"\n### {path} (incumbent mode {mode})\n")
    print("| instance | n | optimum (source) | start incumbent | sym file | test available | nodes off | nodes on | on/off | dropped | final off / on | open off / on | valid |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    tot = {"inst": 0, "valid": 0, "off": 0, "on": 0, "drop": 0, "fired": 0}
    for name, d in by.items():
        if 0 not in d or 1 not in d: continue
        a, b = d[0], d[1]
        ok = a["valid"] and b["valid"]
        tot["inst"] += 1; tot["valid"] += ok
        no, nn = a.get("nodes"), b.get("nodes")
        if no and nn: tot["off"] += no; tot["on"] += nn
        tot["drop"] += b.get("dropped") or 0; tot["fired"] += (b.get("dropped") or 0) > 0
        av = "no: " + b["dom_unavailable"] if b.get("dom_unavailable") else "yes"
        src = a["opt_src"] + (f", {a['n_optima']} optima" if a.get("n_optima") else "")
        ub = a.get("ub_used") if a.get("ub_used") is not None else "none"
        print(f"| {name} | {a['n']} | {a['opt']} ({src}) | {ub} | {'yes' if a['sym_file'] else ''} | {av} | {no} | {nn} | "
              f"{(nn / no if no and nn else float('nan')):.3f} | {b.get('dropped')} | {a.get('final_incumbent')} / {b.get('final_incumbent')} | "
              f"{a.get('open_nodes')} / {b.get('open_nodes')} | {'yes' if ok else 'NO (' + a['status'] + '/' + b['status'] + ')'} |")
    print(f"\n{tot['inst']} instances, {tot['valid']} valid on both arms; test fired on {tot['fired']}; "
          f"nodes {tot['off']} off, {tot['on']} on ({tot['on'] / max(tot['off'], 1):.3f}), {tot['drop']} dropped")
