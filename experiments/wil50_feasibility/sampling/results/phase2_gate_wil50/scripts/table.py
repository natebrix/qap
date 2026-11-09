import json
G = "/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/gate"
rows = json.load(open(G + "/compare.json"))
CLS = {"r.9.39": "first solve fathomed", "r.4.39": "first solve fathomed", "r.3.36": "CP phase fathomed",
       "r.7.35": "CP phase fathomed", "r.6.41": "CP phase stalled", "r.2.26": "CP phase stalled",
       "r.3.15": "CP phase iteration limit", "r.13.4": "CP phase iteration limit", "r.1": "depth 1 (from the root)",
       "r": "root (QPB only)"}
tot = 0
for r in rows:
    tot += r["cpu_replay"]
    t = r["rep_trace"]
    o = r["rec_trace"]
    if t:
        cut = f"{t['cut_status']} {t['cut_best']:.4f} / {t['cut_iters']} it / {t['cut_rounds']} rounds" if t["cut"] else "-"
        tr = f"{t['lb']:.4f} | {t['status']} {t['iters']} | {cut} | fathomed {t['fathomed']} | pool_out {t['pool_out']}"
    else:
        tr = "no DNN"
    kids = r["n_children_replay"]
    kb = r.get("children_bound")
    pm = r.get("children_pred_min_max")
    print(f"{r['side']:4s} {r['task']:7s} {CLS[r['task']]:26s} [{r['spec']}] {tr} | kids {r['n_children_record']}/{kids} "
          f"line {r['line_record']}/{r['line_replay']} inherited {kb} pred {pm} pred_maxdiff {r.get('children_pred_max_absdiff')} "
          f"| trace {r['trace_identical']} kids_vals {r['children_values_identical']} warm {r.get('warm_identical')} "
          f"| cpu rec {r['cpu_record']:.1f} rep {r['cpu_replay']:.1f}")
    if pm and kb and pm[1] >= min(kb):
        print("   NOTE predicted bound reaches the inherited bound")
print(f"replay CPU total {tot:.0f} s over {len(rows)} runs")
