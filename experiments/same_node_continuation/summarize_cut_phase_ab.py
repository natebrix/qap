"""Summarize a completed run_cut_phase_ab.py experiment without hiding failures."""
import argparse
import json
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("directory", type=Path)
args = ap.parse_args()
manifest = json.loads((args.directory / "manifest.json").read_text())
jobs = manifest["jobs"]
if len(jobs) != 4 or any(j["state"] != "completed" for j in jobs):
    raise SystemExit("All four arms must complete successfully before comparison")
rows = []
for node in sorted({j["node"] for j in jobs}):
    pair = {j["arm"]: j for j in jobs if j["node"] == node}
    a, b = pair["Z_only"], pair["YZ"]
    assert a["command"] == [x for x in b["command"] if x != "-warmY"]
    assert a["final"]["first_lb"] == b["final"]["first_lb"], "Uncut starts differ"
    assert a["rounds"][0]["lb"] == b["rounds"][0]["lb"]
    rows.append(dict(node=node, baseline_lb=a["final"]["best_lb"],
                     warmY_lb=b["final"]["best_lb"],
                     bound_gain=b["final"]["best_lb"]-a["final"]["best_lb"],
                     baseline_fathomed=a["fathomed"], warmY_fathomed=b["fathomed"],
                     baseline_iterations=a["iterations"], warmY_iterations=b["iterations"],
                     baseline_process_cpu=a["process_cpu_seconds"],
                     warmY_process_cpu=b["process_cpu_seconds"],
                     cpu_ratio=b["process_cpu_seconds"]/a["process_cpu_seconds"],
                     baseline_wall=a["wall_seconds"], warmY_wall=b["wall_seconds"]))
summary = dict(source_commit=manifest["commit"],
               elapsed_wall_seconds=manifest["elapsed_wall_seconds"],
               timing_note=manifest["timing_note"], pairs=rows)
(args.directory / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
print("| Node | Z-only bound | Y/Z bound | Difference | Fathomed Z / YZ | Iterations Z / YZ | Process CPU ratio YZ/Z |")
print("|---|---:|---:|---:|---|---:|---:|")
for r in rows:
    print(f"| {r['node']} | {r['baseline_lb']:,.6f} | {r['warmY_lb']:,.6f} | {r['bound_gain']:+.6f} | "
          f"{r['baseline_fathomed']} / {r['warmY_fathomed']} | {r['baseline_iterations']} / {r['warmY_iterations']} | {r['cpu_ratio']:.4f} |")
