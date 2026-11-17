#!/usr/bin/env python3
"""Summarize the matched tai60b continuation using its certified checks."""

import argparse
import json
from pathlib import Path

from run_gate import compare, parse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace_dir", type=Path)
    parser.add_argument("early_summary", type=Path)
    args = parser.parse_args()
    early = json.loads(args.early_summary.read_text())
    regional = parse(args.trace_dir / "regional.log")
    coarse = parse(args.trace_dir / "coarse.log")
    comparison = compare(regional, coarse, early["separation_cpu"])
    common = {}
    early_max_diff = None
    for arm, result in (("regional", regional), ("coarse", coarse)):
        if result["status"] != "fathomed" and result["iters"] != 600:
            raise ValueError(f"{arm} stopped at iteration {result['iters']} ({result['status']})")
        common[arm] = {p["it"]: p for p in result["trajectory"]}
        replay_diff = max(abs(common[arm][p["it"]]["best"] - p["best"]) for p in early[arm]["trajectory"])
        early_max_diff = max(early_max_diff or 0, replay_diff)
    if early_max_diff > 1e-5:
        raise ValueError(f"the early certified trajectory changed by {early_max_diff}")

    cutoff = early["ub"] - 0.95
    milestones = []
    for it in (100, 120, 200, 300, 400, 500, 600):
        if it not in common["regional"] or it not in common["coarse"]:
            continue
        a, b = common["regional"][it]["best"], common["coarse"][it]["best"]
        milestones.append({"it": it, "regional": a, "coarse": b, "gain": b - a,
                           "remaining_regional_shortfall": cutoff - a})

    targets = []
    for million in (590, 594, 596, 598, 600, 602, 603, 604, 605, 606, 607):
        target = million * 1_000_000
        reached = {arm: next((p["it"] for p in result["trajectory"] if p["best"] >= target), None)
                   for arm, result in (("regional", regional), ("coarse", coarse))}
        targets.append({"target": target, **reached})

    out = {"case": early["case"], "fixes": early["fixes"], "cutoff": cutoff, "max_iter": 600,
           "capacity": early["capacity"], "n_cuts": early["n_cuts"], "earlier_separation_cpu": early["separation_cpu"],
           "concurrent_arms": True, "early_check_max_difference": early_max_diff,
           "regional": regional, "coarse": coarse, "comparison": comparison,
           "milestones": milestones, "iterations_to_targets": targets}
    (args.trace_dir / "summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"milestones": milestones, "iterations_to_targets": targets}, indent=2))


if __name__ == "__main__":
    main()
