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
    early_max_diff = 0
    for arm, result in (("regional", regional), ("coarse", coarse)):
        if result["status"] != "fathomed" and result["iters"] != 600:
            raise ValueError(f"{arm} stopped at iteration {result['iters']} ({result['status']})")
        common[arm] = {p["it"]: p for p in result["trajectory"]}
        replay_diff = max(abs(common[arm][p["it"]]["best"] - p["best"]) for p in early[arm]["trajectory"])
        early_max_diff = max(early_max_diff, replay_diff)
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
    for target in (590_000_000, 594_000_000, 596_000_000, 598_000_000, 600_000_000, 602_000_000,
                   602_800_000, 603_000_000, 603_100_000, 603_200_000, 604_000_000, 605_000_000, 606_000_000,
                   607_000_000):
        reached = {}
        for arm, result in (("regional", regional), ("coarse", coarse)):
            original_separation = early["separation_cpu"] if arm == "coarse" else 0
            reached[arm] = next(({"iteration": p["it"], "cpu_with_setup": p["cpu"] + result["coarse_setup_cpu"],
                                  "cpu_with_setup_and_original_separation": p["cpu"] + result["coarse_setup_cpu"] + original_separation}
                                 for p in result["trajectory"] if p["best"] >= target), None)
        targets.append({"target": target, **reached})

    out = {"case": early["case"], "fixes": early["fixes"], "cutoff": cutoff, "max_iter": 600,
           "capacity": early["capacity"], "n_cuts": early["n_cuts"], "earlier_separation_cpu": early["separation_cpu"],
           "concurrent_arms": True, "early_check_max_difference": early_max_diff,
           "regional": regional, "coarse": coarse, "comparison": comparison,
           "comparison_replay_only": compare(regional, coarse, 0),
           "milestones": milestones, "iterations_to_targets": targets}
    (args.trace_dir / "summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"milestones": milestones, "iterations_to_targets": targets}, indent=2))


if __name__ == "__main__":
    main()
