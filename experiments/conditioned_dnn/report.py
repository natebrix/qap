#!/usr/bin/env python3
"""Render completed scalar results; never infer a bound from partial runs."""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    text = ["# Conditioned DNN: measured results", "",
            "Experiment off cognac `6293cb7`, 3 November 2026. All bounds below are "
            "certified lower bounds. Child-policy CPU totals include the initial parent solve; "
            "the unsplit control is charged for its full cold solve. "
            "No production search code changed.", "",
            "## Completed cases", ""]
    for name, title in [("nug20", "nug20 root, facility 0"),
                        ("wil50_node9", "wil50 node 9, facility 37")]:
        folder = HERE / "results" / name
        if not (folder / "completed.json").exists():
            raise RuntimeError(f"{name}: incomplete; refusing to publish a final report")
        s = json.loads((folder / "summary.json").read_text())
        m = json.loads((folder / "manifest.json").read_text())
        assert s["complete"] and s["children"] == len(m["locations"])
        text += [f"### {title}", "", f"Complete split: **{s['children']} children**. "
                 f"Cutoff: {s['cutoff']:.2f}. Parent/probe/additional iteration caps: "
                 f"{m['parent']}/{m['probe']}/{m['rest']}.", "",
                 "| Policy | Certified parent bound | Total solver CPU (s) | Children fathomed | Child iterations |",
                 "|---|---:|---:|---:|---:|"]
        text.append(f"| Initial parent | {s['parent']['lb']:.6f} | {s['parent']['cpu']:.3f} | — | — |")
        c = s["matched_parent_control"]
        text.append(f"| Unsplit parent, probe-matched CPU | {c['lb']:.6f} | {c['cpu']:.3f} | — | — |")
        labels = {"probe": "Conditioned: short probes", "reuse": "Probes + cached Y/Z continuation",
                  "discard": "Probes + recomputation from parent Z", "direct": "Direct branching, no probes"}
        for k, label in labels.items():
            r = s[k]
            text.append(f"| {label} | {r['lb']:.6f} | {r['cpu_including_parent']:.3f} | "
                        f"{r['children_fathomed']}/{s['children']} | {r['child_iterations']} |")
        saving = 100 * (1-s['reuse']['cpu_including_parent']/s['discard']['cpu_including_parent'])
        ratio = s['reuse']['cpu_including_parent']/s['direct']['cpu_including_parent']
        disagreements = [r['location'] for r in s['pairs']
                         if (r['reuse'] > s['cutoff']) != (r['fresh'] > s['cutoff'])]
        text += ["", f"The short conditioned bound is **{s['probe']['lb']-c['lb']:+.6f}** "
                 "relative to the CPU-matched unsplit parent. "
                 f"The latter ran {c['iterations']} iterations ({c['status']}); its stopping "
                 f"checkpoint exceeded the target by {c['cpu']-s['probe']['cpu_including_parent']:.3f} CPU seconds.", "",
                 f"Reuse saved **{saving:.2f}%** of total solver CPU versus discarding the probes. "
                 f"Its CPU ratio to direct branching was **{ratio:.4f}**. "
                 f"The final minimum bound differed from direct branching by "
                 f"**{s['reuse']['lb']-s['direct']['lb']:+.6f}**. "
                 f"Maximum absolute per-child bound difference: "
                 f"{max(abs(r['reuse_minus_fresh']) for r in s['pairs']):.6f}. "
                 "Per-child differences include different early-stop iterations and overshoot "
                 "after the cutoff; they are not certificate errors.", "",
                 f"Locations with different reuse/direct pruning decisions: **{disagreements or 'none'}**.", "",
                 f"Weakest locations: short probes {s['probe']['bottleneck_locations']}; "
                 f"reuse {s['reuse']['bottleneck_locations']}; direct {s['direct']['bottleneck_locations']}. "
                 f"Retained child Y/Z payload: {s['saved_child_state_bytes']/1e9:.4f} GB "
                 f"({s['saved_child_state_bytes']/2**30:.4f} GiB). "
                 f"Summed probe snapshot write time: {s['probe_snapshot_write_wall']:.3f} s.", "",
                 f"[Summary JSON](results/{name}/summary.json), "
                 f"[per-child CSV](results/{name}/children.csv), "
                 f"[manifest](results/{name}/manifest.json).", ""]
        with (folder / "children.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(s["pairs"][0]), lineterminator="\n")
            w.writeheader()
            w.writerows(s["pairs"])
    w = json.loads((HERE / 'results/wil50_node9/summary.json').read_text())
    c = w['matched_parent_control']
    text[4:4] = ["## Finding", "",
                 "The fixed-budget conditioning policy is **not supported for promotion on this "
                 "wil50 pilot**. Its short-probe pass gained zero over the initial parent. "
                 f"After the additional child work, reuse gives {w['reuse']['lb']:.6f} "
                 f"at {w['reuse']['cpu_including_parent']:.3f} solver CPU seconds, compared "
                 f"with {c['lb']:.6f} at {c['cpu']:.3f} seconds for the unsplit control.", "",
                 "The cache successfully reuses child computations. Its savings against "
                 "discarding probes should be distinguished from any advantage over direct "
                 "branching: without early stopping, the per-child iteration budgets are "
                 "150 for reuse, 200 for discard, and 150 for direct branching. The direct "
                 "comparison below is the relevant one for judging extra search efficiency.", "",
                 "This is a **finite-budget policy result for one wil50 node and one facility**. "
                 "Exact conditioned DNN dominates the exact parent relaxation; finite-run "
                 "certified bounds can be below those optima. The candidate retains the "
                 "200-iteration parent bound, while the separate control spends more work "
                 "on the parent itself. The experiment does not rule out a stronger parent "
                 "warm start, adaptive refinement of the weakest children, another facility, "
                 "or conditioning with inherited cuts.", ""]
    text += ["## Relation to the existing cut policy", "",
             "The archived [shorter cut-phase study](../same_node_continuation/budget_study/RESULTS.md) "
             "already reports fathoming this same wil50 node 9 with Y12s at 720 total iterations "
             "and retained bound 48816.286927. Its 1050.166 whole-process CPU seconds come "
             "from a different run/build environment and are not a matched timing control here. "
             "Consequently a gain over the uncut parent is not evidence of a gain over "
             "cognac's cut-strengthened policy. This experiment should not displace that policy "
             "on the strength of an uncut comparison.", "",
             "## Interpretation boundaries", "",
             "- The complete disjunction is essential: these minima include every free location, "
             "without symmetry reduction or cherry-picking children.",
             "- Reuse versus discard measures avoidable repeated work. Direct branching is the "
             "relevant baseline for a search-speed claim. No complete subtree search was run.",
             "- This compares uncut DNN. It does not establish superiority over the existing "
             "12-round cut-strengthened DNN policy or measure inherited cut pools.",
             "- The saved state includes Y/Z and the best certified bound, not the tracked "
             "eigensolver workspace. Thus fixed-total-iteration results can differ from direct "
             "solves; no exact continuation claim is made for tracked mode.",
             "- Three one-thread child workers ran concurrently; the parent control overlapped "
             "the second phase. CPU ratios are diagnostic on this machine, not controlled "
             "sequential wall-clock speedups. Disk I/O, cache hashing, scheduling, and file reads "
             "are excluded from solver CPU; raw process wall times are retained.",
             "- Child states are kept on disk, one loaded per active worker. An in-memory "
             "implementation retaining every child would need the reported aggregate payload. "
             "A production cache would need an explicit memory/eviction policy.", "",
             "## Validation and provenance", "",
             "- Ten fully enumerated signed/asymmetric QAPs of orders 3–7: 50 children, "
             "all inherited/probe/reuse/fresh bounds safe; disjunction coverage exact; "
             "50 wrong-node snapshot rejections. Full-spectrum state and bound continuation "
             "errors were both zero.",
             "- Existing continuation checks: 72 identity cases, 648 safe certificates, "
             "72 rescaling cases, 108 changing-cut cases, seven malformed states; passed.",
             "- `make test`: 4,536 brute-force runs, five resume cases, parallel tests, "
             "and seven exact regression matches; passed. Logs are in `results/`.",
             "- The post-run [audit](results/validation.json) checks all raw child logs, "
             "commands/budgets, complete coverage, cache-bound inheritance, input/binary "
             "hashes, summary hashes, and every locally retained child-state hash.",
             "- Measurements used local commit `0ae156a31ae9c2f0af53abe4e2f59c83aef34dd2`. "
             "GitHub publication through the connector produced commit "
             "`b36b6b6cbe71e56f919189e8f4a25b58b00fcf9c` with the **identical tree** "
             "`994a1903d6f002be32c0923c3e3cf9f693077085`. Reporting files were added afterward; "
             "the measured executable was not changed.",
             "- Intel Xeon Platinum 8573C; one OpenBLAS thread per solver; "
             "SciPy LP64 OpenBLAS, compiler `-O3 -march=native -ffp-contract=off`. "
             "AMX was not enabled for these runs. Full command lines and executable/input "
             "SHA-256 hashes are retained in the manifests.", ""]
    (HERE / "RESULTS.md").write_text("\n".join(text))


if __name__ == "__main__":
    main()
