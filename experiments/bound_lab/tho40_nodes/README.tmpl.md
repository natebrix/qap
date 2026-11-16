# tho40 bound lab: the grid cut phase per node (Nov 15 2026)

Agent `tho40_nodes` of `experiments/bound_lab/` (see `../COMMON.md`). Question: is the
regional-triangle cut phase (CP-ADMM form, `latest/dnncut.h`) valid on tho40, and how much
does it add to the certified DNN bound, node by node? Tools as built in the main checkout:
`experiments/dnn_cuts/dnncuts`, `latest/dnntest`. Nothing under `latest/` was changed.

Marks: **[m]** measured here (the log is named), **[e]** estimated. The machine was shared
(load about 6 on 4 cores throughout), so iterations are the metric; CPU seconds are indicative.
This file is assembled by `scripts/make_readme.py` from `README.tmpl.md`; the tables are the
scripts' output on the logs.

## Summary

@@SUMMARY@@

## Files

| | |
|---|---|
| `scripts/path.py` | the best known solution (`QAPSOLN/tho40.sln`, inverted for `tho40t.dat`), its four grid images (all 240,516), node specs along it |
| `scripts/part1_jobs.py`, `scripts/p1_summary.py` | part 1's runs and table (`logs/p1_*.log`) |
| `scripts/sample_nodes.py`, `nodes.tsv` | part 2's sample from `infra/phase1/runs/tho40ab-32/nodes.txt.gz` |
| `scripts/p2_summary.py`, `part2.json` | part 2's table (`logs/p2pl_*.log` plain, `logs/p2cp_*.log` CP phase, `logs/p2u_*.log` and `logs/p2long_*.log` longer CP runs) |
| `scripts/cost_table.py`, `scripts/tree_model.py`, `scripts/level_gain.py` | the derived tables of part 2 |
| `scripts/p3_summary.py` | part 3's tables (`logs/p3r_*.log` rounds form, `logs/p3cpb_*.log` CP form with boxes) |
| `scripts/runq.sh`, `scripts/runseq.sh`, `scripts/*.jobs` | the job runners and lists; each log's first line is its command (the lists were edited as the session went: the log is the record) |
| `logs/validate_cp.out`, `logs/validate_cp/` | brute-force check of the CP path on small instances |
| `cp_d4.log`, `rounds_d4.log` | the gate test of `../COMMON.md` (not mine) |

## Part 1: validity

**Brute force** [m] (`logs/validate_cp.out`): the repository's CP-path validator
(`experiments/wil50_feasibility/roadmap/admm/cpadmm/validate_cp.py --sep --count 1`, 47 s) on
41 nodes of n = 7, 8, 9 random instances (prefix regions) and nug12 with a 3 x 4 grid,
optima by enumeration: 205 CP solves, **10,508 certified values, max(certified - optimum) =
-1e-10**; 3,547 of them from the arms that separate inside the solve (`-cpsep`, the form used
below), same maximum; no fathom with the cutoff at the optimum, 80 correct fathoms at
optimum - 1. (The script ends with "FAILURES (502)": every one is a per-sweep progress line of
self-test (e), which the script's "must end in ok" rule predates; all 123 x 5 self-test
verdict lines (a)-(e) end in `ok`, and no certified-value or fathom failure is listed.)
`tests/dnn_test.py` (`make dnncheck`) covers the plain bound only, so it was not rerun.

**The path to the best known solution** [m] (`logs/p1_*.log`).
`scripts/path.py` inverts `tho40.sln` for `tho40t.dat`, checks that it and its three grid
images evaluate to 240,516, and takes the image with facility 26 at location 2 (path_check's).
A node that fixes d assignments of that solution has optimum <= 240,516 (= 240,516 if the best
known value is optimal), so no certified value may pass 240,516. Each run is
`dnncuts -x SPEC -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 -cpsep 100:1000 -cpiters N -T 1 -v qaplib/tho40t.dat`
(cutoff far above, so nothing stops early): 300 plain iterations cold, then N CP iterations,
certified every 10. Orders: `warm` = the path_check warm dive's order
(26, 30, 2, 16, 35, 24, 34, 21, 33, 37, 39, 15), `cold` = the cold dive's (differs from depth 3
on), `index` = facilities 0, 1, 2, ... (not a branching order: other nodes of the same solution).
N = 2,000 at depths 8 to 12, where the bound comes within reach of the optimum and an error
would show; at depths 1 to 6, where the bound is thousands below, N = 500 (`warm`) and 300
(`index`), cut down from the planned 1,000 for time.

{{p1_summary.py}}

@@P1TEXT@@

## Part 2: what the CP phase gains, by depth and shortfall

**Sample** (`scripts/sample_nodes.py`, `nodes.tsv`): from the production trace of subtree
`26:10,8:0` (22,650 nodes, 644 branched), branched nodes at quantiles of the give-up shortfall:
five at depth 3, five at depth 4 plus the gate node of `../COMMON.md`, four each at depths 5 and
6; the depth-2 node itself, the depth-1 nodes `26:10` (its parent) and `26:2` (the best known
solution's), and the root.

**Runs.** Plain [m]: `dnntest -p experiments/estimates/tho40/tho40.par -x SPEC -u 240516 -i N -H 0 -T 1 -v`
(`logs/p2pl_*.log`), N = 1,500 for the first three, then 1,000 at depths 4 to 6, 1,200 at
depth 3, 1,000 at depths 1 and 2, 600 at the root (cut for time; the "last-300 gain" column
says how converged each is). CP [m]: the brief's command with the cutoff at 241,500 instead
of `-u 240516` and `-cpiters 800` (`logs/p2cp_*.log`); the fathom iteration is the first check
whose certified value passes 240,514.05. That is the same run: at the gate node the
`-u 240516` log (`cp_d4.log`) and the high-cutoff log pass the cutoff at the same CP iteration
(120) with the same value to all digits (240,514.98184145204); checks further from a run's own
cutoff read about 10 lower in the high-cutoff run (a cheaper certificate there), so the
high-cutoff numbers are if anything conservative. `-u 240516` confirmations:
@@P2U@@

Columns: shortfalls are 240,514.05 less the certified bound (negative = past the cutoff).
"production shortfall" is the trace's give-up value (200 iterations, warm-started from the
parent), with the node's subtree size in the production run. "plain@300" is cold, by dnntest
with tho40.par and by dnncuts' own first solve (its DNN defaults; the CP solve starts from that
state). "CP gain" = CP value at its last check less the plain end value; where the CP run
stopped at its 241,500 cutoff or was still rising, it is a lower bound.

{{p2_summary.py}}

@@P2TEXT@@

### What a fathom at the give-up node saves

{{cost_table.py}}

@@COSTTEXT@@

### A threshold model of the subtree [e]

{{tree_model.py 1.2}}

@@MODELTEXT@@

### How much of a level is it worth?

{{level_gain.py}}

@@LEVELTEXT@@

### The root bound with cuts, run long

`logs/p2cp_root.log` [m]: `dnncuts -x - -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 -cpsep 100:1000 -cpiters 20000 -T 1 -v qaplib/tho40t.dat`
under `timeout 2400` (40 minutes of wall on the shared machine). Plain at 300 iterations:
225,681.44. The CP solve passes the plain DNN value (226,522) at CP iteration 120 and reads

| CP iteration | 500 | 1,000 | 2,000 | 3,000 | 4,000 | 5,000 | 5,490 (killed by the timeout) |
|---|---|---|---|---|---|---|---|
| certified | 227,870.2 | 228,603.7 | 229,312.3 | 229,668.5 | 229,884.6 | 230,031.7 | **230,088.31** |

@@ROOTTEXT@@

## Part 3: the other families

@@P3INTRO@@

{{p3_summary.py}}

@@P3TEXT@@

## What it implies, and what I would do next

@@NEXT@@
