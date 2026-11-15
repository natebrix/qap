# tho40 bound lab: the grid cut phase per node (Nov 15 2026)

Agent `tho40_nodes` of `experiments/bound_lab/` (see `../COMMON.md`). Question: is the
regional-triangle cut phase (CP-ADMM form, `latest/dnncut.h`) valid on tho40, and how much
does it add to the certified DNN bound, node by node? Tools as built in the main checkout:
`experiments/dnn_cuts/dnncuts`, `latest/dnntest`. Nothing under `latest/` was changed.

Marks: **[m]** measured here (the log is named), **[e]** estimated. The machine was shared
(load 5 to 6 on 4 cores), so iterations are the metric; CPU seconds are indicative only.

STATUS: work in progress; sections are filled in as the runs finish.

## Files

| | |
|---|---|
| `scripts/path.py` | the best known solution (`QAPSOLN/tho40.sln`, inverted for `tho40t.dat`), its four grid images (all 240,516), node specs along it |
| `scripts/part1_jobs.py`, `scripts/p1_summary.py` | part 1's runs and table (`logs/p1_*.log`) |
| `scripts/sample_nodes.py`, `nodes.tsv` | part 2's sample from `infra/phase1/runs/tho40ab-32/nodes.txt.gz` |
| `scripts/p2_summary.py`, `part2.json` | part 2's table (`logs/p2pl_*.log` plain, `logs/p2cp_*.log` CP phase, `logs/p2u_*.log` CP with `-u 240516`) |
| `scripts/runq.sh`, `scripts/runseq.sh`, `scripts/*.jobs` | the job runners and lists; each log's first line is its command |
| `logs/validate_cp.out`, `logs/validate_cp/` | brute-force check of the CP path on small instances |
| `cp_d4.log`, `rounds_d4.log` | the gate test of `../COMMON.md` (not mine) |

## Part 1: validity. PASS (deep nodes done; depths 1-6 pending)

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

**The path to the best known solution** [m] (`logs/p1_*.log`, table by `scripts/p1_summary.py`).
`scripts/path.py` inverts `tho40.sln` for `tho40t.dat`, checks that it and its three grid
images evaluate to 240,516, and takes the image with facility 26 at location 2 (path_check's).
A node that fixes d assignments of that solution has optimum <= 240,516 (= 240,516 if the best
known value is optimal), so no certified value may pass 240,516. Each run is
`dnncuts -x SPEC -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 -cpsep 100:1000 -cpiters N -T 1 -v`
(cutoff far above, so nothing stops early): 300 plain iterations cold, then N CP iterations,
certified every 10. Orders: `warm` = the path_check warm dive's order
(26, 30, 2, 16, 35, 24, 34, 21, 33, 37, 39, 15), `cold` = the cold dive's (differs from depth 3),
`index` = facilities 0, 1, 2, ... (not a branching order: different nodes of the same solution).

(Interim table, 9 deep runs; the final one replaces it when depths 1-6 finish: `p1_table.md`.)
Three nodes (warm d11, warm d12, index d10) converge to 240,515.99999 and stay below 240,516
for 700 to 1,100 further iterations; four more end 0.6 to 1.7 below. No certified value above
240,516 in 1,800 checks.

## Root bound with cuts, run long (interim; part 2 below when its runs finish)

`logs/p2cp_root.log` [m]: `dnncuts -x - -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 -cpsep 100:1000 -cpiters 20000 -T 1 -v qaplib/tho40t.dat`
under `timeout 2400` (40 minutes of wall on the shared machine). Plain at 300 iterations:
225,681.44. The CP solve passes the plain DNN value (226,522) at CP iteration 120 and reads

| CP iteration | 500 | 1,000 | 2,000 | 3,000 | 4,000 | 5,000 | 5,490 (killed by the timeout) |
|---|---|---|---|---|---|---|---|
| certified | 227,870.2 | 228,603.7 | 229,312.3 | 229,668.5 | 229,884.6 | 230,031.7 | **230,088.31** |

3,568 triangles in the pool at the end, 1,260 CPU s for the CP solve (0.23 s per iteration at m = 40).
Still rising (+57 over the last 490 iterations, the increments shrinking by about a quarter per 1,000).
