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

| run | depth | m | plain@300 | CP iters | CP certified (last) | gain over last 200 it | max certified, all checks | 240,516 - max | checks | |
|---|---|---|---|---|---|---|---|---|---|---|
| p1_cold_d08 | 8 | 32 | 238534.7 | 2000 | 240399.7587 | +19.6 | 240399.7587 | 116.2413 | 200 | PASS |
| p1_cold_d10 | 10 | 30 | 239373.1 | 2000 | 240515.3993 | +0.7 | 240515.3993 | 0.6007 | 200 | PASS |
| p1_index_d08 | 8 | 32 | 239300.2 | 2000 | 240514.3245 | +1.8 | 240514.3245 | 1.6755 | 200 | PASS |
| p1_index_d10 | 10 | 30 | 240237.6 | 2000 | 240516.0000 | +0.0 | 240516.0000 | 0.0000 | 200 | PASS |
| p1_warm_d08 | 8 | 32 | 238458.3 | 2000 | 240411.1237 | +18.8 | 240411.1237 | 104.8763 | 200 | PASS |
| p1_warm_d09 | 9 | 31 | 238604.5 | 2000 | 240440.5206 | +16.5 | 240440.5206 | 75.4794 | 200 | PASS |
| p1_warm_d10 | 10 | 30 | 239311.2 | 2000 | 240514.8094 | +0.8 | 240514.8094 | 1.1906 | 200 | PASS |
| p1_warm_d11 | 11 | 29 | 239802.2 | 2000 | 240516.0000 | -0.0 | 240516.0000 | 0.0000 | 200 | PASS |
| p1_warm_d12 | 12 | 28 | 239971.3 | 2000 | 240516.0000 | -0.0 | 240516.0000 | 0.0000 | 200 | PASS |

worst: max certified - 240516 = 0.000000 over 1800 certified values in 9 runs

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

| node | depth | m | production shortfall (subtree nodes) | plain@300 short: dnntest par / dnncuts | plain end short (its, last-300 gain) | CP: fathom at | CP shortfall at CP it 100/200/400/800 | CP end short (gain/200 its) | CP gain vs plain end | share of plain end shortfall | -u 240516 run |
|---|---|---|---|---|---|---|---|---|---|---|---|
| root `-` | 0 | 40 | - (-) | - / 14833 | - | no | 14057/13530/12871/12151 | 10426 (@5490, +22) | - | - | - |
| d1_path `26:2` | 1 | 39 | - (-) | - / - | - | - | - | - | - | - | - |
| d1_sub `26:10` | 1 | 39 | - (-) | 13457 / 12149 | 11553 (running 620, +1512) | no | 11349/10912/10340/- | 10136 (@500, +454) | +1417 | 0.12 | - |
| d2_sub `26:10,8:0` | 2 | 38 | 8634.2 (22650) | 7079 / 6769 | 6494 (iter_limit 1000, +22) | no | 6104/5719/5217/4653 | 4653 (@800, +232) | +1840 | 0.28 | - |
| d3_a `26:10,8:0,31:32` | 3 | 37 | 73.7 (38) | 287 / 251 | 13 (iter_limit 1200, +4) | it 70 (300+70) | -427/-894/-/- | -989 (@230, +2533) | +1002 | 79.15 | - |
| d3_b `26:10,8:0,28:32` | 3 | 37 | 1906.4 (110) | 1942 / 1890 | 1648 (iter_limit 1200, +5) | it 650 (300+650) | 1238/808/331/-140 | -140 (@800, +189) | +1788 | 1.08 | - |
| d3_c `26:10,8:0,38:32` | 3 | 37 | 2426.9 (287) | 2593 / 2546 | 2168 (iter_limit 1200, +6) | no | 1911/1386/951/509 | 509 (@800, +176) | +1659 | 0.77 | - |
| d3_d `26:10,8:0,3:32` | 3 | 37 | 3679.3 (781) | 3759 / 3698 | 3383 (iter_limit 1200, +5) | no | 2978/2501/2004/1505 | 1505 (@800, +197) | +1879 | 0.56 | - |
| d3_e `26:10,8:0,2:32` | 3 | 37 | 5039.3 (6281) | 5147 / 4940 | 4667 (iter_limit 1200, +8) | no | 4308/3947/3485/2966 | 2966 (@800, +214) | +1701 | 0.36 | - |
| d4_a `26:10,8:0,4:32,39:7` | 4 | 36 | 89.5 (37) | 295 / 305 | 50 (iter_limit 1000, +8) | it 80 (300+80) | -259/-725/-/- | -987 (@306, +728) | +1038 | 20.56 | - |
| d4_b `26:10,8:0,24:32,16:39` | 4 | 36 | 929.5 (37) | 1112 / 1070 | 820 (iter_limit 1000, +10) | it 280 (300+280) | 535/168/-195/-557 | -557 (@800, +146) | +1377 | 1.68 | - |
| d4_c `26:10,8:0,18:32,24:7` | 4 | 36 | 1628.3 (142) | 1812 / 1834 | 1531 (iter_limit 1000, +8) | it 680 (300+680) | 1233/743/321/-91 | -91 (@800, +160) | +1622 | 1.06 | - |
| d4_d `26:10,8:0,39:32,30:7` | 4 | 36 | 2196.2 (212) | 2398 / 2276 | 2050 (iter_limit 1000, +15) | no | 1683/1342/930/505 | 505 (@800, +170) | +1545 | 0.75 | - |
| d4_e `26:10,8:0,2:32,6:39` | 4 | 36 | 3554.5 (1665) | 3844 / 3668 | 3417 (iter_limit 1000, +17) | no | 3076/2749/2339/1873 | 1873 (@800, +193) | +1544 | 0.45 | - |
| d4_gate `26:10,8:0,39:32,38:7` | 4 | 36 | 407.2 (37) | 807 / 786 | 339 (iter_limit 1500, +2) | it 120 (300+120) | 171/-380/-807/- | -988 (@535, +301) | +1327 | 3.91 | - |
| d5_a `26:10,8:0,11:32,6:39,18:7` | 5 | 35 | 101.3 (36) | 332 / 393 | 81 (iter_limit 1000, +7) | it 90 (300+90) | -169/-699/-/- | -987 (@315, +733) | +1068 | 13.18 | - |
| d5_b `26:10,8:0,2:32,16:39,18:24` | 5 | 35 | 364.5 (36) | 654 / 687 | 327 (iter_limit 1000, +12) | it 130 (300+130) | 192/-337/-717/- | -986 (@650, +198) | +1313 | 4.01 | - |
| d5_c `26:10,8:0,39:32,18:7,6:39` | 5 | 35 | 941.6 (36) | 1241 / 1193 | 881 (iter_limit 1000, +10) | it 240 (300+240) | 579/102/-346/-812 | -812 (@800, +190) | +1692 | 1.92 | - |
| d5_d `26:10,8:0,2:32,6:39,33:24` | 5 | 35 | 2359.0 (371) | 2527 / 2519 | 2288 (iter_limit 1000, +8) | no | 1965/1561/1162/734 | 734 (@800, +173) | +1555 | 0.68 | - |
| d6_a `26:10,8:0,18:32,2:7,6:39,33:15` | 6 | 34 | 121.1 (35) | 461 / 542 | 100 (iter_limit 1000, +13) | it 110 (300+110) | 43/-602/-/- | -987 (@332, +760) | +1087 | 10.85 | - |
| d6_b `26:10,8:0,35:32,17:39,2:8,24:24` | 6 | 34 | 342.1 (35) | 518 / 591 | 322 (iter_limit 1000, +6) | it 120 (300+120) | 90/-404/-795/- | -986 (@564, +252) | +1308 | 4.06 | - |
| d6_c `26:10,8:0,2:32,18:39,35:24,7:8` | 6 | 34 | 655.5 (68) | 806 / 793 | 599 (iter_limit 1500, +1) | it 190 (300+190) | 309/-34/-400/-797 | -797 (@800, +158) | +1396 | 2.33 | - |
| d6_d `26:10,8:0,2:32,6:39,33:24,39:8` | 6 | 34 | 1537.8 (134) | 1709 / 1715 | 1495 (iter_limit 1500, +1) | it 740 (300+740) | 1182/753/361/-44 | -44 (@800, +161) | +1538 | 1.03 | - |

@@P2TEXT@@

### What a fathom at the give-up node saves

| node | depth | production shortfall | subtree nodes | CPU below the node (s) | = DNN iterations at depth+1 | CP iterations to fathom | x 1.2 | iterations saved | ratio |
|---|---|---|---|---|---|---|---|---|---|
| d3_a | 3 | 73.7 | 38 | 83 | 541 | 70 | 84 | 457 | 6.4 |
| d3_b | 3 | 1906.4 | 110 | 395 | 2564 | 650 | 780 | 1784 | 3.3 |
| d4_a | 4 | 89.5 | 37 | 70 | 503 | 80 | 96 | 407 | 5.2 |
| d4_b | 4 | 929.5 | 37 | 103 | 747 | 280 | 336 | 411 | 2.2 |
| d4_c | 4 | 1628.3 | 142 | 451 | 3262 | 680 | 816 | 2446 | 4.0 |
| d4_gate | 4 | 407.2 | 37 | 50 | 364 | 120 | 144 | 220 | 2.5 |
| d5_a | 5 | 101.3 | 36 | 70 | 596 | 90 | 108 | 488 | 5.5 |
| d5_b | 5 | 364.5 | 36 | 76 | 645 | 130 | 156 | 489 | 4.1 |
| d5_c | 5 | 941.6 | 36 | 81 | 691 | 240 | 288 | 403 | 2.4 |
| d6_a | 6 | 121.1 | 35 | 55 | 565 | 110 | 132 | 433 | 4.3 |
| d6_b | 6 | 342.1 | 35 | 42 | 436 | 120 | 144 | 292 | 3.0 |
| d6_c | 6 | 655.5 | 68 | 116 | 1194 | 190 | 228 | 966 | 5.2 |
| d6_d | 6 | 1537.8 | 134 | 304 | 3135 | 740 | 888 | 2247 | 3.5 |

@@COSTTEXT@@

### A threshold model of the subtree [e]

sample (production shortfall -> CP iterations, running maximum): 74:70 90:80 101:90 121:110 342:120 364:130 407:130 656:190 930:280 942:280 1538:740 1628:740 1906:740
seconds per DNN iteration by depth (production): 2:0.122 3:0.172 4:0.154 5:0.138 6:0.118 7:0.097 8:0.078
production: 22650 nodes, 78801 CPU s; CP iteration = 1.2 plain iterations
| threshold T | give-up nodes fathomed by CP | nodes left | node ratio | CPU left | CP CPU added | CPU ratio |
|---|---|---|---|---|---|---|
| 250 | 193 | 15904 | 0.702 | 66728 | 3018 | 0.885 |
| 500 | 314 | 11408 | 0.504 | 57743 | 5562 | 0.803 |
| 1000 | 419 | 6417 | 0.283 | 44410 | 10252 | 0.694 |
| 1500 | 409 | 3686 | 0.163 | 33186 | 15561 | 0.619 |
| 1906 | 374 | 2359 | 0.104 | 25778 | 19011 | 0.568 |

@@MODELTEXT@@

### How much of a level is it worth?

| parent depth | children | of them branched (uncensored) | median gain, all children (censored ones at their stop) | median gain of the branched children | quartiles of the branched children's gain | median parent shortfall |
|---|---|---|---|---|---|---|
| 2 | 38 | 31 | 6443 | 6207 | 5027 / 6801 | 8634 |
| 3 | 1147 | 217 | 2770 | 2734 | 2077 / 3469 | 2427 |
| 4 | 7812 | 256 | 1684 | 1576 | 1250 / 1990 | 930 |
| 5 | 8960 | 105 | 1153 | 1081 | 858 / 1372 | 364 |
| 6 | 3570 | 34 | 1123 | 650 | 500 / 777 | 342 |
| 7 | 1122 | 0 | 825 | - | - | 164 |

parents at depth 2 short by > 2,000: 38 children, 18% fathomed, median gain 6443 (fathomed ones censored)
parents at depth 3 short by > 2,000: 777 children, 73% fathomed, median gain 3403 (fathomed ones censored)
parents at depth 4 short by > 2,000: 1116 children, 86% fathomed, median gain 2507 (fathomed ones censored)
parents at depth 5 short by > 2,000: 210 children, 86% fathomed, median gain 2304 (fathomed ones censored)

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


CP form: triangles (p2cp) against triangles + round 0's 200 boxes (p3cpb), certified shortfall to 240,514.05:

| node | CP it 100 | 200 | 400 | 600 | 800 | pool at the end | at the give-up point: violated (max) | at the CP run's last iterate: violated (max) |
|---|---|---|---|---|---|---|---|---|

@@P3TEXT@@

## What it implies, and what I would do next

@@NEXT@@
