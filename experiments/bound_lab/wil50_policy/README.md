# wil50 / sko49: tho40's phase policy on complete subtrees (bound_lab `wil50_policy`, Nov 16 2026)

Question (the coordinator's brief): do the three changes to how the cut phase is *used* that halved
tho40's warm subtrees (`../tho40_search/README.md`, "The stop rule": a wide zone, the hand-over rule
`cut.handover`, early entry) make wil50's and sko49's subtrees cheaper? **[m]** measured here (the run
under `runs/` is named), **[e]** estimated. Iterations are the metric; CPU is this box's (4 cores, four
solver processes at a time, `lscpu.txt`).

## Answer

**A modest gain, not tho40's 2x: about 1.1-1.35x on wil50 and 1.2-1.35x on sko49, from two of the three
changes. The wide zone changes nothing, and the hand-over rule has to be refitted to be *less* patient
far from the cutoff than tho40's, or it loses most of the gain.**

| | subtree (cold `-F`) | base | recommended arm | factor in iterations (CPU) | the tree |
|---|---|---|---|---|---|
| wil50 | a220 `21:22,31:23,40:15` | 232 nodes, 17,191 it | `hoe2`: 12,699 it | **1.35x** (1.36x) [m] | the same 232 nodes, node for node |
| wil50 | b83 `8:9,21:11,45:40,34:0` (a depth-4 child of a251) | 181 nodes, 13,296 it | `hoe2`: 12,230 it | **1.09x** (1.07x) [m] | 181 nodes in both |
| wil50 | a251, a325 (laptop traces, not rerun: 35 and 100 CPU-hours an arm here) | 387,735 and 1,066,452 it | | **1.17x and 1.16x** [e, `extrap.py`] | |
| sko49 | s188 `20:6,40:4,41:16` | 47 nodes, 4,162 it | `hoe2w`: 3,446 it | **1.21x** (1.23x) [m] | the same 47 nodes |
| sko49 | s232 `23:0,40:8,41:10` | 271 nodes, 34,557 it | `hoe2w`: 25,763 it | **1.34x** (1.37x) [m]; 1.56x on the 271 nodes both arms bounded | 404 nodes: 133 more, cheap ones (3,627 it) |

- **wil50: about 1.15-1.2x [e], measured 1.09x and 1.35x.** Both measured factors are inside the 1.4x
  that single subtrees vary by, but they are not noise of that kind: the trees are the same and the
  whole difference sits on the few nodes that run a phase (a220: 7 nodes, 10,649 -> 6,111 it, 1.74x;
  the 225 nodes that fathom in their first solve: 6,542 against 6,588). What limits it is the
  composition: on the large subtrees 46-48% of the iterations are first solves that fathom, which no
  phase policy touches. a220 overstates the gain (44% of its iterations are missed phases, one of
  them a 3,400-iteration phase that ended 0.1 short); a251 and a325 have 25%.
- **sko49: about 1.2-1.35x [m, two subtrees]**, with `wc=1.4` in the arm. Attribution: on s188 the
  arm without `wc` (`hoe2`) is 1.04x and `wc=1.4` alone 1.05x, the two together 1.21x; on s232
  `wc=1.4` alone is 1.08x (1.17x on matched nodes) and `hoe2` alone was stopped at the session's budget after 73 of its nodes: on the 73 nodes it
  shares with the others, base took 10,623 iterations, `wc14` 9,166, `hoe2` 8,769 (1.21x) and
  `hoe2w` 7,726 (1.37x) [m, partial: `runs/s232_hoe2`, incomplete]. So on sko49 each half gives
  about 1.05-1.2x and they add.
- **Recommended blocks** (on top of each instance's current par file; `cut.handover` needs
  `../tho40_search/stop/handover.patch`, which this branch has applied in `latest/`):

  ```
  # wil50 (arm hoe2)                      # sko49 (arm hoe2w)
  dnn.giveup_warm = 100                   dnn.giveup_warm = 100
  cut.handover = 32:0:25                  cut.handover = 27:0:5
                                          dnn.beta_sched = two,f=1.6667,t=60,w=0.25,wc=1.4
  ```
  Not `wc=1.4` on wil50: it made both wil50 subtrees worse than `hoe2` (below).
- **What did not carry over from tho40.**
  1. *The wide zone is already there.* Under the base zones (0.0075: 366 short; 0.018: 421 short)
     every unfathomed node of every subtree here, and of the laptop's a220, a251 and a325, already
     runs a phase (no "no phase" node in any trace). The zone arm is the base arm, bit for bit, at
     depth 3 and below; it could only differ at depths 1-2, outside these subtrees. Not run.
  2. *Patience far from the cutoff is not repaid.* With constants in tho40's spirit (`32:0:250`, arm
     `hoe`) a220's root phase ran 2,660 iterations instead of 980 for 14.7 units more (105.3 -> 90.6
     short), and its children's first solves cost the same (5,405 against 5,314 on the 179 matched
     ones): 1,680 iterations lost, 1.19x instead of 1.35x. The phase's tail is flat here (1 unit per
     100 iterations after 1,500), where tho40's far phases were still closing 20 units per 100.
  3. *The gain is near the cutoff*: phases that the stall rule lets run 1,600-3,400 iterations to end
     0.1-6 short, whose 46 children then fathom in 25-30 iterations each (1,150-1,250 in all). The
     rule stops them at 1,000-1,150.
  4. *Early entry cuts both ways.* It saves 190 iterations at each give-up (first solve 110
     instead of 300), but a warm node that would have fathomed in its first solve between iterations
     110 and 300 now gives up and runs a phase: on b83, 7 such nodes cost 2,701 iterations where
     they cost 2,000 (and the phases fathom them). That is why b83 is 1.09x and not more.
- **`wc` (the coordinator's added arm).** On the far roots `wc=1.4` reaches the base phase's end value
  1.27x sooner on wil50 and 1.44x sooner on sko49; 2.8 is no better and 5.6 worse (table below). In
  the search, alone: sko49 s188 1.05x, s232 1.08x with a larger tree (359 nodes: two phases that
  fathomed after 2,050-2,460 iterations in base stalled 0.7 and 3.5 short); 1.17x on the matched
  nodes. On wil50 on top of `hoe2` it loses: a220 13,451 against 12,699 (277 nodes against 232), b83
  13,645 against 12,230 (268 nodes against 181; 0.97x of base). Near
  misses fathom less often at the higher penalty, as cpfast found on tho40; `cut.cp_wc_gap` (cpfast's
  patch, not in this build) is the thing to try before using `wc` on wil50.
- **Sanity: PASS.** No run printed `New objective value` (wil50 48,816, sko49 23,386 stand) and every
  finished run ends with its subtree fathomed (columns "new obj" and "complete" below). The unpatched
  keys reproduce the documented trees: a220 232 DNN nodes (1, 47, 184), s188 47 nodes and 4,162
  iterations to the digit, s232 271 nodes (1, 46, 180, 44).

## How it was run

- Solver: `latest/bnbdrv` of this branch = `claude/cool-babbage-ilnogj` + `handover.patch` (committed
  here; unset, the key changes nothing; the settings goldens differ by the one new unset key and are
  not regenerated). `pip install scipy-openblas32 scipy`, `make -C latest -j4`.
- `run.sh INST NAME ARM SPEC`: one `bnbdrv -T 1`-style process (one thread) per subtree and arm, cold
  from `-F`, with the instance's symmetry file and incumbent, `QAP_DNN_TRACE` + `QAP_DNN_TRACE_HIST` +
  `QAP_NODE_TRACE` and the full run's `QAP_DNN_POOL_CHECK=20`. `queue.sh` keeps four running from
  `jobs.txt` (`queue.log` has the order). Each `runs/NAME_ARM/` has `cmd.txt`, `log.txt`, `err.txt`,
  `dnn.txt` (one line per DNN node, with each phase's history), `nodes.txt`.
- Arms (`par/INST_ARM.par` = the base file plus the lines at its end):

  | arm | wil50 | sko49 |
  |---|---|---|
  | `base` | `experiments/wil50_feasibility/fullrun/fullrun.par` | `experiments/estimates/sko49/sko49.par` (Cw) |
  | zone | = base (see above) | = base |
  | early | `dnn.giveup_warm = 100` (both instances already have `dnn.gu_zone = 0`): inside `hoe*`, not run alone | the same |
  | `hoe` | early + `cut.handover = 32:0:250` (the first fit, patient far away) | early + `27:0:16` |
  | `hoe2` | early + `cut.handover = 32:0:25` (refit after a220 `hoe`) | early + `27:0:5` |
  | `hoe2w` | `hoe2` + `wc=1.4` | `hoe2` + `wc=1.4` |
  | `wc14`, `wc28`, `wc56` | base with `wc` = 1.4, 2.8, 5.6 in `dnn.beta_sched` | the same |

- **The stall rule "on and off" with the hand-over rule:** in the patch the hand-over rule *replaces*
  the stall rule (`CutProjector::stop()`), so every `ho*` arm is "stall off"; the primal gate
  (`cut.pgate`) is a rounds-form key that the CP form ignores, and neither instance sets `cut.giveup`.
  The base arm is "stall on, no hand-over". A both-rules arm would need a code change; not done.
- Warm roots (the parent's state, as the tho40 study): not built. `experiments/estimates/README.md`
  ("Cold against warm") measured a220 warm and pooled at the same 232 nodes and 4% fewer iterations.
- Subtrees. a220 and sko49's two are the brief's. a251 (387,735 iterations on the laptop: about 35
  CPU-hours an arm on this box, against 17 core-hours in all) could not be run; **b83**, a depth-4
  child of a251 (`21:11,34:0,45:40,8:9`: 82.6 short and 180 nodes below it in the laptop's a251), was
  run cold instead as the second wil50 subtree.

## Every run [m]

`python3 summarize.py`; "stalled" is "handed over" in the `ho*` arms; `a220r`, `s232r` are the
single-phase `wc` probes (`-d 3`). x base is base / arm. `s232_hoe2` reads RUNNING: it was stopped by
hand at 07:00 (its `cmd.txt`), incomplete, and is used only for the 73-node comparison in the answer.

| subtree | arm | DNN nodes (by depth) | first solves: it (give-ups) | phases: n = fathomed + handed over/stalled + limit | phase it: fathomed + stalled + limit | total it | x base | CPU s | x base | new obj | complete |
|---|---|---|---|---|---|---|---|---|---|---|---|
| a220 | base | 232 (3: 1, 4: 47, 5: 184) | 8,692 (7) | 7 = 2 + 4 + 1 | 629 + 4,470 + 3,400 = 8,499 | 17,191 | 1.00x | 5512 | 1.00x | 0 | yes |
| a220 | hoe | 232 (3: 1, 4: 47, 5: 184) | 7,507 (7) | 7 = 2 + 5 + 0 | 630 + 6,350 + 0 = 6,980 | 14,487 | 1.19x | 4795 | 1.15x | 0 | yes |
| a220 | hoe2 | 232 (3: 1, 4: 47, 5: 184) | 7,548 (7) | 7 = 2 + 5 + 0 | 761 + 4,390 + 0 = 5,151 | 12,699 | 1.35x | 4055 | 1.36x | 0 | yes |
| a220 | hoe2w | 277 (3: 1, 4: 47, 5: 184, 6: 45) | 8,619 (7) | 7 = 1 + 6 + 0 | 312 + 4,520 + 0 = 4,832 | 13,451 | 1.28x | 4294 | 1.28x | 0 | yes |
| a220r | wc14 | 1 (3: 1) | 300 (1) | 1 = 0 + 1 + 0 | 0 + 1,010 + 0 = 1,010 | 1,310 | - | 461 | - | 0 | yes |
| a220r | wc28 | 1 (3: 1) | 300 (1) | 1 = 0 + 1 + 0 | 0 + 1,010 + 0 = 1,010 | 1,310 | - | 468 | - | 0 | yes |
| b83 | base | 181 (4: 1, 5: 46, 6: 90, 7: 44) | 8,826 (4) | 4 = 0 + 4 + 0 | 0 + 4,470 + 0 = 4,470 | 13,296 | 1.00x | 3436 | 1.00x | 0 | yes |
| b83 | hoe2 | 181 (4: 1, 5: 46, 6: 90, 7: 44) | 7,019 (11) | 11 = 7 + 4 + 0 | 1,931 + 3,280 + 0 = 5,211 | 12,230 | 1.09x | 3207 | 1.07x | 0 | yes |
| b83 | hoe2w | 268 (4: 1, 5: 46, 6: 90, 7: 88, 8: 43) | 9,420 (12) | 12 = 6 + 6 + 0 | 2,135 + 2,090 + 0 = 4,225 | 13,645 | 0.97x | 3352 | 1.03x | 0 | yes |
| s188 | base | 47 (3: 1, 4: 46) | 2,405 (2) | 2 = 1 + 1 + 0 | 397 + 1,360 + 0 = 1,757 | 4,162 | 1.00x | 1303 | 1.00x | 0 | yes |
| s188 | hoe | 47 (3: 1, 4: 46) | 2,168 (2) | 2 = 1 + 1 + 0 | 350 + 1,610 + 0 = 1,960 | 4,128 | 1.01x | 1254 | 1.04x | 0 | yes |
| s188 | hoe2 | 47 (3: 1, 4: 46) | 2,335 (2) | 2 = 1 + 1 + 0 | 566 + 1,100 + 0 = 1,666 | 4,001 | 1.04x | 1234 | 1.06x | 0 | yes |
| s188 | hoe2w | 47 (3: 1, 4: 46) | 2,116 (2) | 2 = 1 + 1 + 0 | 310 + 1,020 + 0 = 1,330 | 3,446 | 1.21x | 1060 | 1.23x | 0 | yes |
| s188 | wc14 | 47 (3: 1, 4: 46) | 2,314 (2) | 2 = 1 + 1 + 0 | 233 + 1,400 + 0 = 1,633 | 3,947 | 1.05x | 1163 | 1.12x | 0 | yes |
| s232 | base | 271 (3: 1, 4: 46, 5: 180, 6: 44) | 16,908 (20) | 20 = 14 + 6 + 0 | 10,939 + 6,710 + 0 = 17,649 | 34,557 | 1.00x | 9197 | 1.00x | 0 | yes |
| s232 | hoe2 | 73 (3: 1, 4: 40, 5: 32) | 3,910 (9) | 9 = 7 + 2 + 0 | 3,229 + 1,630 + 0 = 4,859 | 8,769 | - | 2559 | - | 0 | RUNNING |
| s232 | hoe2w | 404 (3: 1, 4: 46, 5: 225, 6: 132) | 15,229 (21) | 21 = 12 + 9 + 0 | 4,454 + 6,080 + 0 = 10,534 | 25,763 | 1.34x | 6723 | 1.37x | 0 | yes |
| s232 | wc14 | 359 (3: 1, 4: 46, 5: 180, 6: 132) | 17,743 (17) | 17 = 9 + 8 + 0 | 4,513 + 9,600 + 0 = 14,113 | 31,856 | 1.08x | 8278 | 1.11x | 0 | yes |
| s232r | wc14 | 1 (3: 1) | 300 (1) | 1 = 0 + 1 + 0 | 0 + 950 + 0 = 950 | 1,250 | - | 388 | - | 0 | yes |
| s232r | wc28 | 1 (3: 1) | 300 (1) | 1 = 0 + 1 + 0 | 0 + 990 + 0 = 990 | 1,290 | - | 418 | - | 0 | yes |
| s232r | wc56 | 1 (3: 1) | 300 (1) | 1 = 0 + 1 + 0 | 0 + 990 + 0 = 990 | 1,290 | - | 384 | - | 0 | yes |

Matched nodes (`python3 summarize.py --matched NAME ARM1 ARM2`: the nodes both arms bounded, by how
they ended; where the trees have the same node count but fewer "both", the branching facility
differed under one node, whose 45-46 children then carry other names):

```
a220: base 232 nodes, hoe 232 nodes, both 186
  matched nodes: base 15963 it, hoe 13345 it (1.20x)
    fathomed in the first solve in both     179 nodes:   5314 it against   5405 (0.98x)
    phase in both                             7 nodes:  10649 it against   7940 (1.34x)
  only base: 46 nodes 1228 it; only hoe: 46 nodes 1142 it
a220: base 232 nodes, hoe2 232 nodes, both 232
  matched nodes: base 17191 it, hoe2 12699 it (1.35x)
    fathomed in the first solve in both     225 nodes:   6542 it against   6588 (0.99x)
    phase in both                             7 nodes:  10649 it against   6111 (1.74x)
  only base: 0 nodes 0 it; only hoe2: 0 nodes 0 it
a220: hoe2 232 nodes, hoe2w 277 nodes, both 186
  matched nodes: hoe2 11460 it, hoe2w 11163 it (1.03x)
    fathomed in the first solve in both     179 nodes:   5349 it against   5371 (1.00x)
    phase in both                             7 nodes:   6111 it against   5792 (1.06x)
  only hoe2: 46 nodes 1239 it; only hoe2w: 91 nodes 2288 it
b83: base 181 nodes, hoe2 181 nodes, both 136
  matched nodes: base 11953 it, hoe2 10873 it (1.10x)
    fathomed in the first solve in both     125 nodes:   4243 it against   4252 (1.00x)
    phase in both                             4 nodes:   5710 it against   3920 (1.46x)
    other                                     7 nodes:   2000 it against   2701 (0.74x)
  only base: 45 nodes 1343 it; only hoe2: 45 nodes 1357 it
b83: hoe2 181 nodes, hoe2w 268 nodes, both 181
  matched nodes: hoe2 12230 it, hoe2w 11396 it (1.07x)
    fathomed in the first solve in both     169 nodes:   5527 it against   5651 (0.98x)
    phase in both                            11 nodes:   6621 it against   5485 (1.21x)
    other                                     1 nodes:     82 it against    260 (0.32x)
  only hoe2: 0 nodes 0 it; only hoe2w: 87 nodes 2249 it
s188: base 47 nodes, hoe2w 47 nodes, both 47
  matched nodes: base 4162 it, hoe2w 3446 it (1.21x)
    fathomed in the first solve in both      45 nodes:   1805 it against   1706 (1.06x)
    phase in both                             2 nodes:   2357 it against   1740 (1.35x)
  only base: 0 nodes 0 it; only hoe2w: 0 nodes 0 it
s232: base 271 nodes, wc14 359 nodes, both 271
  matched nodes: base 34557 it, wc14 29437 it (1.17x)
    fathomed in the first solve in both     251 nodes:  10348 it against   9307 (1.11x)
    phase in both                            17 nodes:  22276 it against  19473 (1.14x)
    other                                     3 nodes:   1933 it against    657 (2.94x)
  only base: 0 nodes 0 it; only wc14: 88 nodes 2419 it
s232: base 271 nodes, hoe2w 404 nodes, both 271
  matched nodes: base 34557 it, hoe2w 22136 it (1.56x)
    fathomed in the first solve in both     249 nodes:   9586 it against   8956 (1.07x)
    phase in both                            19 nodes:  23748 it against  12303 (1.93x)
    other                                     3 nodes:   1223 it against    877 (1.39x)
  only base: 0 nodes 0 it; only hoe2w: 133 nodes 3627 it
s232: hoe2 73 nodes, hoe2w 404 nodes, both 73
  matched nodes: hoe2 8769 it, hoe2w 7726 it (1.13x)
    fathomed in the first solve in both      64 nodes:   2730 it against   2589 (1.05x)
    phase in both                             7 nodes:   5579 it against   4811 (1.16x)
    other                                     2 nodes:    460 it against    326 (1.41x)
  only hoe2: 0 nodes 0 it; only hoe2w: 331 nodes 18037 it
s188: base 47 nodes, hoe2 47 nodes, both 47
  matched nodes: base 4162 it, hoe2 4001 it (1.04x)
    fathomed in the first solve in both      45 nodes:   1805 it against   1925 (0.94x)
    phase in both                             2 nodes:   2357 it against   2076 (1.14x)
  only base: 0 nodes 0 it; only hoe2: 0 nodes 0 it
```

## What the phases did [m]

Every node that ran a phase (`python3 summarize.py --phases RUN`): first solve (iterations, shortfall),
the phase (iterations, how it ended, gap left), and the DNN nodes and iterations below it.

a220 base:
```
depth pool | first: it short | phase: it status left | below: nodes it
3 0 |  300   254.0 |   980 stalled      105.3 |  231  15911 | 21:22 31:23 40:15
4 1 |  300    56.8 |   880 stalled       30.7 |   46   2139 | 21:22 31:23 40:15 45:40
5 1 |  330     8.1 |   367 fathomed      -0.0 |    0      0 | 0:0 21:22 31:23 40:15 45:40
4 1 |  300    32.6 |  1620 stalled        4.5 |   46   1228 | 21:22 31:23 40:15 45:30
4 1 |  320    14.1 |   262 fathomed      -0.1 |    0      0 | 21:22 31:23 40:15 45:20
4 1 |  300    30.1 |  3400 iter_limit     0.1 |   46   1144 | 21:22 31:23 40:15 45:10
4 1 |  300    49.7 |   990 stalled       21.9 |   46   1256 | 21:22 31:23 40:15 45:0
```
a220 `hoe` (patient far away: the root's 2,660 iterations):
```
depth pool | first: it short | phase: it status left | below: nodes it
3 0 |  300   254.0 |  2660 stalled       90.6 |  231  11527 | 21:22 31:23 40:15
4 1 |  110    68.8 |   810 stalled       28.5 |   46   1952 | 21:22 31:23 40:15 45:40
5 1 |  110    17.2 |   400 fathomed      -0.2 |    0      0 | 0:0 21:22 31:23 40:15 45:40
4 1 |  110    53.5 |  1020 stalled        3.5 |   46   1142 | 21:22 31:23 40:15 45:30
4 1 |  110    34.2 |   230 fathomed      -0.1 |    0      0 | 21:22 31:23 40:15 45:20
4 1 |  110    46.0 |  1100 stalled        1.4 |   46   1144 | 21:22 31:23 40:15 45:10
4 1 |  110    62.5 |   760 stalled       20.7 |   46   1265 | 21:22 31:23 40:15 45:0
```
a220 `hoe2` (recommended):
```
depth pool | first: it short | phase: it status left | below: nodes it
3 0 |  300   254.0 |  1120 stalled      102.6 |  231  11279 | 21:22 31:23 40:15
4 1 |  110    69.8 |   550 stalled       34.1 |   46   2068 | 21:22 31:23 40:15 45:40
5 1 |  110    20.8 |   491 fathomed      -0.0 |    0      0 | 0:0 21:22 31:23 40:15 45:40
4 1 |  110    54.1 |  1020 stalled        5.7 |   46   1239 | 21:22 31:23 40:15 45:30
4 1 |  110    35.7 |   270 fathomed      -0.0 |    0      0 | 21:22 31:23 40:15 45:20
4 1 |  110    48.0 |  1140 stalled        3.1 |   46   1125 | 21:22 31:23 40:15 45:10
4 1 |  110    63.5 |   560 stalled       25.4 |   46   1273 | 21:22 31:23 40:15 45:0
```
b83 base and `hoe2`:
```
depth pool | first: it short | phase: it status left | below: nodes it
4 0 |  310   199.7 |  1280 stalled       48.9 |  180  11706 | 21:11 34:0 45:40 8:9
5 1 |  330    24.0 |  1540 stalled        4.5 |   45   1343 | 13:10 21:11 34:0 45:40 8:9
5 1 |  300    32.1 |   600 stalled       24.8 |   89   5251 | 0:10 21:11 34:0 45:40 8:9
6 1 |  300    12.0 |  1050 stalled        6.3 |   44   1683 | 0:10 13:20 21:11 34:0 45:40 8:9
depth pool | first: it short | phase: it status left | below: nodes it
4 0 |  310   199.7 |  1010 stalled       53.5 |  180  10910 | 21:11 34:0 45:40 8:9
5 1 |  110    40.4 |   446 fathomed      -0.1 |    0      0 | 21:11 34:0 45:40 47:10 8:9
5 1 |  110    32.6 |   305 fathomed      -0.1 |    0      0 | 21:11 28:10 34:0 45:40 8:9
5 1 |  110    44.6 |  1140 stalled        6.9 |   45   1357 | 13:10 21:11 34:0 45:40 8:9
5 1 |  110    42.9 |   480 stalled       27.7 |   89   5272 | 0:10 21:11 34:0 45:40 8:9
6 1 |  110    22.0 |   220 fathomed      -0.8 |    0      0 | 0:10 21:11 34:0 45:40 47:20 8:9
6 1 |  110    11.2 |   120 fathomed      -1.0 |    0      0 | 0:10 21:11 29:20 34:0 45:40 8:9
6 1 |  110    21.8 |   400 fathomed      -0.1 |    0      0 | 0:10 21:11 28:20 34:0 45:40 8:9
6 1 |  110    24.0 |   650 stalled        9.4 |   44   1874 | 0:10 13:20 21:11 34:0 45:40 8:9
7 1 |  110    14.8 |   260 fathomed      -0.4 |    0      0 | 0:10 13:20 21:11 27:30 34:0 45:40 8:9
6 1 |  110    17.7 |   180 fathomed      -1.2 |    0      0 | 0:10 21:11 34:0 45:40 4:20 8:9
```
s232 base and `hoe2w`:
```
depth pool | first: it short | phase: it status left | below: nodes it
3 0 |  300   245.5 |   910 stalled      120.6 |  270  33347 | 23:0 40:8 41:10
4 1 |  320     4.5 |   141 fathomed      -0.1 |    0      0 | 23:0 40:8 41:10 46:15
4 1 |  310    36.2 |  1800 fathomed      -0.1 |    0      0 | 23:0 40:8 41:10 46:14
4 1 |  300    15.3 |   196 fathomed      -0.2 |    0      0 | 23:0 40:8 41:10 46:13
4 1 |  300    21.9 |   682 fathomed      -0.0 |    0      0 | 23:0 40:8 41:10 46:9
4 1 |  300    73.0 |   790 stalled       50.6 |   45   5786 | 23:0 40:8 41:10 46:7
5 1 |  330     8.3 |   314 fathomed      -0.0 |    0      0 | 13:1 23:0 40:8 41:10 46:7
5 1 |  320    18.0 |  2460 fathomed      -0.2 |    0      0 | 12:1 23:0 40:8 41:10 46:7
5 1 |  300     6.1 |   233 fathomed      -0.1 |    0      0 | 23:0 3:1 40:8 41:10 46:7
4 1 |  300    30.8 |   527 fathomed      -0.1 |    0      0 | 23:0 40:8 41:10 46:6
4 1 |  300    30.3 |   556 fathomed      -0.0 |    0      0 | 23:0 40:8 41:10 46:5
4 1 |  300    32.9 |   970 fathomed      -0.0 |    0      0 | 23:0 40:8 41:10 46:4
4 1 |  300    43.7 |  2320 stalled        6.0 |   45   1371 | 23:0 40:8 41:10 46:3
4 1 |  300    56.4 |   960 stalled       30.5 |   45   2030 | 23:0 40:8 41:10 46:2
5 1 |  330     2.9 |   239 fathomed      -0.0 |    0      0 | 12:1 23:0 40:8 41:10 46:2
4 1 |  300    80.4 |   600 stalled       60.6 |   89   9094 | 23:0 40:8 41:10 46:1
5 1 |  300    21.3 |  2050 fathomed      -0.0 |    0      0 | 23:0 35:2 40:8 41:10 46:1
5 1 |  310    14.2 |   608 fathomed      -0.1 |    0      0 | 13:2 23:0 40:8 41:10 46:1
5 1 |  300    29.9 |  1130 stalled       12.4 |   44   1450 | 12:2 23:0 40:8 41:10 46:1
5 1 |  740     0.3 |   163 fathomed      -0.1 |    0      0 | 0:2 23:0 40:8 41:10 46:1
depth pool | first: it short | phase: it status left | below: nodes it
3 0 |  300   245.5 |  1040 stalled      106.3 |  403  24423 | 23:0 40:8 41:10
4 1 |  110    47.0 |  1140 stalled        2.0 |   45   1189 | 23:0 40:8 41:10 46:14
4 1 |  110    18.1 |   186 fathomed      -0.4 |    0      0 | 23:0 40:8 41:10 46:13
4 1 |  110    21.3 |   364 fathomed      -0.1 |    0      0 | 23:0 40:8 41:10 46:9
4 1 |  110    71.2 |   510 stalled       50.0 |   89   4638 | 23:0 40:8 41:10 46:7
5 1 |  110     8.7 |   201 fathomed      -0.1 |    0      0 | 23:0 35:1 40:8 41:10 46:7
5 1 |  110    16.4 |   410 fathomed      -0.1 |    0      0 | 13:1 23:0 40:8 41:10 46:7
5 1 |  110    26.9 |   690 stalled        7.9 |   44   1212 | 12:1 23:0 40:8 41:10 46:7
5 1 |  110     9.9 |   278 fathomed      -0.1 |    0      0 | 23:0 3:1 40:8 41:10 46:7
4 1 |  110    32.8 |   369 fathomed      -0.0 |    0      0 | 23:0 40:8 41:10 46:6
4 1 |  110    40.3 |   390 fathomed      -0.1 |    0      0 | 23:0 40:8 41:10 46:5
4 1 |  110    38.4 |   610 fathomed      -0.2 |    0      0 | 23:0 40:8 41:10 46:4
4 1 |  110    51.9 |   730 stalled       11.2 |   45   1377 | 23:0 40:8 41:10 46:3
4 1 |  110    61.8 |   520 stalled       31.0 |   45   1851 | 23:0 40:8 41:10 46:2
5 1 |  110     9.0 |   323 fathomed      -0.0 |    0      0 | 12:1 23:0 40:8 41:10 46:2
4 1 |  110    80.3 |   480 stalled       57.6 |  133   7099 | 23:0 40:8 41:10 46:1
5 1 |  110    24.5 |   770 stalled        5.0 |   44   1226 | 23:0 35:2 40:8 41:10 46:1
5 1 |  110    15.0 |   310 fathomed      -0.6 |    0      0 | 23:0 25:2 40:8 41:10 46:1
5 1 |  110    27.7 |   820 fathomed      -1.2 |    0      0 | 13:2 23:0 40:8 41:10 46:1
5 1 |  110    34.1 |   200 stalled       31.0 |   44   1459 | 12:2 23:0 40:8 41:10 46:1
5 1 |  110    18.1 |   193 fathomed      -0.0 |    0      0 | 0:2 23:0 40:8 41:10 46:1
```

- Under `hoe2` no phase runs past 1,140 iterations on wil50; base's ran 1,540, 1,620 and 3,400 to end
  4.5, 4.5 and 0.1 short. The children of the three cost 1,125-1,357 either way.
- On s232 the rule also stops three phases that in base *fathomed* after 1,800, 2,050 and 2,460
  iterations (at 2.0, 5.0 and 7.9 short, after 1,140, 770 and 690): three lost fathoms, each replaced
  by 44-45 children that cost 1,190-1,230 iterations. That is the 133 extra nodes, and it is a gain.
- Phases that fathom quickly are unchanged or shorter (they start from a first solve stopped at 110).

## The `wc` sweep on single phases [m]

`python3 wcsweep.py`: the root of a220 and of s232, cold, first solve 300 iterations, then the phase
under the stall rule; gap to the cutoff by phase iteration (the trace's history), and the iteration at
which each `wc` first reaches the gaps the base phase had at 300, 600 and at its end. Far misses only:
a pooled near miss needs its parent's state, so the near side was measured in the search arms instead
(`wc14`, `hoe2w` above).


wil50 a220 root (254.0 short after the first solve)
| wc | gap at 100 | 200 | 300 | 500 | 700 | 900 | it to base's gap at 300 (152.7) | at 600 (120.2) | at its end, 980 (105.3) |
|---|---|---|---|---|---|---|---|---|---|
| 0.7 | 206.8 | 174.1 | 152.7 | 127.8 | 115.2 | 107.5 | 300 (1.00x) | 600 (1.00x) | 980 (1.00x) |
| 1.4 | 211.4 | 173.4 | 149.3 | 121.7 | 108.5 | 100.6 | 300 (1.00x) | 510 (1.18x) | 770 (1.27x) |
| 2.8 | 228.9 | 184.6 | 156.9 | 126.4 | 112.4 | 104.1 | 320 (0.94x) | 570 (1.05x) | 860 (1.14x) |

sko49 s232 root (245.5 short)
| wc | gap at 100 | 200 | 300 | 500 | 700 | 900 | it to base's gap at 300 (161.7) | at 600 (133.9) | at its end, 910 (120.6) |
|---|---|---|---|---|---|---|---|---|---|
| 0.7 | 209.5 | 179.8 | 161.7 | 140.9 | 128.9 | 120.6 | 300 (1.00x) | 600 (1.00x) | 910 (1.00x) |
| 1.4 | 212.0 | 173.5 | 152.0 | 128.7 | 117.0 | 109.6 | 260 (1.15x) | 440 (1.36x) | 630 (1.44x) |
| 2.8 | 228.3 | 181.0 | 155.2 | 129.1 | 115.4 | 107.5 | 270 (1.11x) | 450 (1.33x) | 610 (1.49x) |
| 5.6 | 245.5 | 207.3 | 177.6 | 139.5 | 126.4 | 117.7 | 360 (0.83x) | 570 (1.05x) | 820 (1.11x) |

## The hand-over cost, refitted [m data, e model]

`python3 fitcost.py 50 48816 ../../estimates/wil50/subtrees/{a251,a325,a220}/trace.txt`: the 364
missed phases of the laptop's complete wil50 subtrees, by the gap they left, and what their children
then cost per child (m = 50 - depth children):

```
missed phases: 364
gap left (units) | phases | phase it mean | children it per child: mean, median | DNN nodes below per child
   0-5    |  39 |  2121 |    32.0    31.9 |  1.00
   5-10   |  47 |  1323 |    34.5    34.0 |  1.00
  10-20   |  79 |  1065 |    40.3    39.8 |  1.01
  20-30   |  55 |   904 |    60.4    49.2 |  1.07
  30-50   |  68 |   783 |   130.9   117.9 |  1.80
  50-75   |  44 |   652 |   300.8   277.7 |  3.57
  75-100  |  14 |   607 |   697.4   664.1 |  7.92
 100-150  |  13 |   629 |  1684.6  1471.1 | 17.60
 150-250  |   4 |   642 |  5607.5  6109.9 | 59.66
 250-999  |   1 |   630 | 22670.7 22670.7 | 253.66
fit (all): N0 380.6 A -1396.66 B 950.205 | rms log10 residual 1.34
```

- 32 iterations a child at the cutoff, doubling about every 22 units from 25 short on: far steeper
  than a quadratic. In the key's unit (g = 0.1% of the cutoff, 48.8 units) the first fit was
  `32:0:250`. sko49 (`fitcost.py 49 23386` on Cs232 and this study's base runs: 13 missed phases)
  is 27-33 a child up to 20 short, 45 at 31, 130 at 50-60, 725 at 121: `27:0:16`.
- **That curve is not what a longer phase buys.** It compares *different nodes* that the same stall
  rule stopped at different gaps: a node left 100 short is a harder node than one left 85 short. Run
  as the rule (`hoe`), a220's root bought 14.7 units for 1,680 iterations and its children were no
  cheaper. `hoe2`'s constants (wil50 `32:0:25`, sko49 `27:0:5`) keep the measured cost at the cutoff
  (N0) and cut the far term so that a far phase stops about where the stall rule stops it (a220's
  root at 1,120, s232's at 1,040-1,160, b83's at 1,010). They are set from that one observation,
  not fitted: the far term is the open constant.
- Offline replay (`replay.py`, after tho40's `replay_rules.py`) on the base runs' histories, written
  before the search arms: `replay_a220_base.out` (7 phases) and `replay_sko49_base.out` (22). Its
  children column uses the curve above and so shares its flaw; its phase-iteration column said what
  the runs then showed on wil50 (8,499 -> 4,829 phase iterations under `32:0:25`; measured 5,151). On
  sko49 it predicted no gain from any stop rule (x1.00-1.02) because it charges a lost fathom at the
  curve's price; the run found the lost fathoms cheap.

## wil50's larger subtrees [e]

`python3 extrap.py`: the laptop's complete subtrees by where their iterations go, and what the two
effects measured on a220 would save there (early: each warm give-up's first solve ends at 110; near
hand-over: a missed pooled phase ends by 1,050 iterations, its children's cost unchanged). It does not
count early entry's cost on slow first-solve fathoms (b83), so read it as an upper estimate.

| subtree | DNN nodes | iterations | first solves that fathom | give-ups' first solves | fathoming phases | missed phases | [e] early saves | [e] hand-over saves | [e] factor |
|---|---|---|---|---|---|---|---|---|---|
| a220 (laptop) | 232 | 18,852 | 7,245 (38%) | 2,100 (11%) | 1,247 (7%) | 8,260 (44%) | 1,140 (6%) | 3,250 (17%) | 1.30x |
| a251 (laptop) | 3935 | 387,735 | 178,459 (46%) | 57,610 (15%) | 53,036 (14%) | 98,630 (25%) | 39,930 (10%) | 17,590 (5%) | 1.17x |
| a325 (laptop) | 11923 | 1,066,452 | 513,233 (48%) | 154,710 (15%) | 124,869 (12%) | 273,640 (26%) | 105,020 (10%) | 40,210 (4%) | 1.16x |

By this model a220 itself is 1.30x (measured 1.35x on this box's own base). For the tree: about
**1.15x, 1.1-1.2x [e]**, i.e. wil50's 7,500-18,600 core-hours less 10-15%.

## Not done, and what I would do next

- a251 complete (too large for this session); sko49 s232 `hoe2` to its end; warm roots; an early-only arm and a hand-over-only arm
  on a220 (the matched nodes above separate the two effects, an arm of each would do it cleanly).
- `dnn.giveup_warm` between 100 and 300 (150, 200): b83 says 100 is too early for the nodes that
  fathom slowly in their first solve.
- The far term of `cut.handover` on a deep subtree. a220 and b83 are shallow; in a251 79% of the
  iterations sit under two depth-4 nodes whose pooled phases stalled at 600 iterations, 108 and 143
  short. Whether a longer (or `wc=1.4`) phase there pays is the largest open question for wil50,
  and needs a251's `0:0` or `47:0` child run both ways (2-20 CPU-hours an arm).
- `wc=1.4` with `cut.cp_wc_gap` (near misses at the old penalty), on sko49 first.
- sko49's estimate rests on probes, not complete subtrees, in its tail: the factor here is for the
  188-232 short class only.

## Files

| | |
|---|---|
| `run.sh`, `queue.sh`, `jobs.txt`, `queue.log` | the runner, the four-at-a-time queue, the jobs in the order run |
| `par/` | the arms |
| `runs/NAME_ARM/` | `cmd.txt` (command, start, end, load), `log.txt`, `err.txt`, `dnn.txt`, `nodes.txt` |
| `tracelib.py`, `summarize.py` | the trace reader; the table, `--phases`, `--matched` |
| `fitcost.py`, `replay.py`, `replay_*.out` | the hand-over cost by gap; the offline replay and its outputs |
| `wcsweep.py`, `extrap.py` | the `wc` table; the estimate on the laptop's subtrees |
| `lscpu.txt` | the box |
| `mkreadme.sh`, `README.tmpl.md` | this file's assembly |
