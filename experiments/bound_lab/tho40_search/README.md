# tho40: the grid cut phase in the search (bound_lab agent `tho40_search`, Nov 15 2026)

DRAFT, runs still in flight: see "Status" at the end. **[m]** measured here (run named), **[e]** estimated.

Question (`../COMMON.md`): tho40 has never run a cut phase. On complete subtrees, end to end, what is
the grid cut phase worth, and which `cut.zone`?

## How it was run

- One `latest/bnbdrv` process per subtree and arm (`run.sh`; two at a time, `queue.sh` over `jobs.txt`),
  `qaplib/tho40t.dat`, `latest/SYM/tho40t.sym`, `-u 240516`, `QAP_NODE_TRACE` and `QAP_DNN_TRACE` on.
  Each run directory `runs/NAME_ARM/` has `cmd.txt` (the command), `log.txt`, `nodes.txt`, `dnn.txt`.
- **Arms** (`par/`): `base` = `experiments/estimates/tho40/tho40.par`, unchanged. `z003`, `z010`, `z030` =
  base plus the block at the end of this file with `cut.zone` = 0.003, 0.01, 0.03.
- **What the zone does** (`latest/bnb.C`, line 884): after a node's first solve ends unfathomed, the phase
  runs when `cutoff - first-solve bound <= cut.zone x |incumbent|`: 722, 2,405 and 7,215 units short for
  the three zones. The first solve is unchanged by the zone: it gives up from iteration `dnn.giveup` = 200
  when its extrapolated bound cannot reach the cutoff, and `dnn.gu_zone` = 0.003 makes a give-up within
  722 of the cutoff wait 100 more iterations. The phase (CP form) then runs at most
  `cut.rounds x cut.iters + cut.freeze` = 3,400 iterations and stops earlier when it fathoms or when,
  after 600 iterations (`cut.stall` 3 x `cut.iters`), the gap is more than 3 x the gain of the last 600.
  An unfathomed phase hands its pool to the children (`pool.on`), whose first solve then carries it.
- **Roots.** `t1`-`t4` and `gate` start cold from `-F`, in every arm, so arms are matched; but a cold root's
  first solve stops at 200 iterations far looser than the same node in a run, which starts from its
  parent (t2: 3,193 short cold, 1,931 in production). `w18`, `w38`, `w30` start **warm**: the A/B root's
  state file (`runs/wroot_z030/root.state`, `-d 2 -o`), cut to one child (`cutstate.py`, as
  `path_check/dive.py --warm`), resumed with `-r`. Their first solves reproduce the production trace's
  bounds to the digit (root 231,879.8596; `18:32` 4,474.4 short).
- **Reference.** The production run of the subtree `-F 26:10,8:0` (`infra/phase1/runs/tho40ab-32/`,
  22,650 nodes, 78,801 c7a CPU s, the base settings). The base arm here gave production's node counts
  exactly on t1, t2, t3 (106, 211, 287) and, where both branched on the same lines, the same outcome at
  every node (`calib.py`); this box's CPU per matched node is 0.70 of production's.
- The box was shared (load 5-7 on 4 cores, my two processes at about 50% of a core each): **iterations and
  nodes are the metrics**; CPU s are the solver's own CPU time, indicative. A phase iteration costs about
  1.25x a first-solve iteration (0.126 against 0.103 s at m = 36).

TABLES_GO_HERE
