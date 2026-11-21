# lsla/grids: the CAP 150 lookahead over every row at the open grids' roots (coordinator, Nov 21)

Task (coordinator, agreed with Nathan): score every row at the roots of tho40, sko49 and wil50 with the
150-iteration lookahead and rank production's root row among them. No whole trees.

**Status: tho40 done; sko49 running (started 17:25 UTC); wil50 to follow.**

## Answer so far

- **tho40: production's root row r26 is 16th of 40 by the lookahead's score, 5.7% behind the best, r21** [m
  `runs/tho40/rank.md`, `analyze.py tho40 12 26`]. **But the ranking decides little at the top: the best three rows
  are within 1.0% (r21, r0, r2), the best ten within 4.2%, and the whole range is 18% (r6 worst, 1.18x).** No
  child fathoms within 150 iterations at any row (all 480 open), so children open says nothing here; only the
  shortfall sum separates the rows.

## Setup

- Build: cognac 0e6c2ecf + `experiments/line_signals/experiment_build.patch` + `../lsla.patch` +
  `../../tail/tail.patch` + `../../root/root.patch` (claude/cool-babbage-ilnogj 4a77c074): with symmetry left,
  a row's children are its J3 columns only.
- Params: each instance's current params with only `dnn.depths` moved to start at 0, so the root gets a DNN bound
  and the lookahead its state:
  - `tho40c_root.par` = origin/cognac `infra/phase1/recipes/tho40c.par`, 1:38 -> 0:38;
  - `sko49_root.par` = `experiments/estimates/sko49/sko49.par`, 1:47 -> 0:47;
  - `fullrun_root.par` = `experiments/wil50_feasibility/fullrun/fullrun.par`, 1:48 -> 0:48.
- `run.sh INST NPROC`: rows interleaved over NPROC processes (each bounds the root itself), one shared
  `runs/INST/children.txt` (QAP_LA_PROGRESS and QAP_LA_RESUME), per part `runs/INST/pN/la.txt`.
- Production's root row: tho40 r26 (`experiments/bound_lab/tho40_dives/logs/root/nodes.txt`), wil50 r21
  (`wil50_feasibility/sampling/README.md`, "Root lines"), sko49 r41 (`estimates/sko49/tune/sampling_C/
  stageA_summary.md`: "root row 41 (10)").

## tho40 (`runs/tho40/`, 16:25-17:24 UTC)

- 40 rows x 12 children, every child solved to the cap: 480 solves, 72,000 iterations, 13,840 CPU s (3.84 core-hours)
  [m: the four parts' la.txt], about 29 s a child; the root bound about a minute per process.
- Checks [m]: every row has exactly 12 children, in columns 0-3, 8-11, 16-19 (the symmetry file's J1, 1-based
  1-4, 9-12, 17-20); `analyze.py`'s sums from `children.txt` match each part's own la.txt to 0.3.
- Best r21 138,880; r0 140,186 (1.009x); r2 140,221 (1.010x); r32 1.021x; r5 1.021x. Production's r26 146,835
  (1.057x, 16th). Worst r6 164,107 (1.182x). Every child's shortfall is about 11,000-14,000 against the cutoff
  240,518 [m].
- So the lookahead would move tho40's root from r26 to r21, a 5.7% smaller sum; whether that is a better tree is
  what the coordinator's tho30 validation (whole trees forced to the best, a middle and the worst row) tests.

