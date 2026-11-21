# lsla/grids: the CAP 150 lookahead over every row at the open grids' roots (coordinator, Nov 21)

Task (coordinator, agreed with Nathan): score every row at the roots of tho40, sko49 and wil50 with the
150-iteration lookahead and rank production's root row among them. No whole trees.

**Status: tho40 running (started 16:25 UTC); sko49 and wil50 to follow.**

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
