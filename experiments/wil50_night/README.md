# wil50_night: screening settings on two complete wil50 subtrees (Nov 24, in progress)

Nathan, Nov 24: "I want to gun for any remaining perf improvements - either in the solver itself, or via parameter
tuning." Each arm is `par/base.par` (a copy of cognac's `infra/phase1/recipes/wil50_d2r10.par`) plus the lines at
its end, run cold from `-F` on the bound lab's two wil50 subtrees (a220 `21:22,31:23,40:15`, b83
`8:9,21:11,45:40,34:0`; `experiments/bound_lab/wil50_policy` in cognac), one process a run, four at a time
(`run.sh`, `queue.sh`, `jobs.txt`), incumbent 48,818, solver as cognac 626230c8. Results to follow here.

Arms so far, each untried on wil50 by the repository's own studies (a subagent's inventory, Nov 24):
`wcg150` (wc = 1.4 with `cut.cp_wc_gap = 150`), `extrap` (`cut.cp_extrap = 100:1:30:2`), `refresh80`
(`dnn.refresh = 80`), `gw150` (`dnn.giveup_warm = 150`), `xl2` (`branch.xbar_lines = 2`).
