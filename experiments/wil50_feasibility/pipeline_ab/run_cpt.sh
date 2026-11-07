#!/bin/bash
# The node-10 A/B in the production configuration (README.md, "Results, Nov 7
# (cloud box): the production configuration"): bestcoret on cores 0-1 and cpcutt on
# cores 2-3, 2 slots each, at the same time, into runs_cpt/ (PAB_RUNS to change it),
# as runs_cp/'s bestcore and cpcut were run. Resumable: run it again.
D=$(cd "$(dirname "$0")" && pwd)
export PAB_RUNS=${PAB_RUNS:-$D/runs_cpt} SLOTS=2
mkdir -p "$PAB_RUNS"
taskset -c 0,1 "$D/run_pair.sh" bestcoret > "$PAB_RUNS/launch_bestcoret.log" 2>&1 &
sleep 3
taskset -c 2,3 "$D/run_pair.sh" cpcutt > "$PAB_RUNS/launch_cpcutt.log" 2>&1 &
wait
python3 "$D/summarize.py" "$PAB_RUNS" | tee "$PAB_RUNS/summary.txt"
