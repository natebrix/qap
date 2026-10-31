#!/bin/bash
# Both tuning experiments, one after the other (give-up first, so its timings
# are not shared with the elimination runs). Resumable.
D=$(cd "$(dirname "$0")" && pwd)
"$D/giveup/run_cloud.sh"
"$D/elim/run_cloud.sh"
