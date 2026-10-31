#!/bin/bash
# Both tuning experiments, one after the other (give-up first, so its timings
# are not shared with the elimination runs). Resumable.
D=$(cd "$(dirname "$0")" && pwd)
"$D/giveup/run_cloud.sh"
"$D/elim/run_cloud.sh"
# then the gap-ratio give-up rule (needs a dnntest with -R: rebuild first)
make -C "$D/../../.." solver dnntest > /dev/null 2>&1
"$D/giveup/run_cloud_R.sh"
