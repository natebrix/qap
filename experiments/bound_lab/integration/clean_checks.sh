#!/bin/bash
# clean_checks.sh : a clean checkout of the merged commit (git worktree /home/user/wt-clean), built from nothing,
# then the repo's checks one after another. Logs in logs/clean_*.log, each ending with its exit code.
L=$(cd "$(dirname "$0")" && pwd)/logs; cd /home/user/wt-clean || exit 1
git rev-parse HEAD > $L/clean_commit.txt; git status --short >> $L/clean_commit.txt
step() { local n=$1; shift; { echo "\$ $*"; date -u +%FT%TZ; "$@"; echo "exit $?"; date -u +%FT%TZ; } > $L/clean_$n.log 2>&1; }
step build make -C latest -j2
step build_dnntest make -C latest dnntest
step build_dnncuts make -C latest dnncuts
step make_test make test
step nodetrace make nodetrace
step shipcheck make shipcheck
step dnncheck make dnncheck
step continuationcheck make continuationcheck
