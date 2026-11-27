#!/bin/bash
# Keep the line study's arms going on this box (README.md) for up to SECS seconds (default 7000): each
# minute, an arm stopped part-way (no "exit 0" in its cmd.txt, no bnbdrv running for it: a container
# restart or the OOM killer) is trimmed (run.sh trim, when its checkpoint holds more than one DNN state)
# and resumed (run.sh resume), when the memory cgroup's RSS leaves room for it under LIMIT bytes; at most
# 3 resumes per arm per call. Returns when every arm has finished. Run it in the background: it also
# keeps the cloud session's container from being reclaimed while the session is idle.
#   supervise.sh [SECS]          ARMS="r1:r1 c26:c26 r45:r45 r8:r8" (NAME:LINE), LIMIT=13500000000
D=$(cd "$(dirname "$0")"; pwd)
ARMS=${ARMS:-"r1:r1 c26:c26 r45:r45 r8:r8"}
LIMIT=${LIMIT:-13500000000}
CG=/sys/fs/cgroup/memory
end=$(( $(date +%s) + ${1:-7000} ))
declare -A tries
running() { for p in $(pgrep -x bnbdrv); do tr '\0' ' ' < /proc/$p/cmdline | grep -q "runs/$1/out.state" && return 0; done; return 1; }
rss() { awk '/^total_rss /{print $2}' $CG/memory.stat; }
log() { echo "$(date -u +%H:%M) $*" | tee -a "$D/runs/supervise.log"; }
while [ "$(date +%s)" -lt $end ]; do
  left=0
  for a in $ARMS; do
    name=${a%%:*}; line=${a#*:}; S=$D/runs/$name/out.state
    grep -q '^exit 0$' "$D/runs/$name/cmd.txt" && continue
    left=$((left+1))
    running "$name" && continue
    [ "${tries[$name]:-0}" -ge 3 ] && continue
    [ -s "$S" ] || { log "$name: no checkpoint (run.sh arm $name $line by hand)"; tries[$name]=3; continue; }
    nw=$(awk '/^warm_starts /{print $2; exit}' "$S")
    [ "${nw:-0}" -gt 1 ] && { log "$name: trim ($nw DNN states)"; bash "$D/run.sh" trim "$name" 2>&1 | tee -a "$D/runs/supervise.log"; }
    need=$(( 300000000 + 7 * $(stat -c %s "$S") ))
    if [ $(( $(rss) + need )) -gt "$LIMIT" ]; then continue; fi
    tries[$name]=$(( ${tries[$name]:-0} + 1 ))
    log "$name: resume (RSS $(( $(rss) / 1000000 )) MB + ~$(( need / 1000000 )) MB)"
    setsid nohup bash "$D/run.sh" resume "$name" "$line" > /dev/null 2>&1 < /dev/null &
    sleep 5
  done
  [ $left -eq 0 ] && { log "all arms finished"; break; }
  sleep 60
done
for a in $ARMS; do name=${a%%:*}
  echo "$name: $(grep -c '^node 3 ' "$D/runs/$name/nodes.txt") children; $(tail -1 "$D/runs/$name/cmd.txt"); checkpoint $(date -u -r "$D/runs/$name/out.state" +%H:%M)"
done
