#!/bin/bash
# Keep the arms going (cognac's tai60b_lines/supervise.sh, adapted): each minute an arm stopped part-way
# (no "exit 0" in its cmd.txt and no bnbdrv for it) is trimmed if its checkpoint holds more than one DNN
# state and resumed when memory allows; at most 3 resumes per arm. Every 15 minutes the logs and traces
# are committed and pushed. ARMS="c54:c54 c47:c47 r31:r31", LIMIT bytes of RSS, SECS (default 25000).
D=$(cd "$(dirname "$0")"; pwd)
ARMS=${ARMS:-"c54:c54 c47:c47 sc54:c54 sc47:c47"}
LIMIT=${LIMIT:-13500000000}
end=$(( $(date +%s) + ${1:-25000} ))
declare -A tries
running() { for p in $(pgrep -x bnbdrv); do tr '\0' ' ' < /proc/$p/cmdline | grep -q "runs/$1/out.state" && return 0; done; return 1; }
rss() { ps -C bnbdrv -o rss= | awk '{s+=$1} END {print s*1024+0}'; }
log() { echo "$(date -u +%H:%M) $*" | tee -a "$D/runs/supervise.log"; }
last=0
while [ "$(date +%s)" -lt $end ]; do
  for a in $ARMS; do
    name=${a%%:*}; line=${a#*:}; S=$D/runs/$name/out.state
    grep -q '^exit 0$' "$D/runs/$name/cmd.txt" 2>/dev/null && continue
    running "$name" && continue
    [ "${tries[$name]:-0}" -ge 3 ] && continue
    [ -s "$S" ] || continue
    nw=$(awk '/^warm_starts /{print $2; exit}' "$S")
    [ "${nw:-0}" -gt 1 ] && { log "$name: trim ($nw DNN states)"; bash "$D/run.sh" trim "$name" 2>&1 | tee -a "$D/runs/supervise.log"; }
    need=$(( 2000000000 + 7 * $(stat -c %s "$S") ))
    [ $(( $(rss) + need )) -gt "$LIMIT" ] && continue
    tries[$name]=$(( ${tries[$name]:-0} + 1 ))
    log "$name: resume"
    if [ "${name#s}" != "$name" ]; then cmd=(sresume "$line"); else cmd=(resume "$name" "$line"); fi
    setsid nohup bash "$D/run.sh" "${cmd[@]}" > /dev/null 2>&1 < /dev/null &
    sleep 5
  done
  if [ $(( $(date +%s) - last )) -ge 900 ]; then
    last=$(date +%s)
    ( cd "$D" && python3 -I score.py c54 c47 sc54 sc47 r31 > out/score.txt 2>&1; cd "$D/../../.." && git add experiments/branching_lab/linetest &&
      git commit -q -m "branching_lab/linetest: logs in progress $(date -u +%H:%M)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014vvMYRPMz8BP8Xc3WigdYm" && git push -q origin claude/cool-babbage-ilnogj-linetest ) > /dev/null 2>&1
  fi
  sleep 60
done
