#!/bin/sh
# queue.sh FILE: run the jobs of FILE one after another (each line: [VAR=..] NAME SPEC options...);
# lines are consumed from the top, so the file can be appended to while the queue runs.
D=/home/user/qap/experiments/bound_lab/cpfast
while true; do
  line=$(grep -v '^#' "$1" | grep -v '^$' | head -1)
  [ -z "$line" ] && break
  # drop the line
  awk -v l="$line" 'done || $0 != l { print; next } { done = 1 }' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
  echo "$(date +%H:%M:%S) start $line" >> $D/logs/queue.log
  sh -c "env $(echo "$line" | awk '{for(i=1;i<=NF;i++){ if ($i ~ /^[A-Z]+=/) printf "%s ", $i; else break } }') $D/scripts/run.sh $(echo "$line" | awk '{s=0; for(i=1;i<=NF;i++){ if (!s && $i ~ /^[A-Z]+=/) continue; s=1; printf "%s ", $i } }')"
  echo "$(date +%H:%M:%S) done  $line" >> $D/logs/queue.log
done
