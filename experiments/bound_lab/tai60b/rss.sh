#!/bin/bash
# rss.sh LOG CMD...: run CMD (output to LOG), then append its wall time and peak resident memory (VmHWM)
L=$1; shift
"$@" > $L 2>&1 &
p=$!; t0=$(date +%s); hwm=0
while kill -0 $p 2>/dev/null; do h=$(awk '/VmHWM/{print $2}' /proc/$p/status 2>/dev/null); [ -n "$h" ] && hwm=$h; sleep 2; done
echo "# wall $(( $(date +%s) - t0 )) s, peak RSS $(( hwm / 1024 )) MB" >> $L
