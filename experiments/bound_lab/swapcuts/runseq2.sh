#!/bin/sh
# runseq.sh JOBFILE: lines "name|arguments of swapdrv"; one process at a time; logs/NAME.log, first line = the command
cd "$(dirname "$0")"
while IFS='|' read -r name args; do
  [ -z "$name" ] && continue
  [ -f "logs/$name.done" ] && continue
  echo "# ./swapdrv2 $args" > "logs/$name.log"
  ./swapdrv2 $args >> "logs/$name.log" 2>&1
  touch "logs/$name.done"
done < "$1"
