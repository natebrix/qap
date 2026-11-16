#!/bin/bash
# A core's queue of dive sets, one after the other:  ./queue.sh CORE "SET SEED STATE LEVEL SECONDS" ...
cd "$(dirname "$0")"; CORE=$1; shift
for job in "$@"; do set -- $job; ./run_dive.sh $1 $2 A $PWD/$3 $CORE $4 $5; done
