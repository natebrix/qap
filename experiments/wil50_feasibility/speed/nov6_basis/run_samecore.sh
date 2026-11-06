#!/bin/bash
# The replicate on one core (README section 5): f9 and f6 of the warm set with the
# switch off and on, one run after the other on CPU (default 0), in the order
# base ship ship base, so that a difference between cores or between the two
# halves of the hour cannot look like a saving. Logs: samecore/<config>_<family>.log;
# then summarize.py --dir samecore sA sB, and the replicates sA sA2, sB sB2.
D=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:?a copy of the dnntest to measure}
export BIN
OUT=$D/samecore FAMILIES=${FAMILIES:-"f9 f6"} CPU=${CPU:-0} \
  CONFIGS="sA@DNN_SHIP_BASIS=0 sB@DNN_SHIP_BASIS=1 sB2@DNN_SHIP_BASIS=1 sA2@DNN_SHIP_BASIS=0" bash "$D/run_warmset.sh"
{ python3 "$D/summarize.py" --dir "$D/samecore" sA sB sB2 --table; python3 "$D/summarize.py" --dir "$D/samecore" sA2 sB sB2;
  python3 "$D/summarize.py" --dir "$D/samecore" sA sA2; python3 "$D/summarize.py" --dir "$D/samecore" sB sB2;
  python3 "$D/summarize.py" --dir "$D/samecore" --profile sA sB sB2 sA2; } > "$D/samecore/summary.txt"
