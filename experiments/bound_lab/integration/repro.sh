#!/bin/bash
# repro.sh : one cheap reproduction per patch with the merged build (/home/user/qap/latest/bnbdrv, copied to bin_local/),
# one process after another. Output repro/NAME/{cmd.txt,log.txt,err.txt,nodes.txt,dnn.txt}. Iterations and node counts
# are the metrics (the box is shared with the identity runs).
set -u
H=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap; L=$R/experiments/bound_lab
BIN=$H/bin_local/bnbdrv_new_$(sha256sum "$R/latest/bnbdrv" | cut -c1-12); [ -x "$BIN" ] || cp "$R/latest/bnbdrv" "$BIN"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
run() {   # run NAME ARGS...
  local NAME=$1; shift; local D=$H/repro/$NAME; mkdir -p "$D"; [ -f "$D/done" ] && return 0
  rm -f "$D/nodes.txt" "$D/dnn.txt"
  ( cd $R/latest; export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1
    { echo "$(basename "$BIN") $*"; date -u +%FT%TZ; } > "$D/cmd.txt"
    "$BIN" "$@" > "$D/log.txt" 2> "$D/err.txt"; echo "exit $?" >> "$D/cmd.txt"; date -u +%FT%TZ >> "$D/cmd.txt" )
  touch "$D/done"
}
# 1. cut.giveup_depth (nug30_phase arm g2): nug30's nine depth-1 nodes; with cut.giveup = 100:5,300:2 and
#    cut.giveup_depth = 2 their phases run to the stall rule (author: 630-1,240 iterations, 6,570 in all, gaveup 0)
run giveup_depth_nug30_g2 -p $H/repro/par/nug30_g2.par -u 6124 -d 2 ../qaplib/nug30.dat
# 2. dnn.dscale = 0.5, dnn.dscale_until = 100 (admm, section 4): the tai30b depth-3 subtree (author: 195 DNN nodes,
#    184 fathomed, 13,965 iterations; 243 nodes)
run dscale_tai30b_F3_ds0.5u100 -p $L/admm/bb_tai30b_ds0.5u100.par -u 637117113 -F 0:3,1:7,2:10 ../qaplib/tai30b.dat
# 3. cut.handover + wc (tho40_search, "With the faster phase", arm zallhw): the parent 26:10,8:0 bounded and saved,
#    then the warm subtree 38:32 from its state (author: 410 iterations, 1 node)
T="-s SYM/tho40t.sym -u 240516"
run handover_wroot_zallhw -p $L/tho40_search/par/zallhw.par $T -F 26:10,8:0 -d 2 -o $H/repro/handover_wroot_zallhw/root.state ../qaplib/tho40t.dat
if [ ! -f $H/repro/handover_w38_zallhw/done ]; then
  mkdir -p $H/repro/handover_w38_zallhw
  python3 $L/tho40_search/cutstate.py $H/repro/handover_wroot_zallhw/root.state $H/repro/handover_w38_zallhw/in.state 3 38:32 || echo "cutstate failed"
fi
run handover_w38_zallhw -p $L/tho40_search/par/zallhw.par $T -r $H/repro/handover_w38_zallhw/in.state ../qaplib/tho40t.dat
# 4. cut.regions (regions, section 10): tai30b to proven optimality with the line-prefix regions (author: 363 nodes;
#    DNN 299 nodes, 250 fathomed; cut phase 49 nodes, 22 fathomed; in three resumed pieces there, one piece here)
run regions_tai30b_cut -p $L/regions/bnb/tai30b_cut.par -u 637117114 ../qaplib/tai30b.dat
