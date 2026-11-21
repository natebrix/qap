#!/bin/bash
# lsla/grids: the CAP 150 lookahead over every row at a grid's root (coordinator, Nov 21; branching_lab/root).
# Build: cognac 0e6c2ecf + experiment_build.patch + ../lsla.patch + ../../tail/tail.patch + ../../root/root.patch
# (in /home/user/wt-cognac/latest_root). With symmetry left, a row's children are its J3 columns only.
#   ./run.sh INST NPROC   start NPROC processes, rows interleaved; all share runs/INST/children.txt (resume)
#   ./run.sh INST NPROC I one part (I = 0..NPROC-1), in the foreground
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
BUILD=/home/user/wt-cognac/latest_root
INST=$1; NP=$2
case $INST in
  tho40) DAT=tho40t; U=240518; N=40 ;;
  sko49) DAT=sko49t; U=23386;  N=49 ;;
  wil50) DAT=wil50t; U=48816;  N=50 ;;
  *) echo "unknown $INST"; exit 2 ;;
esac
PAR=$D/$( case $INST in tho40) echo tho40c_root.par;; sko49) echo sko49_root.par;; wil50) echo fullrun_root.par;; esac )
O=$D/runs/$INST; mkdir -p "$O"
if [ $# -lt 3 ]; then
  for ((i = 0; i < NP; i++)); do setsid nohup "$0" "$INST" "$NP" "$i" > "$O/p$i.out" 2>&1 < /dev/null & done
  exit 0
fi
I=$3; LINES=$(for ((r = I; r < N; r += NP)); do printf "r%d," $r; done); LINES=${LINES%,}
P=$O/p$I; mkdir -p "$P"
{ echo "QAP_LA=1:150:0:0 QAP_LA_LINES=$LINES ./bnbdrv -p $PAR -s SYM/$DAT.sym -u $U -d 0 ../qaplib/$DAT.dat";
  md5sum "$BUILD/bnbdrv"; date -u +%FT%TZ; } > "$P/cmd.txt"
cd "$BUILD"
set +e
QAP_LA=1:150:0:0 QAP_LA_LINES=$LINES QAP_LA_LOG=$P/la.txt QAP_LA_PROGRESS=$O/children.txt QAP_LA_RESUME=$O/children.txt \
  QAP_NODE_TRACE=$P/nodes.txt ./bnbdrv -p "$PAR" -s SYM/$DAT.sym -u $U -d 0 ../qaplib/$DAT.dat > "$P/log.txt" 2> "$P/err.txt"
rc=$?
{ echo "exit $rc"; date -u +%FT%TZ; } >> "$P/cmd.txt"
