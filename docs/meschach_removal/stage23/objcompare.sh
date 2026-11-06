#!/bin/sh
# Machine-code comparison of two builds of latest/ (run from the repo root):
#   sh docs/meschach_removal/stage23/objcompare.sh OLD_LATEST NEW_LATEST
# objdiff.sh over every object, then, for qpbnd.o and qpbnd2.o (where PERM
# was replaced), fndiff.py per function and fpops.py for the functions that
# changed (init_ST, init_w0).
D=docs/meschach_removal/stage23
old=$1 new=$2
echo "# objdiff.sh: every object (same = identical disassembly)"
sh $D/objdiff.sh "$old" "$new"
for o in qpbnd qpbnd2; do
  echo; echo "# fndiff.py $o.o: functions that differ"
  python3 $D/fndiff.py "$old/$o.o" "$new/$o.o"
done
echo; echo "# fpops.py: floating-point instructions and calls of the changed functions"
python3 $D/fpops.py "$old/qpbnd.o" "$new/qpbnd.o" 'init_ST(MAT'
python3 $D/fpops.py "$old/qpbnd2.o" "$new/qpbnd2.o" 'init_w0(MAT'
