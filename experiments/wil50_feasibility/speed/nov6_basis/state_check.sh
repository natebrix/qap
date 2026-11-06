#!/bin/bash
# State files both ways (latest/state.C): the nug20 DNN search of tests/ship_test.py
# stopped after one CPU slice with QAP_DNN_SHIP=32,
#   v6:  switch off            -> version 6, as before
#   v7:  DNN_SHIP_BASIS=1      -> version 7, the warm lines carry the basis (fp32)
#   v7h: ... and QAP_DNN_SHIP_BASIS_BITS=16 (bf16)
# then: the coordinator's State class (bench/parallel.py, coord/common.py) reads a
# version 7 file and writes it back unchanged; the new solver resumes the rewritten
# file, and a version 6 file with the switch on (no basis: the children refresh);
# the old solver (OLD) refuses a version 7 file and still resumes a version 6 one.
# Logs: logs/state/. OLD: the bnbdrv built before this study. CPU: 1.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
OLD=${OLD:?path of the old bnbdrv}
NEW=${NEW:-$R/latest/bnbdrv}
CPU=${CPU:-1}
L=$D/logs/state
T=$(mktemp -d)
mkdir -p "$L"
cd "$R/latest" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 \
  QAP_DNN_REFRESH=40 DNN_CHOLQR_PASSES=1 QAP_DNN_WARM_Y=2 QAP_DNN_EXT=20 QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_GRID=4:5 \
  QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20 QAP_DNN_SHIP=32
FL="-p PARAM/nug_deep_xbar9_all.par -x 0:4:1:300:1"
DAT=$R/qaplib/nug20.dat
run() { local tag=$1 bin=$2; shift 2; env "$@" taskset -c $CPU "$bin" $FL -u 2570 -o "$T/$tag.state" -t 3 "$DAT" > "$L/$tag.log" 2>&1; echo "$tag: exit $? $(head -1 "$T/$tag.state" 2>/dev/null)"; }
run v6 "$NEW" DNN_SHIP_BASIS=0
run v7 "$NEW" DNN_SHIP_BASIS=1
run v7h "$NEW" DNN_SHIP_BASIS=1 QAP_DNN_SHIP_BASIS_BITS=16
python3 - "$T" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), "..", "bench"))
from parallel import State
t = sys.argv[1]
for tag in ("v6", "v7", "v7h"):
    p = os.path.join(t, tag + ".state")
    warm = [l.split() for l in open(p) if l.startswith("warm ")]
    v = int(open(p).readline().split()[1])
    if v >= 7:
        b = [w[-1] for w in warm if w[-1] != "-"]; d = [w[-2] for w in warm if w[-1] != "-"]
        ms = sorted({int(w[1]) for w in warm if w[-1] != "-"})
        print(f"{tag}: version {v}, {len(warm)} warm lines, {len(b)} with a basis (m {ms}); basis token mean "
              f"{sum(map(len, b)) / max(1, len(b)) / 1e3:.1f} kB ({b[0][:8] if b else '-'}...), DNN state token mean "
              f"{sum(map(len, d)) / max(1, len(d)) / 1e3:.1f} kB")
    else:
        print(f"{tag}: version {v}, {len(warm)} warm lines, last field of each: {sorted({w[-1][:5] for w in warm})}")
    s = State(p)
    s.write(p + ".py")
    lines = lambda f: [l for l in open(f) if l.startswith(("QAP-STATE", "warm ", "node "))]
    same = lines(p) == lines(p + ".py")
    print(f"{tag}: State read and write keeps the version, the warm lines and the nodes: {'yes' if same else 'NO'}")
EOF
resume() { local tag=$1 bin=$2 in=$3; shift 3; env "$@" QAP_DNN_TRACE="$T/$tag.trace" taskset -c $CPU "$bin" $FL -r "$in" -o "$T/$tag.out.state" -t 3 "$DAT" > "$L/$tag.log" 2>&1; local rc=$?
  echo "$tag: exit $rc $(grep -h 'read_state' "$L/$tag.log" | head -1)$( [ -f "$T/$tag.trace" ] && echo "; $(grep -c ' seed 1' "$T/$tag.trace") of $(grep -c ' warm 1 ' "$T/$tag.trace") warm-started nodes seeded")"; }
resume new_reads_v7py "$NEW" "$T/v7.state.py" DNN_SHIP_BASIS=1
resume new_reads_v7h "$NEW" "$T/v7h.state" DNN_SHIP_BASIS=1
resume new_reads_v6_on "$NEW" "$T/v6.state" DNN_SHIP_BASIS=1
resume old_reads_v7 "$OLD" "$T/v7.state"
resume old_reads_v6 "$OLD" "$T/v6.state"
rm -rf "$T"
