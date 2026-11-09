#!/bin/bash
# Phase 1 acceptance (infra/phase1/README.md) on the laptop with the local: backend,
# mirroring node.sh's init / coordinator+worker / results / reference / compare steps.
set -uo pipefail
cd ~/dev/qap
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
S=/private/tmp/claude-501/-Users-nathan-brixius-dev-brixius-web/f833c398-2e99-4f2e-a75b-8111d7ca8792/scratchpad/phase1_local
RUN=nug20-local-noub; BK="local:$S/store_c"; OUT=$S/out_c; SLOTS=8
mkdir -p "$OUT" "$S/work"
step() { echo "$(date +%H:%M:%S) $*"; }
INSTANCE= PARAMS= OPT= UB= EXTRA= ENV= COORD=
source infra/phase1/recipes/nug30.sh; INSTANCE=qaplib/nug20.dat; OPT=2570; UB=
step "init $RUN: $INSTANCE $PARAMS ub $UB slots $SLOTS; $COORD"
python3 -m coord.qapctl init "$RUN" --backend "$BK" --instance "$INSTANCE" --params "$PARAMS" --slots $SLOTS "--extra=$EXTRA" --visibility private $COORD > "$OUT/init.log" 2>&1 || { step "init FAILED"; tail -20 "$OUT/init.log"; exit 1; }
step "init: $(grep -o 'frontier at depth .*' "$OUT/init.log")"
python3 -m coord.qapctl coordinator "$RUN" --backend "$BK" > "$OUT/coordinator.log" 2>&1 & CL=$!
python3 -m coord.qapctl worker "$RUN" --backend "$BK" --slots $SLOTS --workdir "$S/work" > "$OUT/worker.log" 2>&1; WR=$?
wait $CL; CR=$?
step "worker exit $WR, coordinator exit $CR, $(grep -c ' -> ' "$OUT/worker.log") units"
python3 infra/phase1/bk.py "$BK" get "runs/$RUN/results.txt" "$OUT/results.txt" || { step "no results.txt"; exit 1; }
step "finished: $(grep -E 'Nodes:|best objective|Total Time' "$OUT/results.txt" | tr -s ' ' | tr '\n' ' ')"
step "reference: bench/parallel.py --jobs $SLOTS"
python3 bench/parallel.py "$INSTANCE" --jobs $SLOTS --params "$PARAMS" "--extra=$EXTRA" > "$OUT/reference.log" 2>&1 || step "reference run failed"
python3 infra/phase1/compare.py "$OUT/results.txt" "$OUT/reference.log" --opt "$OPT" > "$OUT/compare.txt"; RC=$?
cat "$OUT/compare.txt"; step "compare exit $RC"
