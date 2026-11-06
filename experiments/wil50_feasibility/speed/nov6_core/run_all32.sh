#!/bin/bash
# A/B of this study's switches on the 32 wil50 nodes of tuning/giveup/all32.txt,
# cold, in the production `fast` configuration of pipeline_ab/run_pair.sh as far as
# dnntest can express it: degree 3, tight lambda (-L 2), PRSM (-g 1.0 -G 1.0), a full
# eigendecomposition every 40 iterations (-r 40), give-up after 300 (-H 300),
# -i 2500 -u 48816, one CholQR pass, AMX, R's rank term in fp32 off certificate
# iterations, the two-phase beta (cold calls). ONE process at a time on CPU 3.
# The configurations alternate node by node (base, then each variant, next node),
# so a drift of the host's load hits all of them alike. A private copy of the
# binary is taken at the first run, so a rebuild during the run changes nothing.
#
# CONFIGS="name@env ..." (env: assignments joined by '+', or '-' for none);
# NODES: line numbers of all32.txt (default 1..32). Resumable: a node whose log has a
# RESULT line is skipped. Logs: runs/<name>_<id>.log. Tables: summarize.py base lazy ...
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CPU=${CPU:-3}
# REFRESH_R=N in a configuration's env list sets dnntest's -r (default 40); dnntest ignores the variable.
CONFIGS=${CONFIGS:-"base@- lazy@DNN_CORE_LAZYCERT=1 lazysym@DNN_CORE_LAZYCERT=1+DNN_CORE_NOSYMV=1+DNN_NOSYMW=1 bundle@DNN_CORE_LAZYCERT=1+DNN_CORE_NOSYMV=1+DNN_NOSYMW=1+DNN_CORE_REFRESH_PARTIAL=1 bundle_r80@DNN_CORE_LAZYCERT=1+DNN_CORE_NOSYMV=1+DNN_NOSYMW=1+DNN_CORE_REFRESH_PARTIAL=1+REFRESH_R=80"}
NODES=${NODES:-$(seq 1 32)}
PROD="DNN_AMX=1 DNN_CHOLQR_PASSES=1 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7"
mkdir -p "$D/runs" "$D/bin"
if [ ! -x "$D/bin/dnntest" ]; then
  cp "$R/latest/dnntest" "$D/bin/dnntest"
  (cd "$R" && git rev-parse --short HEAD && git diff --quiet HEAD -- latest || echo "latest/ modified") > "$D/bin/commit.txt" 2>/dev/null
fi
BIN=$D/bin/dnntest
echo "$(date -u +%FT%TZ) start: configs [$CONFIGS] nodes [$(echo $NODES | tr '\n' ' ')] load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$D/runs/index.txt"
for k in $NODES; do
  spec=$(sed -n "${k}p" "$D/../../tuning/giveup/all32.txt" | awk '{print $1}')
  for c in $CONFIGS; do
    name=${c%%@*}; envs=${c#*@}
    [ "$envs" = "-" ] && envs=""
    f=$D/runs/${name}_$k.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    r=40; case "$envs" in *REFRESH_R=*) r=$(echo "$envs" | sed -E 's/.*REFRESH_R=([0-9]+).*/\1/');; esac
    echo "$(date -u +%T) $name $k load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$D/runs/index.txt"
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $PROD ${envs//+/ } \
      taskset -c $CPU "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r $r -H 300 -i 2500 -u 48816 -v \
      -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  done
done
echo "$(date -u +%FT%TZ) done" >> "$D/runs/index.txt"
