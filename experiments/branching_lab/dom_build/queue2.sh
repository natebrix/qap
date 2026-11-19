#!/bin/bash
# after queue.sh: the key-on machinery and the repo's own checks on the patched tree
D=/home/user/qap/experiments/branching_lab/dom_build; cd $D
until grep -q "QUEUE DONE" logs/slices.log 2>/dev/null; do sleep 30; done
python3 keyon_mech.py > logs/keyon_mech.log 2>&1
cd /home/user/wt-dombuild
for t in test nodetrace dnncheck continuationcheck shipcheck; do
  echo "== make $t" >> $D/logs/repo_checks.log
  ( time make $t ) >> $D/logs/repo_checks_$t.log 2>&1; echo "make $t: exit $?" >> $D/logs/repo_checks.log
done
echo QUEUE2 DONE >> $D/logs/repo_checks.log
