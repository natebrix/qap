# nug30, the Phase 1 acceptance run (docs/distributed/plan.md section 11): the
# Nov 8 record recipe of experiments/nug30_solve/README.md (rule 9,
# nug_deep_xbar9_all.par, DNN at depths 1-7) with one change: QAP_DNN_SHIP=64
# instead of 32. With 64 a node restored from a state file starts its DNN from
# its parent's exact doubles, so the coordinator's tree, cut into units at
# other places than bench/parallel.py's, gives the same totals node for node
# (RUNNING.md, make shipcheck). 32 rounds that state to fp32, which can change
# a node's DNN iterations and, at the cutoff's edge, whether it is fathomed;
# on nug20 (174 nodes, Nov 8) the counts still agreed with both. 64 costs
# larger state files, not correctness.
INSTANCE=qaplib/nug30.dat
PARAMS=PARAM/nug_deep_xbar9_all.par
OPT=6124          # the known optimum, which the check expects
UB=6124           # the incumbent given to the solver (--no-ub drops it)
EXTRA="-x 1:7:1:200:1"
ENV="QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_BETA=0.21
     QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003
     QAP_DNN_SHIP=64 DNN_R32_NONCERT=1 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1
     DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 DNN_BETA_SCHED=w=0.36"
# units of about a minute (the whole tree is about 1,300 CPU s on the M2 Pro),
# so that units return, split and are stolen as they will at scale
COORD="--budget 60 --target-unit-cpu 60 --steal-after 20 --checkpoint-every 20 --ckpt-upload 60"
