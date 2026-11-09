# nug30, the Phase 1 acceptance run (docs/distributed/plan.md section 11): the
# Nov 8 record recipe of experiments/nug30_solve/README.md (rule 9,
# nug_deep_xbar9_all.par, DNN at depths 1-7) with one change: QAP_DNN_SHIP=64
# instead of 32. With 64 a node restored from a state file starts its DNN from
# its parent's exact doubles, so the coordinator's tree, cut into units at
# other places than bench/parallel.py's, is the same tree: the same nodes,
# fathoms and give-ups at every depth, which compare.py checks. 32 rounds that
# state to fp32, which can change a node's DNN iterations and, at the cutoff's
# edge, whether it is fathomed; on nug20 (174 nodes, Nov 8) the counts still
# agreed with both. 64 costs larger state files, not correctness.
# The DNN iterations themselves can still differ by a few even with 64: nug24
# with this recipe (Nov 9), the coordinator with 5-second units against a
# serial run, 3,165 against 3,175 at depth 3 and 1,406 against 1,404 at depth
# 4, the tree identical. Explained and fixed Nov 9 (state.C): 64 shipped the
# upper triangles of Z and Y, mirrored on read, while with DNN_CORE_NOSYMV the
# in-process matrices are symmetric to rounding only; 64 now ships the whole
# matrices and make shipcheck's ship64 is exact again.
INSTANCE=qaplib/nug30.dat
PARAMS=../infra/phase1/recipes/nug30.par   # relative to latest/: nug_deep_xbar9_all.par plus the [settings] below
OPT=6124          # the known optimum, which the check expects
UB=6124           # the incumbent given to the solver (--no-ub drops it)
EXTRA=""          # the DNN depths (dnn.depths = 1:7) are in the params file
ENV=""            # experiment switches (DNN_*) would go here, as --env K=V
# units of about a minute (the whole tree is about 1,300 CPU s on the M2 Pro),
# so that units return, split and are stolen as they will at scale
COORD="--budget 60 --target-unit-cpu 60 --steal-after 20 --checkpoint-every 20 --ckpt-upload 60"
