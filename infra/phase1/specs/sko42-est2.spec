# sko42's estimate dived again with the estimator's DNN trace (QAP_EST_DNN_TRACE=1): sko42-est-261121b's
# paused frontier (FROM, DIVES: a dives run of a paused run, no new solve phase), the first estimate's settings
# otherwise. Every node the estimator bounds then has its ADMM iterations beside the same node's in
# sko42-261120b's DNN trace (the whole tree, solved): whether a node bounded by the estimator costing 0.79 of
# the same node in the solve (experiments/estimates_nov19, "The estimator against a known tree") is less work
# or the same work timed lower. Also a second draw of the same estimate (the first: 324 CPU h against 539).
# Nov 22: relaunched as sko42-est2-261122. sko42-est2-261121 failed every gate: the coordinating instance
# exports ESTIMATE_ENV, and the solver refused QAP_EST_DNN_TRACE, which it did not register (it does now).
PROBLEM=sko42
RUN=sko42-est2-261122
FROM=sko42-est-261121b
DIVES=900:k1.5
ESTIMATE_ENV="QAP_EST_MIN_DIVES=5 QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03 QAP_EST_DNN_TRACE=1"
VISIBILITY=private
WORKERS=4
HOURS=4
