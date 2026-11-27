# A third estimate of tho40 under tho40c, the second of tho40-est2-261119's frontier: its paused solve's
# open nodes dived again (FROM, DIVES: a dives run of a paused run, no new solve phase) with the
# estimator fixed. Until Nov 19 the estimator dropped the warm start of every stratum past the 8th of a
# unit, cut pool included, so a Knuth dive's first bound from 44% of est2's strata ran cold and was
# counted as the tree's (experiments/estimates_nov19, "Knuth dives from cold strata"): est2 read 239 c7a
# CPU-days, 196 with those strata imputed from the warm ones. Now their states are spilled to disk and
# read back. Everything else as est2 had it, so that the difference is the fix.
PROBLEM=tho40
RUN=tho40-est3-261119
FROM=tho40-est2-261119
DIVES=900:k1.5
ESTIMATE_ENV="QAP_EST_MIN_DIVES=5 QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03"
VISIBILITY=private
WORKERS=2
HOURS=10
