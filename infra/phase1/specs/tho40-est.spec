# A second estimate of tho40 under tho40c (recipes/defaults.txt), with what the first one's diagnosis
# asked for (experiments/tho40_est_diag): its dives starved by setup (68 of 95 estimator CPU hours),
# most strata one dive, sampling by QPB weight at about 20% at every depth. Here: the dive budget
# counts after setup (the solver's default since Nov 19), 5 dives at least per branched open node,
# and sampling by depth (the depth-3 open nodes carry 85% of the CPU: most of them; depths 5-6, which
# hold little, a few percent). First estimate: about 190 c7a CPU-days (140-250).
PROBLEM=tho40
RUN=tho40-est2-261119
ESTIMATE=1:900:k1.5
ESTIMATE_ENV="QAP_EST_MIN_DIVES=5 QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03"
VISIBILITY=private
WORKERS=2
HOURS=10
