# The dives estimator checked against a known tree: sko42 under sko42.par as solved by sko42-261120b (Nov 20:
# 128,431 nodes, 539.3 CPU h; experiments/estimates/sko42/zone_nov19, "Solved"), estimated the way tho40 was
# (tho40-est2/est3: Knuth dives k1.5, 900 s a unit, 5 dives at least per branched open node, sampling by depth)
# and the way sko49 would be. The solve phase is short, 0.3 h of wall from the start (the workers join after
# about 10 minutes), so that dives estimate most of the tree from a shallow frontier, as they would for sko49;
# the decision point off, so that the dives run to the end whatever they cost. The incumbent is start.sh's
# default, 15,814, as sko42-261120b's (which found the optimum, 15,812, on the way). Every Spot vCPU of the
# quota: four 64-slot workers. The first launch, sko42-est-261121, stopped at its coordinator's gates (Nov 21: the
# parallel gate's interrupt-and-resume check hit a race in bench/parallel.py, fixed since).
PROBLEM=sko42
RUN=sko42-est-261121b
PARAMS=experiments/estimates/sko42/sko42.par
ESTIMATE=0.3:900:k1.5
ESTIMATE_ENV="QAP_EST_MIN_DIVES=5 QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03"
ESTIMATE_DECIDE=off
VISIBILITY=private
WORKERS=4
HOURS=8
