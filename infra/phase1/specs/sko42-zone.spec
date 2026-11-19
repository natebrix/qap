# The cut.zone A/B on sko42 from scratch (docs/review_nov18_runs.md: the Nov 18 arms differed in
# tracing and were both reclaimed 30 minutes in): 0.010 (sko42_z10.par, sko42.par until Nov 19) against
# 0.025 (sko42_z25.par: that line changed, rule 9 at the bottom in both), both traced,
# on on-demand workers so that neither arm is interrupted.
PROBLEM=sko42
RUN=sko42-zone-261119
ARMS="z10=experiments/estimates/sko42/sko42_z10.par z25=experiments/estimates/sko42/sko42_z25.par"
MARKET=on-demand
WORKERS=1
HOURS=24
