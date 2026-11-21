# tho40 to the end under tho40c (infra/phase1/recipes/tho40c.par, recipes/defaults.txt's), the settings the
# Nov 19-20 tuning round left as they were (experiments/estimates_nov19, "The arms' result": branch.xbar_lines 4
# neutral, dnn.cross_first 5 1.6x the CPU, cut.zone 0.03 the same trees); dnn.retry = 1 is in it. The tree:
# about 6,700-7,300 CPU h after tho40-est2-261119's solve phase (est3's dives 6,727, the solved sample of the
# tuning arms about 7,250), so about 28 h on 263 slots, a 60 h limit for a tree up to twice that. Public and
# traced: the run of record. The incumbent is start.sh's default, 240,518 (the best known 240,516 + 2).
PROBLEM=tho40
RUN=tho40-261121
PARAMS=infra/phase1/recipes/tho40c.par
WORKERS=4
HOURS=60
