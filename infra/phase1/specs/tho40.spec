# tho40 to the end under tho40c_sb (infra/phase1/recipes/tho40c_sb.par): tho40c, the settings the Nov 19-20
# tuning round left as they were (experiments/estimates_nov19, "The arms' result"; dnn.retry = 1 is in it),
# with DNN strong branching at the root and depth 1 (Nov 22, PR 6; RUNNING.md "DNN strong branching (sb.*)"):
# every row at the root, the 4-row shortlist at the 12 depth-1 nodes, depth 2 left off for this first run, and
# dnn.depths from 0 for the root's DNN state. Its cost first: about 4 core-hours for the root in the frontier
# pass (one process, before the run exists: the worker group stays at 0 until then), then the 12 depth-1 nodes,
# up to 2.6 core-hours each, side by side. The tree under tho40c: about 6,700-7,300 CPU h after tho40-est2-261119's
# solve phase (est3's dives 6,727, the tuning arms' solved sample about 7,250), about 28 h on 263 slots; the
# branching lab's root study predicts 0.68 of that under r2, a bet more than a finding
# (experiments/branching_lab/root on claude/cool-babbage-ilnogj). A 60 h limit for a tree up to twice that.
# Public and traced: the run of record. The incumbent is start.sh's default, 240,518 (the best known 240,516 + 2).
# State files become version 10 once strong branching has chosen the root's line: brixius-web reads them from
# 8782050 (main).
PROBLEM=tho40
RUN=tho40-261122
PARAMS=infra/phase1/recipes/tho40c_sb.par
WORKERS=4
HOURS=60
