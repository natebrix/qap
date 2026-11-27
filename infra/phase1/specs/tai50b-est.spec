# An estimate of tai50b under its current recipe (recipes/tai50b.par: rule 9 at the bottom), whose
# tree is known: tai50b-r9's 6,629 nodes and 419,520 CPU s (117 core-hours; infra/phase1/runs/
# tai50b-r9). The first estimate on worker groups: the root (72 minutes on one core) on the
# on-demand coordinator, then 1 hour of solving on 2 workers, then 600 s of survivor dives a unit.
PROBLEM=tai50b
RUN=tai50b-est-261119
ESTIMATE=1:600
VISIBILITY=private
WORKERS=2
HOURS=8
