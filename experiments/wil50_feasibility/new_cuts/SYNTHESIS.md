# New cut families for wil50: synthesis of the three-agent study (Nov 4)

Three research agents, about five to seven hours each (05:26-12:05 UTC),
looked for valid inequalities or stronger relaxations that would add another
level of certified bound at depths 3-5 (about +130, a factor of about 47 in
the tree) on top of the production regional-triangle phase. Reports:
`literature/REPORT.md`, `grid/REPORT.md`, `facility/REPORT.md`; running logs
in each `NOTES.md`. The coordinator re-read the headline numbers from the
logs (split bounds, box-cut certificates, node-10 disjunction, validity
checks). The bar set in advance: a median of +40 on top of the converged
triangle phase for at most about twice its CPU, or a credible path to +130.

## Verdict

**No family meets the bar, and none of the three angles finds a credible path
to +130 through linear inequalities on the DNN's matrix.** The three agents
reach the same diagnosis independently:

- **After the triangle phase, every single grid threshold is saturated.**
  Each of the 13 threshold terms, and the row and column problems as wholes,
  sit at or above their own optimum in the objective's direction on nodes 13
  and 10 (grid). That is why every single-region family has failed: 5-point
  counts, local-7, and by the same argument the per-region moment block.
- **The facility side carries nothing extra.** wil50's flows are statistically
  random, and the objective is a sum over the 13 thresholds of flow-weighted
  cuts, so every regional triangle is already an inequality of the flow
  graph's bisections; no facility-side inequality is violated at converged
  iterates (facility).
- **RLT-type strength does not exist at this size.** RLT2 is weaker than the
  DNN at every published n >= 24 (-0.6 gap units at nug30, -1.25 at sko42);
  RLT3 and SDP+RLT2 are 1e3-1e4 times beyond the per-node budget; RLT2's three
  local projections (facility, location and mixed triples) are satisfied at
  our iterates to within 0.5-4% of probability mass and add 0 to +0.5 in a
  pool (literature, facility).
- **What remains is global coupling of the middle thresholds.** Against the
  exact optimum of the 12-region problem built from the five middle
  thresholds, the triangle-phase iterate is 117 short at node 13 and 143 at
  node 10, more for more thresholds (grid). No few-facility or single-region
  family sees it. Harvesting it needs that coarse problem's optimum proven,
  which costs about five levels of a coarse branch-and-bound (+40 per level):
  branching in disguise. Disjunctive and lift-and-project cuts are capped by
  their worst branch, the child the bound already favours (node 10: line
  disjunctions at most +42 to +59, one N+ level +19 to +24), so ordinary
  branching dominates them (facility).

## The one candidate left

**Large-Q 1-box count inequalities on the threshold regions** (Junger-Kaibel),
from the literature agent: on nug12 they close the DNN's gap where the
triangles add nothing (an approximate conic-solver value, not a certified
bound); on wil50 they act on the edge regions where the triangles are weak.
At a fixed multiplier on top of the triangle phase the Lagrangian rises +17.0
on node 13 and +16.8 on d2, but certified it is +14.3 and -4.0: the dense
term slows the certificate. With optimised multipliers inside the pool and a
1,500-iteration quiet round the estimate is +10 to +25, about a fifth of a
level (roughly 1.3-2x on the tree at best), for about a day of work. The
literature agent's kill rule: a median below +15 over the triangle reference
on nodes 10, 13 and d2.

A second, weaker candidate (RLT2 location-side equations on coarse
super-cells, 1-2 weeks, 0 to +30 estimated) is listed by the literature agent
and not recommended.

## Closed by this study

Facility-side families (clusters, facility regions, boxes over facility sets,
quadratic-semi-assignment clique/cut inequalities, stars, parity, dispersion,
flow-weighted sub-QAP cuts); location-side geometric families (parity and
checkerboard regions, clique and star dispersion, free-region, rectangle,
quadrant-XOR and interior row/column triangles, l1/hypermetric/negative-type
conditions); the row/column decomposition; objective splitting with a coarse
optimum (killed by its certification cost); the three local RLT2 projections;
RLT2, RLT3, SDP+RLT2, Lasserre level 2; disjunctive, lift-and-project and
count-split cuts; the certificate's child bounds (0 of 2,209 children
fathomed at nodes 10 and 13, 14 of 2,116 at d2); the per-region moment block
(the bound memo's item 3), by the saturation argument.

## By-products

- **Restart from Y as well as Z.** Restarting a converged node from Z alone
  lost 32 certified units at 400 iterations, against +1.2 with Y and Z
  (facility). This fits the parent-Y child warm start (`../warm_y/`) and
  argues for carrying Y wherever a converged state is resumed: checkpoints,
  state files, children.
- **The host changes across container restarts** (see `README.md`): CPU
  seconds from different restarts are not comparable unless the host is
  logged.
- **Branching.** The best line at the lines-study parents is the lightest
  facility's (facility 45), which rules 8 and 9 already find; the grid agent
  suggests a branching rule preferring the facilities whose coarse marginals
  are most spread over the middle super-cells (untested).

## What it means for the wil50 plan

The bound is not where another factor of 47 will come from. The roadmap's
remaining levers stand as they were: pool inheritance in the search (the
largest unmeasured one), the coordinator, per-node speed and hardware (GPUs),
and above all the sampling plan, which decides whether the tree sits at the
small end of its range. The box-inequality pool test is a cheap optional
side-bet worth at most a fifth of a level.
