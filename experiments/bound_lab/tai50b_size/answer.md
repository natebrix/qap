**tai50b with the cluster-region cut phase is a run of about 100 core-hours on this machine: roughly 5,000 to
13,000 nodes and 2e5 to 6e5 CPU-seconds. Hundreds of core-hours at the outside, not thousands.** [e]

| estimate | nodes | CPU-s | core-hours here | what it is |
|---|---|---|---|---|
| E1 | 4,573 (2,473 to 7,135) | 228,101 | **63** (35 to 105) | depths 0-1 measured, 18 of the 50 depth-1 nodes; below them 275 Knuth dives, expanded by the sampling design alone |
| E2 | 12,023 (5,786 to 19,493) | 420,631 | **117** (75 to 164) | the same runs, with the two depth-1 children that hold a near-optimal solution counted as copies of the on-path child's subtree |
| E3 | 13,189 | 601,790 | **167** | a model of which nodes cannot fathom (those holding a solution within the gap measured on the best known solution's path), priced with measured node costs; 60 core-hours at half that gap |
| far case | 66,159 | 3,615,595 | 1,004 | E3 at twice the measured gap, which the nodes bounded here contradict (1 of 199 nodes that far out branched) |

- Intervals are 5th to 95th bootstrap percentiles and are too narrow: they do not know that the levels that
  carry the cost have few dives (see "What could make it wrong"). **My range is 60 to 250 core-hours**, with E2
  as the number I would plan on. On a c7a core [e]: about 1.4 times these figures, if this box's 0.70 of a c7a
  per matched tho40 node carries over (section 1).
- The old estimate (no cut phase) was a floor of 1e11 CPU-s. This one is 5 to 6 orders of magnitude lower.
- **Where the cost is** (E2 by depth, section 7): depth 5 36%, depth 2 20%, depth 6 17%, depth 4 12%, depth 1
  10%, the root 1%; everything from depth 7 down 3%. The DNN iterations are the whole cost (0.3 to 0.7 s each).
- **The shape of the tree.** The root ends 5.27M (1.15%) under the best known value. 11 of 17 sampled
  depth-1 children fathom there. Under the ones that branch, all 93 depth-2 node bounds of the dives (79 distinct nodes) fathomed.
  What does not fathom is the path of the best known solution and the few nodes that hold a solution nearly
  as good: the bound on that path is 3.2M short at depth 1, 0.63M at depth 5, 0.11M at depth 8, 5k at depth
  13, and a node holding a solution within about that gap has to be branched. The landscape is flat (150
  2-swaps of the best known solution are within 1M of it), so this "web" is a few dozen open nodes wide at
  depths 4 and 5 and thin again by depth 7.

**Settings** (`par/A.par`: `experiments/estimates/tai50b/tai50b.par` with these changed or added; the root
separately with `cut.stall = 0` and `cut.cp_iters = 6000`, `par/root_long.par`):

```
dnn.depths = 0:48        # the root too: it carries the cut pool and the state its children start from
dnn.giveup = 100         # a pooled first solve hands over to the phase early
pool.on = 1
cut.regions = experiments/bound_lab/regions/regions/tai50b.comb.txt     # 61 location-side regions (absolute path in the .par)
cut.side = 0
cut.cp = 100:1000
cut.warmy = 1
cut.zone = 1             # a phase at every DNN miss
cut.rounds = 12          # x cut.iters = 2,400 phase iterations at most; I would raise it to 20 (section 9) [e]
cut.iters = 200
cut.freeze = 0
cut.stall = 3
```

Needs `../regions/regions.patch`. These are the settings every number here was measured under, chosen once from
the tho40 block and the regions README, not tuned: no A/B was run (section 9 says what the runs suggest).

**What it rests on** [m]

- The root: **453,552,676 certified** after 470 + 6,000 iterations (4,341 CPU-s), a new best root bound, and
  not yet the plateau (section 3).
- 18 of the 50 depth-1 nodes, warm from the root with its pool: the on-path child and a random 17 of the other
  49 (section 4).
- The path of the best known solution bounded node by node to depth 13, and the path of a second near-optimal
  solution (under `31:32`) to depth 4, which behaves the same (section 5).
- 275 dives from saved states, stratified by the path's levels, 355 node bounds in all (section 6); a model of
  the open nodes checked against the 235 distinct nodes bounded (section 7).
- Validity: 19 node bounds on the path of the best known solution, none fathomed, none above 458,821,517;
  tai30b solved from an incumbent one above its optimum in 363 nodes, as in the regions session (sections 2, 8).
- 49,000 CPU-s of node bounds in total, on 4 cores.

**What could make it wrong**

1. **The web is thinly sampled.** Levels 4 and 5 carry about half of E2 and have 8 and 5 dives, each stopped
   by its wall limit with a long dive in flight; the dive estimator weights a rare child by 100 or more, so one
   dive moves those levels by tens of thousands of CPU-s. E1 and E2 can be off by 2x either way from this alone.
2. **The twins are assumed, not measured, below depth 4.** E2 gives `31:0` (never bounded) and `31:32` the
   on-path child's subtree. E3 assumes a twin's gap is the spine's, that 2-swap descent finds the best solution
   inside a node, and that every node branches on the spine's line. More near-optimal solutions than the
   descent finds would widen the web; the nodes bounded so far agree with the rule at scale 1 (section 7).
3. **32 depth-1 nodes were not bounded.** By the landscape only two of them (`31:0`, `31:42`) are near, and E2
   and E3 count them; a far node that stalls anyway costs 49 quick children (4 of the 15 sampled far nodes did).
4. **Below depth 13 on the path nothing was run.** The estimates extrapolate the path to depth 45 at 13,000
   CPU-s per path (E2 counts three such paths: 9% of it), on the pattern that the bound there is within 5k and its children fathom in seconds.
   tai30b's path behaved that way. If the open nodes near the bottom multiply (many solutions within a few
   thousand of the best known: the nearest 2-swap is 1,534 away), the bottom of the tree grows; those nodes
   are small (m under 37) and cheap, so I expect nodes more than CPU.
5. **The cut code's validity is the regions session's**, plus the path checks here. A wrong cut that fathoms
   off the path would make every number here too small and the run's answer wrong. tai35b and tai40b solved
   end to end are the test, not yet run.
6. If 458,821,517 is not optimal, the search finds the better value and the tree is smaller.
7. CPU-seconds are this machine's (263 ms per plain iteration at m = 48); nodes and iterations carry over.
8. In the other direction: a longer root, a longer phase limit and any tuning would shrink it (section 9).
