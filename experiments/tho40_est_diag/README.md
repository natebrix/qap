# tho40-est: what its estimate says, and how far to trust it (Nov 18)

The estimate under diagnosis: `infra/phase1/runs/tho40-est-dives/`.
- **Setup:** tho40c.par, incumbent 240,518, 1-hour solve on 128 c7a slots.
- **Dives:** 20% of the open nodes sampled by QPB weight, Knuth dives `k1.5`, 600 s per unit.
- **Coverage:** read part-way, with 524 of 552 units merged.
- **Result:** 1.87M nodes and 4,504 c7a CPU h. `docs/review_nov18_runs.md` read the 4,504 CPU h as a floor.

Scripts here, all reading the committed record (`python3 -I SCRIPT`):

| file | what |
|---|---|
| `load.py` | units, strata and dives, each tied to its unit (units.tsv's line ranges) |
| `strata_dives.py` → `strata_dives.txt` | strata and dives by depth, dive lengths, per-depth contributions, the 28 never-leased units |
| `tail.py` → `tail.txt` | the largest dives, the estimate without them, a delete-a-group jackknife over units |
| `sim_full.py` → `sim_full_comb.txt` | Knuth against "exact below depth D" at equal CPU on the comb arm's complete tho40ab subtree |

## 1. Not a truncation: depth 6 is where this tree ends

**The code does not cut dives short.**
- A Knuth dive in `latest/estimate.C` (`knuth_dive`) runs until the node it picked is not branched.
- There is no depth cap, and the per-unit budget is checked only between dives.

**The dives end early because the nodes stop branching.** Share of dives that ended one level below their stratum:

| stratum depth | dives ending one level below | also, by `strata_dives.txt` |
|---|---|---|
| 3 | 95.6% | 4.1% ended at +2, 0.2% at +3 |
| 4 | 99.3% | |
| 5 | 100% | |

- One dive reached depth 7, from a depth-6 stratum.

**Independent evidence agrees that tho40c's tree ends at depth 6:**
- **The bound lab's comb arm** solved every subtree below a sampled depth-4 node exactly. It saw nothing below depth 6 (`../bound_lab/tho40_dives/README.md`).
- **The comb arm's complete tho40ab subtree** has 877 nodes, the deepest at depth 6, with 0 of its 70 depth-6 nodes branched.
- **The pooled branch rate of a depth-6 node** is 4 in 1,464, or 0.27%. That pools:
  - the dives' depth-6 visits: 0 of 169;
  - tho40-est's depth-6 strata: 1 of 53;
  - its solve's done part: about 3 of 927;
  - the bound lab: 0 of 315.

**What is missing below depth 6 is small.** At the pooled rate, depths 7 and below add:

| | central | range (branch rate 0.1-0.9%) | share of the estimate |
|---|---|---|---|
| nodes | 94k | 34k-309k | +5% (2-17%) |
| CPU | 94 CPU h | 34-309 CPU h | +2% (1-7%) |

- The +2 incumbent adds the path to 240,516: about 36 branched levels, on the order of a thousand nodes.
- **So 4,504 CPU h is not a floor.** Its error is sampling error (section 3), mostly from depths 5 and 6.

## 2. The 28 never-leased units: negligible here

- **Size:** 585 open nodes, 3.3% of the 17,474 open, all at depths 5-6, all in the small lane.
- **Value:** merged strata put an open node at depth 5 at 1.6 subtree nodes and 12 CPU s, and at depth 6 at 1.6 and 8.
- **Effect:** they would add about 950 nodes and 1.7 CPU h, +0.05% and +0.04%.
- **Direction:** low, since they are left out entirely, not even counted as done.
- **Elsewhere it would matter:** they were left out by the queue starvation, not by design. In a tree whose small lane holds heavy nodes, the same starvation would not be negligible.

## 3. The heavy tail

**Where the estimate comes from.**
- **Depth-3 open nodes** (7,721 of the 16,889) carry 85% of the CPU: 1,738 s and 216 subtree nodes each.
- **Depth 4:** 210 s each; **depth 5:** 12 s; **depth 6:** 8 s.

**Seven dives carry half the nodes.** Of 7,165 dives:
- The seven largest are 52% of the nodes and about 30% of the CPU (`tail.txt`).
- All seven came from depth-3 strata and went to depth 6.
- Six of them were the only dive in their stratum.
- **Where they come from:** 4.4% of the depth-3 dives found a branched depth-4 child, and only 8 of those 168 found a branched depth-5 grandchild. Each of those 8 dives estimates about 40,000 nodes for its stratum.
- **By depth:** the depth-6 estimate (1.01M nodes) is 97% from ten dives. The depth-5 estimate (652k) is spread over 3,180 dives, with the top ten at 17%.
- **Dropping the largest dives:**

| dives dropped | nodes | CPU h |
|---|---|---|
| none | 1.87M | 4,504 |
| largest 1 | 1.68M | 4,326 |
| largest 3 | 1.36M | 3,874 |
| largest 10 | 0.80M | 3,015 |

**Why so few dives.** The setup bounding the sampled open nodes took 68.3 of the estimator's 94.6 CPU h. The dives themselves cost 13 s each on average and took 26.3 h.
- Setup used up each unit's 600 s, so the second pilot round was cut.
- 828 of the 1,108 dived strata got one dive. Unit wall times had a median of 607 s and a maximum of 3,263 s; 102 units ran over 900 s.

**The range.** `combine()` bootstraps dives within strata, and sampled strata as if drawn with replacement.

| | combine()'s bootstrap, 5-95% | delete-a-group jackknife over units, 30 groups |
|---|---|---|
| total CPU | 3,502-5,553 h | se 699 h (16%) |
| nodes | 1.25M-2.60M | se 0.44M (24%) |

The two agree.

**Standard error by depth:**

| depth | nodes se | CPU se |
|---|---|---|
| 4 | 8% | 8% |
| 5 | 12% | 19% |
| 6 | 40% | 37% |

- **Credibility:** the range is credible as the sampling error of these dives. The depth-6 part rests on about 8 events (Poisson ±35%).
- **What it can't show:** events rarer than about 1 in 7,000 dives. Section 1 bounds those to a few percent.
- **The depth-6 count may be high, not low.** Where branched depth-5 nodes were counted directly, the rate came out lower:

| source | branched depth-5 nodes |
|---|---|
| depth-4 strata's dives | 0.73% (22 of 3,011) |
| depth-5 open nodes | 1.7% |
| solve's done part | about 0.5% |
| implied by the estimate's depth-6 count | 4.6% |
| depth-3 strata's dives (the population that matters) | 4.8% (8 of 168), confidence interval 2-9% |

  The 8-of-168 figure is consistent with the estimate, so the estimate stays unbiased, with ±40% at depth 6.

## 4. The whole tho40c tree, and the next estimate's design

**Most likely: about 190 c7a CPU-days (4,600 CPU h), defensible range 140-250 (3,400-6,000).** Nodes are about 1.9M, range 1.0M-2.9M.

- **tho40-est:** 4,504 CPU h, or 188 CPU-days, with a jackknife se of 29 days.
- **Plus:** the deep tail (+2%) and the pending units (+0.04%).
- **The bound lab's comb arm agrees, independently:** 215 c7a CPU-days (se 82) and 1.87M nodes (se 0.74M), measured with full enumeration below depth 4 and incumbent 240,516.
- **Combined:** weighting the two by their inverse variances gives 191 ± 27.
- **The bound lab's same-path route (about 110) is the outlier.**
- **Wall time:**
  - about 19 h on four c7a.16xlarge (256 slots at about 95% busy), 14-25 h over the range;
  - about 38 h on two.

**For the next tho40 estimate, in order of value:**

1. **Spend the budget on dives, not setup.**
   - A dive costs 13 s, while bounding a sampled open node costs 30-100 s.
   - Give the setup its own allowance, and then guarantee every branched stratum ≥5 dives.
   - Today, 75% of the dived strata got one.
2. **Stratify the sample by depth.** QPB's g² weights barely separate open nodes under the DNN bound: the sampling probabilities were 0.17-0.27 at every depth.
   - Sample in proportion to each depth's CPU per open node: depth 3 heavily (50-100%: it holds 85% of the CPU), depth 4 at about 20%, depths 5-6 at a few percent (8-12 s each).
   - At 50% of depth 3 alone, the setup is about 3,860 nodes × 96 s ≈ 100 CPU h.
3. **Solve exactly below depth 4, with a node cap.** In tho40c a depth-4 subtree averages 13 nodes and 210 s.
   - Completing the subtree of the depth-4 node a dive reaches removes the depth-5/6 variance, which is 8 events today.
   - Fall back to Knuth when a subtree passes a cap: the bound lab saw single depth-4 subtrees of 5.7-13.3M CPU s, where phases stall.
   - **Cost:** about +10 s per depth-3 dive on average. Only 4.4% of the depth-4 picks branch, at about 35 × 9.7 s plus depth 6 each.
   - **On the comb subtree** (`sim_full_comb.txt`), with 50% of its CPU for each rule:
     - Knuth reads a median of 0.82 of the CPU, typically low;
     - exact below depth 3 reads 1.02;
     - its 10-90% band is about as wide, on a single 877-node subtree.
4. **A global budget, allocated by depth.** Today each unit gets 600 s whatever it holds. Allocate by each depth's share of the estimate's variance instead, Neyman on CPU per depth.
5. **Readout diagnostics, computed by the estimate itself:**
   - the top-1 and top-10 dives' shares by depth;
   - the share of one-dive strata;
   - a jackknife se by depth;
   - the deepest dive against the depth holding the most mass;
   - the pooled branch rate at the deepest level, giving a below-the-deepest-level tail as here.
