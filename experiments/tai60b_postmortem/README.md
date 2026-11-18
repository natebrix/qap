# tai60b-261117: where the tree went, and why the estimate missed it (Nov 18 2026)

The run closed the tree at 608,215,054 in 331,149 nodes and 3,597 m7a core-hours, against an estimate
of 29,357 nodes and 279 cloud-box core-hours (`docs/next_runs.md`). From the website's public tree
shards (`chains.py`, `chains.txt`). The collected node trace, committed Nov 18
(`infra/phase1/runs/tai60b-261117/nodes.txt.gz`), has the same 331,149 nodes and 3,597 CPU h.

## The tree is two chains and the subtrees hanging off them [m]

The root branches on r44. Two of its 60 children hold 99.95% of the tree: `44:45` (the best known
solution's path, 155,549 nodes) and `44:58` (its mirror, the facilities of locations 43 and 44
exchanged, 175,311 nodes; no solution under the cutoff, so it dies out by depth 46). Below each, one
child per depth carries on (the chain) and the chain node's other children are fathomed or hold
subtrees. The two chains take the same lines through depth 20.

| chain depth | rule | line | off-chain children branched (both chains) | their subtrees: nodes | % of the tree | CPU h | % of the CPU | deepest |
|---|---|---|---|---|---|---|---|---|
| 2 | 9 | r1 | 10 | 10,659 | 3.2 | 226 | 6.3 | 12 |
| 3 | 9 | r6 | 14 | 83,111 | 25.1 | 1,588 | **44.1** | 15 |
| 4 | 9 | r7 | 12 | 3,310 | 1.0 | 68 | 1.9 | 12 |
| 6 | 9 | r10 | 4 | 13,211 | 4.0 | 187 | 5.2 | 18 |
| 7 | **2** | c54 | **18** | 114,669 | **34.6** | 873 | **24.3** | 24 |
| 8 | 2 | c2 | 12 | 49,429 | 14.9 | 332 | 9.2 | 24 |
| 9 | 2 | r8 | 6 | 13,461 | 4.1 | 96 | 2.7 | 23 |
| 11 | 2 | r2 | 6 | 15,046 | 4.5 | 74 | 2.1 | 33 |

(Every chain depth with at least 0.2% of the CPU is in `chains.txt`.)

- **Two branchings decide two thirds of the cost**: the chain node's line at depth 3 (rule 9, r6) leaves
  7 children per chain that branch, whose subtrees take 44% of the CPU; at depth 7 (the first depth of
  rule 2, the bottom block) c54 leaves 9 per chain, 35% of the nodes and 24% of the CPU. c54 is the
  U-sum tie-break's location 54 (`../rule9_bottom/README.md`: location 54 at all 97 depth-7 branchings).
- The per-depth humps of the run page (depths 9-11 and 20-22) are these subtrees: the depth-3 ones
  reach depth 15, the depth-7 and depth-8 ones depth 24.
- Few nodes have such subtrees: the 16 chain nodes at depths 1-8 (49,000 to 175,000 nodes each) and
  the 32 off-chain children that branch at depths 3 and 7 (about 6,000 nodes each). A branching rule that spends minutes or hours deciding a line there
  (`../line_signals/README.md`: a DNN lookahead over 8 lines at 150 iterations, about 18 CPU-hours at
  n = 60) costs a fraction of a percent of what the line decides.

## Why the estimate missed [m]

The estimate (`experiments/bound_lab/tai60b/README.md` on `claude/cool-babbage-ilnogj`) measured the
path node by node and estimated the subtrees off it by Knuth dives from the path's off-path children:
- **284 dives in all, 2 to 8 per set, and none went more than two levels below its start.** The subtrees
  that hold the cost are 12 (depth 3) and 17 (depth 7) levels deep. A dive reaches their bottom only if
  it takes a surviving child at every level: 1 in 6 to 8 of the chain node's children at depths 3 and 7,
  and 1 to 7% of the nodes below (the run page's fathom rates by depth), and nothing in the dives biased
  them toward the survivors. Knuth's warning in thesis section 4.3: "almost always low, except
  for rare occasions when they will be much too high."
- **It estimated another tree.** Its path took `12:0` at depth 3 (its line at depth 2 was r12); the run
  took r1 there, `1:33` at depth 3 and `12:0` at depth 15. Its dives ran the bound lab's build and
  `kidW.par` (the stall rule on), not `recipes/tai60b.par`. Shares by depth do not compare one to one.
- E1 (the path alone) was 15,207 nodes; each chain of the run is 10 to 11 times that. E2 counted the
  mirror as a copy of the path, which it nearly is (175,311 against 155,549).

## What this says for the next estimates

1. Estimate the production tree: the same binary, params file and incumbent as the run.
2. The dives must reach the bottom of the near-miss subtrees: many more dives, a strong bias toward
   surviving children, and a check that the deepest levels are sampled (`../../docs/estimates_plan.md`).
3. The line at the few heavy nodes near the top is where settings can save most: `docs/estimates_plan.md`, 1.
