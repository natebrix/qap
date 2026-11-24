# wil50_pair: rule 10 against rule 9 at wil50's depth 2, on the laptop (Nov 24)

wil50's next-run file, `infra/phase1/recipes/wil50_d2r10.par`, branches depth-2 nodes by rule 10 (every line, the
size score on the children's QPB bounds). That has been measured on sko42 only (five depth-2 subtrees: 0.51 of rule
9's nodes, 0.70 of its CPU). This is the same test on wil50 itself, sized for a night on the laptop: a depth-2
node's whole subtree, solved under both arms side by side.

    caffeinate -i experiments/wil50_pair/run.sh

needs `latest/bnbdrv` built. It runs the two nodes below in turn, each under both arms at once, half the
performance cores an arm. Ctrl-C drains (a minute or two); the same command goes on from where it stopped and
skips what is finished. After each node it prints `summary.py`'s table (also `summary.txt`): each arm's line,
nodes and CPU, and rule 10's share of rule 9's. Each arm's record is in `r10/subtrees/NODE/` and
`r9/subtrees/NODE/` (`results.txt`, `trace.txt`, `summary.txt`), through `experiments/estimates/subtree.sh`.

## The arms

- `r10`: `r10/r10.par`, a copy of `infra/phase1/recipes/wil50_d2r10.par` as of Nov 24 (a copy so that a later
  edit of the recipe does not change this arm).
- `r9`: `r9/r9.par`, the same `[settings]` under `experiments/estimates/wil50/wil50.par`'s rule blocks: rule 9 at
  depth 2. The two files give every other depth the same block (checked when `wil50_d2r10.par` was written), so
  the arms differ in the line the depth-2 node is branched on and in nothing below it.

Both from an incumbent of 48,818, the next run's.

## The nodes

Four depth-2 nodes of the Nov 7 frontier (`experiments/wil50_feasibility/sampling`) were bounded and branched
alone under each arm (`lines.txt`; about 510 to 620 CPU s each on a cloud container):

| node | first solve, short of the cutoff | after its cut phase | rule 10's line | rule 9's line |
|---|---|---|---|---|
| `21:0,45:13` | 274 | 133 | r31 | c49 |
| `21:3,34:0` | 271 | 131 | r31 | c9 |
| `21:22,31:16` | 289 | 122 | r45 | r45 |
| `21:2,24:0` | 266 | 114 | c9 | c9 |

Where the two rules take the same line the arms are the same run, so the pair is the first two. Rule 9's lines
there are the ones the Nov 7 run took at these nodes, warm (`frontier_stats.tsv`: col 49, col 9).

## What it can and cannot say

- Where these nodes stand among wil50's depth-2 nodes under the new settings is not known. On Nov 7 (a narrower
  cut zone) the first three were 133 to 147 short after a phase and `21:2,24:0` 247 short without one; with
  `cut.zone = 0.012` every depth-2 node runs a phase, and all four end 114 to 133 short. The nodes that carried
  most of the Nov 7 frontier's cost were the ones without a phase, so that picture has changed and has not been
  measured again. The estimate run (`infra/phase1/specs/wil50-est.spec`) will show it.
- Two nodes of 304, chosen from four tried for the rules disagreeing. On the other two rule 10 takes rule 9's
  line: at such nodes the recipe's change does nothing.
- Each arm starts cold at its node (`-F`): no state from the depth-1 parent, so the node's own bound is not the
  full run's, though here rule 9 takes the same lines either way.
- Size: not measured. From the sampling study's subtrees, a few hundred to a few thousand DNN nodes an arm, at
  about 13 s a node on the laptop: one to ten CPU hours an arm [e], within a factor of three.

A result near sko42's (rule 10 well under rule 9 in nodes and CPU at both nodes) supports the recipe as it is;
rule 10 clearly worse at both is a reason to put depth 2 back on rule 9 before wil50's run, which
`experiments/estimates/wil50/wil50.par`'s blocks with `wil50_d2r10.par`'s settings would be (`r9/r9.par`).
