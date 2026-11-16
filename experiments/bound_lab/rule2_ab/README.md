# tho40: rule 2 against rule 9 at the bottom of the tree (Nov 16 2026)

Question (Nathan): `tho40.par` uses rule 9 through depth 6 and rule 2 below. Is that still right?

Two depth-5 subtrees of the A/B subtree (`26:10,8:0,2:32,18:39,35:24` = a, `26:10,8:0,2:32,6:39,35:24` = b;
404 nodes each in production, six branchings at depth 7 apiece), cold from `-F`, `-u 240516`,
`base.par` = `experiments/estimates/tho40/tho40.par`, `r9all.par` = the same with the third block
(rule 9) through depth 38, so the bottom block is never used. `run.sh`; four processes at once on a
shared box, so iterations are the metric. [m] all.

| subtree | arm | nodes | branchings at depth 7 (params level) | DNN iterations | CPU s |
|---|---|---|---|---|---|
| a | base (rule 2 below 6) | 404 | 6 (level 3) | 10,636 | 748 |
| a | rule 9 everywhere | 404 | 6 (level 2) | 10,583 | 745 |
| b | base | 404 | 6 (level 3) | 11,127 | 774 |
| b | rule 9 everywhere | 404 | 6 (level 2) | 10,925 | 762 |

Same trees; 0.5% and 1.8% fewer iterations with rule 9. Every depth-8 child fathoms in both arms, so
the rule can only change the children's iterations here. Depth 8 and below is 1.9% of the A/B
subtree's CPU and 7.5% of tho40-b's merged tree at 25 h (`experiments/tho40b_analysis/README.md`).
