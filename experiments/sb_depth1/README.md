# DNN strong branching at depth 1 under sko42's recipe: a laptop probe (Nov 22)

[m] = measured, [e] = estimate.

**The question.** DNN strong branching (`sb.*`, PR 6) was measured on nug30 and tho30 under the nug30 record
recipe: no cut phase, no pools, no `dnn.retry` (`experiments/branching_lab/port` on `claude/cool-babbage-ilnogj`).
The open runs' recipes have all three (tho40c, `sko42.par`, `sko49.par`): a depth-1 node's cut phase hands its pool
to its children, so strong branching's child solves run pooled, and with `dnn.retry = 1` a failed one is solved again.
The PR left two things unmeasured there: how often a pooled child solve fails or is retried at `sb.iters = 150`, and
what a firing costs at this size. tho40-261122 runs it at depth 1 now, but without `QAP_SB_TRACE`, so its counts are
not visible.

**The probe.** sko42 is solved (sko42-261120b, `sko42.par`, 128,431 nodes, 539.3 CPU h), its node trace committed.
`probe.sh` runs one firing at each of its 12 depth-1 nodes: `bnbdrv -F 36:k -d 1` under `sko42_sb1.par`
(`sko42.par` plus `sb.depths = 1:1`, every other `sb.*` at its default): the node bounded as in the run, its line
chosen by strong branching, then stopped. `summary.py` sets each beside sko42-261120b: rule 9's line there and the
subtree it grew, strong branching's line, where rule 9's line ranks among the candidates, and the firing's child
solves (pooled or not), retries, failures and CPU.

It does not measure whether strong branching's line grows a smaller subtree: that takes solving the subtree under
both lines (14 to 103 c7a CPU h each but one, `36:0`).

**Its cost [e].** The depth-1 nodes' own bounds took 502 to 841 CPU s on c7a in sko42-261120b. A firing where no
symmetry is left (9 nodes: rule 9 chose c0 at 8, r22 at 36:0, 41 children each) bounds the children of 4 rows and 4
columns, about 300 solves; where symmetry is left (36:3, 36:10, 36:17: r22, 23 children) 4 rows, about 92. At the laptop's ~2.5x a c7a core (`../bound_lab/laptop`),
about 10 s a child: up to an hour a node, run side by side.

**A smoke test here (36:3, 15 minutes on the cloud container) [m].** The node's own bound as in the run: its first
solve gives up at 15,468 after 300 iterations (44 s), its cut phase reaches 15,615 (1,900 iterations, 342 s) and
hands on a pool (67,699 bytes). The firing's child solves run under that pool, 22 s each here, ending 15,656 to
15,693, above the node's own bound: 12 of them, none retried.

```
caffeinate -i experiments/sb_depth1/probe.sh
```

## Results (Nov 23, the laptop: 10 at a time, stopped once and resumed) [m]

`runs/summary.txt` (`summary.py`; a firing's cost from its progress file, which keeps every child solve across the
stop: the trace line of a resumed firing counts only its new solves):

| node | rule 9 | its subtree: nodes, CPU h | strong branching | rule 9's line: rank of candidates, score / best | the runner-up / best | child solves | pooled | retried | failed | firing CPU s |
|---|---|---|---|---|---|---|---|---|---|---|
| 36:0 | r22 | 242, 1.6 | c6 | 3 of 8, 1.07x | 1.022x | 312 (172 fathomed) | yes | 0 | 0 | 4,887 |
| 36:1 | c0 | 3,571, 15.1 | c0 (the same) | 1 of 8, 1.00x | 1.161x | 312 (47 fathomed) | yes | 1 | 0 | 6,730 |
| 36:2 | c0 | 6,898, 28.9 | c0 (the same) | 1 of 8, 1.00x | 1.115x | 312 (21 fathomed) | yes | 1 | 0 | 6,973 |
| 36:3 | r22 | 3,757, 17.1 | r22 (the same) | 1 of 4, 1.00x | 1.315x | 92 (0 fathomed) | yes | 0 | 0 | 1,913 |
| 36:7 | c0 | 3,110, 13.8 | c0 (the same) | 1 of 8, 1.00x | 1.129x | 312 (48 fathomed) | yes | 1 | 0 | 6,663 |
| 36:8 | c0 | 16,486, 69.0 | c0 (the same) | 1 of 8, 1.00x | 1.145x | 312 (10 fathomed) | yes | 1 | 0 | 7,009 |
| 36:9 | c0 | 25,340, 103.3 | c0 (the same) | 1 of 8, 1.00x | 1.079x | 312 (6 fathomed) | yes | 0 | 0 | 6,973 |
| 36:10 | r22 | 13,175, 48.8 | r22 (the same) | 1 of 4, 1.00x | 1.170x | 92 (0 fathomed) | yes | 0 | 0 | 2,145 |
| 36:14 | c0 | 4,994, 22.3 | c0 (the same) | 1 of 8, 1.00x | 1.063x | 312 (30 fathomed) | yes | 0 | 0 | 6,847 |
| 36:15 | c0 | 21,425, 86.8 | c0 (the same) | 1 of 8, 1.00x | 1.050x | 312 (6 fathomed) | yes | 0 | 0 | 6,979 |
| 36:16 | c0 | 20,131, 90.9 | c0 (the same) | 1 of 8, 1.00x | 1.001x | 312 (1 fathomed) | yes | 2 | 0 | 6,982 |
| 36:17 | r22 | 9,410, 41.8 | r22 (the same) | 1 of 4, 1.00x | 1.161x | 92 (0 fathomed) | yes | 2 | 0 | 2,225 |

- **Strong branching takes rule 9's line at 11 of the 12 nodes.** The one other line is at 36:0, the smallest
  subtree (242 nodes, 1.6 CPU h of 539), where rule 9's r22 scores 1.07x strong branching's c6. Elsewhere the
  runner-up is 0.1% (36:16: c0 against c41, a coin flip) to 31% behind. On sko42 rule 9's marginal score already
  picks the line the children's bounds pick.
- **Every firing ran under the node's pool** (each depth-1 node's cut phase handed one on, 62 to 70 kB), and the
  pooled children mostly end above the node's own bound.
- **The retry without the pool fired, and worked** (the path PR 6 had not seen run): 8 of 3,084 child solves were
  solved again without the pool (`retry 1` in the progress files), each 300 iterations in all and 84 to 104 CPU s,
  about 4x a child; none still failed after it.
- **The cost**: 312 child solves at each of the 9 nodes with no symmetry left (6,700 to 7,000 CPU s; 36:0 4,887,
  with 172 children fathomed inside their 150 iterations), 92 at the 3 with symmetry (about 2,000 s): 66,300 CPU s,
  18.4 CPU h for the 12 firings, about 22 s a child. The laptop under ten such processes was no faster than the cloud
  container (22 s a child there), not the ~2.5x a c7a core expected from `../bound_lab/laptop`: the firings took
  about 2.5 h, not one.
- **What it says [e].** On sko42 under `sko42.par`, strong branching at depth 1 would change one line in twelve, at
  the smallest subtree, for about 3% of the tree's CPU (18 laptop CPU h against 539 c7a CPU h; not the same machine).
  On tho30 the depth-1 lookahead changed lines and took 25% off the tree (`../branching_lab/grid` on
  `claude/cool-babbage-ilnogj`), on nug30 7 to 10%: its value is the instance's. For sko49, under the same recipe
  family, depth 1 is unlikely to pay; the root may (`../branching_lab/root`: r45 against production's r41,
  uncalibrated), at about 6 h of root at n = 49 with `dnn.depths` from 0.
- Not checked: the scores are the plain sum at depth 1 (`sb.power`'s default); whether the 8th power would agree.
