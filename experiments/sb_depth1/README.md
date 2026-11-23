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

## Results

(to come)
