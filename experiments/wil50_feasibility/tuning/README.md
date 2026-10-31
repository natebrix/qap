# wil50: cheap CPU levers before any cluster run (started Oct 30 2026, M2 Pro)

Two experiments from the discussion after the 100-core estimate
(../NOTES.md section 9). Both were cut short when the laptop session ended;
scripts and partial results are here so the cloud session can finish them.
All runs use the solver's mode (subspace tracker) with the Peaceman-Rachford
step and degree-4 filter: `dnntest -T 1 -d 4 -g 1.0 -G 1.0`.

## 1. Eliminated assignments as a bound-strength lever (`elim/`)

Aim: the DNN bound with eliminated assignments (DNNParams::elim, commit
9862021) is a strictly stronger relaxation of a node. If cheap eliminations
exist at wil50 depth 3 to 5, they move fathoming up and divide the tree.

Learned so far:

- **QPB's U matrix eliminates nothing at wil50 depth 3.** Console (`b`, `u`)
  on the 5 near-miss and 2 fathomed depth-3 sample nodes: parametric QPB
  46,928 to 47,081 against the cutoff 48,814 (gap about 1,870), largest U
  entry 222 to 282. Zero cells of 2,209 satisfy bound + U >= cutoff on any
  node (`elim/console_out.txt`, `elim/elims.txt`). The QPB route is closed.
- **The DNN certificate's own child bounds are the candidate instead.**
  `dnntest -C 2` prints, from the certificate's outer assignment problem,
  a certified bound for every child (lb + reduced cost), with no extra SDP
  solve. Children above the cutoff are exactly the assignments that can be
  eliminated. `elim/round.sh <id> <spec>` bounds a node with `-C 2`, forbids
  the children the certificate fathoms (`-E`), bounds again, and reports both
  bounds and both fathomed-children counts. It was started on the 7 nodes
  (`elim/order.txt`) with `-H 0 -i 2500` and killed before the first round
  converged; `elim/r1_*.log` are partial.

Open: how much does round 2 gain on the near misses (59 to 116 short), how
many more children fall, and does it iterate? A gain of 60 on a near miss is
half a level; the division-by-47 table in ../NOTES.md section 5 prices that.

## 2. Give-up threshold under PRSM (`giveup/`)

Aim: clear misses are most of the nodes and each costs the full give-up
budget (`-H 300`). PRSM's bound at 300 iterations is 60 better than ADMM's,
so the hopeless rule might fire at 150 or 200 with the same decisions, a
1.3 to 1.5x cut on the dominant node class.

Run: all 32 sample nodes (`giveup/all32.txt`: 16 depth-3, 16 depth-4
children of nodes 9 and 6) with `-H 150`, `-H 200`, `-H 300`, 8 at a time
(so CPU per node is about 2x the alone figures). Results in
`giveup/H150_*.txt`, `giveup/H200_*.txt`; the H=300 pass and 6 of the H=200
nodes did not run.

Learned so far, and it matters:

- **H=150 is unsafe.** It abandons 9 nodes that fathom: depth-3 node 12
  (given up at 170 iterations, lb 48,703; fathoms at 280 iterations with
  H=200) and depth-4 nodes 17, 18, 22, 24, 27, 29, 31, 32 (given up at 150
  with lb 48,644 to 48,732; 17 and 18 fathom at 210 and 270 iterations with
  H=200, and all 16 depth-4 children fathomed in the original -f/ADMM sample
  with H=300). A false give-up costs n-d children, about 46 nodes, so one
  in four is ruinous.
- **H=200 is unproven.** It fathoms everything that finished (node 12, and
  depth-4 nodes 17 to 21, 23, 25, 26, 28, 30), but the six slowest depth-4
  fathoms (22, 24, 27, 29, 31, 32) were still running when the session ended.
- The hopeless extrapolation was tuned on ADMM trajectories. Under PRSM the
  certified bound at 150 iterations sits 100 to 170 below its final value on
  nodes that fathom, so the rule's "recent gains" extrapolation is too
  pessimistic early. Either keep H=300 with PRSM, or re-fit the extrapolation.

## Suggested completion (cloud)

1. Finish `giveup/` at H=200 and H=300 for the 32 nodes (alone or 4 at a
   time); report false give-ups and CPU per class. If H=200 has any false
   give-up, look at the hopeless rule's extrapolation under PRSM rather
   than at H.
2. Run `elim/round.sh` on the 7 nodes to convergence; report round-1 and
   round-2 bounds, children fathomed by the certificate in each round, and
   whether a third round still gains. If it gains, the next step is wiring
   `child_lb`-based eliminations into bnb.C's DNN call (children the
   certificate fathoms need no bound at all, and the surviving children can
   inherit the parent's elimination set).
3. Only then the line-comparison experiment on wil50 depth-2 nodes
   (../NOTES.md discussion; port of experiments/dnn_bound/predictor/lines).


# Results (cloud, Oct 31 2026: Xeon 2.1 GHz, 4 cores; solver mode with PRSM and degree 4)

All runs `dnntest -T 1 -d 4 -g 1.0 -G 1.0 ... -u 48816` (cutoff 48,814.05),
at most 4 at a time. CPU figures are for this machine and are not comparable
with the laptop's.

## 1. Give-up threshold (`giveup/cloud/`, `run_cloud.sh`)

All 32 nodes of `all32.txt` (id = line number) with `-H 300` and `-H 200`,
`-i 2500`. Logs `cloud/H<H>_<id>.log`.

| | H = 300 | H = 200 |
|---|---|---|
| fathomed | 18 (depth-3 nodes 11, 12; all 16 depth-4 nodes) | the same 18, same iterations and bounds |
| CPU per fathomed node (mean) | 105 s | 105 s |
| given up | 14 (all other depth-3 nodes) | the same 14 |
| CPU per given-up node (mean) | 175 s | 95 s |
| certified lb at the give-up point | 48,359-48,752 (62-455 short; the near misses 6, 9, 14, 15, 16 at 360-480 iterations, 62-130 short) | 48,279-48,648 (166-535 short, all at 200 iterations) |
| false give-ups | - | **0** |
| total CPU, 32 nodes | 4,329 s | 3,205 s (-26%) |

**H = 200 makes no false give-up on this sample, but its margin is thin.**
`-v` traces of the four slowest fathoms (nodes 22, 24, 27, 29: 340-400
iterations; `cloud/trace/`) replayed through the rule (`replay.py`): it would
fire at iterations 130-180 on all four, so they are safe only from H = 180-190,
10-20 iterations below 200. And at H = 200 the near misses happen to give up
exactly at 200; at H = 220 or 250 they run on to 230-440 (`rules.py`). Varying
the rule's factor (5 = current, 10, 20, 40) does not help: anything that
protects the fathoms earlier than 180 stops giving up on near misses at all.

**Why: under PRSM every node's bound is still rising steeply before about 180
iterations** (the far misses gain about 1,300 between iterations 100 and 200),
so an extrapolation of recent gains cannot tell fathoms from misses that early.

**The gap to the cutoff separates them.** The gap of the running best shrinks
by a factor of 0.34-0.55 every 50 iterations on all four fathoms, right up to
the fathom; on the misses the 50-iteration ratio is 0.58-0.77 at iteration 200
and rises steadily toward 0.9-0.98 (`gaprule.py`). A rule "from iteration H,
give up when g(t)/g(t-50) > rho" with H = 150, rho = 0.65 never fires on the
fathoms (worst ratio 0.55) and gives up on the four traced misses at 190-230.
It is in the solver as an opt-in (`DNNParams::hopeless_gap_ratio`, `dnntest
-R rho`; off by default), and was run on the 32 nodes (`run_cloud_R.sh`,
`cloud/H150R065_<id>.log`):

%%RRESULT%%

## 2. Eliminations from the DNN certificate (`elim/`, `run_cloud.sh`)

`round.sh` on the 7 nodes of `order.txt`, `-H 0 -i 2500 -C 2`
(`r1_<id>.log`, `summary_<id>.txt`):

| node | round-1 certified lb | short of the cutoff | children the certificate fathoms | round 2 |
|---|---|---|---|---|
| 6 | 48,758.2 | 55.8 | 0 of 2,209 | nothing to eliminate |
| 9 | 48,770.6 | 43.4 | 0 of 2,209 | nothing to eliminate |
| 14 | 48,730.5 | 83.5 | 0 of 2,209 | nothing to eliminate |
| 15 | 48,753.7 | 60.4 | 0 of 2,209 | nothing to eliminate |
| 16 | 48,758.8 | 55.3 | 0 of 2,209 | nothing to eliminate |
| 11 | 48,815.4 | fathomed (290 iterations) | - | - |
| 12 | 48,814.9 | fathomed (280 iterations) | - | - |

The five near misses run to the 2,500-iteration limit with the primal and the
certified value within 2 of each other, so these are the converged DNN values:
43 to 84 short (the earlier ADMM sample had them 59 to 116 short at 720-800
iterations; PRSM to convergence gains 15-30). The certificate's child bounds
(lb + reduced cost of the outer assignment problem) are equal to lb for 90%
of the children; the largest gain is 2.8 to 6.5, against the 43 to 84 a child
needs; no child is within 38 of the cutoff. At a converged multiplier the
outer assignment problem is almost flat, as on nug30
(`experiments/dnn_bound/predictor`): the certificate cannot eliminate
assignments, and round 2 never ran. The answer to the question is no: this
adds nothing to the near misses.

## 3. Wiring child_lb eliminations into bnb.C

Not worth doing, given part 2: no child of any of these nodes is eliminated.
For the record, it would be: pass `child_bounds = 1` in bnb.C's DNN call
when the node is not fathomed; drop the children whose `child_lb` exceeds the
cutoff (no bound needed); and give each surviving child the parent's
elimination set (stored on the node, as the warm start's Z is) plus its own,
passed to the DNN call as `DNNParams::elim` (commit 9862021). The machinery
exists; the eliminations do not.
