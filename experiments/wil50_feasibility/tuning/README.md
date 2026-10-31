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
