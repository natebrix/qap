# Can wil50 (or sko49) be solved with the certified DNN bound? Notes from the laptop (Oct 29 2026)

Research notes for the cloud session. Written on the M2 Pro laptop (Accelerate
build, `VECLIB_MAXIMUM_THREADS=1`, see RUNNING.md) after the nug30 replication.
The depth-3 sample below is the one already quoted in docs/dnn_bound (commit
7f106c6); the depth-4 follow-up, the tree estimate, the "lever" arithmetic and
the literature summary are new. Raw logs are in `depth3/` and `depth4/`
(`nodes.txt` lists the node specs, `run.sh` the exact commands).

## 1. Summary

- wil50 at depth 3 is a bound-strength problem, not a speed problem: the DNN
  relaxation is 60 to 500 short of the cutoff on 14 of 16 random depth-3 nodes.
- One more assignment buys 60 to 80: all 16 depth-4 children of two near-miss
  depth-3 nodes are fathomed, with certified margins of 0.1 to 5.4 above the cutoff.
- Tree estimate: about 4e7 DNN nodes, 1e10 CPU-s on this laptop. This agrees with
  Fujii et al. (ZIB 21-02, Table 5: 2.1e7 to 3.9e7 for wil50). About 95% of the
  work is the depth-5 descendants of the depth-3 "clear misses"; the near misses
  are about 1%.
- The lever is bound strength at depths 3 to 5, worth about 0.27% per level
  (about 130 on a cutoff of 48,814). Each level by which fathoming moves up divides
  the tree by about 47. Two levels (about 0.55%) would make wil50 a one-week laptop
  run; one level is about nine months.
- Nothing in the literature supplies that 0.5%: DNN is the strongest relaxation
  anyone has computed on these instances, triangle inequalities do not help (your
  `experiments/dnn_cuts` result), RLT3 is 10 TB at n = 50, and symmetry is already
  spent (order 4). Section 5 has the details and two cheap experiments that would
  pin the estimate down.

## 2. Depth-3 sample (16 random nodes, seed 20261029)

`dnntest -f -T 1 -H 300 -i 800 -u 48816 -x <spec> qaplib/wil50.dat`, 8 at a
time, cutoff 48,814.05 (best known 48,816, even objective, tolerance 1.95).
Root DNN iteration: 530 to 600 ms (20x nug30's m = 28 node); depth 3 (m = 47):
about 0.9 s per iteration alone, about 0.55 to 0.75 s when 8 run concurrently
on the shared AMX units. Summary (full lines in `depth3/node*.log`):

| outcome | nodes | certified lb below the cutoff by | iterations | CPU s |
|---|---|---|---|---|
| fathomed | 2 (nodes 11, 12) | (above by 1.1, 1.8) | 440 | 391, 393 |
| converged, missed | 5 (nodes 6, 9, 14, 15, 16) | 59, 74, 79, 83, 116 | 720 to 800 | 439 to 696 |
| given up at 300 iterations | 9 | 300 to 515 (at 300 iterations, not converged) | 300 | 275 to 283 |

The five converged misses have primal and dual values within a few units of
each other, so they are not iteration-limited: the DNN bound at depth 3 is
simply 0.1 to 0.2% below the cutoff on them.

## 3. Depth-4 follow-up (16 children of the two closest misses, seed 20261030)

Parents: node 9 (`26:27,8:38,34:28`, 59 short) and node 6
(`20:12,40:28,17:25`, 74 short). One random free facility per parent (37 and
43), eight random free locations each. Same `dnntest` settings, cold start.

| | iterations | CPU s | certified margin above the cutoff |
|---|---|---|---|
| 16 of 16 fathomed | 170 to 690 (median about 300) | 145 to 534 (mean about 270) | 0.1 to 5.4 |

The small margins are the early stop working (the certificate clears the cutoff
and the node stops), so the gain from one assignment is at least 60 to 80 here,
and from the root's gap of about 680 (root DNN 1.4% below the optimum) down to the
depth-3 median, the gap closes by roughly 130 per level.

## 4. Tree estimate and where the work is

5x10 grid, symmetry of order 4, branching on a facility with n - d children:
about 3e4 nodes at depth 3, 1.4e6 at depth 4, 6.4e7 at depth 5. Applying the
depth-3 sample's proportions:

| depth-3 class | share | what happens below | DNN nodes |
|---|---|---|---|
| fathomed at depth 3 | 12% | nothing | 3e4 (the depth-3 nodes themselves) |
| near miss (converged, up to about 120 short) | 31% | all children fathom at depth 4 | about 4.5e5 |
| clear miss (given up, 300 to 515 short at 300 iterations) | 56% | children still 170 to 370 short at depth 4; grandchildren mostly fathom at depth 5, a tail to depth 6 | about 8e5 at depth 4, about 3.7e7 at depth 5 |

Total about 4e7 DNN nodes. At about 250 CPU s per node (depths 3 to 5 cost
about the same per iteration, 300 to 700 iterations) that is 1e10 CPU-s: about
300 CPU-years, 40 years on this laptop at 8 workers, four months on 1,000 cores,
three weeks on the 5,184 Fujii et al. used. Fujii et al. estimate 2.1e7 to 3.9e7
nodes (Table 5, their three branching rules); two independent methods agree, and
the agreement is the reason they wrote "could be solved" and have not.

A correction to what I said after the depth-3 sample: the near misses are not
the lever. Their subtrees end at depth 4 and are about 1% of the work; 95% is
below the clear misses. So a stronger bound only helps if it moves every node's
fathoming depth up by a level.

The biggest uncertainty is the clear misses' converged values: 300 to 515 short
is at 300 iterations with `-H 300`, and the near misses gained a lot between
iteration 300 and 800. If the clear misses converge to 150 to 300 short, the
estimate above holds; if they stay above 400, more of the tree goes to depth 6
and the total is closer to 1e8 to 1e9 (the sko49 end of Fujii's table).

## 5. The lever: bound strength per level

Fathoming depth moves up one level for every 130 (0.27%) of extra bound at
depths 3 to 5, and each level divides the tree by the branching factor, about 47:

| extra bound strength | DNN nodes | CPU-s | laptop, 8 workers |
|---|---|---|---|
| none (today) | 4e7 | 1e10 | about 40 years |
| about 0.27% (one level) | about 8e5 | 2e8 | about 9 months |
| about 0.55% (two levels) | about 2e4 | 5e6 | about 1 week |

A 3 to 5x faster node (the eigensolver ideas in section 6) multiplies into any
row but does not change the row. For reference, on nug30 the plain SDP to DNN
step (adding Y >= 0) closed about half the gap; a single constraint family has
mattered that much before, which is why cuts were the natural thing to try.

### What the literature says about getting 0.5%

- DNN is the strongest relaxation anyone has computed on these instances.
  Povh and Rendl (2009) show their copositive-derived model is equivalent to the
  strongest of the earlier SDP models (Zhao, Karisch, Rendl, Wolkowicz; Burer and
  Vandenbussche); nothing computed since is stronger at n = 50.
- Triangle inequalities: your `experiments/dnn_cuts` says +0.2 and +0.4 on nug30
  depth-2 misses with gaps 76 and 29, nothing at the root, violations at most
  0.02. That matches what I would expect from the literature: with the gangster
  constraints and Y >= 0 the DNN relaxation nearly implies them.
- Level-2 RLT (Adams, Guignard, Hahn, Hightower 2007): the literature-survey
  agent read its dual-ascent root bounds on the nug instances as slightly below
  the DNN values (nug30 about 5934 against DNN about 5950). I have not verified
  that number against the paper; please check before spending effort on RLT2.
  docs/dnn_bound currently names level-2 RLT as the next thing to try, so this
  matters.
- Level-3 RLT (Hahn, Roth, Saltzman, Guignard 2013) is much tighter but the
  coefficient matrix is O(S^8): they quote 144 GB at S = 30, 1.6 TB at S = 40 and
  10.1 TB at S = 50, and 3 to 4 times that for branch and bound. Not a path.
- DNN plus cuts at n >= 30 is otherwise untested in print (Kim and Kojima's
  successive cutting inequalities, arXiv 2107.08665, are for general binary QPs,
  not QAP at this size).
- Symmetry (de Klerk and Sotirov 2009) helps when a data matrix has a large,
  ideally transitive, automorphism group. The 5x10 grid's group has order 4 and
  is already used by the tree. No extra reduction available for wil50 or sko49.
- Region branching (branch on which block of the grid a facility goes to, then
  refine) has no precedent I could find for QAP bounds. It would change the
  tree's shape at the top; whether the DNN bound is stronger on such nodes is
  unknown.

## 6. Eigenproblem: what is left per node

The literature-survey agent on the ADMM eigenproblem found no fundamental fix
(the projection onto the PSD cone of a (K x K) matrix, K = (m-1)^2 + 1, is the
cost and has to happen every iteration). Its shortlist, each a constant factor:

1. Split `prof[1]` to see how much of the eigen step is the subspace update
   against the periodic `ssyevd` refresh (measure first; `dnntest -v`).
2. Inertia-gated refresh: a cheap `ssytrf` inertia count tells when the tracked
   negative subspace has changed dimension, so the full refresh runs only then.
3. Side-switching subspace tracking: track whichever of the positive or negative
   eigenspace is smaller (near convergence the positive side is tiny), plus a
   parent-basis warm start of the subspace when the parent's DNN bound missed
   (your idea 4a).
4. Peaceman-Rachford (symmetric dual update) in place of ADMM, as in Graham et al.
   2022; fewer iterations if it converges.
5. Elimination-aware Y-step: exploit that eliminated rows and columns are zero.

Compounded, 3 to 5x per node at best. Worth doing for nug30-class instances
and for the sample experiments below, but it does not change the wil50 verdict.

## 7. Two cheap experiments that would sharpen this

1. Converge three of the depth-3 clear misses (nodes 1, 2, 3 say) with
   `-H 0 -i 2500`, about 40 minutes each on the laptop, to replace "300 to 515
   short at 300 iterations" with converged values. This decides whether the
   total is 4e7 or 1e8 or more.
2. Score any bound-strengthening prototype (cuts, a different dual update, a
   stronger relaxation) on the 16 depth-3 nodes in `depth3/nodes.txt` rather
   than on nug30: the question is whether it adds 130 at wil50 depth 3, and
   nug30 depth 2 cannot answer that.

## 8. Not done, available on request

- H=200 against H=300 A/B on nug30 (the laptop's 20-minute run makes this cheap).
- A region-branching probe: bound a few wil50 nodes that fix a facility to a
  2x5 block rather than a cell and compare the DNN values with cell-level nodes.
