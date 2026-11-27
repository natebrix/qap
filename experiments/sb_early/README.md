# Strong branching: leaving losing lines early, and using the chosen line's solves; the cut projector's stores

Nov 27, 2026. Three changes from an outside review of the solver (its items 1b, 1a and 2a), on
`claude/cool-babbage-ilnogj` only, each behind a default that changes nothing in a run that does not ask for it
(the third changes no result at all). `make test` passes on this container, where the coordinator's two
AWS-mock scenarios are skipped (no boto3); `tests/sb_test.py` has cases for both keys.

## `sb.early = 1 | 2`: a candidate line is left once it cannot win

A line's score is a sum of nonnegative terms (shortfall^P over its children), so once its partial sum passes the
best complete line's it is out, and its other children are not solved. The line chosen is the one every solve
would choose. Lines are taken in decreasing marginal score; within a line the cells already known first, then
the children by decreasing x_ij of the node's marginals (1) or by increasing QPB bound (2).

On wil50, replayed from the child solves of `../wil50_sb` (`cached.sh`: the firing's progress file under the new
settings hash, so no child is solved again and the count of children needed is exact):

| node | every child (Nov 25-26) | `sb.early = 1` | `sb.early = 2` | line |
|---|---|---|---|---|
| root (50 rows, P = 8) | 750 solves, 14.7 CPU h | 192 | 194 | r45, as before |
| 21:10 (4 rows, 4 columns, P = 1) | 376 solves, 6.7 CPU h | 226 | 213 | r45, as before |

So about a quarter of the root's firing and 57 to 60% of a depth-1 firing. On the test cases (`tests/sb_test.py`):
tai15b 322 child solves to 107, nug18 448 to 418, the same lines and the same trees.

What it costs: the trace no longer has every line's full score (a line left early shows `>` and its partial
score), and scores are over the largest shortfall a child of the node can have, not the largest found. The
analysis of `../wil50_sb` (which line a tail-weighted score would take) needs the full table: `sb.early = 0`.

## `sb.reuse = 1`: the chosen line's children are not bounded twice from nothing

A child of the chosen line that the firing fathomed is not created (it counts with the children eliminated at
branching); the others start from the bound they reached in the firing. Both are certified bounds. The tree can
differ where a child's higher starting bound puts it in another block of the params file. On the test cases:
tai15b proves its optimum on 16 nodes for 52, nug18 on 20 for 138 (most children of those small nodes fathom
inside the firing's 150 iterations). Not measured on a large instance: at wil50's depth 1 the chosen line has
0 to 22 fathomed children of 49 (`../wil50_sb`), each worth at most its own first solve.

## The cut projector stores Y once an iteration (`latest/dnncut.C`, side 1)

With `cut.cp_sweeps = 3:1` the projector's clipping step ran four times a CP-ADMM iteration and stored the
full matrix each time (36 MB at m = 47); between two of them only the group sums are read. Now only the last
one stores. wil50's depth-1 node 21:0 under `wil50_d2r10_g150xl2.par` (`2a/`: first solve of 300 iterations,
then a cut phase of 1,860), cognac's binary against this one, side by side on this container (4 cores):

| at once | cut phase CPU s, before | after | ratio |
|---|---|---|---|
| 2 processes | 943.3 | 926.1 | 0.982 |
| 4 processes | 955.3, 957.6 | 937.9, 936.1 | 0.980 |

The same bound (48,434.7559), rounds, iterations and pool to the byte in all six runs; the first solve's CPU is
unchanged (118 to 130 s in both). About 2% of the cut phases' CPU, so about 1% of a wil50 run. Not measured on
c7a with 64 processes, where memory bandwidth is scarcer.

## Follow-ups

- **The children's states (the rest of the review's 1a), to measure separately.** A surviving child of the
  chosen line could start from its state after the firing's `sb.iters` iterations instead of the node's. Not
  low risk: it changes every such child's solve (its give-up test, its pool, its cut phase) and so the tree,
  and it needs the states of a line kept until the winner is known (about 72 MB a child at m = 47) and handed
  through `WarmStart`. To measure: paired subtrees with and without, as the tuning rounds do.
- `sb.early` makes strong branching at depth 2 (`sb.gap`) and at the root cheaper than when they were left
  off for their cost; neither is re-examined here.
- Not on cognac. The recipes do not set either key.
