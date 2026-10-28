# Experiment: a certified DNN bound at the top of the tree

Offline studies, in Python, of a doubly nonnegative (DNN) SDP bound: the
facially reduced ZKRW / Oliveira-Wolkowicz-Xu relaxation, solved by ADMM. The
bound is certified, i.e. computed from a dual-feasible point, and every
reported value is valid. The code and full reports for each study are in its
directory. The scripts were run from a scratch directory, so their paths
need adjusting before they will run here.

- `nug20/`: the DNN bound fathoms 96 of the 98 depth-2 nodes and all 200
  sampled depth-3 nodes that QPB leaves open. But at 8 s and 4 s per node it
  is about 5x (depth 2) to 30x (depth 3) slower than letting QPB solve the
  subtree. Not worth it at n <= 24.
- `nug30/`: the DNN root bound is 5949 against QPB's 5332 (optimum 6124). At
  depth 2 it fathoms an estimated 82% of the 225 QPB-open nodes, at about
  130 s per node. A node's QPB value predicts which ones it misses, and
  those are the heavy subtrees. At depth 3 under the survivors it fathoms
  19 of 20. The solver's estimator puts QPB-only nug30 at 30-120 CPU-days
  on this machine. A DNN-at-depth-2-and-3 + QPB hybrid is estimated at
  2.5-5 CPU-days in Python, or 1-2 with a C/LAPACK implementation.
  The estimates are labelled in the report as measured, estimated or
  recalled.
