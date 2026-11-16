- **The answers.** Both searches start at the optimum + 1, report `New objective value` = the optimum
  (283315445 in `runs/bfs35/d14_c*/log.txt`, 637250948 in `runs/bfs40/d16_c*/log.txt`), end with every process at
  exit 0 and no `node` line in the last level's state files. The values of the `QAPSOLN` permutations were
  recomputed from the `.dat` files and equal them.
- **The path check** (`collect.py`, in `logs/solve_bfs35.txt` and `logs/solve_bfs40.txt` with every node): a node
  all of whose fixes agree with the optimal permutation must not be fathomed, and no DNN or cut bound there may
  exceed the optimum. tai35b: 32 nodes on the path (depths 1 to 32), all branched or enumerated, largest
  certified bound 283,315,444.9994. tai40b: 37 nodes, largest 637,250,947.9995. **0 violations.** The roots:
  280,789,962 and 636,110,765.
- The check is sharp where it matters: from depth 13 (tai35b) and 20 (tai40b) the bounds on the path are within
  10 of the optimum, most within 0.1, and never above it, through 30 more solves each with an inherited pool.
@FILE:validity_second.md@
- Not a check of every cut: a cut that removed a solution worse than the optimum only would change node
  counts, not answers, and would not show here. What these runs exclude is a cut that removes the optimum on
  these two instances with these regions.
