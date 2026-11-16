Three estimates, from the same runs:

- **E1, the sampling design alone**: depths 0 and 1 measured (the root, the on-path child, 17 random others
  expanded 49/17), everything below from the dives: the on-path child's subtree by spine level, the sampled
  branched nodes' subtrees from `LB*`.
- **E2, post-stratified by the landscape**: the two depth-1 children that hold a solution within 1M of the best
  known value (`31:0`, not sampled, and `31:32`) are counted as copies of the on-path child's subtree (the twin
  probe of section 5 supports it for `31:32`); `31:43` and `31:42` get `31:43`'s dives; the other 45 are
  expanded 45/15 from the 15 sampled. E2 is above E1 because the dives under `31:32` (few, and blind to its
  path) saw none of what the spine-stratified dives see under `31:7`.
- **E3, the web model priced with measured costs** (`web.py`, `webcost.py`, `logs/web.txt`, `logs/webcost.txt`,
  section 5's rule: a node is open iff it holds a solution within g(depth) of the best known value). The rule
  against the @NCHECK@ distinct nodes bounded here (`webcheck.py`, `logs/webcheck.txt`):

@WEBCHECK@

  Below half the spine's gap every node branched; between 0.5 and 1.5 about half do, most of them near misses
  at the phase limit; above twice the gap one node in @ABOVE2@ did (a depth-1 near miss). So gap scale 1.0 is the
  model to read, 0.5 a floor, and 2.0 is contradicted by these nodes (it is in the table as the far case).

@WEBCOST@

`estimate.py`'s output (E1, E2, the on-path child's subtree T_A, and E2 by depth):
