# Region branching probe (Oct 30 2026)

Question (experiments/wil50_feasibility/NOTES.md, section 5): instead of
fixing a facility to one location, fix it to a block of the grid. Is the DNN
bound of such a region node strong enough to fathom whole blocks?

`region_probe.py` bounds, for one node and one free facility f, the region
node of every block (f forbidden at every free location outside the block,
via the DNN bound's eliminated assignments, `dnntest` node-file third field)
and every single-location child, with the solver's settings (`-H 300`).

## nug30, depth-2 far miss `0:13,6:14`, facility 21, 2x3 blocks (`nug30_0-13_6-14.txt`)

| block | cells | region lb | lowest cell lb | cells fathomed | region CPU | cells CPU |
|---|---|---|---|---|---|---|
| (0,0) | 6 | 6049.0 (gave up) | 6058.1 | 4/6 | 22 | 140 |
| (0,1) | 6 | 6097.1 (gave up) | 6097.0 | 5/6 | 26 | 58 |
| (1,0) | 4 | 6047.7 (gave up) | 6058.1 | 2/4 | 22 | 66 |
| (1,1) | 6 | 6072.2 (gave up) | 6081.0 | 4/6 | 22 | 71 |
| (2,0) | 3 | 6122.1, fathomed | 6122.1 | 3/3 | 93 | 87 |
| (2,1) | 3 | 6122.3, fathomed | 6122.6 | 3/3 | 12 | 20 |

(cutoff 6122.05; "gave up": the 300-iteration hopeless rule.)

The region bound is within about 10 of the lowest cell bound, and equal to it
where it converges: a region node is about as strong as its weakest cell,
which is the most a relaxation of "f in the block" can be. So region branching
does not move fathoming up a level; a block fathoms only when all its cells
would, and the gain is at most one DNN solve instead of one per cell for
blocks that fathom whole. Here that does not pay: region branching would cost
about 530 CPU s at this node against 440 s for cell branching. One node of
nug30; running the probe on wil50's depth-3 near misses would confirm it there.
