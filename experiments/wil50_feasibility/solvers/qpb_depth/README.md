# QPB down the wil50 tree: a spot check (laptop, Nov 6)

Question: could the search hand deep wil50 nodes to the QPB bound, as the
nug30 search does from depth 7 or 8? Known before this check: QPB is about
1,870 short at the depth-3 sample nodes (`../../tuning/README.md` section
1) and the Nov 2 run_ab2 arm, which limited the DNN to four levels and left
the rest to QPB, ran past 10^7 QPB nodes under each depth-7 miss
(`../../stronger_bounds/README.md`, step 9). Not known: how fast QPB's gap
closes with depth.

`cmds.txt` through the console (`latest/bnbdrv -c -u 48816 -p
PARAM/nug_deep.par ../qaplib/wil50.dat < cmds.txt`, output
`console_out.txt`; the console's fix prompts read stdin, so `-a` does not
work for this) on one path of the node-10 subtree: the nodes the Nov 6
`fast` arm bounded (`../../pipeline_ab/runs/fast/trace.txt`) down to depth
6, then arbitrary legal extensions of the depth-6 node (`dnn_d*.log`:
`dnntest` fathoms them in 60-110 iterations). Cutoff 48,814.05.

| depth | node | free | parametric QPB | short | DNN at the same node |
|---|---|---|---|---|---|
| 3 | node 10, `39:27,30:48,25:49` | 47 | 46,833.8 | 1,980 | 48,600.3 at 300 it (given up), 48,655.8 converged |
| 4 | `8:45,...` | 46 | 46,996.8 | 1,817 | fathomed, 48,815.5 at 140 it |
| 4 | `49:45,...` | 46 | 46,927.0 | 1,887 | hopeless at 300 it, 48,763.6 |
| 5 | `9:45,19:0,...` | 45 | 47,101.4 | 1,713 | fathomed, 48,814.4 at 150 it |
| 6 | `0:26,9:47,...,49:45` | 44 | 47,268.5 | 1,546 | fathomed, 48,814.1 at 70 it |
| 7 | + `1:0` | 43 | 47,570.9 | 1,243 | fathomed, 48,830.2 at 110 it |
| 8 | + `2:1` | 42 | 47,725.8 | 1,088 | fathomed, 48,835.5 at 80 it |
| 10 | + `3:2,4:3` | 40 | 48,096.4 | 718 | fathomed, 48,931.9 at 60 it |

(QPB without the parametric step is 20-110 lower at each node; GLB is
39,000-42,400.)

Reading:

- **QPB closes 100-300 units per level and is still 718 short at depth 10,
  where the DNN fathoms in 60 iterations.** At the depth-6 node, where the
  DNN fathoms at once, QPB is 1,546 short; on this path's rate it would
  fathom around depth 13 or 14, seven or eight levels and 40-47 children
  per level below the DNN's fathom: 10^5-10^7 QPB nodes under one depth-6
  DNN fathom, which is the Nov 2 observation again (10^7 and more under a
  depth-7 miss), now with the gap profile behind it. And this is an easy
  path: the extensions are arbitrary assignments, which make a subproblem
  easier than the search's marginal lines do, so the gaps along real
  search paths will be at least these.
- **The difference from nug30 is the scale of the gap, not the shape.** At
  nug30 the QPB root is 12% below the optimum and the DNN 3%; the search
  hands nodes to QPB from depth 7 or 8 because by then QPB's gap is a few
  units and its nodes cost milliseconds. At wil50 QPB is 4% short at depth
  3 and 1.5% at depth 10; the DNN is 0.3% at depth 3 and fathoms from depth
  4 on. One wil50 unit is also worth less: 130 units is a level here
  (`../../NOTES.md`), so 718 short at depth 10 is still about five levels.
- **So there is no depth in the wil50 tree where handing over to QPB
  pays**, and the rule that every wil50 node is a DNN node (`QAP_DNN_EXT`,
  step 9) stands. The one QPB use that survives is the lookahead inside
  the branching rules, which costs a few QPB solves per node and does not
  need the bound to be tight.
