# Does a dive from tho40-b's incumbent recover the best known solution? (Nov 15)

tho40-b was given the best known value + 2 (240,518; tho40's objective values are
even), as the 2000 nug30 solve was given the optimum + 2. Its result means something
only if it finds a solution of 240,516 itself, and it can only if no node on the way
to one is fathomed: a bound above a subtree's optimum anywhere on that path loses the
solution. The cutoff is 240,516.05 (the incumbent less the even-objective tolerance
of 1.95), so a bound 0.05 too high on that path would do it, and the run would show it only
at its end.

`dive.py` follows the path to QAPLIB's `tho40.sln` (`QAPSOLN/`, added Nov 15) one node at
a time with tho40-b's recipe (`infra/phase1/recipes/tho40.par`, `SYM/tho40t.sym`,
`-u 240518`): each node is one `bnbdrv` process run to its own depth, so the solver
bounds it and chooses its line as the run would, and the next node is the child on
that line that agrees with the solution.
- **Warm** (`--warm`): the parent's state file, cut to the path child alone, is resumed
  (`-r`), so each node starts from its parent's shipped DNN state, as a coordinator unit does.
- **Cold**: each node is started by `-F` with the path so far.

On `tho40t.dat` the solution is `tho40.sln` inverted. The symmetry file keeps the grid's
180° turn of it, which puts facility 26 at location 2: tho40-b's depth-1 node `26:2`.

| | cold (`cold.tsv`) | warm, as the run (`warm.tsv`) |
|---|---|---|
| nodes on the path | 38, depths 0 to 37 | 38 |
| every node branched, path child kept | yes | yes |
| the leaf (depth 37, 3 facilities free, enumerated) | 240,516, the solution | 240,516, the solution |
| CPU on this cloud box (Xeon 2.1 GHz) | 634 s | 577 s |

**Two paths.** The root and `26:2` match tho40-b's own nodes: `26:2`'s DNN bound,
226,246.9108, and its line, location 0. At `30:0`, depth 2, the cold node branches on
location 39 and the warm one on location 32. Below that the two dives follow different
paths to the same solution.

**Comparison still to make.** tho40-b had not reached `30:0` by 17:05 UTC Nov 15. When
it does, its bound should be the warm dive's 231,440.5089, and its line location 32. A
second run of `dive.py --warm` reproduced every bound and line.

**The bound along the path** (240,516 less the bound):

| depth | 1 | 2 | 4 | 6 | 8 | 9 | 10 | 11 to 36 |
|---|---|---|---|---|---|---|---|---|
| cold | 14,269 | 10,052 | 7,009 | 3,978 | 2,301 | 1,125 | 57.7 | 0.00 to 0.10 |
| warm | 14,269 | 9,075 | 5,238 | 3,220 | 755 | 525 | 0.89 | 0.00 |

From depth 11 down, the DNN bound on both paths is 240,516.0000 to four places at
almost every node: the subtree's optimum exactly, 0.05 below the cutoff. Such a node
cannot be fathomed until the incumbent is 240,516 itself, after which the cutoff is
240,514.05 and it fathoms at once. That is the extra tree an incumbent of the best known + 2
costs, which the recipe put at about 0.1%. From depth 22 down, branching drops children
whose bounds clear the cutoff (2 of 18 at depth 22, 9 of 10 at depth 30), and the
path child was always kept.

**What this shows and what it doesn't.** The DNN bound is certified, so in exact arithmetic it cannot
exceed a subtree's optimum. This checks the certificate, the fathoming tolerance, the
branching filters and the symmetry restriction on 2 x 37 nodes where an error of 0.05
would have shown, and that the run's settings, given time, reach 240,516. It covers
two paths, not every path to every optimal solution.

```
python3 dive.py OUT --warm      # about 10 min on one core; prints PASS or FAIL, writes OUT/path.tsv
```

For the next instance, run it before the launch: `QAPSOLN/` now has `sko42.sln` and
`wil50.sln` too, and `dive.py` takes the instance, solution, recipe, symmetry file and
grid as options (a grid of r x c locations; a symmetry axis through the root's
location is not handled).

## With the next run's settings (Nov 16): PASS, both dives

`infra/phase1/recipes/tho40c.par` (tho40-b's settings with the bound lab's full block: the grid cut phase
at zone 0.06, pools handed to children, the hand-over rule, the phase penalties, `dnn.retry = 1`), on
cognac 52adc4d, from the same incumbent 240,518:

| | cold (`cold_tho40c.tsv`) | warm, as the run (`warm_tho40c.tsv`) |
|---|---|---|
| nodes on the path | 38, depths 0 to 37 | 38 |
| every node branched, path child kept | yes | yes |
| the leaf (depth 37, enumerated) | 240,516, the solution | 240,516, the solution |
| highest bound on the path | 240,516.0000 | 240,516.0000 |
| CPU on this cloud box | 2,798 s | 1,742 s |

**The bound along the warm path** (240,516 less the bound) against tho40-b's settings:

| depth | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| tho40-b | 14,269 | 9,075 | 6,871 | 5,238 | 3,891 | 3,220 | 1,716 | 755 | 525 | 0.89 |
| tho40-c | 7,426 | 5,644 | 4,177 | 3,374 | 2,263 | 1,380 | 1,100 | 327 | 0.00 | 0.00 |

- The root and `26:2` are tho40-b's nodes: the same QPB root, the same plain DNN bound at `26:2`
  (226,246.9108), the cut phase then 233,089.6108, and the same line (location 0). At `30:0` the warm
  node branches on location 32 as before; the cold one on 39, and the two paths part from depth 3.
- From depth 9 (warm) or 10 (cold) the bound is the subtree's optimum to four places, as with tho40-b's
  settings from depth 11: the cuts close the gap two levels sooner, and never pass the optimum.
- The path costs about three times tho40-b's CPU: every node on it runs a cut phase that cannot fathom
  it. Off the path the phases are what fathom nodes.

```
python3 dive.py OUT --warm --params ../infra/phase1/recipes/tho40c.par --resume
```

`--resume` (Nov 16) keeps the nodes OUT already holds, so an interrupted dive goes on where it stopped;
a depth-1 node takes about 11 minutes here with these settings.
