# The DNN give-up and Intel AMX, on 48 nodes tho40-b bounded (Nov 15-16)

Two questions on nodes from tho40-b's part-way export (21:19 UTC Nov 15, 1.8e6 nodes):
- **The give-up:** would a node that gave up fathom if its DNN bound ran on?
- **Per-core speed:** what a DNN iteration costs here, with and without Intel AMX
  (`dnn.amx = 1`), against tho40-b's c7a.

`probe.py sample` drew 16 give-up nodes (DNN status hopeless, branched) at each of depths 4, 5
and 6 (`sample.json`, seed 1). `probe.py run` rebuilt each as the run bounded it: a
chain of `bnbdrv` processes from the root down its branching path, each run to its own
depth and resumed from its parent's state file cut to the next node, so each node starts
from its parent's shipped DNN state. Then it ran three arms on the node itself:

| arm | settings |
|---|---|
| A0 | the run's (`infra/phase1/recipes/tho40.par`) |
| A1 | A0 with `dnn.amx = 1` |
| G | A0 with `dnn.giveup = 0`: never give up; the bound stops at its fathom or its 4,000 iterations |

**The copies are faithful.** A0 reproduces the run's DNN line: the same status on 48 of 48 nodes and
the same iteration count on 44. Bounds differ by 0.08 units at the median (0.9 at the
90th percentile), from the run's fp32 parent states.

## The give-up: running on does not pay

| depth | nodes | fathom by 400 iterations | by 1,000 | by 4,000 | net per give-up at L = 400 / 1,000 / 4,000 (c7a s) |
|---|---|---|---|---|---|
| 4 | 16 | 0 | 0 | 0 | -22 / -104 / -514 |
| 5 | 16 | 0 | 1 | 1 | -15 / -80 / -422 |
| 6 | 16 | 0 | 2 | 2 | -11 / -62 / -327 |

Only 3 of the 48 fathom within 4,000 iterations (at 465, 570 and 930), and their subtrees
were small. A limit L in place of the give-up loses CPU at every depth, even L = 400
(`probe.py report`: savings are the run's subtree CPU below a node that fathoms by L,
costs the extra iterations at the node's own c7a rate).

So the 40% of CPU under give-ups whose children all fathomed at once (`../giveups.py`) is
not a give-up called too soon. The parents' DNN bounds converge below the cutoff, and the
children fathom because the branching fixes an assignment.

**The other way may pay.** In the export, give-up nodes are 3% of the nodes but 30% of the
DNN CPU (1.92e6 of 6.41e6 c7a s), because they run 200 iterations, or about 330 within
`dnn.gu_zone` of the cutoff. Their iterations past the 100th are 19% of all the CPU. The nodes
that fathom rarely need that long: past depth 4, at most 1.6% fathom after iteration 100.
A give-up decided at about 100 iterations, at a cost of those late fathoms branching,
is the next test: an A/B of `dnn.giveup` (100 against 200) and `dnn.gu_zone` on complete
tho40-b subtrees.

## Per-core speed

Median ms per DNN iteration, same nodes and same iterations:

| depth | c7a in tho40-b (64 workers) | this box (Sapphire Rapids, 4 cores), amx 0 | amx 1 |
|---|---|---|---|
| 4 | 158.6 | 114.3 | 89.5 |
| 5 | 140.3 | 99.0 | 81.0 |
| 6 | 117.0 | 86.5 | 69.1 |

- **This box is faster than tho40-b's c7a:** 1.38x (10%-90%: 1.29-1.49) without AMX.
- **AMX adds 1.26x** per iteration (1.14-1.33). It leaves the status unchanged on 48 of
  48 nodes and the iteration count on 45; bounds move by 0.09 units at the median.
- **On EC2** (`infra/phase0/runs/20261115-2317-eu-north-1`, 32 wil50 nodes, every core
  loaded), a c7i (AMX) core is 1.81x a c7a core. But a c7i.4xlarge has 8 cores to the
  c7a.4xlarge's 16 at the same Spot price, so the c7a is 12% cheaper per node.
- **Load is not the cause:** loading every core slows the c7a by only 1.14x.
- **Nor is the BLAS:** the same OpenBLAS (scipy-openblas32 0.3.34) gives the Zen 4 AVX-512
  kernels too; the Zen 4 runs AVX-512 at half width.

## Running it

```
python3 probe.py sample EXPORT sample.json --per-depth 16 --depths 4,5,6 --seed 1
python3 probe.py run sample.json WORK --workers 4
python3 probe.py run sample.json WORK --workers 4 --arms G --max-jobs 4   # one batch, within a 10-minute command
python3 probe.py report sample.json WORK
```

`run` keeps finished work and skips it, and deletes each 12 MB state once it is spent.
In a cloud session, run the arms in batches in the foreground: the container is
reclaimed when the session is idle, and background jobs go with it.
