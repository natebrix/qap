# sko42, first attempt (November 8 2026): stopped after 5.5 hours, not finished

Skorin-Kapov's 6 x 7 grid instance (QAPLIB optimum 15,812, found by heuristics
and never proven), tried with the nug30 record recipe unchanged
(`../nug30_solve/README.md`), 8 workers on the M2 Pro laptop, binary at f60f77e.

## Files

`qaplib/sko42.dat` has the grid distances as the first matrix. The symmetry
files act on the second matrix (locations), so the search ran on
`qaplib/sko42t.dat`, the same instance with the two matrices swapped (B then
A), and a new `latest/SYM/sko42t.sym` built the way `nug30.sym` is for the
6 x 7 grid: the 12 cells of the middle column and the 24 cells of the three
left columns as the fixed sets, one pair of reflections.

```
12
1 2 3 4 8 9 10 11 15 16 17 18
1
6 24
4 11 18 25 32 39
1 2 3 4 8 9 10 11 15 16 17 18 22 23 24 25 29 30 31 32 36 37 38 39
```

Checked before the run: `bnbdrv -p PARAM/nug_deep.par -d 1 -u 15812
../qaplib/sko42t.dat` reads the file by name and gives 13 nodes, the root plus
the 12 Burnside orbits of the 42 locations under the grid's 4 symmetries.

## The run

```sh
# the nug30 record environment, then
python3 bench/parallel.py qaplib/sko42t.dat --jobs 8 --ub 15812 \
    --params PARAM/nug_deep_xbar_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint sko42t.ckpt --checkpoint-every 300
```

Started 01:31, stopped by hand at 07:03 (signal 15; the checkpoint written at
06:58 is intact, 2.4 GB, resume with `--resume`). Progress at the stop:

| | |
|---|---|
| nodes bounded | 9,651 (9,650 DNN, 9,318 DNN-fathomed, 96.6%) |
| open nodes | 3,168, still growing (1,539 at 3 h, 2,551 at 4.5 h) |
| tasks done | 124 |
| CPU | 32,988 s, summed over 8 workers (9.2 CPU-hours) |
| per node | 3.4 CPU s (0.35 s on tho30, 0.35 s on nug30) |
| log | `sko42t_j8_m2pro.log` |

The frontier was widening throughout, so no finish estimate can be made from
this run; what it says is that sko42 is not a laptop-night problem on the
nug30 recipe. Each node costs 10x a tho30 node (n = 42 against 30: the DNN
works on n^2 x n^2 matrices, and the nodes that survive to the give-up
threshold cost the full 200 iterations). The fathoming rate per node, 96.6%,
is tho30's (96.1%) and nug30's (96.3%); what differs is the depth the tree has
to reach and the price of each node.

Nothing here is tuned for n = 42: the DNN depths 1-7, beta 0.21 m/3, the 200
iteration give-up and rule 8 at every node are nug30's. If sko42 is attempted
again: rule 9 (26% on tho30, `../tho30_solve/README.md`), a root-line check
of the kind in `../wil50_feasibility/tuning/lines/root/`, and a machine with
more than 8 workers, resumed from the checkpoint rather than restarted.
