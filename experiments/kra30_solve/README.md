# kra30a and kra30b, proven optimal (November 7 2026)

Krarup's hospital-layout instances, the nug30 record recipe unchanged
(`../nug30_solve/README.md`), 8 workers on the M2 Pro laptop, binary at b326767.
The data files are the scaled variants in `qaplib/` (flow and distance swapped,
distances divided by 5; see `RUNNING.md`, "Three data files are scaled
variants"), so the objectives are QAPLIB's 88,900 and 91,420 divided by 5, and
every bound here is one fifth of the literature's.

| | kra30a | kra30b |
|---|---|---|
| optimum (scaled / QAPLIB) | 17,780 / 88,900 | 18,284 / 91,420 |
| upper bound given | 17,780 | 18,284 |
| root QPB (scaled) | 13,595 (x5 = 67,976) | 13,674 (x5 = 68,370) |
| DNN root (scaled, from `../dnn_bound/`) | 17,367 (x5 = 86,837; gap 2.3%) | 17,572 (x5 = 87,858; gap 3.9%) |
| nodes | 198 (DNN: 197, 190 fathomed) | 3,446 (DNN: 3,445, 3,315 fathomed) |
| DNN CPU | 100.5 s | 1,251.6 s |
| CPU, summed over 8 workers | 101.4 s | 1,260.8 s |
| wall | 17.8 s (158 tasks, 10 steals, 5.69x) | 187.2 s, 3 min 7 s (333 tasks, 15 steals, 6.73x) |
| log | `kra30a_nug30recipe_j8_m2pro.log` | `kra30b_nug30recipe_j8_m2pro.log` |

Both searches closed with no open node and the given upper bound as the best
objective, so the QAPLIB values are optimal. Against nug30 on the same recipe the
same day (3,934 nodes, 1,370 CPU s, 196 s): kra30a is 20x smaller, kra30b 0.9x,
although kra30b's DNN root gap (3.9%) is the widest of the three. The DNN closes
kra's gap fast below the root, as the GLB/QPB history did not (QPB's root gap is
23-25% here against 5% on nug30; kra30b took the 2002 grid run years).

Recipe (from `RUNNING.md`, "The nug30 record setting"), with the instance and
upper bound changed:

```sh
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 \
       QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 \
       QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32 DNN_R32_NONCERT=1 \
       DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 \
       DNN_CORE_REFRESH_PARTIAL=1 DNN_BETA_SCHED=w=0.36
python3 bench/parallel.py qaplib/kra30a.dat --jobs 8 --ub 17780 \
    --params PARAM/nug_deep_xbar_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint kra30a.ckpt --checkpoint-every 300
# kra30b: qaplib/kra30b.dat --ub 18284
```

Not tuned: rule 8 at every node, the DNN depths 1-7, beta 0.21 m/3 and the rest
are nug30's values. The symmetry files `latest/SYM/kra30a.sym` and
`kra30b.sym` (11 and 6 root orbits) were picked up by name.

## The unscaled files and kra32 (November 8)

Run overnight on the laptop, the same recipe and binary (f60f77e), 8 workers.
None of the three has a symmetry file (the `.sym` lookup is by basename, so
`kra30a-old.dat` finds nothing), so these trees branch on all 30 or 32 lines
at the root and are not node-for-node with the scaled runs above; they prove
the same optima.

| | kra30a-old | kra30b-old | kra32 (scaled) |
|---|---|---|---|
| optimum proved | 88,900 | 91,420 | 17,740 (x5 = 88,700, QAPLIB's) |
| upper bound given | 88,900 | 91,420 | 17,740 |
| nodes | 873 (842 fathomed) | 28,112 (27,072 fathomed) | 10,852 (10,484 fathomed) |
| DNN CPU | 369 s | 17,929 s | 7,368 s |
| CPU, summed over 8 workers | 373 s | 18,063 s | 7,406 s |
| wall | 56 s (144 tasks, 6.67x) | 2,844 s, 47 min 24 s (811 tasks, 6.35x) | 1,065 s, 17 min 45 s (709 tasks, 6.95x) |
| vs the scaled run with symmetry | 4.4x nodes, 3.7x CPU | 8.2x nodes, 14x CPU | - |
| log | `kra30a-old_j8_m2pro.log` | `kra30b-old_j8_m2pro.log` | `kra32_j8_m2pro.log` |

- All three searches closed with no open node and the given upper bound as the
  best objective. The unscaled kra30a and kra30b give QAPLIB's 88,900 and 91,420,
  so the scaled files and their `/5` optima are consistent with the originals
  (the kra30a root anomaly of `../dnn_bound/cpp_report.md` was the scaling, not
  the solver).
- Symmetry is worth more than the orbit count suggests on kra30b: 8x the nodes
  and 14x the CPU without it, against 4.4x and 3.7x on kra30a. Anyone running the
  `-old` files should copy the `.sym` file to the `-old` basename first
  (`latest/SYM/kra30a-old.sym`), or pass `-s`.
- kra32 is the first 32-location instance proven here: 10,852 nodes and 2.1
  CPU-hours, 3x kra30b's scaled tree, with no symmetry file. `QAPSOLN/kra32.sln`
  said 17,780; fixed on Nov 8 to 17,740 (88,700/5) with an optimal permutation from a
  tabu search, checked against both data files.
