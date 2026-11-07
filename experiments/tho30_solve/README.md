# tho30, proven optimal (November 7 2026)

Thonemann and Bölte's 3 x 10 instance (QAPLIB optimum 149,936), the nug30
record recipe unchanged (`../nug30_solve/README.md`), 8 workers on the M2 Pro
laptop, binary at afe298b, `qaplib/tho30.dat` as distributed (flow first, the
grid second; `latest/SYM/tho30.sym`, 10 root orbits, picked up by name).

| | tho30 | kra30b (same day) | nug30 (same day) |
|---|---|---|---|
| optimum proved | 149,936 | 18,284 | 6,124 |
| root QPB | 123,702 (17.5% gap) | 13,674 | - |
| DNN root (from `../dnn_bound/`) | 143,575 (4.2% gap) | 17,572 (3.9%) | 5,950 (2.8%) |
| nodes | 18,410 (DNN: 18,403, 17,680 fathomed) | 3,446 | 3,934 |
| DNN CPU | 6,354 s | 1,252 s | 1,358 s |
| CPU, summed over 8 workers | 6,401 s | 1,261 s | 1,370 s |
| wall | 937 s, 15 min 37 s (525 tasks, 11 steals, 6.83x) | 3 min 7 s | 3 min 16 s |
| log | `tho30_nug30recipe_j8_m2pro.log` | | |

The search closed with no open node and the given upper bound as the best
objective, so 149,936 is optimal. 4.7x nug30's tree and CPU for a root gap
1.5x as wide in relative terms, in line with the exponential tree model;
kra30b (3.9% gap, 0.9x nug30) is the outlier on the other side. tho30 took
the 2002 grid run about 15 CPU-years at ~34e9 QPB nodes [unverified from
memory; see Anstreicher, Brixius, Goux, Linderoth 2002]; this is 1.8 CPU-hours.

Recipe (from `RUNNING.md`, "The nug30 record setting"), with the instance and
upper bound changed:

```sh
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 \
       QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 \
       QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32 DNN_R32_NONCERT=1 \
       DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 \
       DNN_CORE_REFRESH_PARTIAL=1 DNN_BETA_SCHED=w=0.36
python3 bench/parallel.py qaplib/tho30.dat --jobs 8 --ub 149936 \
    --params PARAM/nug_deep_xbar_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint tho30.ckpt --checkpoint-every 300
```

Not tuned: rule 8 at every node, the DNN depths 1-7, beta 0.21 m/3 and the rest
are nug30's values. With a tree this size the cut phase at depth 2 (break-even on
nug30) and rule 9 are the first things to try if a faster tho30 is wanted.
