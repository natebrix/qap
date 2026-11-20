# lpbound: the LP bound without the PSD condition (Kurt Anstreicher's suggestion, Nov 20 2026)

**In progress.** Early results [m]: the LP (the DNN relaxation without R >= 0) is exactly the Adams-Johnson RLT-1
LP (checked numerically, `t0check.py`); it reproduces the published RRD95 values (nug12 523, nug15 1041). At tai60b's
depth-17 node 22:31 (m = 43) PDLP solves it to 1e-4 in 520 s wall, 0.9 GB: 541.9M, against QPB 527.7M, DNN 608.0M,
cutoff 608.2M, i.e. far closer to QPB than to DNN there. Line ranking tests running.
