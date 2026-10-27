# Experiment: per-node diagonal shift A + Diag(d)

`diag_shift.patch` (against 437ccae) adds environment-switched code to try idea 1
of the "bigger ideas" report: before a node's bound, replace the reduced A by
A + Diag(d) and C by C - d diag(B)'. Every permutation keeps its objective, but
the QPB changes. d comes from a few supergradient steps on the projected
eigenvalue bound (PB), which is concave in d. Not merged: it trades about 10%
fewer nodes for 12-20% more time at nug20.

Apply with `git apply experiments/diag_shift/diag_shift.patch`, then
`make clean && make`, and run with e.g.
`QAP_SHIFT_DEPTH=99 QAP_SHIFT_STEPS=2 ./bnbdrv -p PARAM/default.par -u 2570 ../qaplib/nug20.dat`.
The switches are listed at the top of `latest/shiftexp.C`. The default is
off, which matches 437ccae.

## Results (ub = optimum, PARAM/default.par, idea 3 included)

The shift raises QPB a lot at the root (PB proxy, 5 steps):

| instance | root QPB | shifted |
|---|---|---|
| nug12 | 478.5 | 498.5 |
| nug15 | 986.6 | 1006.4 |
| nug16a | 1441.7 | 1458.8 |

Below the root the gain is much smaller. With QAP_SHIFT_GUARD=1, each
unfathomed node is bounded with and without the shift and keeps the better
bound. The second FW run alone (STEPS=0) is the control:

| instance, nodes shifted | mean QPB gain, shift | control | shift effect |
|---|---|---|---|
| nug15, all (m <= 15) | +3.0 | +0.8 | about +2 |
| nug16a, all | +1.6 | +0.8 | about +1 |
| nug20, depth <= 4 (m >= 16) | +6.4 | +1.5 | about +5 |

Full trees, shift applied before the bound (no guard), 2 ascent steps:

| run | nodes | CPU (3-4 jobs on 4 cores) |
|---|---|---|
| nug20, no shift | 565,399 | ~158 s |
| nug20, all nodes | 502,028 (-11%) | 191 s |
| nug20, m >= 14 | 511,428 (-9.5%) | 177 s |
| nug20, m >= 17 | 550,659 (-2.6%) | 166 s |
| nug16a, all | 9,749 (+2%) | |
| nug17, all | 48,667 (+2%) | |

Other variants:
- 5 steps instead of 2: same nodes, slower.
- 1 step: worse than no shift (nug16a 10,576, nug17 52,570).
- Ascent on the eigenvalue term alone (QAP_SHIFT_MODE=0): about the same as PB.
- One step along diag(X B X') at the FW solution (MODE=2): +2 at nug15, worse.
- Children inherit the branched parent's d without re-optimizing
  (QAP_SHIFT_INHERIT=1): no gain (nodes +-3%).
- Shifting B instead of A (QAP_SHIFT_WHICH=1): +0.8, no gain.

Where the time goes: each ascent step costs an eigendecomposition of V'(A+D)V
and a LAP, and the B eigendecomposition is done twice. Reusing those in init_ST
would remove about half of the overhead, which is roughly break-even at nug20.
The per-node gain grows with m (about +2 at m = 12, +5 at m = 16, +20 at
m = 20 at the root), so it may pay at n >= 24. That has not been measured.
