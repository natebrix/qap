# eigcuts: eigenvalue cuts on the LP (Steve Wright's question, Nov 25 2026)

Steve's question: between the LP and the SDP, add curvature cuts v'Yv >= 0 to the LP for fixed directions v,
adaptively (the most negative eigenvectors of the LP's Y, by Lanczos or the power method), re-solving with the
dual simplex. Could that be the bound, a fast-fathom filter, or a guide for branching?

## Answer so far

**Not at our sizes.** The LP is the Adams-Johnson level-1 RLT LP (`../branching_lab/lpbound/`, Kurt Anstreicher's
question of Nov 20): from n = 18 it is below QPB, and PDLP needs 20-31 minutes for it at m = 52-53, as much as the
node's DNN bound. The cuts then tail off. On nug12's root (LP 522.9, DNN 568.0, optimum 578):

| cuts (10 per round) | bound | share of LP -> DNN | CPU s of the last re-solve |
|---|---|---|---|
| 0 | 522.9 | 0% | 108 (the LP itself, dual simplex) |
| 100 | 535.6 | 28% | 12 |
| 200 | 538.0 | 33% | 24 |
| 400 | 540.4 | 39% | 69 |
| 800 | 543.2 | 45% | 142 |

After 800 cuts and 5,400 CPU s, Y still has 38 negative eigenvalues (the least -0.12, from -0.37) and 529 of the
cuts are tight. 40 cuts per round is no better per cut: 1,536 cuts reach 542.2 (43%). The LP has 8,856 variables;
wil50's root LP has about 3.0 million, and a cut is a dense row over all of them. Times are from a loaded
container (four solver processes beside these two) and are rough; the cut counts are not affected.

Why it tails off here: the DNN optimum's dual PSD matrix has high rank (about 600 at wil50's size, R* about 1,500 of
2,117: `../wil50_feasibility/roadmap/admm/alternatives.md`), so at least that many rank-one cuts are needed even
with the right directions; the same fact closed the spectral bundle there.

## Files

- `kelley.py INSTANCE [K] [ROUNDS] [SECONDS]`: the loop, on `ajlp.py`'s LP with HiGHS (`pip install highspy`).
- `out/nug12_k10.txt`, `out/nug12_k40.txt`: a line per round.

Not run: nug15, nug20 (where the LP starts below QPB), sparsified cuts, dropping slack cuts.
