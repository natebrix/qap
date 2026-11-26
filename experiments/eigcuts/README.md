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

## Cuts taken once from the SDP's dual, kept under branching (`dualcuts.py`, Nov 26)

Jeff Linderoth's variant (Gunluk, Junger, Linderoth, Lodi, Luedtke, "Sparse Cuts for the Positive Semidefinite
Cone", arXiv 2603.09864; its section 4.3 calls these "dense dual SDP cuts"): solve the SDP once, take a cut
v'Yv >= 0 for each eigenvector of the dual PSD matrix S* with a positive eigenvalue, put them in the LP, branch,
and hope the geometry does not change much. nug12 (`out/nug12_dualcuts.txt`; the root's dual from Clarabel,
`out/nug12_rootdual.npz`; children's DNN bounds from `latest/dnntest`, run to convergence):

| root, LP 522.9, DNN 568.0 | bound | share of the gap |
|---|---|---|
| 10 cuts (largest eigenvalues of S*) | 544.8 | 49% |
| 20 | 556.6 | 75% |
| 40 | 567.9 | 99.8% |
| 63 (all above 1e-6 of the largest) | 568.0 | 100% |

So at the node the cuts come from, a modest number recovers the SDP: 20 of them do better than 800 of the
adaptive cuts above. Under branching, the 63 root cuts kept and the LP solved again at each child of row 0:

| | LP | LP + root cuts | DNN | share of the child's gap closed |
|---|---|---|---|---|
| depth 1, 12 children | 551.9 to 560.3 | 571.1 to 575.8 | 578.0 to 593.7 | 53 to 73%, mean 59% |
| depth 2 under 0:7, 11 children | 577.9 to 593.1 | 578.0 to 602.3 | 578.0 to 616.0 | 47 to 100%, mean 65% |

With the optimum 578 as incumbent the DNN fathoms all twelve depth-1 children; the LP with the root's cuts
fathoms none (it needs depth 2, where the plain LP already fathoms nine of the eleven). CPU a child: LP 1.1 s
(HiGHS IPM), LP + 63 cuts 24 to 29 s (Clarabel), DNN 0.5 to 2.7 s (our ADMM, to convergence). HiGHS could not
solve the LP with 40 or more of these dense rows: its IPM stopped short of optimal and its dual simplex did not
finish in 50 minutes.

At wil50's size the dual has rank about 600 (`../wil50_feasibility/roadmap/admm/alternatives.md`), each cut is
a dense row over 3.7 million variables (some 2 billion nonzeros against the LP's 10 to 15 million), the instance
is dense so the paper's sparsification has nothing to drop, and the LP without any cut already takes PDLP 20 to
31 minutes (`../branching_lab/lpbound`) against about 70 CPU s for a warm child DNN bound. Not pursued.
