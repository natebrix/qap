# ADMM research: synthesis (Nov 3)

Three memos, written by research agents from the general ADMM and SDP
literature (Boyd et al., OSQP, SCS, COSMO, SDPAD, SDPNAL+, spectral bundle,
Burer-Monteiro, the QAP-DNN papers of Wolkowicz and coauthors), checked
against our code and logs. Web access to arxiv and most publisher pages was
blocked for the agents; citations rest on abstracts, excerpts and, for SCS
and OSQP, their source code. Figures marked [m] are measured, [e] estimated.

- `alternatives.md`: replace the PRSM? **No.** Every published solver for
  this relaxation at n = 49-50 runs cold to 1e-5 in 4,700-15,300 iterations
  of 1.3-1.7 s [m]; ours runs warm-started at 0.33 s and stops at 150-1,500.
  Spectral bundle, Burer-Monteiro, SDPLR and CGAL are ruled out by the
  ranks we measured (primal about 1,500 of 2,117, dual about 600), and
  SDPNAL+ costs 5-25 of our iterations per Newton step [e].
- `projection.md`: the PSD projection is about 74% of an iteration and
  within 1.25-1.45x of the cheapest projection at today's accuracy [e];
  looser accuracy has always cost iterations in our runs.
- `splitting.md`: (gamma1, gamma) = (1, 1) is already the largest dual step
  the convergence theory allows, and equilibration and chordal splitting do
  not apply. What we lack is state management: the penalty beta is fixed,
  the child warm start carries only Z, and the eigensolver accuracy is the
  same all run.

## Candidates worth a test, ranked

| # | idea | from | expected [e] | effort | first test, and what kills it |
|---|---|---|---|---|---|
| 1 | Warm-start children with the parent's Y as well as Z (today `dnn_warm_child` maps Z only; checked) | splitting | **measured Nov 3 (`../../warm_y/README.md`): -24% iterations on the 7 shared fathoms (1.29x CPU), +67 at iteration 100 on the 9 misses, no outcome changes; passes** | done; in the search as `QAP_DNN_WARM_Y=2` | in the stacked node-10 A/B (`../../stronger_bounds/README.md` step 13: bundle D and parent Y together 1.67x in the search); laptop nug30 search 25% with one-pass CholQR |
| 2 | Graham et al.'s restricted PRSM pieces we never built: trace-constrained R-step, doubly stochastic Y-step row (listed "not done" in `experiments/eigensolver/README.md`; checked) | alternatives | **measured Nov 5 (`../../speed/rprsm/README.md`), killed**: trace R-step -1% iterations; DS x-row + arrow (`DNN_RPRSM=10`) all32 cold -10.2% iterations in total, -9.1% median, same fathoms, lb@200 +33 on misses, but -4% CPU; warm children -6.5% iterations, +1 to +4 at 200 | done | opt-in `DNN_RPRSM` bits in `dnn.C`, off by default |
| 3 | R's rank term in fp32 except on certificate iterations (the certificate reads only its own iteration's R; checked) | projection | 4% per iteration [m: 409 -> 395 ms] | 0.5 day | margin back at 5e-6 on certificate iterations and identical fathoms on all32 |
| 4 | Two-phase or residual-balanced beta, changed only at refreshes, frozen after 200 | splitting | 0-15% on wil50; pres/(beta dY) is 0.32-0.45 early and 0.12 late [m] | 2 days | all32: lb at 200 up by under 20 and fathoms under 10 iterations earlier kills it |
| 5 | Complementary-subspace Rayleigh-Ritz and per-column filter degree | projection | 25-55 ms per iteration | 3-5 days | microbenchmark first: under 10 ms saved kills it |

Closed by the memos with reasons: spectral bundle, Burer-Monteiro, SDPLR,
CGAL and SketchyCGAL, SDPNAL+, LOBPCG, Davidson, FEAST, randomized
projections, equilibration, chordal and multi-block splitting, further
over-relaxation, and variants of Anderson, Halpern and momentum. Items 1-4
are a few days together; each is a 32-node run of about 15 minutes, as for
the AMX bundles.
