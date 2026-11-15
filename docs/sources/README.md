# Primary sources

The papers and the thesis this code descends from. Read these before reasoning from
memory about the solver's history, its bound, its branching or its estimator.

| source | where | what it holds |
|---|---|---|
| N. W. Brixius, *Solving Large-Scale Quadratic Assignment Problems*, PhD thesis, University of Iowa, 2000 | `nb_thesis.pdf` here, `nb_thesis.txt` (pdftotext, for grep); also https://nathanbrixius.com/papers/nb_thesis.pdf | the branch-and-bound code's design (ch. 3: the bounds, branching rules, gap-dependent strategies), the tree estimator (ch. 4: Knuth's procedure, importance sampling by g^q, breadth-first top levels, tuning by estimate), the nug30 solve's statistics (ch. 5, Table 5.12, the source of `experiments/nug30_solve/nug30_2000_levels.txt`) |
| K. M. Anstreicher, N. W. Brixius, J.-P. Goux, J. Linderoth, "Solving large quadratic assignment problems on computational grids", *Mathematical Programming* 91 (2002) 563-588, doi:10.1007/s101070100255 | not in this repository (publisher's copyright) | the 2000 nug30 solve: the grid (Condor, MW), the estimator's use for tuning, the incumbent at the optimum + 2 |
| K. M. Anstreicher, N. W. Brixius, "A new bound for the quadratic assignment problem based on convex quadratic programming", *Mathematical Programming* 89 (2001) 341-357, doi:10.1007/s101070100214 | not in this repository (publisher's copyright) | QPB, the bound every node gets first |
| D. E. Knuth, "Estimating the efficiency of backtrack programs", *Mathematics of Computation* 29 (1975) 121-136, doi:10.1090/S0025-5718-1975-0373371-6 | not in this repository | the random-dive estimator the thesis's ch. 4 builds on |
| QAPLIB instances and solutions | `../../qaplib/`, `../../QAPSOLN/` | the instances and the published best solutions (no tho40 solution file) |
| this project's DNN-bound paper | `../dnn_bound/dnn_bound.tex`, `.pdf` | the DNN bound as implemented in `latest/`, and the 2000 vs 2026 nug30 comparison |
| the literature the cut work read | `../../experiments/wil50_feasibility/new_cuts/literature/papers/` (paper.pdf + text.txt each) | RLT, symmetric RLT, box and boolean-quadric cuts, recent QAP bounds |
