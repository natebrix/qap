/*
============================================================
 DNN bound: a certified lower bound for a QAP subproblem

     min  tr(A X B X') + <C, X> + shift   over m x m permutations X

 from the doubly nonnegative (DNN) SDP relaxation of Zhao, Karisch,
 Rendl and Wolkowicz in the facially reduced form of Oliveira,
 Wolkowicz and Xu (OWX 2018), solved by the OWX ADMM:

     min <L, Y>  s.t.  Y = V R V',  R psd,  Y_00 = 1,
                       Y_J = 0 (gangster),  0 <= Y <= 1

 with Y of size m^2+1 on [1; vec(X)] and V the (m^2+1) x ((m-1)^2+1)
 orthonormal basis of the row/column-sum feasible subspace. This is a
 port of experiments/dnn_bound/nug20/dnn.py (see that file and the
 study reports for the derivation); dnn.C lists what differs (faster
 products, eigensolver and certificate, smaller step). With
 prototype_cert = 1, poly = BOX, eig_mode = DOUBLE, beta = 1 and
 check_every = 50 it reproduces dnn.py's certified values (dnntest -P).

 Only the certified value `lb` may be used as a bound. It is the
 Lagrangian dual value of a multiplier Z' whose projection V'Z'V is
 negative semidefinite, checked by a Cholesky factorization (see
 dnn.C, "Certificate"). The ADMM primal value `primal` is not a bound.

 Typical use, for a node reduced with reduceQAP (assign.h):

     DNNParams par;                    // defaults: see below
     DNNResult r = dnn_bound(qr->A, qr->B, qr->C, qr->shift,
                             incumbent - DELTA, par);
     if (r.status == DNN_FATHOMED) ...  // r.lb > cutoff: prune
     node.bound = max(node.bound, r.lb) // r.lb is valid in every case

 Memory: about 100 m^4 bytes (75 MB at m = 28), allocated per call.

 dnn_bound has no global state and allocates its own workspace, so it
 is reentrant (one call per thread; set the BLAS to one thread, see
 dnn_set_blas_threads). For fixed input and parameters the result is
 bitwise reproducible on the same machine and BLAS, unless the CPU
 time limit stops the run.
============================================================
*/

#ifndef DNN_H
#define DNN_H

#include <vector>
#include "meschach.h"
#undef catch

enum DNNStatus {
  DNN_FATHOMED = 0,   // certified lb > cutoff
  DNN_CONVERGED,      // ADMM residuals below tol
  DNN_STALLED,        // certified lb stopped improving and is near the primal value
  DNN_HOPELESS,       // certified lb extrapolates to below the cutoff (heuristic)
  DNN_ITER_LIMIT,
  DNN_TIME_LIMIT,
  DNN_ERROR           // LAPACK failure or bad input; lb is still valid (maybe -inf)
};

const char *dnn_status_name(int status);

// Which dual-bound evaluation is used for the certified value.
enum DNNPoly {
  DNN_POLY_BOX = 0,   // prototype: (L+Z)_00 + sum of negative non-gangster entries
  DNN_POLY_GL  = 1    // Gilmore-Lawler-type evaluation (see dnn.C); always >= BOX
};

// How the inner (m-1) x (m-1) assignment problems of DNN_POLY_GL are solved.
enum DNNInnerLAP {
  DNN_LAP_REDUCE = 0, // row then column reduction: O(m^2) each
  DNN_LAP_EXACT  = 1, // Hungarian method: O(m^3) each
  DNN_LAP_AUTO   = 2  // REDUCE, plus EXACT when the result is close to the cutoff
};

// Eigensolver of the R-step (the projection onto the psd cone).
enum DNNEig {
  DNN_EIG_DOUBLE = 0,   // dsyevd, as dnn.py
  DNN_EIG_SINGLE = 1,   // ssyevd: about 2x faster, same iteration counts in tests
  DNN_EIG_SUBSPACE = 2  // track the negative eigenspace by a filtered subspace
                        // iteration, full ssyevd every `refresh` iterations
};

struct DNNParams {
  int max_iter;        // ADMM iteration limit (default 4000)
  double time_limit;   // CPU seconds, <= 0 for none (default 0)
  int check_every;     // certify (and test convergence) every this many iterations (10)
  double tol;          // convergence: primal residual and step below tol (1e-5)
  double gamma;        // dual step length (1.618)
  double beta;         // ADMM penalty, as a multiple of dnn.py's m/3 for the problem
                       // scaled to max |L| = 1. dnn.py uses 1; 0.3 needs 2-3x fewer
                       // iterations to certify on nug20/nug30 nodes (0.3)
  int poly;            // DNNPoly (DNN_POLY_GL)
  int inner_lap;       // DNNInnerLAP (DNN_LAP_AUTO)
  double auto_gap;     // DNN_LAP_AUTO: solve exactly when lb > cutoff - auto_gap;
                       // <= 0 means 1e-3 |cutoff| (0)
  int prototype_cert;  // 1: certify like dnn.py (project V'ZV with an extra
                       //    eigendecomposition); 0: use the R-step multiplier (0)
  int eig_mode;        // R-step eigensolver, DNNEig (DNN_EIG_SUBSPACE). Certificates
                       // are computed in double in every mode, and re-projected in
                       // double when the difference matters.
  int refresh;         // DNN_EIG_SUBSPACE: full eigendecomposition every this many
                       // iterations (20)
  int cheb_degree;     // DNN_EIG_SUBSPACE: Chebyshev filter degree (2)
  int stall_window;    // stop as STALLED when over the last stall_window checks the
                       // certified lb moved less than stall_lb and is within
                       // stall_gap of the primal value, as dnn.py; 0 = off (0)
  double stall_lb, stall_gap;  // (0.05, 0.5), in objective units as in dnn.py
  int hopeless_after;  // after this many iterations, stop as HOPELESS when the
                       // certified lb, extrapolated generously from its recent
                       // gains, stays below the cutoff. A heuristic (it may give
                       // up on a fathomable node); lb stays valid. 0 = off (200)
  const double *L_add; // optional (m^2+1)^2 symmetric matrix added to the lifted
                       // objective L (unscaled, row-major), e.g. multiplier-weighted
                       // valid inequalities; the caller adjusts the shift for their
                       // constants. NULL: none (default)
  double L_add_abs;    // sum of the magnitudes of the terms in L_add, for the
                       // rounding margin (0)
  int child_bounds;    // 1: also return child bounds (DNNResult::child_lb) from the
                       // best certificate's outer assignment problem (0)
  int verbose;         // 1: print the certified bound at every check, 2: also a
                       // CPU profile at the end

  DNNParams();
};

struct DNNResult {
  double lb;           // certified lower bound, shift included (-HUGE_VAL if none)
  double primal;       // ADMM primal value <L,Y> + shift (NOT a bound)
  int status;          // DNNStatus
  int iterations;
  double cpu;          // CPU seconds
  int certificates;    // certified evaluations made
  int cert_failures;   // Cholesky checks that failed (bound not updated then)
  int refinements;     // double-precision re-projections (eig_mode != DOUBLE only)
  double margin;       // largest safety margin subtracted from a certified value
  // With par.child_bounds, unless fathomed: child_lb[i * m + j] is a certified
  // lower bound (shift included) on the node with its facility i at location j
  // (local indices, as reduceQAP numbers them), from the certificate that gave
  // lb: lb + r_ij for the reduced costs r of the outer assignment problem of
  // the Gilmore-Lawler evaluation, less a rounding allowance. No new SDP solve.
  std::vector<double> child_lb;
};

// Final ADMM multiplier of a bounded node, to warm-start its children
// (optional; 2-3x fewer iterations for a child in the tests).
struct DNNWarm {
  int m;
  double scale;                 // the node's objective scale
  std::vector<double> Z;        // lifted (m^2+1)^2, row-major, scaled units
  std::vector<double> Y;        // the final primal iterate (not a solution; for
                                // separating cuts), (m^2+1)^2 row-major
  DNNWarm() : m(0), scale(0) {}
};

// start: a warm start for this node (from dnn_warm_child), or NULL.
// save: if not NULL, receives this node's final ADMM state.
DNNResult dnn_bound(const MAT *A, const MAT *B, const MAT *C, double shift,
                    double cutoff, const DNNParams &par, const DNNWarm *start = NULL,
                    DNNWarm *save = NULL);

// Warm start for the child of a node that fixes its facility i to location
// j (indices local to the node's reduced problem, as reduceQAP numbers them).
// The warm start only changes where the ADMM starts; certificates are
// computed as usual. Returns false on bad input.
bool dnn_warm_child(const DNNWarm &parent, int i, int j, DNNWarm &child);

// Set the number of BLAS threads if the linked BLAS supports it (OpenBLAS);
// returns false otherwise. The bound is meant to run on one core per node.
bool dnn_set_blas_threads(int n);

// Name of the LAPACK/BLAS the module was built against.
const char *dnn_lapack_name();

#endif
