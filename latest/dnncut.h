/*
 dnncut.h: Boolean triangle cuts on top of the certified DNN bound, in
 Lagrangian form (experiments/dnn_cuts, experiments/wil50_feasibility/
 stronger_bounds). A library for the search (bnb.C) and the dnncuts driver.

 Two families, both valid for every permutation and folded into the lifted
 objective with multipliers mu >= 0, so that every certified value of the
 modified problem (shift lowered by sum mu_t h_t) is a valid bound:

 Regional (the default): for a set S of original indices on one side (the
 locations, side 0, or the facilities, side 1) and three free elements
 i, j, k of the other side, the membership variables z_i = sum_{a in S} x_ai
 are 0/1, their first and second moments p_i = sum_{a in S} Y_(ai),(ai) and
 W_ij = sum_{a,b in S} Y_(ai),(bj) are entries of the DNN matrix, and the
 four Boolean triangle inequalities hold for (p, W):
   W_ij + W_ik - W_jk - p_i <= 0  (and the apexes j, k),
   p_i + p_j + p_k - W_ij - W_ik - W_jk <= 1.
 They are not implied by the DNN relaxation (an exact witness is in the
 memo). The regions are the row and column prefixes of a rows x cols grid
 (wil50: 5:10 on the locations, the first matrix; nug30: 5:6 on the second),
 which decompose a Manhattan distance exactly, mapped to the node's free
 elements; or the prefixes of the free elements (validation).

 Atomic (cp.atomic): the same inequalities on single lifted variables.

 Each round solves the DNN relaxation with the current multipliers
 (warm-started), separates the kcuts most violated cuts in the last primal
 iterate, and updates every multiplier by a Polyak subgradient step toward
 target = min(cutoff, best + delta), halving the step factor after a round
 that lowered the bound. The run stops at the first certified bound above
 the cutoff.
*/
#ifndef DNNCUT_H
#define DNNCUT_H

#include "dnn.h"
#include <vector>

struct DNNCutParams {
  int rounds;         // cut rounds after the uncut solve (12)
  int iters;          // iterations per cut round (300)
  int first_iters;    // > 0: run the uncut round here with this many iterations;
                      // 0: the caller's warm start (Y, Z) is the uncut solve
                      // and first_lb is its certified bound (0)
  int kcuts;          // new cuts per round (1000)
  double delta;       // Polyak target = min(cutoff, best + delta) (50)
  double theta;       // Polyak step factor (1.0); 0: multipliers stay at mu0
  double mu0;         // multiplier of a newly added cut (0)
  int grid_rows, grid_cols;  // regions: the prefixes of this grid over the
                      // original indices of `side`; 0: none
  int side;           // 0: regions over the first index (A), 1: the second (B)
  bool prefix;        // regions: prefixes of the free elements (validation)
  bool atomic;        // atomic triangle cuts instead of regional ones
  int stall_rounds;   // > 0: after this many cut rounds, stop when the gap to the
                      // cutoff exceeds stall_factor times the gain of the last
                      // stall_rounds rounds (the gain comes early: on wil50 depth-4
                      // nodes +53 in three rounds, then +3 a round). 0 = off
  double stall_factor;  // (3.0)
  int verbose;        // >= 1: a line per round, >= 2: per region
  DNNCutParams();
};

struct DNNCutResult {
  double first_lb;    // the uncut certified bound
  double best_lb;     // best certified bound over the rounds: valid
  int status;         // DNN_FATHOMED when best_lb > cutoff, else the last round's status
  int rounds;         // cut rounds run
  int iterations;     // ADMM iterations of the cut rounds
  double cpu;         // CPU seconds of the cut rounds (and the uncut one when run here)
  int active;         // cuts with a positive multiplier at the end
  int nregions;
  int stalled;        // 1 when the stall rule stopped the rounds
};

// A, B, C, shift: the node's reduced QAP (reduceQAP). free_orig: for each local
// index t of `side`, its original index (reduceQAP keeps the original order);
// n_orig: the instance size. warm: in, the uncut solve when cp.first_iters is
// 0 (its Y is separated first); out, the last iterate. first_lb: the uncut
// certified bound when cp.first_iters is 0.
DNNCutResult dnn_cut_rounds(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                            const DNNParams &base, const std::vector<int> &free_orig,
                            int n_orig, const DNNCutParams &cp, DNNWarm *warm,
                            double first_lb = -HUGE_VAL);

#endif
