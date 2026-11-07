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

 Five-point count inequalities (cp.kcuts5 > 0, regional only): for five free
 elements T, a flip set F of them (z'_i = 1 - z_i on F) and alpha = 1, 2, the
 count q' = sum_{i in T} z'_i is an integer, so (q' - alpha)(q' - alpha - 1)
 >= 0, which in the moments (p, W) is the clique/cut inequality
   alpha sum_i p'_i - sum_{i<j} W'_ij <= alpha (alpha + 1) / 2
 (Padberg 1989, facets of the Boolean quadric polytope; the triangles are
 the three-element case, and the pure pentagonal inequality is F empty,
 alpha 2). Not implied by the triangles nor by any degree-2 moment matrix
 (research_bound_quality.md section 2.2); separated exhaustively, 48
 patterns per 5-subset per region.

 Local-7 inequalities (cp.kcuts7 > 0, regional only): for seven free elements
 T of a region S, the moments (p_T, W_TT) are the first and second moments of
 the 0/1 vector z_T, so they lie in the convex hull of the admissible
 configurations c in {0,1}^7 (sum c <= |S|, 7 - sum c <= m - |S|). A point
 outside it is separated by projecting onto that hull (Frank-Wolfe over the
 at most 128 vertices): the direction y - x, scaled to max |coefficient| = 1,
 with the right-hand side beta = max over the admissible c, is a valid
 inequality a.p_T + <B, W_TT> <= beta for every permutation. The triangles are
 the 3-element case. Separation samples cp.sub7 random 7-subsets per region
 plus subsets of the most fractional elements (stronger_bounds/local7: the
 violations are 1.7-2.6x the triangles' at wil50's converged iterates).

 Box (count) inequalities (cp.kbox > 0, regional only; Junger-Kaibel,
 Math. Prog. 91, 2001, on the regions' membership vectors): for a set Q of
 free elements of any size and an integer b >= 1, the count q = sum_{i in Q}
 z_i is an integer, so (q - b)(q - b - 1) >= 0, which in the moments reads
   b sum_{i in Q} p_i - sum_{i<j in Q} W_ij <= b (b + 1) / 2
 (the sum triangle is |Q| = 3, b = 1; the 5-point family with no flips is
 |Q| = 5). Separated per region by a local search over Q (add/drop moves
 from cp.box_starts starts per b = 1..cp.box_bmax, sizes 3..cp.box_qmax;
 new_cuts/literature/wil50/box_gate.py): on wil50 the violated ones sit on
 the edge regions with |Q| about 16, where the triangles are weak. The pool's
 box terms are summed per region before they are written into the objective.

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
#include <string>

// A node's cut pool in original indices, handed to its children (pool
// inheritance, bnb.C QAP_DNN_POOL; experiments/pool_inheritance): the cuts with
// a positive multiplier in the solve whose state the children start from, with
// those multipliers. regions: sorted original indices of the region side;
// mem: original indices of the other side (a < b < c for a triangle of type
// 0-3, the sorted set Q for a box, type 8, with its b). Every cut of these
// families is valid for any region, members and b >= 1, so a child keeps a
// cut on its own free part of the region: a triangle loses nothing unless one
// of its members was fixed (then dropped); a box member fixed inside the region
// counts one (b - 1), outside zero.
struct DNNPoolCut { int type, region, b; double mu; std::vector<int> mem; };
struct DNNPool {
  int side;
  std::vector<std::vector<int>> regions;
  std::vector<DNNPoolCut> cuts;
  DNNPool() : side(0) {}
};
// one whitespace-free token for state files ("-" for an empty pool)
std::string dnn_pool_encode(const DNNPool &p);
bool dnn_pool_decode(const std::string &s, DNNPool &p);

struct DNNCutParams {
  int warm_primal;    // opt-in same-node Y/Z restart between rounds (0); see DNNParams
  int rounds;         // cut rounds after the uncut solve (12)
  int iters;          // iterations per cut round (300)
  int first_iters;    // > 0: run the uncut round here with this many iterations;
                      // 0: the caller's warm start (Y, Z) is the uncut solve
                      // and first_lb is its certified bound (0)
  int kcuts;          // new cuts per round (1000)
  int kcuts5;         // > 0: also separate the 5-point count inequalities and add
                      // this many per round, on top of kcuts triangles (0 = off)
  int kcuts7;         // > 0: also separate the local-7 inequalities and add this many
                      // per round (0 = off)
  int sub7;           // random 7-subsets sampled per region and round (1000); a third as
                      // many again are drawn from the 15 most fractional elements
  int kbox;           // > 0: also separate the box (count) inequalities and add this many
                      // per round (0 = off)
  int box_qmax;       // largest |Q| of a box cut (24)
  int box_starts;     // local-search starts per region and b (40)
  int box_bmax;       // b = 1..box_bmax (7)
  int add_until;      // > 0: add new cuts only in rounds <= add_until; later rounds
                      // move the multipliers of the frozen pool (the dual-method test,
                      // stronger_bounds/dualtest). 0 = every round (default)
  double delta;       // Polyak target = min(cutoff, best + delta) (50)
  double theta;       // Polyak step factor (1.0); 0: multipliers stay at mu0
  double damp1;       // factor on the first Polyak step of a phase that starts from an
                      // inherited pool (1.0: none). That step starts from near-optimal
                      // multipliers and overshoots (-88 to -116 at 49:45,
                      // roadmap/convergence2 E2); QAP_DNN_CUT_DAMP1 in bnb.C
  double damp_rest;   // the same factor on the later steps of such a phase (1.0: none;
                      // QAP_DNN_CUT_DAMP1=F:G sets damp1 F and damp_rest G)
  double mu0;         // multiplier of a newly added cut (0)
  int grid_rows, grid_cols;  // regions: the prefixes of this grid over the
                      // original indices of `side`; 0: none
  int side;           // 0: regions over the first index (A), 1: the second (B)
  bool prefix;        // regions: prefixes of the free elements (validation)
  bool atomic;        // atomic triangle cuts instead of regional ones
  int cap_at_cutoff;  // 1: Polyak target = min(cutoff, best + delta) (default); 0: best +
                      // delta, so the step does not vanish within 0.01 of the cutoff
                      // (nodes crawled there for 4-6 rounds; research_bound_quality.md)
  int freeze_after;   // > 0: run only this many cut rounds before the quiet round
  int freeze_iters;   // > 0: a quiet round at the end when not fathomed (multipliers
                      // frozen, no new cuts, this many iterations): the certificate
                      // trails the Lagrangian primal by 25-90 at 300-iteration rounds
  int freeze_if_primal_above;  // 1: the quiet round only when the last primal value is
                      // above the cutoff (the search's use); 0: always (0)
  int freeze_best;    // 1 (default): the quiet round runs from the multipliers (and
                      // iterate) of the best round, not the last one: a last-round
                      // overshoot otherwise freezes a bad multiplier (laptop stepA_2,
                      // stepA_14, short4_d1 in stronger_bounds/laptop/dual). 0: the last
  int stall_rounds;   // > 0: after this many cut rounds, stop when the gap to the
                      // cutoff exceeds stall_factor times the gain of the last
                      // stall_rounds rounds (the gain comes early: on wil50 depth-4
                      // nodes +53 in three rounds, then +3 a round). 0 = off
  double stall_factor;  // (3.0)
  int pgate_rounds;   // > 0: the primal gate (wil50_feasibility/trims): after this many
                      // cut rounds, stop when the best Lagrangian primal value of the
                      // rounds so far is more than pgate_margin below the cutoff: the
                      // certified bound never passes the primal, so the phase cannot
                      // fathom soon. Replayed on 88 logged phases: 4 rounds, margin 20
                      // stops 44 of 63 misses with no false stop (closest fathom 12.4
                      // short, closest caught miss 26.6). 0 = off
  double pgate_margin;  // (20)
  int verbose;        // >= 1: a line per round, >= 2: per region
  const double *L_ext;  // optional (m^2+1)^2 row-major term added to every round's
                      // objective on top of the cuts (unscaled; e.g. externally
                      // separated cuts at fixed multipliers, stronger_bounds/triples);
                      // NULL: none
  double L_ext_abs;   // sum of |entries| of L_ext (rounding margin)
  double ext_const;   // constant subtracted from the shift for L_ext's right-hand sides
  // pool inheritance (DNNPool above). pool_in: a parent's pool, remapped to
  // this node; with first_iters 0 the caller's warm start (first_lb,
  // first_primal) was solved with it at its multipliers (dnn_pool_bound), else
  // round 0 solves with it. The phase then runs as usual from that start, the
  // best round so far being round 0, and its quiet-round test uses the best
  // round's primal (the start of an inherited phase is usually its best round
  // and the rounds after it overshoot). pool_out: out, the pool for the
  // children, at the multipliers of the solve left in *warm (the quiet round's,
  // else the best round's: *warm is set to that round's state); empty when a
  // phase that did not inherit never beat the uncut bound (first_lb).
  const DNNPool *pool_in;
  DNNPool *pool_out;
  const std::vector<int> *other_orig;  // with a pool: the original indices of the other
                      // side's free elements, in local order
  const std::vector<int> *fixed_at;    // with pool_in: for each original index of the other
                      // side, the original region-side index it is fixed to (-1: free)
  int pool_check;     // > 0: evaluate the pool's term on this many random permutations of
                      // the node and abort if it is positive (a guard on the assembly)
  int pool_all_round; // >= 0 with pool_all_out: the whole pool after this round's new cuts,
  DNNPool *pool_all_out;  // zero multipliers included (the CP-ADMM study's frozen pool;
                      // needs other_orig). -1 / NULL: off (default)
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
  int gated;          // 1 when the primal gate stopped the rounds
};

// A, B, C, shift: the node's reduced QAP (reduceQAP). free_orig: for each local
// index t of `side`, its original index (reduceQAP keeps the original order);
// n_orig: the instance size. warm: in, the uncut solve when cp.first_iters is
// 0 (its Y is separated first); out, the last iterate. first_lb: the uncut
// certified bound when cp.first_iters is 0.
DNNCutResult dnn_cut_rounds(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                            const DNNParams &base, const std::vector<int> &free_orig,
                            int n_orig, const DNNCutParams &cp, DNNWarm *warm,
                            double first_lb = -HUGE_VAL, double first_primal = HUGE_VAL);

// The node's bound with cp.pool_in remapped to it at the inherited multipliers
// (frozen): one dnn_bound(par, start, save) on the modified problem; every
// value is certified. cp.side, cp.other_orig and cp.fixed_at describe the node;
// cp.pool_out (optional) gets the remapped pool. kept: cuts kept.
DNNResult dnn_pool_bound(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                         const DNNParams &par, const std::vector<int> &free_orig, int n_orig,
                         const DNNCutParams &cp, const DNNWarm *start, DNNWarm *save,
                         int *kept = NULL, int *dropped = NULL);

// CP-ADMM (experiments/wil50_feasibility/roadmap/admm/cpadmm): the node's bound
// with cp.pool_in (remapped, zero multipliers kept) projected inside the Y-step by
// cyclic Dykstra sweeps (DNNYCuts in dnn.h), certified at every check with the
// sweeps' multipliers. One dnn_bound(par, start, save); every value is certified.
// Regional cuts on side 0 only.
struct DNNCPParams {
  int sweeps;         // Dykstra sweeps (P-step, half-space pass) per Y-step (1)
  int passes;         // half-space passes per sweep (1)
  int diag_every;     // > 0: every this many Y-steps, the distance of the output to a
  int diag_sweeps;    // run of diag_sweeps sweeps from the same start (0; 30)
  int selftest;       // 1: check the cut extraction against assemble() first (0)
  DNNCPParams() : sweeps(1), passes(1), diag_every(0), diag_sweeps(30), selftest(0) {}
};
DNNResult dnn_cp_bound(MAT *A, MAT *B, MAT *C, double shift, double cutoff, const DNNParams &par,
                       const std::vector<int> &free_orig, int n_orig, const DNNCutParams &cp,
                       const DNNCPParams &cpp, const DNNWarm *start, DNNWarm *save, int *kept = NULL);

#endif
