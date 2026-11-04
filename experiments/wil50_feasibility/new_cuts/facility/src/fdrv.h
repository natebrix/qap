// fdrv.h: export hook of dnncut_f.C (a renamed copy of latest/dnncut.C) for the
// facility-side experiments (new_cuts/facility). Nothing here changes a bound.
#ifndef FDRV_H
#define FDRV_H
#include "dnncut.h"
#include <vector>
struct CutExport {
  int stop_before_quiet = 0;      // return right before the quiet round
  int have = 0;                   // 1: the fields below were filled at the quiet round
  std::vector<double> Y0;         // the iterate separated in round 0 (the give-up point)
  std::vector<double> Ltri;       // sum mu_t G_t over the pool at the quiet round's multipliers
  double hsum_tri = 0, abs_tri = 0;   // sum mu_t h_t, the rounding magnitude
  int active = 0;
  double best_mu_lb = 0;
  DNNWarm state;                  // the ADMM state the quiet round starts from
  DNNWarm final_state;            // after the quiet round (if run)
  double quiet_lb = 0, quiet_primal = 0;
};
#include <string>
// extra regions for dnncut_f.C: original indices on a side (0 = first index,
// 1 = second), mapped to the node's free elements; nogrid drops the grid regions
struct FdrvRegion { int side; std::string name; std::vector<int> orig; };
struct FdrvRegions { int nogrid = 0; std::vector<FdrvRegion> extra; std::vector<int> free_orig[2]; };
extern FdrvRegions g_fdrv_regions;
DNNCutResult dnn_cut_rounds_f(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                              const DNNParams &base, const std::vector<int> &free_orig,
                              int n_orig, const DNNCutParams &cp, DNNWarm *warm,
                              double first_lb, CutExport *ex);
#endif
