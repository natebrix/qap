/*

  double precision LAP solver.

  a minor modification of J-V code in lap.h/lap.C

 */

#ifndef DLAP_H
#define DLAP_H

#include <math.h>

// "no minimum yet"; a finite value (it was 1e8) breaks on costs above it
const double DLAP_BIG=HUGE_VAL;

double dlap(int dim, double **assigndouble,
	    int *rowsol, int *colsol, double *u, double *v,
	    double augment_tol = 10.0, bool warm_start = false);
//double augment_tol = 1e-5);


#endif

