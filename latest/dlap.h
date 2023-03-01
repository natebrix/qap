/*

  double precision LAP solver.

  a minor modification of J-V code in lap.h/lap.C

 */

#ifndef DLAP_H
#define DLAP_H

const double DLAP_BIG=100000000.0;

double dlap(int dim, double **assigndouble,
	    int *rowsol, int *colsol, double *u, double *v,
	    double augment_tol = 10.0);
//double augment_tol = 1e-5);


#endif

