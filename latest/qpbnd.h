/*
============================================================

Lower bounds for QAP:

qpbnd():            plain ol' QP bound.

qpbnd_parametric(): QP bound where S and T are recalculated
                    every 'step' iterations.  Appears to 
                    perform better than qpbnd in practice.

glbnd(A,B,C,shift): gilmore-lawler bound.  Not best 
                    implementation, but not bad.

pbnd(A,B,C,shift):  projected eigenvalue bound (PB).  
                    QPB >= PB.

============================================================
*/

#ifndef QPBND_H
#define QPBND_H

#include "qap.h"
#include "fw.h"

class QPBParameters 
{
 public:
  bool compute_lapU;       // whether or not to solve LAP(U)
  int maxier;              // maximum # Frank-Wolfe iterations
  int maxier_nofathom;     // maximum # Frank-Wolfe iterations, when
                           //  we know we can't fathom current node
  double shift;            // constant to add to bound
  double inc;              // current upper bound (incumbent)
  double lapU;             // storage for solution of LAP(U)
  int step;                // update S and T every step iterations
  bool feasible;           // is given X doubly stochastic?
  double lap_scale;        //
  bool improve;            // try bound improvement procedure
  bool save_best;          // keep the best lower bound over all
                           // FW iers, or just the last lb?
  int bound_used;          // which bound was used?

  int fw_iter;             // total FW iterations performed

  bool debug;

  QPBParameters() {
    debug = false;
  }
};

const double INIT_LOBND = -1.0e99;

void check_for_improvement(MAT *A, MAT *B, MAT *C, double shift,
			   MAT *U_qpb,double z_qpb);

void init_gradient(MAT *A,MAT *B,MAT *C,MAT *G,double &C_X,double &G_X,int n);

double init_ST(MAT *A,MAT *B,MAT *C,
	       VEC *&w,VEC *&u,MAT *&W,MAT *&U,VEC *&s,VEC *&t,
	       MAT *&S,MAT *&T,bool need_WU);

void update_STG(VEC *w,VEC *u,MAT *W, MAT *U,VEC *s, VEC *t, 
		MAT *G,MAT *X,MAT *S, MAT *T,double factor);

VEC	*v_ins_sort_descending(VEC *x, PERM *order);
VEC	*v_ins_sort_ascending(VEC *x, PERM *order);
void vtav(MAT *A, MAT *OUT);
void vavt(MAT *A, MAT *OUT);

void find_basis(VEC *w,VEC *u,VEC *s,VEC *t);

double qpbnd(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1, QPBParameters *p);

double qpbnd_parametric(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1, QPBParameters *p);

double qp_glb_bnd(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
		  QPBParameters *p);

double pbnd(MAT *A,MAT *B,MAT *C,double shift);

void leader_matrix(MAT *A,MAT *B,MAT *C,MAT *L);
double glbnd(MAT *A,MAT *B,MAT *C,MAT *U,double shift);


double qpbnd0(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
	     QPBParameters *p,double lam);

int *qpb_heuristic(MAT *X);

#endif
