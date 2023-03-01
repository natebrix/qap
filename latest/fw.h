/*
============================================================
Contains code to solve the quadratic program

    min f(X) = tr (AXB - SX - TX + C)X'
    s.t.
           sum_i( X_ij ) = 1,    j = 1..n
(*)        sum_j( X_ij ) = 1,    i = 1..n
           X_ij >= 0,          i,j = 1..n

  We write X \in P to denote that X satisfies the constraints
  of (*).

  We use the Frank-Wolfe method, 
  (aka 'conditional gradient method)

     1) compute X_0 \in P, set i = 0.
     2) solve
         min grad(X_i)' D
           D \in P
     3) compute steplength k minimizing
         f( (1-k) X_i + k D)
     4) update: 
         X_(i+1) = (1-k) X_i + k D
     5) i = i + 1, go back to 2.

  The problem in 2) is a linear assignment problem, which
  is solved easily.

  -----

  fw() takes *lots* of parameters:

  fw(A,B,C,S,T,G,X,U,xstar,C_X,G_X,n,ier,startit,maxit,maxit_nofathom,
  offset, cutoff, feasible, status, up_bnd, bnd, mult, save_best, 
  improve, bnd_used, need_U);

  (A,B,C,S,T) :  define QP to solve
  X:             solution of the QP
  U:             dual matrix of QP: provides bound estimates for 
                 child subproblems.
  C_X:           value of C.X, modified by fw().
  G_X:           value of G.X, modified by fw().
  n:             size of A,B,C,S,T
  ier:           (output) number of FW iterations performed
  maxit:         perform at most maxit iterations
  maxit_nofathom if up_bound is < cutoff, there is no way
                 we will be able to fathom at this node, so if
                 ier > maxit_nofathom, return.  we may not want
                 to return right away, since we want accurate
                 U for branching purposes.
  offset:        constant term added on to bound
  cutoff:        terminate if lower bound is greater than this value.
  feasible:      does the provided X have row/col sums equal to 1?
  status:        reason why algorithm terminated -- if status==0
                 then there is no point in trying to improve the bound.
  bnd/up_bnd:    the solution value of the QP is in interval 
                 [bound,up_bound].

  mult:          scale LAPs by this factor before solving (but only if
                 INTEGER_LAP is defined)

  save_best:     if true, return the bound,X,U corresponding the largest
                 lower bound, otherwise just return the last bnd.
  improve:       if true, solve a few additional LAPs at the end to
                 try to improve the bound.
  bnd_used:      records whether the bound improvement procedure was
                 actually used or not; for bookkeeping
  need_U:        if true, U will be computed, otherwise it will not.

============================================================
*/

#ifndef FW_H
#define FW_H

// exit status for fw procedure:
const int FATHOM    = -1;
const int NO_FATHOM =  0;
const int CONTINUE  =  1;
const int MAX_ITER  =  2;

#include "meschach.h"
#undef catch
#include "lap.h"
#include "dlap.h"

void gradient(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,
	      MAT *G,MAT *X,int n);

void fw(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,MAT *G,MAT *X,MAT *U,
	int *xstar,double &C_X, double &G_X,
	int n,int &ier,int startit,int maxit,int maxit_nofathom,double offset,
	double cutoff,
	bool &feasible,int &status,double &up_bnd,double &bnd,double mult,
	bool save_best,bool improve,int &bnd_used,bool need_U=false);

double compute_U(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,MAT *X,MAT *U,int n,
	       bool compute_lapU, double &lapU,double mult,double offset);

void init_X(MAT *X,int n);

int try_to_improve(MAT *U_qpb,int *x_U1,double &z_qpb,
		   MAT *U_glb,int *x_U2,double &z_glb,
		   double incumbent,const int steps);

#endif
