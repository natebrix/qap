#include "fw.h"
#include "qpbnd2.h"
#include <math.h>
#include "linalg.h"
#include "util.h"

#include <iostream>

/*
  [please see comments in qpbnd2.h]
 */

extern int accurate_ier;
extern double accurate_scaling;

//extern stopwatch s_init,s_fw,s_update,s_lap;

void init_w0(MAT *A,VEC *&w,MAT *&W,bool ascending)
{
  const int n = A->m;
  MAT *A1=MNULL;

  int *order_w = new int[n-1];
  for (int k = 0; k < n-1; k++) order_w[k] = k;

  /* compute a1 = v'av */
  A1 = m_get(n-1,n-1);
  vtav(A,A1);

  my_svd(A1,w);

  /* now oldA1 = A1'*diag(u)*A1 */
  if(ascending)
    ins_sort_ascending(w->ve,  order_w, n-1);
  else
    ins_sort_descending(w->ve, order_w, n-1);

  // Callers pass W sized (n x n) but it has to be (n-1 x n-1). Meschach's
  // px_rows() allocated a new matrix then, and the caller's old one leaked:
  // ~5 GB by 20M nodes on nug24. Free it first and allocate the right size.
  if (W != MNULL) M_FREE(W);
  W  = m_get(n-1,n-1);  permrows(A1,order_w,W);   // W[k] = A1[order_w[k]]

  /* free memory */
  M_FREE(A1);
  delete [] order_w;
}

// why the hell does this cause a core dump at the end of
// qpbnd_parametric2()?
void init_w(MAT *A,VEC *&w,MAT *&W,bool ascending)
{
  const int n = A->m;
  MAT *A1 = MNULL;

  int *order_w0 = NULL;

  order_w0 = new int[n-1];
  for(int i=0;i<n-1;i++)
    order_w0[i] = i;

  A1 = m_get(n-1,n-1);

  /* compute a1 = v'av */
  vtav(A,A1);

  my_svd(A1,w);

  if(ascending)
    ins_sort_ascending(w->ve, order_w0, n-1);
  else
    ins_sort_descending(w->ve, order_w0, n-1);

  permrows(A1,order_w0,W);

  /* free memory */
  M_FREE(A1);
  delete [] order_w0;
}

double qpbnd2(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
	      VEC *w,VEC *u,MAT *W,MAT *U,QPBParameters *qpbp)
{
  bool improve = false;
  const int maxier = qpbp->maxier;
  const int maxier_nofathom = qpbp->maxier_nofathom;
  int status,bnd_used;
  int iter, total_iter;
  double temp;
  double eig;
  double lo_bound, up_bound;
  MAT *S,*T;
  MAT *G;
  VEC *s,*t;

  int *xstar;

  double C_X, G_X;

  // s_init.go(); 
  eig = init_ST(A,B,C,w,u,W,U,s,t,S,T,false);
  // s_init.stop(); 
  // s_fw.go();
  total_iter = 0;
  up_bound =  1.0e99;
  lo_bound = -1.0e99;

  QPBWork &qw = qpb_work(X->n);
  xstar = qw.xstar;
  G = qw.G;   // the gradient, zeroed

  if(is_uniform(X)) {
    init_gradient(A,B,C,G,C_X,G_X,X->n);
  }
  else {
    gradient(A,B,C,S,T,G,X,X->n);
    C_X = mm_trace(C,X);
    G_X = mm_trace(G,X);
  }

  fw(A,B,C,S,T,G,X,U1,xstar,C_X,G_X,A->m,iter,total_iter,
     maxier-accurate_ier,maxier_nofathom-total_iter,
     qpbp->shift+eig,qpbp->inc - qpbp->cutoff_slack,qpbp->feasible,
     status, up_bound, lo_bound, qpbp->lap_scale, false, false, bnd_used);
  total_iter += iter;

  if(lo_bound <= qpbp->inc - qpbp->cutoff_slack) {
    if(accurate_ier > 0) {
      lo_bound = INIT_LOBND;     // (was -1e-99, i.e. 0: too high when bounds are negative)
      fw(A,B,C,S,T,G,X,U1,xstar,C_X,G_X,A->m,iter,total_iter,accurate_ier,1000,
	 qpbp->shift+eig,qpbp->inc - qpbp->cutoff_slack,qpbp->feasible,
	 status, up_bound, lo_bound, accurate_scaling,true,improve,bnd_used);
      total_iter += iter;
    }
    else 
      lo_bound = compute_U(A,B,C,S,T,X,U1,A->m,false,temp,
			   accurate_scaling,qpbp->shift+eig);
  }
  qpbp->fw_iter += total_iter;
  //  s_fw.stop();
  /*
  printf("init %.10f   fw %.10f    up %.10f   lap %.10f\n",
	 s_init.elapsed(),s_fw.elapsed(),
	 s_update.elapsed(),s_lap.elapsed());
  */
  M_FREE(S);  M_FREE(T);
  V_FREE(s);  V_FREE(t);

  return lo_bound;
}

double qpbnd_parametric2(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
			 VEC *w,VEC *u,MAT *W,MAT *U,
			 QPBParameters *qpbp)
{
  int status, bnd_used;
  int iter, total_iter;
  double temp;
  double eig;
  double lo_bound, up_bound;
  MAT *S=MNULL,*T=MNULL;
  MAT *G;
  VEC *s=VNULL,*t=VNULL;

  int *xstar;

  double C_X, G_X;

  eig = init_ST(A,B,C,w,u,W,U,s,t,S,T,false);
  total_iter = 0;

  up_bound =  1.0e99;
  lo_bound =  INIT_LOBND;

  QPBWork &qw = qpb_work(X->n);
  xstar = qw.xstar;
  G = qw.G;   // the gradient, zeroed

  if(is_uniform(X)) {
    init_gradient(A,B,C,G,C_X,G_X,X->n);
  }
  else {
    gradient(A,B,C,S,T,G,X,X->n);
    C_X = mm_trace(C,X);
    G_X = mm_trace(G,X);
  }

  do {
    fw(A,B,C,S,T,G,X,U1,xstar,C_X,G_X,A->m,iter,total_iter,
       std::min(qpbp->step,qpbp->maxier - total_iter-accurate_ier),
       qpbp->maxier_nofathom-total_iter,
       qpbp->shift+eig,qpbp->inc - qpbp->cutoff_slack,qpbp->feasible,
       status, up_bound, lo_bound, qpbp->lap_scale, false, false, bnd_used);
    total_iter += iter;

    if(status != CONTINUE) break;

    if(total_iter < qpbp->maxier - accurate_ier) {
      // update S, T matrices.  The gradient can be updated
      // at the same time.
      update_STG(w,u,W,U,s,t,G,X,S,T,1.0);
      G_X = mm_trace(G,X);
    }
  }  while(total_iter < qpbp->maxier - accurate_ier);

  if(status != FATHOM) {
    if(accurate_ier > 0) {
      lo_bound = INIT_LOBND;
      fw(A,B,C,S,T,G,X,U1,xstar,C_X,G_X,A->m,iter,total_iter,accurate_ier,1000,
	 qpbp->shift+eig, qpbp->inc - qpbp->cutoff_slack,qpbp->feasible,
	 status, up_bound, lo_bound, accurate_scaling,true,qpbp->improve,
	 bnd_used);
      total_iter += iter;
    }
    else 
      lo_bound = compute_U(A,B,C,S,T,X,U1,A->m,false,temp,
			   accurate_scaling,qpbp->shift+eig);
  }
  qpbp->fw_iter += total_iter;

  V_FREE(t);
  V_FREE(s);
  M_FREE(S);
  M_FREE(T);


  return lo_bound;
}
