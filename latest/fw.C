/*
  bugs fixed: 
  12/5/99: once in a great while, bound was returned as 1e-99.
  this was due to a zero steplength after several calls to fw().
  basically, the QP was solved exactly.  reordered the test for
  zero steplength, so that the bound is updated before returning.

  12/17/99: basically scrapped the old fw code and rewrote it
  to elimiate unnecessary O(n^3) work.

  4/5/00: can either use integer or fp LAP solver...#define
          INTEGER_LAP in util.h to use integer LAP solves.

 */

#include <stdio.h>
#include <math.h>

#include "fw.h"
#include "util.h"
#include "linalg.h"

//#define PRINT

//extern stopwatch s_lap, s_dir;

/*****************************************************************/
/*
============================================================
 init_X(X,m);

 Initialize the matrix X, of size m x m, to E/m, where
 E is the matrix of ones.
============================================================
 */
void init_X(MAT *X,int m)
{
  register int i,j;
  double one_over_m = 1.0/m;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      X->me[i][j] = one_over_m;
}

/*
============================================================
obj = qp_obj(A,B,C,S,T,X,n);

compute objective value 
  f(X) = tr AXBX' - SX - TX + CX' 
  f(X) = tr AXBX' - SXX' - XTX' + CX' (or is it this?!?)

  [...this function is quite expensive to call, and there's
      usually no need to call it...]
============================================================
*/
double qp_obj(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,MAT *X,int n)
{
  double obj;
  int i,j;
  MAT *temp1=MNULL,*temp2 = MNULL;

  /*  A*X*B */
  temp1 = my_m_mlt(A,X,MNULL);
  temp2 = my_m_mlt(temp1,B,MNULL);

  /*  (A*X*B - S*X - X*T)*X' + C*X'  */
  M_FREE(temp1);
  temp1 = my_m_mlt(S,X,MNULL);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      temp2->me[i][j] -= temp1->me[i][j];

  M_FREE(temp1);
  temp1 = my_m_mlt(X,T,MNULL);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      temp2->me[i][j] = temp2->me[i][j] - temp1->me[i][j] + C->me[i][j];

  obj =  mm_trace(temp2, X);

  M_FREE(temp1);
  M_FREE(temp2);
  return obj;
}

/*
============================================================
 compute gradient of objective function...that is:
  G =  2(A*X*B - S*X - X*T) + C

 [we avoid calling this function as much as possible, 
  since it is quite expensive: in practice it should only
  be called once, for the initial iterate X_0, and even
  then if X is not equal to E/n.  If X_0=E/n, then call
  init_gradient instead.]
============================================================
 */
void gradient(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,
	      MAT *G,MAT *X,int n)
{
  int i,j;
  MAT *temp = MNULL;
  MAT *temp1 = MNULL;

  temp1 = my_m_mlt(A,X,temp1);
  temp  = my_m_mlt(temp1,B,temp);
  
  temp1 = my_m_mlt(S,X,temp1);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      temp->me[i][j] -= temp1->me[i][j];

  temp1 = my_m_mlt(X,T,temp1);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++) {
      temp->me[i][j] -= temp1->me[i][j];
      G->me[i][j]     = 2.0*temp->me[i][j] + C->me[i][j];
    }
  M_FREE(temp);
  M_FREE(temp1);
}

/*
============================================================
   compute D1 = D - X
   [currently unused]
============================================================
 */
void search_dir(int *d, MAT *X, MAT *D1, int n)
{
  register int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      D1->me[i][j] = -X->me[i][j]; 

  for(i=0;i<n;i++)
      D1->me[i][d[i]] += 1.0;
}

/*
============================================================
  dG = 2(AX*_kB-SX*_k-X*_kT) + C - G_k

  temp1, temp2 are global variables.

  (below are several other functions that do the same thing,
   only faster.  Eventually only one of these deltaG functions
   will remain.  The efficiency of the code depends quite
   heavily on this function (and the LAP solver) )
============================================================
void deltaG(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,
	    MAT *G,MAT *dG,MAT *X,int *xstar,int n)
{
  int i,j;

  // temp2 = A*diag(xstar)*B
  permcols(A,xstar,temp1);
  temp2 = my_m_mlt(temp1,B,temp2);

  permcols(S,xstar,temp1);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      temp2->me[i][j] -= temp1->me[i][j];
  
  permrows(T,xstar,temp1);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      dG->me[i][j] = 2.0*(temp2->me[i][j] - temp1->me[i][j]) +
	C->me[i][j] - G->me[i][j];

}
*/

// no temporary storage required.
inline void deltaG(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,
		   MAT *G,MAT *dG,MAT *X,int *p,int *ip,
		   double &dG_X,int n)
{
  register int i,j,k;
  register double v;
  int *ipk;
  double *Ai, *Si, *Tpi;
  //  double *Ci, *Gi;

  dG_X = 0.0;
  for(i=0;i<n;i++) {
    Ai = A->me[i];    Si = S->me[i];
    Tpi = T->me[p[i]];
    for(j=0;j<n;j++) {
      v = 0.0; ipk = ip;
      for(k=0;k<n;k++) {
	//	v += A->me[i][ip[k]]*B->me[k][j];
	v += Ai[*ipk++]*B->me[k][j];
      }
      dG->me[i][j] = 2.0*(v - Si[ip[j]] - Tpi[j]) 
	+ C->me[i][j] - G->me[i][j];
      dG_X += X->me[i][j]*dG->me[i][j];
    }
  }
}

// faster? 
inline void deltaG_alt(double **a,double **b,double **C,double **S,double **T,
		       double **G,double **dG,double **X,int *p,int *ip,
		       double &dG_X,int n)
{
  register int i,j,k;
  register double aik, *ai, *bkj;
  register double *dgi;

  for(i=0;i<n;i++) {
    dgi = dG[i];  ai = a[i];
    
    bkj = b[0];
    aik = ai[ip[0]];

    for(j=0;j<n;j++) {
      dgi[j] = aik*(*bkj++);
    }

    for(k=1;k<n;k++) {
      bkj = b[k];
      aik = ai[ip[k]];
      
      for(j=0;j<n;j++) {
	dgi[j] += aik*(*bkj++);
      }
    }
  }

  register double *Si, *Tpi;
  dG_X = 0.0;
  for(i=0;i<n;i++) {
    Si = S[i];    Tpi = T[p[i]];
    for(j=0;j<n;j++) {
      dG[i][j] = 2.0*(dG[i][j] - Si[ip[j]] - Tpi[j]) 
	+ C[i][j] - G[i][j];
      dG_X += X[i][j]*dG[i][j];
    }
  }
}

inline void deltaG_alt2(double **a,double **b,double **C,double **S,double **T,
			double **G,double **dG,double **X,int *p,int *ip,
			double &dG_X,double &dG_D,int n)
{
#define DG_UNROLL  (8)
#define DG_UNROLL2 (4)
  register int i,j,k;
  register double aik, *ai, *bkj;
  register double *dgi;

  for(i=0;i<n;i++) {
    dgi = dG[i];  ai = a[i];

    bkj = b[0];
    aik = ai[ip[0]];
    for(j=0;j<n;j++) {
      dgi[j] = aik*(*bkj++);
    }

    for(k=1;k<n;k++) {
      bkj = b[k];
      aik = ai[ip[k]];
      
      // if you change the value of DG_UNROLL, change the body of this loop:
      for(j=0;j<n-(DG_UNROLL-1);j+=DG_UNROLL) {
	//	v += A->me[i][ip[k]]*B->me[k][j];
	dgi[j  ] += aik*(*bkj++);
	dgi[j+1] += aik*(*bkj++);
	dgi[j+2] += aik*(*bkj++);
	dgi[j+3] += aik*(*bkj++);
	dgi[j+4] += aik*(*bkj++);
	dgi[j+5] += aik*(*bkj++);
	dgi[j+6] += aik*(*bkj++);
	dgi[j+7] += aik*(*bkj++);
      }

#if (DG_UNROLL > DG_UNROLL2)
      // if you change the value of DG_UNROLL2, change the body of this loop:
      for(;j<n-(DG_UNROLL2-1);j+=DG_UNROLL2) {
	dgi[j  ] += aik*(*bkj++);
	dgi[j+1] += aik*(*bkj++);
	dgi[j+2] += aik*(*bkj++);
	dgi[j+3] += aik*(*bkj++);
      }
#endif

      // finish up
      for(;j<n;j++) {
	dgi[j]   += aik*(*bkj++);
      }


    }
  }

  register double *Si, *Tpi;
  dG_X = 0.0;
  dG_D = 0.0;
  for(i=0;i<n;i++) {
    Si = S[i];    Tpi = T[p[i]];
    for(j=0;j<n;j++) {
      dG[i][j] = 2.0*(dG[i][j] - Si[ip[j]] - Tpi[j])
      	+ C[i][j] - G[i][j];
      dG_X += X[i][j]*dG[i][j];
    }
    dG_D += dG[i][p[i]];
  }
  dG_D = dG_D - dG_X;
}

/*
============================================================
(the documentation here is incomplete, see fw.h)

bound = fw(A,B,C,S,T,X,n,maxit,offset,cutoff,status);

solve the following quadratic program using the Frank-Wolfe
(conditional gradient) algorithm:

min trace (AXB - SX - XT + C)X' = f(X)
    X doubly stochastic and nonnegative

xstar          : solution of last LAP
maxit          : maximum # of F-W iterations
maxit_nofathom : maximum # of F-W iterations, if the upper
                 bound on f(X) falls below cutoff.
offset: constant term added on to bound
cutoff: terminate if lower bound is greater than this value.
feasible: is the provided X double stochastic and nonnegative?
status: reason why algorithm terminated
upbnd:  upper bound on f(X)
lobnd:  lower bound on f(X)
scaling: scale LAPs by this amount (if necessary...we only scale
         if the LAP solver requires integer input)
save_best: keep the best lower bound, which might not be the last
           lower bound computed.  better bnd quality, more expensive
improve: call bound improvement procedure?

Upper and lower bounds on f(X) are stored in upbnd and lobnd.

At each step, we form G(x), the gradient of f(x), and solve

      min G(X) . D 
         D

Then choose a steplength to 
      minimize trace f((1-k)*X + k*D)
(which is easy since function is quadratic)

============================================================
 */
void fw(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,MAT *G,MAT *X,MAT *U,
	int *xstar,double &C_X, double &G_X,
	int n,int &ier,int startit,int maxit,int maxit_nofathom,double offset,
	double cutoff,bool &feasible,int &status,
	double &upbnd, double &lobnd, double scaling,bool save_best,
	bool improve,int &bnd_used,bool need_U)
{
  const int imp_steps = 3;
  const double augment_tol = 1.0;        // for LAP code, see dlap.h

  int j,k;
  bool done, have_U_matrices;
  MAT *dG;
  double obj, alpha, one_minus_alpha;
  double bnd;
  double G_D, dG_D, dG_X;
  int *inv_xstar, *xstar_old = NULL;
  double *ys, *yt;                       // dual variables of LAP
#ifdef INTEGER_LAP
  int **D0;
  int *ys_int, *yt_int; 

  ys_int = new int[n];  yt_int = new int[n];
  
  D0 = new int*[n];
  for(j=0;j<n;j++)
    D0[j] = new int[n];
#endif

  bnd_used = -100;

  double lobnd_old;
  MAT *best_X=MNULL, *U_old=MNULL;

  if(save_best) {
    best_X = m_copy(X,MNULL);
  }

  dG = m_get(n,n);
  //  xstar = new int[n];   
  inv_xstar = new int[n];
  ys = new double[n];   yt = new double[n];

  status = CONTINUE;  
  done = false;
  have_U_matrices = false;

  for(ier=1;(ier<=maxit)&&(!done);ier++)
    {
      //      s_lap.go();
#ifdef INTEGER_LAP
      obj = lap_jv(G->me,D0,xstar,inv_xstar,ys,ys_int,yt,yt_int,
		   n,scaling); 
#else
      obj = dlap(n,G->me,xstar,inv_xstar,ys,yt,augment_tol);
#endif
      //      s_lap.stop();

      // (the search direction is D = X - diag(xstar))

      // compute change in gradient
      //      s_dir.go();
      deltaG_alt2(A->me,B->me,C->me,S->me,T->me,G->me,dG->me,X->me,
		 xstar,inv_xstar,dG_X,dG_D,n);
      //      s_dir.stop();
      // printf("%f %f %f %f\n", dG_X, dG_D, S->me[0][0], S->me[0][1]);
      // [G_X was updated at the end of the last iteration]
      // [dG_X was computed by deltaG]
      // [dG_D was computed by deltaG]
#ifdef INTEGER_LAP
      G_D  = mm_trace_alt(G,xstar)  -  G_X;
      obj = obj - G_X;
#else
      G_D  = obj  -  G_X;
      obj  = obj  -  G_X;
#endif

      alpha = -G_D/dG_D;
      alpha = std::max(alpha,0.0);
      alpha = std::min(alpha,1.0);
      one_minus_alpha = 1.0 - alpha; // duh.. :)

      if(!feasible) {alpha = 1.0; feasible = true; } 

      for(j=0;j<n;j++) {
	for(k=0;k<n;k++) 
	  X->me[j][k] *= one_minus_alpha;
      }
      for(j=0;j<n;j++)
	X->me[j][xstar[j]] += alpha;

      C_X = one_minus_alpha * C_X + alpha * mm_trace_alt(C, xstar);

      upbnd = 0.5 * (G_X + alpha * G_D + C_X) + offset;
      bnd = upbnd + obj; 

      // wait to see if we can terminate before updating
      // gradient

#ifdef PRINT
      printf("%3d %f %f %f   %f\n",ier+startit,bnd,upbnd,alpha,G_X);
#endif

      if((ier+1==maxit)&&(improve)) {
	// save U into U_old, bnd into oldbnd
	U_old = m_get(U->m, U->n);
	reducedCosts(G->me, U_old->me, ys, yt, n);
	xstar_old = copy_array(xstar,n);
	lobnd_old = bnd;
      }
      else if((ier==maxit)&&(improve)) {
	// save U into U_old, bnd into oldbnd
	reducedCosts(G->me,U->me,ys,yt,n);
	need_U = false;
	lobnd = bnd;
	have_U_matrices = true;
      }
      else if(lobnd < bnd) {

	lobnd = std::min(bnd,upbnd);//	lobnd = bnd;
	if(save_best&&(lobnd < cutoff)) {
	  m_copy(X, best_X);
	  reducedCosts(G->me,U->me,ys,yt,n);
	  need_U = false;
	}
	
	if(lobnd > cutoff) {
#ifdef PRINT
	  printf("terminated fw early due to bound: ier %d\n",ier);
#endif
	  status = FATHOM;	    done = true;
	  break;
	}
      }
           
      if ((upbnd < cutoff)&&(ier>=maxit_nofathom))
	{
#ifdef PRINT
	  printf("terminated early -- cannot fathom: ier %d\n",ier);
#endif
	  status = NO_FATHOM;    done = true;
	  reducedCosts(G->me,U->me,ys,yt,n);
	  need_U = false;
	  break;
	}
      if (upbnd == bnd) {
#ifdef PRINT
	printf("terminated early -- closed bound gap: ier %d\n",ier);
#endif
	status = NO_FATHOM;	done = true;
	reducedCosts(G->me,U->me,ys,yt,n);
	need_U = false;
	break;
      }
      
      if (alpha <= 0.0)
	{
#ifdef PRINT
	  printf("terminated early -- solved QP exactly at ier %d\n",ier);
#endif
	  alpha = 0.0; 
	  status = NO_FATHOM;	done = true;
	  reducedCosts(G->me,U->me,ys,yt,n);
	  need_U = false;
	  break; 
	}

      if ((ier==maxit)&& need_U && (status == CONTINUE)) {
	reducedCosts(G->me,U->me,ys,yt,n);
      }

      // update gradient (needed for next iteration)
      for(j=0;j<n;j++)
	for(k=0;k<n;k++) {
	  G->me[j][k] += alpha*dG->me[j][k];
	}
      G_X += alpha*(G_D + dG_X + alpha*dG_D);
          }

  if(ier > maxit) ier = maxit;
 
  if((status != FATHOM) && improve && have_U_matrices)
  {
    bnd_used = try_to_improve(U,xstar,lobnd,U_old,xstar_old,lobnd_old,cutoff, imp_steps);
  }

  if(save_best) {    M_FREE(best_X);  }
  if(improve)   {    M_FREE(U_old); delete [] xstar_old;  } 

  M_FREE(dG); 
  delete [] inv_xstar;
  delete [] ys; delete [] yt;

#ifdef INTEGER_LAP
  delete [] ys_int; delete [] yt_int;
  for(j=0;j<n;j++)
    delete [] D0[j];
  delete [] D0;
#endif
}


// may no longer be necessary.
double compute_U(MAT *A,MAT *B,MAT *C,MAT *S,MAT *T,MAT *X,MAT *U1,int n,
	       bool compute_lapU, double &lapU,double scaling,double offset)
{
  MAT *G;
  int **D0;
  int i,j;
  int *p, *ip;
  double *ys, *yt; 
  int *ys_int, *yt_int; 
  double lower_bound,obj,bound;

  G = m_get(n,n);  
  p = new int[n];  ip = new int[n];
  ys_int = new int[n]; yt_int = new int[n];
  ys = new double[n]; yt = new double[n];

  D0 = new int*[n];
  for(i=0;i<n;i++)
    D0[i] = new int[n];

  //  scaling = 50.0;
  gradient(A,B,C,S,T,G,X,n);

  obj = lap_jv(G->me,D0,p,ip,ys,ys_int,yt,yt_int,n,scaling);

  lower_bound =  obj - mm_trace(G,X);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++) {
      U1->me[i][j] = G->me[i][j] - (yt[j] + ys[i]);
      G->me[i][j] = G->me[i][j] - C->me[i][j];
    }

  bound = mm_trace(G,X)/2.0 + mm_trace(C,X) + offset;
  lower_bound = lower_bound + bound;

  M_FREE(G);
  delete [] ip; delete [] p;
  delete [] ys; delete [] yt;
  delete [] ys_int; delete [] yt_int;

  for(i=0;i<n;i++)
    delete [] D0[i];
  delete [] D0;

  return lower_bound;
}


// U = lam*U1 + (1-lam)U2
inline void convex_comb(MAT *U1,MAT* U2,MAT *U,double lam)
{
  for(int i=0;i<U->m;i++)
    for(int j=0;j<U->n;j++)
      U->me[i][j] = lam*U1->me[i][j] + (1-lam)*U2->me[i][j];
}

//#define PRINT_IMPROVE
int try_to_improve(MAT *U_qpb,int *x_U1,double &z_qpb,
		   MAT *U_glb,int *x_U2,double &z_glb,
		   double incumbent,const int steps)
{
  int bound_used;

  //  const int steps = 3;  // make this a parameter
  const int n = U_qpb->m;

  double z1, z2, zl, zbest;

  double Uqpb_x1, Uglb_x1, Uqpb_x2, Uglb_x2, 
    Uqpb_xl, Uglb_xl;
  double d1, d2, lam1, lam2, lam, dlam;

  //  z1 = dlap(n,U1->me,x_U1,temp,s1,t1);
  z1 = z_qpb;// + z1;

  Uqpb_x1 = mm_trace_alt(U_qpb,x_U1);  
  Uglb_x1 = mm_trace_alt(U_glb,x_U1);
  d1 = z_glb - z_qpb + Uglb_x1 - Uqpb_x1;

  // check to see if QPB is best.
  if((d1 < 0)||(z1 > incumbent)) {
    bound_used = QP_BND;

#ifdef PRINT_IMPROVE
    printf("QPB is best: %f + %f < 0\n",z_glb-z_qpb,Uglb_x1-Uqpb_x1);
#endif

    z_qpb = z1;
  }
  else { // check to see if GLB is best

    //z2 = dlap(n,U2->me,x_U2,temp,s2,t2);

    z2 = z_glb;// + z2;
    Uqpb_x2 = mm_trace_alt(U_qpb,x_U2);  
    Uglb_x2 = mm_trace_alt(U_glb,x_U2);
    d2 = z_glb - z_qpb + Uglb_x2 - Uqpb_x2;

    if((d2 > 0)||(z2 > incumbent)) {
#ifdef PRINT_IMPROVE
      printf("GLB is best: %f + %f > 0\n",z_glb-z_qpb,Uglb_x2-Uqpb_x2);
#endif
      bound_used = GLB_BND;

      z_qpb = z2;
      m_copy(U_glb,U_qpb);  // copy new U matrix
    }
    else {

      MAT *Ulam = MNULL, *Ubest = MNULL;
      int *x_Ulam, *temp;
      double *slam, *tlam;

      Ulam   = m_get(n,n);
      Ubest  = m_get(n,n);
      x_Ulam = new int[n];
      temp   = new int[n];
      slam   = new double[n];
      tlam   = new double[n];
      
#ifdef PRINT_IMPROVE
      printf("possibility for improvement: QPB: %f  GLB %f\n",z_qpb,z_glb);
#endif

      bound_used = QP_GLB_BND;
      
      lam1 = 0.0; lam2 = 1.0;
      if(z1 >= z2) {
	//	Ubest = m_copy(U_qpb,Ubest);
	zbest = z1;
      }
      else {
	//	Ubest = m_copy(U_glb,Ubest);
	zbest = z2;
      }
      
      for(int i=0; i<steps; i++) {
	// get a trial lam \in [lam1,lam2]
	lam = (z2 - z1 + lam1*d1 - lam2*d2)/(d1 - d2);
	
	convex_comb(U_glb,U_qpb,Ulam,lam);
	
	// zl = lap_jv(Ulam->me,x_Ulam,slam,tlam,n,scaling);
	zl = dlap(n,Ulam->me,x_Ulam,temp,slam,tlam);
	zl = lam*z_glb + (1-lam)*z_qpb + zl;

	if(zl > zbest) {
	  zbest = zl;
	  reducedCosts(Ulam->me,Ubest->me,slam,tlam,n);
	  if(zbest > incumbent) break;
	}
	
	Uqpb_xl = mm_trace_alt(U_qpb,x_Ulam);  
	Uglb_xl = mm_trace_alt(U_glb,x_Ulam);
	
#ifdef PRINT_IMPROVE
	printf("[%.3f,%.3f] lam = %.3f, bnd = %f\n",lam1,lam2,lam, zl);
#endif
	
	dlam = z_glb - z_qpb + Uglb_xl - Uqpb_xl;
	
	if(dlam > 0) { // move lam1
	  lam1 = lam;
	  z1 = zl;
	  d1 = dlam;
	  copy_array(x_Ulam,x_U1,n);
	}
	else if(dlam < 0) { // move lam2
	  lam2 = lam;
	  z2 = zl;
	  d2 = dlam;
	  copy_array(x_Ulam,x_U2,n);
	}
	if(zl > incumbent) break;
      }

      // copy best of z1, z2 into z_qpb.  also generate appropriate U matrix.
      if((zbest > z_qpb)&&(zbest > z_glb)) {
	z_qpb = zbest;
	m_copy(Ubest,U_qpb);
      }
      else if(z_qpb < z_glb) {
	z_qpb = z_glb;
	m_copy(U_glb,U_qpb);
      }
#ifdef PRINT_IMPROVE
      printf("final bound: %f\n",z_qpb);
#endif

      M_FREE(Ulam);        M_FREE(Ubest);  
      delete [] x_Ulam;    delete [] temp;
      delete [] slam;      delete [] tlam;

    } // end of improvement phase

  }

  return bound_used;
}
