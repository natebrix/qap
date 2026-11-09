#include "fw.h"
#include "qpbnd.h"
#include <math.h>
#include "util.h"
#include "linalg.h"
#include <fstream>
#include <vector>
#include <cstring>

#ifdef TIME_STUFF
stopwatch s_init,s_fw,s_update,s_lap,s_dir;
#endif

// a couple functions only for internal use...
inline double _glb(MAT *A,MAT *B,MAT *C,MAT *U,int *p,double shift);
inline double _qpbnd_parametric(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,int *xstar,
				QPBParameters *p);
int accurate_ier = 5;
double accurate_scaling = 10.0;

int *qpb_heuristic(MAT *X) {
  if (X->m != X->n) {
    printf("ERROR: matrix should be square.\n");
    return NULL;
  }
  int n = X->m;
  int *p, *ip;
  double *ys, *yt;

  MAT *X1 = m_get(X->m, X->n);
  for (int i = 0; i < X1->m; i++) {
    for (int j = 0; j < X1->n; j++) {
      X1->me[i][j] = -X->me[i][j];
    }
  }
  
  p = new int[n];     ip = new int[n];
  ys = new double[n]; yt = new double[n];

  double obj = dlap(n, X1->me, p, ip, ys, yt, 1.0);
  delete[] ip; delete[] ys; delete[] yt;
  return p;
}

/*
============================================================
 compute initial value for gradient G, given that 
 X = E/n, where E is the matrix of ones.

 From the second paper:
 init G: (2/n) A*E*B + C 
   G_ij = (2/n) sum_k r(A)_i B_kj + C_ij
   G_ij = ((2/n)*r(A)_i*c(B)_j) + C_ij

 r(.) is row sum, c(.) is column sum

 this is much more efficient than simply calling 
 gradient() -- this is O(n^2), not O(n^3)
============================================================
*/
void init_gradient(MAT *A,MAT *B,MAT *C,MAT *G,double &C_X,double &G_X,
	       int n)
{
  int i,j,k;
  double *ra, *ai;
  double *cb;
  register double v;

  C_X = 0.0;  G_X = 0.0;

  // ra: row sum of A
  ra = new double[n];
  for(i=0;i<n;i++) {
    ai = A->me[i]; 
    ra[i] = 0.0;
    for(j=0;j<n;j++) {
      ra[i] += ai[j];
    }
  }

  // cb: col sum of B
  // (since B is symmetric could compute row sum of B)
  cb = new double[n];
  for(j=0;j<n;j++) {
    cb[j] = 0.0;
    for(i=0;i<n;i++) {
      cb[j] += B->me[i][j];
    }
  }

  for(i=0;i<n;i++)
    for(j=0;j<n;j++) {

      v = 2.0*ra[i]*cb[j]/n + C->me[i][j];

      G->me[i][j] = v;      G_X += v;
      C_X += C->me[i][j];
    }

  G_X *= 1.0/n;     C_X *= 1.0/n;

  delete [] ra;     delete [] cb;
}

/*
========================================================================

 find_basis(g,h,s,t);

 solve the linear program

  max e's + e't
 subject to
      t_i + s_j <= g_i h_j

  Input: Vectors g and h
 Output: Vectors s and t

 All vectors should be allocated to be of same size.

========================================================================
*/
void find_basis(VEC *w,VEC *u,VEC *s,VEC *t)
{
  int i,m;

  m = w->dim;

  s->ve[0] = 0.0;
  t->ve[0] = w->ve[0]*u->ve[0] - s->ve[0];
  i = 1;
  while(i < m)
    {
      t->ve[i] = w->ve[i-1]* u->ve[i] - s->ve[i-1];
      s->ve[i] =  w->ve[i] * u->ve[i] - t->ve[i];
      if(i+1 < m)
	{
	  s->ve[i+1] = w->ve[i+1]* u->ve[i]  - t->ve[i];  
	  t->ve[i+1] = w->ve[i+1]*u->ve[i+1] - s->ve[i+1];
	}
      i = i + 2;
    }
}

/*
========================================================================
========================================================================
 */
void update_basis(VEC *w,VEC *u,VEC *g_in,VEC *h_in,VEC *s,VEC *t)
{
  VEC *g, *h;
  int i,m;
  m = u->dim;
  
  g = v_copy(g_in,VNULL);
  h = v_copy(h_in,VNULL);

  t->ve[0] = 0;
  s->ve[0] = u->ve[0]*w->ve[0] - t->ve[0];
  i = 0;
  while(i < m-1)
    {
      if(h->ve[i] < g->ve[i])  // then solve for s->ve[2], then t->ve[2]
	{
	  s->ve[i+1] = u->ve[i]*w->ve[i+1] - t->ve[i];
	  i = i+1;
	  t->ve[i] = u->ve[i]*w->ve[i] - s->ve[i];
	 
	  g->ve[i] = g->ve[i] + (g->ve[i-1] - h->ve[i-1]);
	}
      else if(g->ve[i] < h->ve[i])
	{
	  t->ve[i+1] = w->ve[i]*u->ve[i+1] - s->ve[i];
	  i = i + 1;
	  s->ve[i] = w->ve[i]*u->ve[i] - t->ve[i];
	  
	  h->ve[i] = h->ve[i] + (h->ve[i-1] - g->ve[i-1]);
	}
      else
	{
	  t->ve[i+1] = w->ve[i]*u->ve[i+1] - s->ve[i];
	  i = i + 1;
	  s->ve[i] = w->ve[i]*u->ve[i] - t->ve[i];
	}
    }

  V_FREE(g);  V_FREE(h);

}


/*
======================================================================
  Given A and B.

  let A1 = V'*A*V,
      B1 = V'*B*V, where V' means the transpose of matrix V.

  let W,U,r,s satisfy
      W*diag(r)*W' = A1,
      U*diag(s)*U' = B1, 

  where r is sorted in ascending order, s sorted in descending order.

  using r and s, compute g and h.

  Then 
      S = V*W*diag(g)*W'*V'
      T = V*U*diag(h)*U'*V'

  *** NOTE ***
  W and U may be precomputed.  For instance, when using branching
  rules that compute bounds at some of the subproblems, some subproblems
  share the same A or B matrices, and therefore the same W or U matrices.
  In these instances, set need_WU to false to indicate that the values
  of W and U are already known.
  
====================================================================== 
*/

// The eigendecomposition of V'AV (and V'BV) behind W and U (review E13, Nov 9).
// It depends on A alone, and the children of a node that branches on a facility
// all have the parent's A less the same row and column, so a depth-first search
// bounds them with the same A one after another, their subtrees (smaller n) in
// between: one entry per size n, the last A seen at that size, serves every
// sibling. Column branching does the same for B. The cached result is copied
// out (the callers free W, U, w, u), bitwise what the computation gives.
namespace {
struct EigEntry { MAT *key = MNULL; VEC *vals = VNULL; MAT *vecs = MNULL; };
std::vector<EigEntry> eig_cache_A, eig_cache_B;
long eig_hits = 0, eig_misses = 0;

bool same_matrix(const MAT *a, const MAT *b)
{
  if (!a || !b || a->m != b->m || a->n != b->n) return false;
  for (unsigned i = 0; i < a->m; i++)
    if (memcmp(a->me[i], b->me[i], a->n * sizeof(double)) != 0) return false;
  return true;
}

// w (sorted) and W = rows of the eigenvectors of V'AV in that order
void eig_of(MAT *A, bool descending, std::vector<EigEntry> &cache, VEC *&w, MAT *&W)
{
  const int n = A->m;
  if ((int)cache.size() <= n) cache.resize(n + 1);
  EigEntry &e = cache[n];
  if (same_matrix(e.key, A)) {
    eig_hits++;
    w = v_get(n-1);  memcpy(w->ve, e.vals->ve, (n-1) * sizeof(double));
    W = m_get(n-1,n-1);
    for (int i = 0; i < n-1; i++) memcpy(W->me[i], e.vecs->me[i], (n-1) * sizeof(double));
    return;
  }
  eig_misses++;
  std::vector<int> order(n-1);
  for (int k = 0; k < n-1; k++) order[k] = k;
  MAT *A1 = m_get(n-1,n-1);
  vtav(A,A1);
  w = v_get(n-1);
  my_svd(A1,w);
  if (descending) ins_sort_descending(w->ve, order.data(), n-1);
  else ins_sort_ascending(w->ve, order.data(), n-1);
  W = m_get(n-1,n-1);  permrows(A1,order.data(),W);   // W[k] = A1[order[k]]
  M_FREE(A1);
  if (!e.key || (int)e.key->m != n) {
    M_FREE(e.key); V_FREE(e.vals); M_FREE(e.vecs);
    e.key = m_get(n,n); e.vals = v_get(n-1); e.vecs = m_get(n-1,n-1);
  }
  for (int i = 0; i < n; i++) memcpy(e.key->me[i], A->me[i], n * sizeof(double));
  memcpy(e.vals->ve, w->ve, (n-1) * sizeof(double));
  for (int i = 0; i < n-1; i++) memcpy(e.vecs->me[i], W->me[i], (n-1) * sizeof(double));
}
} // namespace

void init_ST_cache_stats(long &hits, long &misses) { hits = eig_hits; misses = eig_misses; }

double init_ST(MAT *A,MAT *B,MAT *C,
	       VEC *&w,VEC *&u,MAT *&W,MAT *&U,VEC *&s,VEC *&t,
	       MAT *&S,MAT *&T,bool need_WU)
{
  int n = A->m;

  double shift;
  MAT *S1=MNULL, *T1=MNULL;

  if (need_WU) {
    /* w descending with W, u ascending with U (the eigenvalue bound pairs them) */
    eig_of(A, true, eig_cache_A, w, W);
    eig_of(B, false, eig_cache_B, u, U);
  }

  s = v_get(n-1);  t = v_get(n-1);
  find_basis(w,u,s,t);

  S1 = mtvm_mlt(W,s,MNULL);
  T1 = mtvm_mlt(U,t,MNULL);

  /* compute S = V*S1*V' and T = V*T1*V' */
  S = m_get(S1->m+1,S1->n+1);
  vavt(S1,S);

  T = m_get(T1->m+1,T1->n+1);
  vavt(T1,T);

  shift = in_prod(w,u);

  /* free memory */
  M_FREE(S1);  M_FREE(T1);  
  //printf("%f\n", w->ve[0]);
  //printf("%f\n", shift);
  //nice_array_print(w->ve, n-1);
  //nice_array_print(u->ve, n-1); // todo
  return shift;
}

/*
============================================================
update S and T.  Also, update G, the gradient of 

   f(X) = tr AXBX' - SX - TX + CX'

as follows:  (where S',T',G' are the new values of
S,T and the gradient)

   S' = S + dS
   T' = T + dT
   G' = 2(A*X*B - S'*X - X*T') + C 
      = 2(A*X*B - S*X - X*T) + C - 2(dS*X + X*dT)
      = G - 2(dS*X + X*dT)

 (that's 2 matrix mults instead of the 4 required
  by gradient() )
============================================================
*/
void update_STG(VEC *w,VEC *u,MAT *W, MAT *U,VEC *s, VEC *t, 
		MAT *G,MAT *X,MAT *S, MAT *T,
		double factor)
{
  double val;
  unsigned int i,j;
  int m,n;
  MAT *M1 = MNULL, *M2=MNULL, *X1;
  VEC *g, *h;
  VEC *s_new, *t_new;

  double shift;

  n = X->m;
  m = W->m;

  g = v_get(m);
  h = v_get(m);

  s_new = v_get(m);
  t_new = v_get(m);

  //(X1 = V'*X*V;)
  //    X1 = m_get(n-1,n-1);
  X1 = m_get(n,n);
  vtav(X,X1);

  //  g = diag(W'*X1*X1'*W);
  M1 = my_m_mlt(W,X1,M1);
  for(i=0;i<m;i++){
    val = 0.0;
    for(int j=0;j<M1->m;j++)
      val += M1->me[i][j]*M1->me[i][j];
    g->ve[i] = val;
  }

  //  h_i = coefficient of t_i in tr(X1 T1 X1'), T1 = U' diag(t) U with the
  //  eigenvectors u_i as the ROWS of U (my_svd's convention, as mtvm_mlt
  //  below uses it): h_i = ||X1 u_i||^2, the squared norm of row i of U X1'.
  //  Until Nov 8 2026 this took the row norms of U' X1 instead, i.e. the
  //  wrong orientation of U and the transpose of X1, so update_basis chose
  //  (s, t) for the wrong functional (the bound stayed valid: any feasible
  //  (s, t) is; see docs/review_nov8.md B1). g above is right: tr(X1' S1 X1)
  //  = sum_i s_i ||X1' w_i||^2 and my_m_mlt(W, X1) has row i = w_i' X1.
  for(i=0;i<m;i++){
    val = 0.0;
    for(int k=0;k<m;k++){
      double r = 0.0;
      for(j=0;j<(unsigned)m;j++)
	r += X1->me[k][j]*U->me[i][j];
      val += r*r;
    }
    h->ve[i] = val;
  }

  shift = in_prod(g,s) + in_prod(h,t);

  update_basis(w,u,g,h,s_new,t_new);

  for(i=0;i<s_new->dim;i++) {
    s->ve[i] = (1-factor)*s->ve[i] + factor*s_new->ve[i];
    t->ve[i] = (1-factor)*t->ve[i] + factor*t_new->ve[i];
  }

  // G - 2(dS*X + X*dT)
  // save old S for later use
  M2 = m_copy(S,M2);

  M1 = mtvm_mlt(W,s,M1);
  // (S = V*M1*V')
  vavt(M1,S);

  // update G
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      X1->me[i][j] = (S->me[i][j] - M2->me[i][j]);
  M2 = my_m_mlt(X1,X,M2); // M2 = dS*X
  for(i=0;i<n;i++)
    for(j=0;j<n;j++) {
      G->me[i][j] -= 2*M2->me[i][j];
    }

  // save old T for later use
  M2 = m_copy(T,M2);

  M1 = mtvm_mlt(U,t,M1);
  vavt(M1,T);

  // update G
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      X1->me[i][j] = (T->me[i][j] - M2->me[i][j]);
  M2 = my_m_mlt(X,X1,M2); // M2 = X*dT
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      G->me[i][j] -= 2*M2->me[i][j];

  M_FREE(M1);  M_FREE(M2);
  V_FREE(g);  V_FREE(h);
  M_FREE(X1);

  V_FREE(s_new);  V_FREE(t_new);
}

/*
========================================================================

 bound = qpbnd(A,B,C,X,U,...);

 Compute a lower bound on the QAP (A,B,C) using the quadratic
 programming bound given in "...".

 The bound is computed by solving the quadratic program:
  z = min (AXB - SX - XT)X' + CX'
        Xe  = e
        X'e = e
        x_ij >= 0
 where S and T are computed using (A,B,C).

 The returned value is z + shift.  The QP is solved iteratively
 using the Frank-Wolfe algorithm (see fw3.c for details).  Maxier
 controls the maximum number of iterations performed.  If at any
 time, the computed bound is larger than inc, then computation
 is halted, and the current value of the bound is returned.

 [This function usually doesn't get called...qpbnd_parametric()
  usually gives better performance, see below]

========================================================================
*/

QPBWork &qpb_work(int n)
{
  static std::vector<QPBWork> ws;
  if ((int)ws.size() <= n) ws.resize(n + 1, QPBWork{MNULL, 0});
  QPBWork &w = ws[n];
  if (!w.G) { w.G = m_get(n, n); w.xstar = new int[n]; }
  else m_zero(w.G);
  for (int i = 0; i < n; i++) w.xstar[i] = 0;
  return w;
}

double qpbnd(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
	     QPBParameters *p)

{
  int status;
  int iter,total_iter;
  double eig;
  double lo_bound, up_bound;
  MAT *S,*T;
  MAT *G;

  MAT *W,*U;
  VEC *s,*t,*w,*u;

  int *xstar;

  double C_X, G_X;

  up_bound =  1.0e99;
  lo_bound = -1.0e99;

  S = MNULL;  T = MNULL;
  W = MNULL;  U = MNULL;
  w = VNULL;  u = VNULL;  s = VNULL;  t = VNULL;

  eig = init_ST(A,B,C,w,u,W,U,s,t,S,T,true);
  total_iter = 0;

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
     p->maxier,p->maxier_nofathom,p->shift + eig,p->inc - p->cutoff_slack,p->feasible,
     status,up_bound,lo_bound,p->lap_scale,p->save_best,p->improve,
     p->bound_used);
  total_iter += iter;

  if(lo_bound <= p->inc - p->cutoff_slack)
    lo_bound=compute_U(A,B,C,S,T,X,U1,A->m,p->compute_lapU,p->lapU,p->lap_scale,
		       p->shift+eig);
  
  p->fw_iter += total_iter;

  M_FREE(W);  M_FREE(U);
  M_FREE(S);  M_FREE(T);
  V_FREE(w);  V_FREE(u);
  V_FREE(s);  V_FREE(t);

  return lo_bound;
}

void ampl_mat_param(std::ofstream &f, const char *name,  MAT *A)
{
  matrix_to_ampl_file(f, name,  A->me, A->m, A->n);
  f << std::endl;
}

void ampl_const_param(std::ofstream &f, const char *name, double d)
{
  f << "param " << name << " := " << d << ";" << std::endl;
}

void text_mat_param(std::ofstream &f, const char *name,  MAT *A)
{
  matrix_to_file(f, A->me, A->m, A->n);
  f << std::endl;
}

void text_const_param(std::ofstream &f, const char *name, double d)
{
  f << d << std::endl;
}


void write_ampl_qpb(MAT *A, MAT *B, MAT *C, MAT *S, MAT *T, double d)
{
  std::ofstream f;
  std::string filename = "qpb_dat.dat";
  std::cout << "dumping QPB info to file " << filename << std::endl;
  f.open(filename);
  //  dump_qpb(A, B, C, S, T);
  ampl_const_param(f, "n", A->n);
  ampl_mat_param(f, "a", A);
  ampl_mat_param(f, "b", B);
  ampl_mat_param(f, "c", C);
  ampl_mat_param(f, "s", S);
  ampl_mat_param(f, "t", T);
  ampl_const_param(f, "d", d);
  f.close();
}

void write_text_qpb(MAT *A, MAT *B, MAT *C, MAT *S, MAT *T, double d)
{
  std::ofstream f;
  std::string filename = "qpb_dat.txt";
  std::cout << "dumping QPB info to file " << filename << std::endl;
  f.open(filename);
  text_const_param(f, "n", A->n);
  text_mat_param(f, "a", A);
  text_mat_param(f, "b", B);
  text_mat_param(f, "c", C);
  text_mat_param(f, "s", S);
  text_mat_param(f, "t", T);
  text_const_param(f, "d", d);
  f.close();
}

/*
========================================================================

 bound = qpbnd_parametric(A,B,C,X,U,p);

 Same as qpbnd(), except that the matrices S and T are updated every
 p->step iterations.  Doing so usually results in a tighter bound.  A
 decent value for step is 30.

======================================================================== 
*/
inline double _qpbnd_parametric(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,int *xstar,
				QPBParameters *p)
{
  MAT *S=MNULL, *T=MNULL;
  MAT *G=MNULL;
  MAT *W=MNULL, *U=MNULL;
  VEC *s=VNULL, *t=VNULL, *w=VNULL, *u=VNULL;

  double C_X, G_X;
  double eig, lo_bound, up_bound;

  bool feasible = p->feasible;
  int status;
  int total_iter, iter;

  //  s_init.go();
  eig = init_ST(A,B,C,w,u,W,U,s,t,S,T,true);
  //  s_init.stop();
  //  s_fw.go();

  total_iter = 0;

  up_bound =  1.0e99;
  lo_bound =  INIT_LOBND;

  // init gradient G (the per-size workspace's, zeroed)
  G = qpb_work(X->n).G;

  // init_gradient() is the gradient at the uniform X only. (This used to
  // test X[0][0] == X[0][1], which a warm start from the parent's X can
  // also satisfy.)
  if(is_uniform(X)) {
    init_gradient(A, B, C, G, C_X, G_X, X->n);
    if (p->debug) {
      write_ampl_qpb(A, B, C, S, T, p->shift+eig);
      write_text_qpb(A, B, C, S, T, p->shift+eig);
    }
  }
  else {
    gradient(A, B, C, S, T, G, X, X->n);
    C_X = mm_trace(C, X);
    G_X = mm_trace(G, X);
  }

  status = CONTINUE;

  do {
    fw(A,B,C,S,T,G,X,U1,xstar,C_X,G_X,A->m,iter,total_iter,
       std::min(p->step,p->maxier - total_iter - accurate_ier),
       p->maxier_nofathom - total_iter,
       p->shift+eig, p->inc - p->cutoff_slack, feasible,
       status, up_bound, lo_bound, p->lap_scale, false, false, 
       p->bound_used);
    total_iter += iter;

    //    s_fw.stop();
    //    s_update.go();

    // no need to update S,T,G if we are about to quit.
    if (status == CONTINUE) {
      if (total_iter < p->maxier - accurate_ier) {
	update_STG(w,u,W,U,s,t,G,X,S,T,1.0);
	G_X = mm_trace(G, X);
      }
    }

    //    s_update.stop();
    //    s_fw.go();
  }  while((total_iter < p->maxier - accurate_ier) && (status==CONTINUE));

  if (status == CONTINUE) {
    if (accurate_ier > 0) {
      lo_bound = INIT_LOBND;
      fw(A,B,C,S,T,G,X,U1,xstar,C_X,G_X,A->m,iter,total_iter,accurate_ier,
	 1000,p->shift + eig, p->inc - p->cutoff_slack, feasible,
	 status, up_bound, lo_bound, accurate_scaling,p->save_best,
	 p->improve,p->bound_used,true);
      total_iter += iter;
    }
    else {
      lo_bound = compute_U(A,B,C,S,T,X,U1,A->m,p->compute_lapU,p->lapU,
			   accurate_scaling,p->shift+eig);
    }
  }
  p->fw_iter += total_iter;
  
  M_FREE(W);  M_FREE(U);
  M_FREE(S);  M_FREE(T);
  V_FREE(w);  V_FREE(u);
  V_FREE(s);  V_FREE(t);
  //  s_fw.stop();
  /*
  printf("init %.4f fw %.4f up %.4f   lap %.4f  dir %.4f\n",
	 s_init.elapsed(),s_fw.elapsed(),
	 s_update.elapsed(),s_lap.elapsed(),s_dir.elapsed());
  */

  return lo_bound;
}

double qpbnd_parametric(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
			QPBParameters *p)
{
  double bnd;
  int *xstar = new int[X->n];

  // _qpbnd_parametric returns the best bound of its last accurate_ier
  // iterations; an earlier S,T phase often had a better one. Return the
  // best (bound, U) pair over all iterations instead (see fw_best_begin).
  // Not with the improvement procedure, which builds its own U.
  MAT *U_best = m_get(U1->m, U1->n);

  fw_patience_begin(p->nofathom_patience);
  fw_best_begin(U_best);
  bnd = _qpbnd_parametric(A,B,C,X,U1,xstar,p);
  const double best = fw_best_end();
  fw_patience_end();

  if (!p->improve && (best > bnd)) {
    bnd = best;
    m_copy(U_best, U1);
  }

  M_FREE(U_best);
  delete [] xstar;

  return bnd;
}

double qp_glb_bnd(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
		  QPBParameters *p)
{
  const int imp_steps = 3;
  double z_glb, z;
  int *x_qpb = new int[X->n];
  int *x_glb = new int[X->n];
  MAT *U_glb = m_copy(U1,MNULL);


  p->bound_used = GLB_BND;
  z_glb = _glb(A,B,C,U_glb,x_glb,p->shift);
  z = z_glb;

  if(z < p->inc) {
    p->bound_used = QP_BND;
    bool save_best = p->save_best;
    p->save_best = false;
    z = _qpbnd_parametric(A,B,C,X,U1,x_qpb,p);
    p->save_best = save_best;

    if(z < p->inc) 
      p->bound_used = try_to_improve(U1,x_qpb,z,U_glb,x_glb,z_glb,
				     p->inc,imp_steps);
  }
    

  delete [] x_qpb;
  delete [] x_glb;
  M_FREE(U_glb);

  return z;
}

/*
========================================================================
 compute the projection bound as described in ... 
========================================================================
*/
double pbnd(MAT *A,MAT *B,MAT *C,double shift)
{
  int i,j,n;
  MAT *T1;
  MAT *A1,*B1,*V;
  VEC *gsort,*hsort;
  VEC *ra, *rb;
  double sa,sb,bnd;
  
  bnd = 0.0;
  sa = 0.0;
  sb = 0.0;
  n = A->m;

  // get eigs of a1, b1
  V = projectionMatrix(n);

  /* compute a1 = v'av and b1 = v'bv */
  T1 = my_m_mlt(A,V,MNULL); 
  A1 = my_mtrm_mlt(V,T1,MNULL); 

  M_FREE(T1);
  T1 = my_m_mlt(B,V,MNULL); 
  B1 = my_mtrm_mlt(V,T1,MNULL); 

  gsort = v_get(n); gsort->dim--;
  hsort = v_get(n); hsort->dim--;
  my_svd(A1,gsort);
  my_svd(B1,hsort);

  ins_sort_descending(gsort->ve, NULL, n);
  ins_sort_ascending(hsort->ve, NULL, n);

  // get row sums of a and b
  ra = v_get(n+1);  ra->dim--;
  rb = v_get(n+1);  rb->dim--;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	ra->ve[i] += A->me[i][j];
	rb->ve[i] += B->me[i][j];
      }
  //  ra = v_ins_sort_descending(ra, PNULL);
  //  rb = v_ins_sort_ascending(rb, PNULL);

  // form this lap and solve it...
  // D = (2/n)*ra*rb' + c;  
  T1 = m_resize(T1,n,n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      T1->me[i][j] = 2.0*ra->ve[i]*rb->ve[j]/n + C->me[i][j];
  {
    std::vector<int> p(n), ip(n);
    std::vector<double> ys(n), yt(n);
    bnd = dlap(n, T1->me, p.data(), ip.data(), ys.data(), yt.data(), 1.0);
  }

  // sum up a and b.

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      sa += A->me[i][j];
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      sb += B->me[i][j];

  // printf("%f %f %f\n",in_prod(gsort,hsort),bnd,-sa*sb/(n*n));
  bnd = -sa*sb/(n*n) + in_prod(gsort,hsort) + bnd + shift;

  //printf("%f %f %f\n",-sa*sb/(n*n),in_prod(gsort,hsort),2*in_prod(ra,rb)/n);
  //bnd = -sa*sb/(n*n) + in_prod(gsort,hsort) + 2*in_prod(ra,rb)/n + shift;
  
  M_FREE(A1);  M_FREE(B1);
  M_FREE(T1);  M_FREE(V);
  V_FREE(ra);
  V_FREE(rb);
  V_FREE(gsort);
  V_FREE(hsort);

  return bnd;
}

// gilmore lawler bound.  does not modify (A,B,C).
inline double _glb(MAT *A,MAT *B,MAT *C,MAT *U,int *p,double shift)
{
  double temp;
  int i,j,k;
  int n = A->m;
  MAT *L, *br;
  double *ar_i;
  int *ip;
  double *ys, *yt;


  L = m_get(n,n);
  ar_i = new double[n-1];
  br = m_get(n,n-1);
  
  ip = new int[n];
  ys = new double[n];  yt = new double[n];

  // copy all elements of B_i into br_i except B_ii
  for(i=0;i<n;i++)
    {
      for(j=0;j<i;j++)
	br->me[i][j] = B->me[i][j];
      for(j++;j<n;j++)
	br->me[i][j-1] = B->me[i][j];
      ins_sort_descending(br->me[i],NULL,n-1);
    }

  // l_ij = a_ii b_jj + <a'_i,b'_j>_ + c_ij

  for(i=0;i<n;i++)
    {
      for(j=0;j<i;j++)
	ar_i[j] = A->me[i][j];
      for(j++;j<n;j++)
	ar_i[j-1] = A->me[i][j];
      ins_sort_ascending(ar_i,NULL,n-1);

      for(j=0;j<n;j++)
	{
	  temp = A->me[i][i]*B->me[j][j] + C->me[i][j];
	  for(k=0;k<n-1;k++)
	    {
	      temp += ar_i[k]*br->me[j][k];
	    }
	  L->me[i][j] = temp;
	}

    }

  // solve LAP(L) to get bound.
  temp = dlap(n, L->me,p,ip,ys,yt) + shift;

  reducedCosts(L->me, U->me, ys, yt, n);

  M_FREE(L);     M_FREE(br);
  delete [] ip;  delete [] ar_i;
  delete [] ys;  delete [] yt;
  return temp;
}

double glbnd(MAT *A,MAT *B,MAT *C,MAT *U,double shift)
{
  double bnd;
  int *xstar = new int[U->n];
  
  bnd = _glb(A,B,C,U,xstar,shift);
  
  delete [] xstar;
  
  return bnd;
}

