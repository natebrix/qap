#include <algorithm>
#include <vector>
/*
============================================================
              (see branch.h for details)
============================================================
*/

// 12/5/2022

// This is honestly the weakest part of the implementation because it
// is poorly factored and overly complex.
//
// Branch functions use criteria to select from rows and
// columns. There are two functions per branching technique - one that
// selects from rows only, and another that looks at both rows and
// columns.
//
// The branch functions accept a U matrix. Regrettably, the U matrix
// often serves two purposes 1) as the thing to examine to make the
// branching decision, 2) as a matrix that is later used for
// fathoming.

#include "branch.h"
#include "util.h"
#include "linalg.h"
#include "qpbnd.h"
#include "qpbnd2.h"
#include "lgb_train.h"
#include <fstream>
#include <string>
#include <math.h>

#define INCORRECT_WEIGHT
#define BRANCH_USES_QPBND

void set_branch_parameters(BranchParameters *bp, BNBParameters *param,
			   bool use_J, bool reorder_children) {
  bp->use_best = param->use_best;
  bp->how = param->rule;
  bp->use_J = use_J;
  bp->reorder_children = reorder_children;
  bp->which_bound = param->which_bound;
}


void determine_J3(const QAPAssignment &assign, int *&myJ3, int &myJ3_dim);

//#define WARMSTART

/*
============================================================
          GLOBAL VARIABLES USED IN BRANCH
*/

int *J1 = NULL;
int J1_dim = 0;

int J_dim = 0;     // number of J2/J3 sets

int **J2 = NULL;   // list of J2 sets
int *J2_dim = 0;   // length of each J2 set

int **J3 = NULL;   // list of J3 sets
int *J3_dim = 0;   // length of each J3 set

int *J3_default = NULL;
int J3_default_dim = 0;

// for boost
int n_feat = 16; // number of features
double *ui;
double *ufill;
MAT *lgb_U;
double *lgb_coeff;
int n_lgb_coeff;

double *read_array(std::ifstream &f, int &n) {
  f >> n;
  double *d = new double[n];
  for (int i = 0; i < n; i++) {
    f >> d[i];
  }
  return d;
}

void read_lgb_coeffs() {
  int n;
  std::ifstream f;
  f.open("lgb_init.txt");
  lgb_coeff = read_array(f, n_lgb_coeff);
  f.close();
}

/*
============================================================
*/


/*
============================================================
branch_init(n);
branch_shutdown();

initialize/deallocate global variables associated with branching
routines.  In particular, initialize the default J3 set to 1..n, that
is, the set of columns that will be looked at when no more symmetry
exists in the problem.

============================================================ 
*/
void branch_init(int n)
{
  if(!use_weighted_branching) 
    std::cout << "[not using weight for branching] " << std::endl;
  else
    std::cout << "branching_weight = " << branching_weight << std::endl;
  J_dim = 0; // only one J2 set

  J1 = NULL;
  J1_dim = 0;

  J2 = NULL;
  J2_dim = 0;

  J3 = NULL;
  J3_dim = 0;
  
  J3_default_dim = n;
  J3_default = new int[J3_default_dim];
  for(int i=0;i<J3_default_dim;i++)
    J3_default[i] = i;

  ui = new double[n];
  ufill = new double[n_feat];
  lgb_U = m_get(n, n);

  read_lgb_coeffs();
}

void branch_shutdown()
{
  for(int i=0;i<J_dim;i++) {
    delete [] J2[i];
    delete [] J3[i];
  }

  delete [] J1;
  delete [] J2;
  delete [] J2_dim;
  delete [] J3_default;
  delete [] J3;
  delete [] J3_dim;
  delete [] ui;
  delete [] ufill;
  M_FREE(lgb_U);
  delete [] lgb_coeff;
}

/*
============================================================
read_symmetry("nug12.sym", n);

File format:                  |
                              |
|J1|                          | (J1 = set of columns to consider
J1                            |       at root of tree.)
number of J2/J3 sets          |
|J2|, |J3|                    |
J2                            | (if all assigned locations are in
J3                            |  the J2 set, then only need to 
|J2|, |J3|                    |  create children for elements of J3)
J2                            |
J3                            |
.                             |
.                             |
.                             |
============================================================
*/
// A malformed file used to be accepted silently (nug22.sym lacked its
// "|J2| |J3|" line until Nov 8 2026, so J3 held garbage and the search
// branched on a fraction of the locations: a wrong proof). Every count and
// entry is now checked, each list must be strictly increasing (branch() and
// row_sum() merge-scan them), entries must lie in 1..n, and the file must
// hold nothing else. Any failure stops the run: a symmetry file that cannot
// be trusted must not be worked around by guessing.
static void sym_fail(const char *filename, const std::string &what)
{
  std::cerr << "symmetry file [" << filename << "]: " << what << std::endl;
  exit(1);
}

static void read_sym_list(std::ifstream &f, const char *filename, const char *name,
			  int *v, int dim, int n)
{
  for (int j = 0; j < dim; j++) {
    if (!(f >> v[j])) sym_fail(filename, std::string("short file while reading ") + name);
    if (v[j] < 1 || v[j] > n)
      sym_fail(filename, std::string("entry ") + std::to_string(v[j]) + " of " + name +
	       " is outside 1.." + std::to_string(n));
    v[j]--;
    if (j > 0 && v[j] <= v[j-1])
      sym_fail(filename, std::string(name) + " is not strictly increasing");
  }
}

static int read_sym_count(std::ifstream &f, const char *filename, const char *name, int n)
{
  int d;
  if (!(f >> d)) sym_fail(filename, std::string("short file while reading ") + name);
  if (d < 0 || d > n)
    sym_fail(filename, std::string(name) + " = " + std::to_string(d) + " is outside 0.." + std::to_string(n));
  return d;
}

bool read_symmetry(const char *filename, int n)
{
  std::ifstream f;
  f.open(filename);
  if (!f) return false;   // no file: no symmetry, the caller reports it

  J1_dim = read_sym_count(f, filename, "|J1|", n);
  J1 = new int[J1_dim];
  read_sym_list(f, filename, "J1", J1, J1_dim, n);

  J_dim = read_sym_count(f, filename, "the number of J2/J3 sets", n);
  J2 = new int *[J_dim];
  J2_dim = new int[J_dim];
  J3 = new int *[J_dim];
  J3_dim = new int[J_dim];

  for (int i = 0; i < J_dim; i++) {
    J2_dim[i] = read_sym_count(f, filename, "|J2|", n);
    J3_dim[i] = read_sym_count(f, filename, "|J3|", n);
    J2[i] = new int[J2_dim[i]];
    read_sym_list(f, filename, "J2", J2[i], J2_dim[i], n);
    J3[i] = new int[J3_dim[i]];
    read_sym_list(f, filename, "J3", J3[i], J3_dim[i], n);
  }
  std::string extra;
  if (f >> extra) sym_fail(filename, "unexpected token [" + extra + "] after the last J3 list");
  f.close();
  return true;
}

/*
============================================================
============================================================
 */
void reduced_B(QAP *qap,QAPAssignment &assign,int *J3,int J3_dim,
	       MAT **&Bs, MAT **&Wbs, VEC **&wbs,int &wb_sz,int m)
{
  int j, dj, k, J3_count;

    // build up a table of wb,Wb
    // wb_sz = # of columns we will actually look at.
    
    wb_sz = 0;  J3_count = 0;
    for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
      if(assign.pinv[j] <= NONE) dj++;
      if(j != J3[J3_count]) continue;
      J3_count++;
      if(assign.pinv[j] > NONE) continue;
      
      wb_sz++;
    }
    
    Bs = new MAT *[wb_sz];
    Wbs = new MAT *[wb_sz];
    wbs = new VEC *[wb_sz];
    for(k=0;k<wb_sz;k++) {
      Bs[k] = m_get(m-1,m-1);
      Wbs[k] = m_get(m-1,m-1);
      wbs[k] = v_get(m-1);
    }
    
    k = 0;  J3_count = 0;
    for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
      if(assign.pinv[j] <= NONE) dj++;
      if(j != J3[J3_count]) continue;
      J3_count++;
      if(assign.pinv[j] > NONE) continue;


      // reduced B is needed for this index.

      assign.pinv[j] = 0;
      reduceB(qap->B,Bs[k],assign.pinv,m-1,assign.n);
      assign.pinv[j] = NONE;
      
      init_w0(Bs[k],wbs[k],Wbs[k],true);
      k++;
    }
}

/*
============================================================
  merge_U(U1,bound1,U2,bound2);

============================================================
*/
void merge_U(MAT *U1, double bound1, MAT *U2, double bound2)
{
  int i,j;

  for(i=0;i<U1->m;i++)
    for(j=0;j<U1->n;j++)
      if(bound2 + U2->me[i][j] > bound1 + U1->me[i][j])
	U1->me[i][j] = bound2 + U2->me[i][j] - bound1;
}

/*
============================================================
  get_row_in_assign(assign,row_in_U,row_in_assign);

  We have decided that 'row' is the row of D that
  is best to branch on.
  But we also need to know which row of the original
  QAP this corresponds to.  This will be stored in row_in_assign.

  Example: suppose n = 5, and the assignment 2->3 has been made.
  Then U will be of size m = 4.  Suppose we have determined that
  we want to branch on row 3 of U.  get_row_in_assign() will 
  determine that this corresponds to row 4 of the original QAP.
============================================================
*/
inline void get_row_in_assign(const QAPAssignment &assign,
			      int &row_in_U, int &row_in_assign)
{
  int j;
  j=0; row_in_assign = 0;
  do{
    if(assign.p[row_in_assign]<=NONE)
      j++;
    row_in_assign++;
  } while((j<=row_in_U)&&(row_in_assign<assign.n));
  row_in_assign--;
}

inline void get_col_in_assign(const QAPAssignment &assign,
			      int &col_in_U, int &col_in_assign)
{
  int j;
  j=0; col_in_assign = 0;
  do{
    if(assign.pinv[col_in_assign]<=NONE)
      j++;
    col_in_assign++;
  } while((j<=col_in_U)&&(col_in_assign<assign.n));
  col_in_assign--;
}

/*
============================================================ 
  find largest element in array v, and store it in best.  
  v[index]==best.
============================================================ 
 */
inline void find_max(double *v, double &best, int &index, int m)
{
  best = v[0];
  index = 0;
  for (int i = 1; i < m; i++) {
    if (v[i] > best) {
      best = v[i];
      index = i;
    }
  }
}

/*
============================================================
Find permutation most closely matching matrix X, by solving
a linear assignment problem.

The contents of X are modified.
============================================================
*/
void find_best_match(MAT *X, int *p)
{
  // might want to use FP LAP solve...
  lap_jv(X->me,p,X->m,-100.0);
}

/*
============================================================ 
  find largest element in array f, and return index of this element.  
  in case of tie, look at v.
============================================================ 
 */
inline void find_most_fathoms(int *f,double *v, int &best_f, double &best_v,
			      int &index, int m)
{
  index = 0;
  best_f = f[0];
  best_v = v[0];
  for(int i=1;i<m;i++) {
    if((f[i] > best_f) ||
       ((v[i] > best_v)&&(f[i]==best_f))) 
      {
	best_f = f[i];
	best_v = v[i];
	index = i;
      }
  }
}

/*
============================================================ 
 row_sum(assign,J3,J3_dim,U,u_row_sum);

  compute row sums of U, storing result in previously 
  allocated u_row_sum
============================================================ 
 */
inline void row_sum(const QAPAssignment &assign,int *J3,int J3_dim,
			MAT *U,double *u_row_sum)
{
  int i,di,j,dj,J3_count;
  di = -1;
  for (i = 0; i < assign.n; i++) {
    if(assign.p[i] > NONE) continue;
    di++; // working on row di of U.

    dj = -1; J3_count = 0;
    for(j=0; (j<assign.n) && (J3_count < J3_dim); j++) {
	
      // skip over columns not in J3.
      if(assign.pinv[j] <= NONE) dj++;
      if(j != J3[J3_count]) continue;
      J3_count++;
      
      // skip over columns that have already been assigned
      if(assign.pinv[j] > NONE) continue;
      
      // ignore disallowed assignments.
      if(assign.disallowed(i,j)) continue;
      //      printf("u set di=%d, dj = %d, Usize = %d %d\n", di, dj, U->m, U->n);
      u_row_sum[di] += U->me[di][dj];
    }
  }
}

void mean_center(MAT *U, MAT *Uc) {
  int di,dj;
  double u_mean = 0.0;
  for (di = 0; di < U->m; di++) {
    for (dj = 0; dj < U->n; dj++) {
      u_mean += U->me[di][dj];
    }
  }
  u_mean /= (U->m * U->n);
  if (abs(u_mean) <= 1e-6) {
    u_mean = 1e-6;
  }
  for (di = 0; di < U->m; di++) {
    for (dj = 0; dj < U->n; dj++) {
      Uc->me[di][dj] = U->me[di][dj] / u_mean;
    }
  }
}

void stretch(double *u_small, int n_small, double *u_big, int n_big)
{
  // can I stretch in one pass?
  // I don't think I can because I am overwriting.
  // I need to warrays
  // ii = list(map(int, (n_big/n_small) * np.array(range(n_small+1))))

  for (int i = 0; i < n_big; i++) {
    u_big[i] = u_small[(int)(1.0 * i * n_small / n_big)];
  }
}

//     return np.concatenate([a_big[:n_small-n_small//2], a_big[-n_small//2:]])
void squeeze(double *u_big, int n_big, double *u_small, int n_small) {
  int half = n_small / 2;
  for (int i = 0; i < half; i++) {
    u_small[i] = u_big[i];
  }
  int offset = n_big - half;
  for (int i = half; i < n_small; i++) {
    u_small[i] = u_big[i + offset];
  }
}

void force_to_size(double *u_from, int n_from, double *u_to, int n_to)
{
  if (n_from < n_to) {
    stretch(u_from, n_from, u_to, n_to);
  }
  else if (n_from > n_to) {
    squeeze(u_from, n_from, u_to, n_to);
  }
  else {
    for (int i = 0; i < n_from; i++) {
      u_to[i] = u_from[i];
    }
  }
}

// things to try:
// - max fathom trick (look at this code)
// - weight observations
// - try to really debug why the performance is worse
// - better training data
// - just do a damn regression

double row_variance(MAT *U, int i) {
  int n = U->m;
  double s = 0;
  for (int j = 0; j < n; j++) {
    s += U->me[i][j];
  }
  double m = s / n;
  s = 0;
  for (int j = 0; j < n; j++) {
    s += (U->me[i][j] - s) * (U->me[i][j] - s);
  }
  s = s / (n - 1);
  return s;
}

double col_variance(MAT *U, int j) {
  int n = U->m;
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += U->me[i][j];
  }
  double m = s / n;
  s = 0;
  for (int i = 0; i < n; i++) {
    s += (U->me[i][j] - s) * (U->me[i][j] - s);
  }
  s = s / (n - 1);
  return s;
}


void select_index_lgb(const QAPAssignment &assign, bool can_branch_on_cols,
		      int *J3, int J3_dim,
		      MAT *U, double bound, QPBParameters *qpbp,
		      int &ind_in_assign, int &ind_in_U, bool &chose_row)
{
  //nice_matrix_print(U->me, U->m, U->n);
  int J3_count;
  int n = U->n;
  mean_center(U, lgb_U);
  //bound_ratio(U, lgb_U, bound, qpbp->inc); // todo are we doing this right???
  /*for (i = 0 ; i < U->m ; i++) {
    for (j = 0; j < U->n; j++) {
      lgb_U->me[i][j] = U->me[i][j];
    }
    }*/
  int *f_row_count = new_zero<int>(n);
  int *f_col_count = new_zero<int>(n);
  double *u_row_sum = new_zero<double>(n);
  double *u_col_sum = new_zero<double>(n);
  ind_in_U = ind_in_assign = -1;
  int di = -1;
  chose_row = true;
  for (int i = 0; i < assign.n; i++) {
    if (assign.p[i] > NONE) continue;
    di++; // working on row di of U.

    for (int dj = 0; dj < U->n; dj++) {
      ui[dj] = lgb_U->me[di][dj];
      if (U->me[di][dj] + bound >= qpbp->inc - DELTA_BRANCH) {
	f_row_count[di]++;
      }
    }

    //nice_array_print(ui, n, ',');
    ins_sort_ascending(ui, NULL, n);
    //nice_array_print(ui, n, ',');
    force_to_size(ui, n, ufill, n_feat); // todo: is this right???
    //nice_array_print(ufill, n_feat, ',');
    //u_row_sum[di] = boost_(ufill); // a good score is SMALL, so negating
    u_row_sum[di] = row_variance(U, di); // a good score is SMALL, so negating
    //printf("Score = %f\n", u_row_sum[di]);
  }

  int best_row_f;
  double best_u_row_sum;
  chose_row = true;
  // The first criteria is number of fathoms as stored in f_row_count.
  // The tiebreaker criteria is in u_row_sum. Higher is better
  //nice_int_array_print(f_row_count, n);
  //nice_array_print(u_row_sum, n, ',');
  find_most_fathoms(f_row_count, u_row_sum, best_row_f, best_u_row_sum, ind_in_U, n);
  get_row_in_assign(assign, ind_in_U, ind_in_assign);
  int dj = -1;
  if (can_branch_on_cols) {
    for (int j = 0; j < assign.n; j++) {
      if (assign.pinv[j] > NONE) continue;
      dj++;

      for (di = 0; di < U->m; di++) {
	ui[di] = lgb_U->me[di][dj];
	if (U->me[di][dj] + bound >= qpbp->inc - DELTA_BRANCH) {
	  f_col_count[dj]++;
	}
      }
      force_to_size(ui, n, ufill, n_feat);
      ins_sort_ascending(ufill, NULL, n_feat);
      //u_col_sum[dj] = boost_(ufill);
      u_col_sum[dj] = col_variance(U, dj);
      //printf("Score = %f\n", u_col_sum[dj]);
    }
    int best_col_f, col_in_U;
    double best_u_col_sum;
    find_most_fathoms(f_col_count,u_col_sum,best_col_f,best_u_col_sum, col_in_U,n);
    if ((best_col_f > best_row_f) ||
	((best_u_col_sum > best_u_row_sum) && (best_col_f==best_row_f))) {
      chose_row = false;
      ind_in_U = col_in_U;
      get_col_in_assign(assign,ind_in_U,ind_in_assign);
    }
  }

  //printf("choice %d %d rc=%d\n", ind_in_U, ind_in_assign, chose_row);
  delete [] f_row_count;
  delete [] u_row_sum;
  delete [] f_col_count;
  delete [] u_col_sum;
}

/*
============================================================ 
 row_sum_with_sym(assign,J3,J3_dim,U,u_row_sum);

  compute row sums of U, storing result in previously 
  allocated u_row_sum
============================================================ 
 */
inline void row_sum_with_sym(const QAPAssignment &assign,
				 int *J3,int J3_dim,double opt,
				 MAT *U,double *u_row_sum)
{
  int i,di,j,dj,J3_count;

  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of U.

      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// count columns not in J3 as a fathom.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) {
	  u_row_sum[di] += opt - assign.bound;
	  continue;
	}
	J3_count++;

	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;
	
	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) continue;
	
	u_row_sum[di] += U->me[di][dj];
      }
    }
}

/*
============================================================ 
 row_col_sums(assign,J3,J3_dim,U,u_row_sum,u_col_sum);

  compute row and column sums of U, storing results in 
  previously allocated u_row_sum and u_col_sum.
============================================================ 
 */
void row_col_sums(const QAPAssignment &assign,int *J3,int J3_dim,
		  MAT *U,double *u_row_sum,double *u_col_sum)
{
  int i,di,j,dj,J3_count;

  di = -1;
  for (i=0;i<assign.n;i++) {
    if (assign.p[i] > NONE) continue;
    di++; // working on row di of U.
    
    dj = -1; J3_count = 0;
    for (j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
      if (assign.pinv[j] <= NONE) dj++;
      if (j != J3[J3_count]) continue; // skip over columns not in J3.
      J3_count++;
      
      // skip over columns that have already been assigned
      if(assign.pinv[j] > NONE) continue;
      
      // ignore disallowed assignments.
      if(assign.disallowed(i,j)) continue;
      
      u_row_sum[di] += U->me[di][dj];
      u_col_sum[dj] += U->me[di][dj];
    }
  }
}

/*
============================================================ 
  weighted_row_sum(assign,J3,J3_dim,U,u_row_sum);

  compute row sums of U, storing result in previously 
  allocated u_row_sum
============================================================ 
 */
inline void weighted_row_sum(const QAPAssignment &assign,
				 int *J3,int J3_dim,double opt,
				 MAT *U,double *u_row_sum)
{
  int i,di,j,dj,J3_count;
  double grf;

  //  cout << "Using weighted row sum" << std::endl;
  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of U.

      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// skip over columns not in J3.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) continue;
	J3_count++;

	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;

	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) continue;
	
	grf = std::min(1.0, (U->me[di][dj])/(opt - assign.bound) );
#ifdef INCORRECT_WEIGHT
	u_row_sum[di] += ((1 - branching_weight)*grf + branching_weight) * grf;
#else
	u_row_sum[di] += (1 - (1 - branching_weight)*grf) * grf;
#endif

      }
    }
}

/*
============================================================ 
 row_col_sums(assign,J3,J3_dim,U,u_row_sum,u_col_sum);

  compute row and column sums of U, storing results in 
  previously allocated u_row_sum and u_col_sum.
============================================================ 
 */
void weighted_row_col_sums(const QAPAssignment &assign,int *J3,int J3_dim,
			       double opt,MAT *U,
			       double *u_row_sum,double *u_col_sum)
{
  int i,di,j,dj,J3_count;
  double grf;
  //  cout << "Using weighted row-col sum" << std::endl;

  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of U.

      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// skip over columns not in J3.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) continue;
	J3_count++;

	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;
	
	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) continue;
	
	grf = std::min(1.0, (U->me[di][dj])/(opt - assign.bound) );

#ifdef INCORRECT_WEIGHT
	u_row_sum[di] += ((1 - branching_weight)*grf + branching_weight) * grf;
	u_col_sum[dj] += ((1 - branching_weight)*grf + branching_weight) * grf;
#else
	u_row_sum[di] += (1 - (1 - branching_weight)*grf) * grf;
	u_col_sum[dj] += (1 - (1 - branching_weight)*grf) * grf;
#endif
      }
    }

}

// rank[d] = position of line d when v is sorted descending (0 = largest).
// v is left as it is.
inline void rank_descending(const double *v, int *rank, int m)
{
  double *w = new double[m];
  int *order = new_perm<int>(m);
  for (int k = 0; k < m; k++) w[k] = v[k];
  ins_sort_descending(w, order, m);   // order[k] = line in position k
  for (int k = 0; k < m; k++) rank[order[k]] = k;
  delete [] w;
  delete [] order;
}

// !! Allocates memory
// u_row_sum[d], u_col_sum[d]: U sums of row/column d; row_order[d],
// col_order[d]: its rank (0 = largest sum). Callers keep the lines with
// rank < use_best. (These used to be the sorted sums and the inverse
// permutation, so callers tested the wrong lines and started each line's
// score from another line's sum.)
inline void row_col_rank_sum(const QAPAssignment &assign, QPBParameters *qpbp,
			 int *J3,int J3_dim,MAT *U,bool weighted,
			 double *&u_row_sum,double *&u_col_sum,
			 int *&row_order,int *&col_order,
			 int m)
{
  u_row_sum = new_zero<double>(m);
  u_col_sum = new_zero<double>(m);

  if(weighted)
    weighted_row_col_sums(assign,J3,J3_dim,qpbp->inc,U,u_row_sum,u_col_sum);
  else
    row_col_sums(assign,J3,J3_dim,U,u_row_sum,u_col_sum);

  row_order = new int[m];
  col_order = new int[m];

  rank_descending(u_row_sum,row_order,m);
  rank_descending(u_col_sum,col_order,m);
}

// as row_col_rank_sum, rows only
inline void row_rank_sum(const QAPAssignment &assign,QPBParameters *qpbp,
		     int *J3,int J3_dim,MAT *U,bool weighted,
		     double *&u_row_sum,int *&row_order,int m)
{
  u_row_sum = new_zero<double>(m);
  if (weighted) {
    weighted_row_sum(assign,J3,J3_dim,qpbp->inc,U,u_row_sum);
  }
  else {
    row_sum(assign,J3,J3_dim,U,u_row_sum);
  }

  row_order = new int[m];

  rank_descending(u_row_sum,row_order,m);
}

/*
============================================================ 

The row selection functions all take basically the same parameters:

void select_row(QAP *qap,const QAPAssignment &assign,
                int *J3, int J3_dim,MAT *X,MAT *U,
		double bound,QPBParameters *qpbp,
                int &row_in_U,int &row_in_assign)

n = the size of the original QAP.
m = the number of unassigned facilities/locations.  m<=n.

qap:     QAP of size m.

assign:  basically, a list of the assignments made to this point,
         see assign.h.  assign contains lists of which facilities
         have been assigned, and which locations have been
         assigned.  By looking at assign, we can determine
         how to map unassigned location i to a row of U.

J3:      We will only create children for facilities in J3.
J3_dim:  size of J3

X:       Solution of the QP used to obtain the bound. (unused)
U:       slack vars.

bound:   the computed QPB for this problem.
qpbp:    parameters for QP bound

row_in_assign:
         (output) we will create subproblems by fixing
         each available X(row_in_assign,j)=1.
row_in_U:
         the row in U containing information about row_in_assign.

Also, we have functions that choose either a row OR a 
column for branching.  These functions have an extra 
parameter, chose_row.  If chose_row==false, then 
row_in_assign and row_in_U are really col_in_assign and
col_in_U.

============================================================ 
*/

/*
============================================================ 
 select_row_max_sum()

 Choose the row i with largest \sum_{j=1}^n U_{ij}. The best
 value will be returned. The index will be stored in row_in_X.
============================================================ 
 */
double select_row_max_sum(QAP *qap,const QAPAssignment &assign,
			  int *J3, int J3_dim,MAT *U,
			  double bound,QPBParameters *qpbp,
			  bool weighted,int &row_in_U,int &row_in_assign)
{
  int m = U->n;
  double *u_sum,best_sum;

  u_sum = new_zero<double>(m);

  if (weighted) {
    weighted_row_sum(assign,J3,J3_dim,qpbp->inc,U,u_sum);
  }
  else {
    row_sum(assign,J3,J3_dim,U,u_sum);
  }

  find_max(u_sum, best_sum, row_in_U, m);
  get_row_in_assign(assign, row_in_U, row_in_assign);
  
  delete [] u_sum;
  return best_sum;
}


/*
============================================================ 
 select_ind_max_sum()

 Choose the row i with largest \sum_{j=1}^n U_{ij},
 or choose row j with largest \sum_{i=1}^n U_{ij},
 whichever is best.
============================================================ 
 */
void select_ind_max_sum(QAP *qap,const QAPAssignment &assign,
			int *J3, int J3_dim,MAT *X,MAT *U,
			double bound,QPBParameters *qpbp,
			bool weighted,int &ind_in_U,int &ind_in_assign,
			bool &chose_row)
{
  int row_in_U, col_in_U;
  int m = U->n;
  double *u_row_sum,best_row_sum;
  double *u_col_sum,best_col_sum;

  u_row_sum = new_zero<double>(m);
  u_col_sum = new_zero<double>(m);

  if(weighted)
    weighted_row_col_sums(assign,J3,J3_dim,qpbp->inc,U,u_row_sum,u_col_sum);
  else
    row_col_sums(assign,J3,J3_dim,U,u_row_sum,u_col_sum);

  find_max(u_row_sum,best_row_sum,row_in_U,m);
  find_max(u_col_sum,best_col_sum,col_in_U,m);

  if(best_row_sum >= best_col_sum) {
    chose_row = true;
    ind_in_U = row_in_U;
    get_row_in_assign(assign,row_in_U,ind_in_assign);
  }
  else {
    chose_row = false;
    ind_in_U = col_in_U;
    get_col_in_assign(assign,col_in_U,ind_in_assign);
  }

  delete [] u_row_sum;
  delete [] u_col_sum;
}

/*
============================================================ 
 Choose the row i with largest {j | U_{ij} > opt}.
 In case of tie, pick row with bigger row sum.
============================================================ 
 */
void select_row_max_fathoms(QAP *qap,const QAPAssignment &assign,
			    int *J3, int J3_dim,MAT *X,MAT *U,
			    double bound,QPBParameters *qpbp,     
			    int &row_in_U,int &row_in_assign)
{
  int i,j,di,dj,J3_count;
  int m = U->n;
  int *f_count,best_f;
  double *u_sum,best_u_sum;
  
  u_sum = new_zero<double>(m);
  f_count = new_zero<int>(m);

  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of U.

      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// skip over columns not in J3.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) continue;
	J3_count++;
	if(assign.pinv[j] > NONE) continue;

	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) {
	  f_count[di]++;
	  continue;
	}
	
	if(U->me[di][dj] + bound >= qpbp->inc - DELTA_BRANCH)
	  f_count[di]++;
	else
	  u_sum[di] += U->me[di][dj];
      }
    }

  find_most_fathoms(f_count,u_sum,best_f,best_u_sum,row_in_U,m);
  get_row_in_assign(assign,row_in_U,row_in_assign);

  delete [] f_count;
  delete [] u_sum;
}

/*
============================================================ 
 Choose the row i with largest {j | U_{ij} > opt}, or
 col j with largest {i | U_{ij} > opt}.
 In case of tie, pick row with bigger row/col sum.
============================================================ 
 */
/*
============================================================
select_col_eccentric(...)

Branching rule 6: branch on the free location (column) that is
farthest from the other free locations, i.e. the largest
max_k B[j][k] over free k != j. On grid instances (nug, sko, tho, ...)
these are the corners and edges: fixing a facility there gives the
children the most spread-out linear terms, so more of them fathom.
It costs nothing per node, unlike lookahead, which bounds every child.
Ties go to the column with more children fixed out by bound + U, then
to the larger U column sum, as in rule 2. Only for column branching
(no symmetry left); see select(). Returns false when every free
location is equally eccentric (e.g. esc16a's hypercube), where the
rule says nothing and the caller uses rule 2 instead.
============================================================
*/
// largest B entry from free location j to another free location
static double location_eccentricity(QAP *qap,const QAPAssignment &assign,int j)
{
  double ecc = -1.0e100;
  for(int k=0;k<assign.n;k++)
    if((k != j) && (assign.pinv[k] <= NONE) && (qap->B->me[j][k] > ecc))
      ecc = qap->B->me[j][k];
  return ecc;
}

bool select_col_eccentric(QAP *qap,const QAPAssignment &assign,MAT *U,
			  double bound,QPBParameters *qpbp,
			  int &ind_in_U,int &ind_in_assign)
{
  const int n = assign.n;
  double best_ecc = -1.0e100, best_sum = -1.0e100;
  double low_ecc = 1.0e100;
  int best_f = -1;
  int dj = -1;

  ind_in_U = -1;  ind_in_assign = -1;
  for(int j=0;j<n;j++) {
    if(assign.pinv[j] > NONE) continue;
    dj++;

    const double ecc = location_eccentricity(qap,assign,j);
    if (ecc < low_ecc) low_ecc = ecc;

    int f = 0;  double sum = 0.0;
    for(int i=0,di=-1;i<n;i++) {
      if(assign.p[i] > NONE) continue;
      di++;
      if(assign.disallowed(i,j) ||
	 (U->me[di][dj] + bound >= qpbp->inc - DELTA_BRANCH))
	f++;
      else
	sum += U->me[di][dj];
    }

    if((ecc > best_ecc) ||
       ((ecc == best_ecc) && ((f > best_f) || ((f == best_f) && (sum > best_sum))))) {
      best_ecc = ecc;  best_f = f;  best_sum = sum;
      ind_in_U = dj;  ind_in_assign = j;
    }
  }
  return best_ecc > low_ecc;
}

/*
============================================================
select_ind_dnn_xbar(...)

Branching rule 8: branch on the free row (facility) or column (location)
whose DNN primal marginals x_ij are most concentrated, i.e. the largest
sum_j x_ij^2 (row) or sum_i x_ij^2 (column). The DNN relaxation puts most
of a concentrated line's mass on a few assignments, so most of the line's
children are assignments the relaxation considers near-zero, and they
fathom sooner. Found on wil50 depth-3 misses and a nug30 depth-2 miss
(experiments/wil50_feasibility/tuning/lines): the line with the largest
sum x^2 was the best or second-best of about 10 candidates at every node,
eccentric (rule 6) 2 to 4.7x worse, lookahead (rule 4) 1 to 1.6x worse.
It costs nothing per node once the DNN bound has run. Columns only when
no symmetry is left (rows_only otherwise, as the other rules). Returns
false when no marginals were set for this node (no DNN bound ran); the
caller then uses lookahead.
============================================================
*/
static const double *dnn_xbar = NULL;
static int dnn_xbar_m = 0;

void branch_set_dnn_marginals(const double *x, int m)
{
  dnn_xbar = x;
  dnn_xbar_m = m;
}

// The free row and column (local indices) with the largest sum of squared
// DNN marginals, and their scores; false when no marginals are set for a
// node of this size.
static bool dnn_xbar_lines(const QAPAssignment &assign, int &row_in_U, double &row_score,
			   int &col_in_U, double &col_score)
{
  const int m = assign.n - assign.nfix;
  if (dnn_xbar == NULL || dnn_xbar_m != m) return false;
  const double *x = dnn_xbar;
  row_in_U = col_in_U = -1;  row_score = col_score = -1.0;
  for (int di = 0; di < m; di++) {
    double s = 0.0;
    for (int dj = 0; dj < m; dj++) s += x[di * m + dj] * x[di * m + dj];
    if (s > row_score) { row_score = s;  row_in_U = di; }
  }
  for (int dj = 0; dj < m; dj++) {
    double s = 0.0;
    for (int di = 0; di < m; di++) s += x[di * m + dj] * x[di * m + dj];
    if (s > col_score) { col_score = s;  col_in_U = dj; }
  }
  return row_in_U >= 0 && col_in_U >= 0;
}

static bool select_ind_dnn_xbar(const QAPAssignment &assign, bool rows_only,
				int &ind_in_U, int &ind_in_assign, bool &chose_row)
{
  int dr, dc;  double sr, sc;
  if (!dnn_xbar_lines(assign, dr, sr, dc, sc)) return false;
  double best;
  if (rows_only || sr >= sc) { chose_row = true;  ind_in_U = dr;  best = sr;  get_row_in_assign(assign, dr, ind_in_assign); }
  else                       { chose_row = false; ind_in_U = dc;  best = sc;  get_col_in_assign(assign, dc, ind_in_assign); }
  if (getenv("QAP_XBAR_DEBUG") && ind_in_U >= 0)
    fprintf(stderr, "xbar: depth %d %s %d sumsq %.4f\n", assign.nfix, chose_row ? "row" : "col", ind_in_assign, best);
  return ind_in_U >= 0;
}

void select_ind_max_fathoms(QAP *qap,const QAPAssignment &assign,
			    int *J3, int J3_dim,MAT *X,MAT *U,
			    double bound,QPBParameters *qpbp,
			    int &ind_in_U,int &ind_in_assign,
			    bool &chose_row)
{
  int i,j,di,dj,J3_count;
  int row_in_U,row_in_assign;
  int col_in_U,col_in_assign;
  int m = U->n;
  int *f_row_count,best_row_f;
  int *f_col_count,best_col_f;
  double *u_row_sum,best_u_row_sum;
  double *u_col_sum,best_u_col_sum;
  
  u_row_sum = new_zero<double>(m);
  u_col_sum = new_zero<double>(m);
  f_row_count = new_zero<int>(m);
  f_col_count = new_zero<int>(m);

  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of D.

      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	
	// skip over columns not in J3.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) continue;
	J3_count++;

	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;

	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) {
	  f_row_count[di]++;
	  f_col_count[dj]++;
	  continue;
	}
	
	if(U->me[di][dj] + bound >= qpbp->inc - DELTA_BRANCH) {
	  f_row_count[di]++;
	  f_col_count[dj]++;
	}
	else {
	  u_row_sum[di] += U->me[di][dj];
	  u_col_sum[dj] += U->me[di][dj];
	}
      }
    }

  find_most_fathoms(f_col_count,u_col_sum,best_col_f,best_u_col_sum,
		    col_in_U,m);
  find_most_fathoms(f_row_count,u_row_sum,best_row_f,best_u_row_sum,
		    row_in_U,m);

  if((best_row_f > best_col_f) ||
     ((best_u_row_sum > best_u_col_sum)&&(best_row_f==best_col_f))) 
     {
       chose_row = true;
       ind_in_U = row_in_U;
       get_row_in_assign(assign,ind_in_U,ind_in_assign);
     }
   else 
     {
       chose_row = false;
       ind_in_U = col_in_U;
       get_col_in_assign(assign,ind_in_U,ind_in_assign);
     }

  delete [] f_row_count;
  delete [] u_row_sum;
  delete [] f_col_count;
  delete [] u_col_sum;
}

/*
============================================================ 
 Choose the row i with largest \sum_{j=1}^n B_{ij}, where
 B_ij is the QP bound obtained by making the assignment i-->j.
============================================================ 
 */
void select_ind_max_sum_using_B(QAP *qap,QAPAssignment &assign,
				int *J3, int J3_dim,MAT *X,MAT *U,       
				double bound,QPBParameters *qpbp,int use_best,
				bool weighted,int which_bound,
				int &ind_in_U,int &ind_in_assign,
				bool &chose_row)
{
  int i,j,di,dj,J3_count;
  int m = U->n, n = assign.n;
  double *u_row_sum,best_row_sum;
  double *u_col_sum,best_col_sum;
  int *row_order, *col_order;
  
  // qpb only
  MAT *Wa=MNULL, *Wb=MNULL;
  VEC *wa=VNULL, *wb=VNULL;

  QAPAssignment child;
  MAT *X1 = MNULL, *U1 = MNULL;
  QAP *qr = new QAP(n);

  X1 = m_resize(X1,m-1,m-1);
  U1 = m_resize(U1,m-1,m-1);
  //printf("in B: m = %d, n = %d\n", m, n);
  row_col_rank_sum(assign,qpbp,J3,J3_dim,U,weighted,
		   u_row_sum,u_col_sum,row_order,col_order,m);

  if(is_qpbnd(which_bound)) {
    wa = v_get(m-1);      wb = v_get(m-1);
    Wa = m_get(m-1,m-1);  Wb = m_get(m-1,m-1);
  }

  // Do a few iterations of qpbnd() for the best use_best rows
  // of the current QAP, as determined by U.

  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of U.

      // get qr->A  
      assign.p[i] = 0;
      reduceA(qap->A,qr->A,assign.p,m-1,assign.n);
      assign.p[i] = NONE;

      if(is_qpbnd(which_bound)) 
	init_w0(qr->A,wa,Wa,false);
 
      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// skip over columns not in J3.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) continue;
	J3_count++;

	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;
	if(assign.disallowed(i,j)) continue;
	if((row_order[di] >= use_best)&&(col_order[dj] >= use_best)) continue;

	child = assign;  	child.fix(i,j);

	reduceB(qap->B,qr->B,child.pinv,m-1,assign.n);
	reduceC(qap,qr,child);

	init_X_from_parent(X1,X,di,dj);  // start from the parent's point

	qpbp->shift = qr->shift;
	qpbp->feasible = true;
	qpbp->lap_scale = 1.0;
	qpbp->improve = false;

	if(is_qpbnd(which_bound)) {
	  init_w0(qr->B,wb,Wb,true);
#ifdef BRANCH_USES_QPBND
	  U->me[di][dj] = qpbnd2(qr->A,qr->B,qr->C,X1,U1,wa,wb,Wa,Wb,
				 qpbp) - bound;
#else
	  U->me[di][dj] = qpbnd_parametric2(qr->A,qr->B,qr->C,X1,U1,
					    wa,wb,Wa,Wb,qpbp) - bound;	
#endif
	}
	else {
	  U->me[di][dj] = glbnd(qr->A,qr->B,qr->C,U1,qpbp->shift) - bound;
	}
      }
    }

  //  select_index_lgb(assign, true, J3, J3_dim, U, bound, qpbp, ind_in_assign, ind_in_U, chose_row);
  select_ind_max_sum(qap,assign,J3,J3_dim,X,U,bound,qpbp,weighted, ind_in_U,ind_in_assign,chose_row);

  delete qr;
  delete [] u_row_sum;  delete [] u_col_sum;
  delete [] row_order;  delete [] col_order;

  M_FREE(X1);  M_FREE(U1);

  if(is_qpbnd(which_bound)) {
    M_FREE(Wa); M_FREE(Wb);
    V_FREE(wa); V_FREE(wb);
  }
}


// returns the max row sum. The row that attains it is in row_in_X.
void select_row_max_sum_using_B(QAP *qap, QAPAssignment &assign,
				int *J3, int J3_dim, MAT *U,
				double bound, QPBParameters *qpbp,int use_best,
				bool weighted, int which_bound,
				int &row_in_U, int &row_in_assign)
{
  int i,j,di,dj,J3_count;
  int m = U->n, n = assign.n;
  double *u_sum,best_sum;
  int *order;

  row_rank_sum(assign, qpbp, J3, J3_dim, U, weighted, u_sum, order, m);
  
  QAPAssignment child;
  MAT *X1 = MNULL, *U1 = MNULL;
  MAT *Wa, *Wb;
  VEC *wa, *wb;
  QAP *qr = new QAP(m-1);

  X1 = m_resize(X1, m-1, m-1);
  U1 = m_resize(U1, m-1, m-1);

  if (is_qpbnd(which_bound)) {
    wa = v_get(m-1);      wb = v_get(m-1);
    Wa = m_get(m-1,m-1);  Wb = m_get(m-1,m-1);
  }

  // Do a few iterations of qpbnd() for the best use_best rows
  // of the current QAP, as determined by U.

  di = -1;
  for (i=0; i<assign.n; i++) {
    if (assign.p[i] > NONE) continue;
    di++; // working on row di of D.
    
    if(order[di] >= use_best) continue;
    
    // get qr->A
    assign.p[i] = 0; assign.nfix++; 
    reduceA(qap->A, qr->A, assign.p, m-1, assign.n);
    assign.p[i] = NONE; assign.nfix--;

    if(is_qpbnd(which_bound)) 
      init_w0(qr->A,wa,Wa,0);
    
    dj = -1; J3_count = 0;
    for(j=0; (j<assign.n) && (J3_count < J3_dim); j++) {
      // skip over columns not in J3.
      if(assign.pinv[j] <= NONE) dj++;
      if(j != J3[J3_count]) continue;
      J3_count++;
      
      // skip over columns that have already been assigned
      if (assign.pinv[j] > NONE) continue;
      
      // ignore disallowed assignments.
      if (assign.disallowed(i,j)) continue;

      child = assign; // this copies
      child.fix(i,j);	  

      // complete QAP reduction
      reduceB(qap->B, qr->B, child.pinv, m-1, assign.n);
      reduceC(qap, qr, child);

      init_X(X1, qr->A->m);
      
      qpbp->shift = qr->shift;
      qpbp->feasible = true;
      qpbp->lap_scale = 1.0;
      qpbp->improve = false;

      if (is_qpbnd(which_bound)) {
	init_w0(qr->B, wb, Wb, 1);
#ifdef BRANCH_USES_QPBND
	U->me[di][dj] = qpbnd2(qr->A, qr->B, qr->C, X1, U1,
			       wa, wb, Wa, Wb, qpbp) - bound;
#else
	U->me[di][dj] = qpbnd_parametric2(qr->A,qr->B,qr->C,X1,U1,
					  wa,wb,Wa,Wb,qpbp) - bound;
#endif
      }
      else {
	U->me[di][dj] = glbnd(qr->A,qr->B,qr->C,U1,qpbp->shift) - bound;
      }
      //printf("bound was %f\n", U->me[di][dj]);
    }
  }
  
  select_row_max_sum(qap,assign,J3,J3_dim,U,bound,qpbp,
		     weighted, row_in_U, row_in_assign);
  
  // Free memory
  delete qr;
  delete [] u_sum;  delete [] order;
  
  M_FREE(X1);  M_FREE(U1);
  if (is_qpbnd(which_bound)) {
    M_FREE(Wa);  M_FREE(Wb);
    V_FREE(wa);  V_FREE(wb);
  }
}

/*
============================================================ 
  "lookahead strategy"
  
  Compute row sum of U:
  for each row i, for each j \in J3,
     add contribution to u_sum
     
  sort row sum

  for each row i, where rank(i) <= k, for each j \in J3,
     compute B_ij, also yielding U1, X1.
     compute J3' set using current QAPAssignment
     compute row sum of U1 over J3'
     keep best result.
============================================================ 
 */
void select_row_lookahead(QAP *qap,QAPAssignment &assign,
			  int *J3, int J3_dim,MAT *X,MAT *U,
			  double bound,QPBParameters *qpbp,int use_best, 
			  bool weighted, int which_bound,
			  int &row_in_U,int &row_in_assign)
{
  int n = assign.n;
  int i,j,di,dj,k,l,bi;
  int m = U->n;
  double *u_sum,best_sum;
  int *order;
  int J3_count;
  double temp;

  MAT *Wa = MNULL;  VEC *wa = VNULL;
  MAT **Wbs = NULL;  VEC **wbs = NULL;  MAT **Bs = NULL;
  int wb_sz;

  if(is_qpbnd(which_bound)) {
    wa = v_get(m-1);      
    Wa = m_get(m-1,m-1);  

    reduced_B(qap,assign,J3,J3_dim,Bs,Wbs,wbs,wb_sz,m);
  }

  // rank row sums
  row_rank_sum(assign,qpbp,J3,J3_dim,U,weighted,u_sum,order,m);

  MAT *X1 = MNULL, *U1 = MNULL;
  int *J3_child, J3_child_dim,J3_child_count;
  double *u1_sum, best_u1_sum;
  QAP *qr = new QAP(m-1);

  u1_sum = new double[m-1];

  // we may set X1 and U1 to a fixed size because we never explicitly
  // disallow assignments.
  X1 = m_resize(X1, m-1, m-1);
  U1 = m_resize(U1, m-1, m-1);

  di = -1;
  for (i=0;i<assign.n;i++) {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di 

      /*
      for (dj = 0; dj < m; dj++) {
	lgb_U->me[di][dj] = 0;
	}*/

      if (order[di] >= use_best) continue;
      
      // get qr->A  
      assign.p[i] = 0;
      reduceA(qap->A,qr->A,assign.p,m-1,assign.n);
      assign.p[i] = NONE;

      if(is_qpbnd(which_bound)) 
	init_w0(qr->A,wa,Wa,false);

      dj = -1; J3_count = 0; bi = -1;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// skip over columns not in J3.
	if(assign.pinv[j] <= NONE) dj++;
	if(j != J3[J3_count]) continue;
	J3_count++;

	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;

	bi++;  // next entry of cached B's
	
	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) continue;
	
	assign.fix(i,j);

	// (already have qr->B)
	// get qr->C
	reduceC(qap,qr,assign);

	init_X_from_parent(X1,X,di,dj);  // start from the parent's point

	qpbp->shift = qr->shift;
	qpbp->feasible = true;
	qpbp->lap_scale = 1.0;
	qpbp->improve = false;

	if(is_qpbnd(which_bound)) {
#ifdef BRANCH_USES_QPBND
	  temp = qpbnd2(qr->A,Bs[bi],qr->C,X1,U1,wa,wbs[bi],Wa,Wbs[bi],
			qpbp) - bound;
#else
	  temp = qpbnd_parametric2(qr->A,Bs[bi],qr->C,X1,U1,
				   wa,wbs[bi],Wa,Wbs[bi],qpbp) - bound;
#endif
	}
	else {
	  temp = glbnd(qr->A,qr->B,qr->C,U1,qpbp->shift) - bound;
	}

	// find the J3 set for the child.
	determine_J3(assign, J3_child, J3_child_dim);

	U->me[di][dj] = std::max(U->me[di][dj], temp);
	// lgb_U->me[di][dj] = U->me[di][dj];
 
	// find the best row sum in U1.
	for(k=0; k<U1->m; k++) {
	  u1_sum[k] = 0.0;
	}

	if(0)
	  weighted_row_sum(assign,J3_child,J3_child_dim,qpbp->inc,U1,u1_sum);
	else
	  row_sum(assign,J3_child,J3_child_dim,U1,u1_sum);

	// find the largest row sum, and add it to total for row i.
	find_max(u1_sum, best_u1_sum, l, U1->m);

	// the (di, dj) entry. So an alternative would be to:
	//  a) log all these in a matrix
	//  b) compute the row variance of the matrix.
	double u_di = J3_child_dim * (U->me[di][dj] + bound) + best_u1_sum;
	u_sum[di] += u_di;
	// lgb_U->me[di][dj] += u_di;
	assign.unfix(i, j);
      }
  }

  //for (di = 0; di < m; di++) {
  //  u_sum[di] = row_variance(lgb_U, di);
  //}
  //nice_array_print(u_sum, m, ',');
  
  // todo: if I want to log, log right here:
  // log everything in u_sum.
  // above, before I modify U, log U.
  
  find_max(u_sum, best_sum, row_in_U, m);
  get_row_in_assign(assign,row_in_U,row_in_assign);

  // print bestk
  //printf("%2d %2d %2d %2d 0\n",m,use_best,row_in_U,order[row_in_U]);

  delete [] u_sum;
  delete [] u1_sum;
  delete [] order;

  delete qr;
  M_FREE(X1);  M_FREE(U1);

  if(is_qpbnd(which_bound)) {
    M_FREE(Wa);    V_FREE(wa);
    for(k=0;k<wb_sz;k++) {
      M_FREE(Bs[k]);
      M_FREE(Wbs[k]);
      V_FREE(wbs[k]);
    }
    delete [] Bs;
    delete [] Wbs;
    delete [] wbs;
  }
}

/*
 Lookahead (rules 4 and 7): bound every child of the candidate lines and
 branch on the line whose children score best. Rule 4's candidates are
 the use_best rows and the use_best columns with the largest U sums.
 Rule 7 (eccentric_cols) takes the use_best most eccentric free columns
 (ties by U sum, see select_col_eccentric) and the use_best-1 rows with
 the largest U sums: as good a choice for grid-like instances at a
 fraction of the cost, and the rows keep it working where locations
 don't differ (esc16a). See experiments/restricted_lookahead.
*/
void select_ind_lookahead(QAP *qap,QAPAssignment &assign,
			  int *J3, int J3_dim,MAT *X,MAT *U,
			  double bound,QPBParameters *qpbp,int use_best,
			  bool weighted, int which_bound,
			  int &ind_in_U,int &ind_in_assign,
			  bool &chose_row, bool eccentric_cols = false,
			  bool xbar_lines = false)
{
  int i,j,di,dj,k,l,bi;
  int m = U->n;
  double *u_row_sum,best_row_sum;
  double *u_col_sum,best_col_sum;
  int row_in_U, col_in_U, row_in_assign, col_in_assign;
  int *row_count, *col_count, J3_count;
  int *row_order, *col_order;
  double temp;
  
  // used by qpb only:
  MAT *Wa = MNULL;  VEC *wa = VNULL;
  MAT **Wbs = NULL;  VEC **wbs = NULL;  MAT **Bs = NULL;
  int wb_sz;

  if (is_qpbnd(which_bound)) {
    wa = v_get(m-1);      
    Wa = m_get(m-1,m-1);  
    reduced_B(qap,assign,J3,J3_dim,Bs,Wbs,wbs,wb_sz,m);
  }

  row_col_rank_sum(assign,qpbp,J3,J3_dim,U,weighted,
		   u_row_sum,u_col_sum,row_order,col_order,m);

  int kr = use_best, kc = use_best;   // candidate rows, columns
  if (eccentric_cols) {
    std::vector<double> ecc(m);
    std::vector<int> by(m);
    for (int jj = 0, d = -1; jj < assign.n; jj++) {
      if (assign.pinv[jj] > NONE) continue;
      d++;
      ecc[d] = location_eccentricity(qap, assign, jj);
    }
    for (int d = 0; d < m; d++) by[d] = d;
    std::sort(by.begin(), by.end(), [&](int a, int b) {
	if (ecc[a] != ecc[b]) return ecc[a] > ecc[b];
	if (u_col_sum[a] != u_col_sum[b]) return u_col_sum[a] > u_col_sum[b];
	return a < b; });
    for (int k = 0; k < m; k++) col_order[by[k]] = k;
    kr = std::max(use_best - 1, 0);
  }
  bool xbar_used = false;
  if (xbar_lines) {
    // rule 9: the only candidates are the marginal row and the marginal
    // column (rule 8's two lines); the children's bounds choose between them
    int dr, dc;  double sr, sc;
    if (dnn_xbar_lines(assign, dr, sr, dc, sc)) {
      for (int d = 0; d < m; d++) { row_order[d] = m;  col_order[d] = m; }
      row_order[dr] = 0;  col_order[dc] = 0;
      kr = kc = 1;
      xbar_used = true;
    }
  }

  int n = assign.n;
  MAT *X1 = MNULL, *U1 = MNULL;
  int *J3_child, J3_child_dim,J3_child_count;
  double *u1_row_sum, *u1_col_sum, best_u1_sum;
  QAP *qr = new QAP(m-1);

  u1_row_sum = new_zero<double>(m-1);
  u1_col_sum = new_zero<double>(m-1);

  // we may set X1 and U1 to a fixed size because we never explicitly
  // disallow assignments.
  X1 = m_resize(X1,m-1,m-1);
  U1 = m_resize(U1,m-1,m-1);

  di = -1;
  for (i=0; i<assign.n; i++) {
    if(assign.p[i] > NONE) continue;
    di++; // working on row di of D.
    
    // get qr->A  
    assign.p[i] = 0;
    reduceA(qap->A,qr->A,assign.p,m-1,assign.n);
    assign.p[i] = NONE;
    
    if (is_qpbnd(which_bound)) {
      init_w0(qr->A,wa,Wa,false);
    }
    
    dj = -1; J3_count = 0; bi = -1;
    for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
      
      // skip over columns not in J3.
      if(assign.pinv[j] <= NONE) dj++;
      if(j != J3[J3_count]) continue;
      J3_count++;
      
      // skip over columns that have already been assigned
      if(assign.pinv[j] > NONE) continue;
      
      bi++;  // next entry of cached B's
      
      // ignore disallowed assignments.
      if (assign.disallowed(i,j)) continue;

      // skip over entries that are are not 'good enough' as determined
      // by row_ / col_order.
      if((row_order[di] >= kr) && (col_order[dj] >= kc)) continue;
      
      assign.fix(i,j);
      
      // find the J3 set for the child.
      determine_J3(assign, J3_child, J3_child_dim);
      
      reduceC(qap,qr,assign);
      init_X_from_parent(X1,X,di,dj);  // start from the parent's point
      
      qpbp->shift = qr->shift;
      qpbp->feasible = true;
      qpbp->lap_scale = 1.0;
      qpbp->improve = false;

      // This is where we determine the 'score' for the cell.
      if (is_qpbnd(which_bound)) {
#ifdef BRANCH_USES_QPBND
	temp = qpbnd2(qr->A,Bs[bi],qr->C,X1,U1,wa,wbs[bi],Wa,Wbs[bi],
		      qpbp) - bound;
#else
	temp = qpbnd_parametric2(qr->A,Bs[bi],qr->C,X1,U1,wa,wbs[bi],Wa,
				 Wbs[bi],qpbp) - bound;
#endif
      }
      else {
	temp = glbnd(qr->A,qr->B,qr->C,U1,qpbp->shift) - bound;
      }
      
      U->me[di][dj] = std::max(U->me[di][dj],temp);
      
      // find the best row sum in U1.
      for(k=0; k<U1->m; k++) {
	u1_row_sum[k] = 0.0;
	u1_col_sum[k] = 0.0;
      }
      
      if(0)
	weighted_row_col_sums(assign,J3_child,J3_child_dim,qpbp->inc,
			      U1,u1_row_sum,u1_col_sum);
      else
	row_col_sums(assign,J3_child,J3_child_dim,U1,
		     u1_row_sum,u1_col_sum);
      
      // find the largest row sum, and add it to total for row i.
      
      find_max(u1_row_sum,best_row_sum,row_in_U,U1->m);
      find_max(u1_col_sum,best_col_sum,col_in_U,U1->m);
      best_u1_sum = std::max(best_row_sum, best_col_sum);
      
      u_row_sum[di] += J3_child_dim*(U->me[di][dj]+bound) + best_u1_sum;
      u_col_sum[dj] += J3_child_dim*(U->me[di][dj]+bound) + best_u1_sum;
      
      assign.unfix(i,j);
    }
  }
  
  // todo: if I want to log, log right here:
  // log everything in u_sum. (row and column)
  // above, before I modify U, log U.

  
  // This is where we are actually selecting the best rows and cols.
  find_max(u_row_sum, best_row_sum, row_in_U, m);
  find_max(u_col_sum, best_col_sum, col_in_U, m);

  if(best_row_sum >= best_col_sum) {
    chose_row = true;
    ind_in_U = row_in_U;
    get_row_in_assign(assign,row_in_U,ind_in_assign);
    //printf("%2d %2d %2d %2d 1\n",m,use_best,ind_in_U,row_order[ind_in_U]);
  }
  else {
    chose_row = false;
    ind_in_U = col_in_U;
    get_col_in_assign(assign,col_in_U,ind_in_assign);
    //printf("%2d %2d %2d %2d 1\n",m,use_best,ind_in_U,col_order[ind_in_U]);
  }
  if (xbar_used && getenv("QAP_XBAR_DEBUG"))
    fprintf(stderr, "xbar9: depth %d %s %d (row %.1f col %.1f)\n", assign.nfix,
	    chose_row ? "row" : "col", ind_in_assign, best_row_sum, best_col_sum);

  delete [] u_row_sum;  delete [] u_col_sum;
  delete [] u1_row_sum; delete [] u1_col_sum;
  delete [] row_order;  delete [] col_order;

  delete qr;
  M_FREE(X1);  M_FREE(U1);
  if(is_qpbnd(which_bound)) {
    M_FREE(Wa);    V_FREE(wa);
    for(k=0;k<wb_sz;k++) {
      M_FREE(Bs[k]);
      M_FREE(Wbs[k]);
      V_FREE(wbs[k]);
    }
    delete [] Bs;
    delete [] Wbs;
    delete [] wbs;
  }
}

/*
============================================================ 
  determine if there is symmetry to be exploited at this node.
  if we are at the root, then we need only consider locations (columns)
  in J1.
  
  if all assigned locations are in the set J2, then we need only
  consider locations in J3.

  (when I say 'consider' I mean we only need look at those columns
   to determine which row to branch on, and we need only create
   children corresponding to those columns)
============================================================ 
 */
void determine_J3(const QAPAssignment &assign, int *&myJ3, int &myJ3_dim)
{
  // check J1 -- if in there, set J3 and get out
  // check J2 -- if all in there, set J3 and get out
  // otherwise, set J3 = 1..n

  int i=0,n=assign.n,k;

  bool use_J3 = false;

  if(assign.depth == 0) // at root, set J3 = J1 and quit.
    {
      if(J1_dim) {

	myJ3 = J1;
	myJ3_dim = J1_dim;
      }
      else {
	myJ3 = J3_default;
	myJ3_dim = J3_default_dim;
      }
      return;
    }

  for(k=0; k<J_dim; k++) {

    //    if(assign.nfix < J2_dim[k])
    if(assign.nfix <= J2_dim[k])
      {
	use_J3 = true;
	for(i=0;i<assign.n;i++)
	  {
	    // if some vars in this row are fixed to 0, then there is
            // no symmetry to exploit.
	    if(
	       (assign.p[i]==PARTIAL) || 
	       ( (assign.p[i]>NONE) && !find(J2[k],J2_dim[k],assign.p[i]) )
	       )
	      {
		use_J3 = false;
		i = assign.n;
		break;
	      }
	  }
	if(use_J3) {
	  break;
	}
      }
  }

  if(use_J3) {
    myJ3_dim = J3_dim[k];
    myJ3 = J3[k];

  }
  else {
    myJ3_dim = J3_default_dim;
    myJ3 = J3_default;
  }

}

#ifdef STORE_J_MATRIX
void scan_and_disallow(QAPAssignment &assign,int *J3,int J3_dim,
		       MAT *U,double bound,double opt)
{
  int i,di,j,dj,J3_count;

  di = -1;
  for(i=0;i<assign.n;i++)
    {
      if(assign.p[i] > NONE) continue;
      di++; // working on row di of U.

      dj = -1; J3_count = 0;
      for(j=0; (j<assign.n)&&(J3_count < J3_dim); j++) {
	
	// skip over columns not in J3.
	if(j != J3[J3_count]) continue;
	J3_count++;
	
	// skip over columns that have already been assigned
	if(assign.pinv[j] > NONE) continue;
	dj++;
	
	// ignore disallowed assignments.
	if(assign.disallowed(i,j)) continue;
	
	if(U->me[di][dj] + bound > opt - DELTA_BRANCH) {
	  assign.disallow(i,j);
	}
      }
    }
}
#endif

inline void reorder(QAPAssignment *child_tmp,QAPAssignment *child,
		    double d_vals[],int kids)
{
 // Now we have determined which children to create and have them
  // stored in kids[].  Sort d_vals[] to determine how to reorder kids[].

  int *d_order;
  int j,jj;
  double save; int save_idx;
  /*
  d_order = new int[kids];
  for(jj=0;jj<kids;jj++)
    d_order[jj] = jj;
  */
  d_order = new_perm<int>(kids);
  
  for(jj = kids-1; jj>0; jj--)
    {
      j = jj;
      save = d_vals[jj-1];
      save_idx = d_order[jj-1];
      d_vals[kids] = save; 
      while(save < d_vals[j])
	{
	  d_vals[j-1]  = d_vals[j];
	  d_order[j-1] = d_order[j];
	  j = j + 1;
	}
      d_vals[j-1] = save;
      d_order[j-1] = save_idx;
    }
  
  // now d_order holds the correct ordering ... children with
  // large D_ij should go to the front of kids[] since they're easier.
  // now kids[] needs to be reordered to match the order specified
  // by d_order.
  
  // for (jj == 0 ; jj < kids ; jj++){
  //   cout << jj << " : "  << d_vals[jj] << " ; " ;
  // }
  
  for (jj = 0 ; jj < kids ; jj++){
    
    child[jj] = child_tmp[d_order[kids-1-jj]];
      
  }
  delete [] d_order;
}

// Select a row or column to branch on. chose_row will indicate which.
inline void select(QAP *qap, QAPAssignment &assign, MAT *X, MAT *U,
		   double bound, QPBParameters *qpbp,
		   BranchParameters *bp,
		   bool &chose_row,
		   int &row_in_U, int &row_in_assign,
		   int *&J3, int &J3_dim, int &J3_count)
{
  int buff_depth = 2;
  chose_row = true;
  int z = 0;

  // only need to look at columns in J3.
  determine_J3(assign,J3,J3_dim);

#ifdef STORE_J_MATRIX
  if (use_J) {
    scan_and_disallow(assign,J3,J3_dim,U,bound,opt);
  }
#endif

  bool can_branch_on_cols = (J3_dim == J3_default_dim);
  /*
    if there is no symmetry left, then we may do either row- or
    column-branching.  Otherwise, we can only do row-branching.
  */
  switch (bp->how) {
  case BRANCH_FIXED:
    chose_row = bp->choose_row;
    row_in_U = bp->fixed_index_in_U;
    row_in_assign = bp->fixed_index_in_assign;
    //printf("fixed: %d %d %d\n", chose_row, row_in_U, row_in_assign);
    break;

  case BRANCH_U_SUM:
    if (can_branch_on_cols) {
      select_ind_max_sum(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			   use_weighted_branching,row_in_U,row_in_assign,chose_row);
    }
    else
      select_row_max_sum(qap,assign,J3,J3_dim,U,bound,qpbp,
			 use_weighted_branching,row_in_U,row_in_assign);
    break;

  case BRANCH_U_FATHOM:
    if(can_branch_on_cols) {
      select_ind_max_fathoms(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			     row_in_U,row_in_assign,chose_row);
    }
    else
      select_row_max_fathoms(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			     row_in_U,row_in_assign);
    break;

  case BRANCH_STRONG: 
    if(J3_dim == J3_default_dim) {
      select_ind_max_sum_using_B(qap,assign,J3,J3_dim,X,U,bound,qpbp,
				 bp->use_best,
				 use_weighted_branching, bp->which_bound,
				 row_in_U,row_in_assign,chose_row);
    }
    else
      select_row_max_sum_using_B(qap,assign,J3,J3_dim,U,bound,qpbp,
				 bp->use_best,
				 use_weighted_branching, bp->which_bound,
				 row_in_U,row_in_assign);
    break;

  case BRANCH_DNN_XBAR:
    if (select_ind_dnn_xbar(assign, !can_branch_on_cols,
			    row_in_U, row_in_assign, chose_row))
      break;
    // no DNN marginals at this node: lookahead, as rule 4
  case BRANCH_LOOKAHEAD:
  case BRANCH_XBAR_LOOKAHEAD:
    if(can_branch_on_cols) {
      select_ind_lookahead(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			   bp->use_best,
			   use_weighted_branching, bp->which_bound,
			   row_in_U,row_in_assign,chose_row,false,
			   bp->how == BRANCH_XBAR_LOOKAHEAD);
    }
    else
      select_row_lookahead(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			   bp->use_best,
			   use_weighted_branching, bp->which_bound,
			   row_in_U,row_in_assign);
    break;

  case BRANCH_LOOKAHEAD_ECC:
    if(can_branch_on_cols) {
      select_ind_lookahead(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			   bp->use_best,
			   use_weighted_branching, bp->which_bound,
			   row_in_U,row_in_assign,chose_row,true);
    }
    else
      select_row_lookahead(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			   bp->use_best,
			   use_weighted_branching, bp->which_bound,
			   row_in_U,row_in_assign);
    break;

  case BRANCH_ECCENTRIC:
    if(can_branch_on_cols) {
      chose_row = false;
      if(!select_col_eccentric(qap,assign,U,bound,qpbp,row_in_U,row_in_assign))
	select_ind_max_fathoms(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			       row_in_U,row_in_assign,chose_row);
    }
    else  // symmetry left: rows only, as for rule 2
      select_row_max_fathoms(qap,assign,J3,J3_dim,X,U,bound,qpbp,
			     row_in_U,row_in_assign);
    break;

  case BRANCH_BUFF:
    select_index_lgb(assign, can_branch_on_cols, J3, J3_dim,
		     U, bound, qpbp, row_in_assign, row_in_U, chose_row);

    break;
    
  default:
    std::cout << "branch::Unrecognized branching rule, aborting." << std::endl;
    exit(-1);
    break;
  }
}

void branch(QAP *qap, QAPAssignment &assign, MAT *X, MAT *U,
	    QAPAssignment *child, int &kids,
	    double bound, QPBParameters *qpbp, BranchParameters *bp,
	    int &loc_in_row, int &loc_kept)
{
  int m = assign.n - assign.nfix; // or U->m
  bool chose_row = true;
  int col_in_assign,col_in_U,row_in_U,row_in_assign;
  int *J3, J3_dim, J3_count;

  //  will get row_in_assign, col_in_assign as selection
  select(qap, assign, X, U, bound, qpbp, bp,
	 chose_row, row_in_U, row_in_assign, J3, J3_dim, J3_count);

  if (!chose_row) {
    col_in_assign = row_in_assign;
    col_in_U = row_in_U;
  }

  kids = 0;
  QAPAssignment child_tmp[m]; 
  //  QAPAssignment child_tmp[U->m]; 
  double *d_vals;
  if (bp->reorder_children)
    d_vals = new double[m+1];  // there can be at most D->m children.
                                  // extra slot is for insertion sort

  // we have now decided to fix the assignments in row row_in_assign to 1.
  // for each permitted assignment in this row, create a new child and
  // make the assignment.
  
  if (chose_row) {
    /* 
       for each potential subproblem in the chosen row, create it
       if we are not sure it can be fathomed.
     */
    col_in_U = -1;     J3_count = 0;
    loc_in_row = 0;    loc_kept = 0;
    for (col_in_assign=0;
	(col_in_assign < assign.n) && (J3_count < J3_dim);
	col_in_assign++)
      {
	// skip over columns not in J3.
	if(assign.pinv[col_in_assign] <= NONE) col_in_U++;
	if(col_in_assign != J3[J3_count]) continue;
	J3_count++;

	// skip over columns that have already been assigned
	if (assign.pinv[col_in_assign] > NONE) continue;
	
	loc_in_row++;
	// if(ceil(U->me[row_in_U][col_in_U] + bound) < qpbp->inc) 
	if (U->me[row_in_U][col_in_U] + bound < qpbp->inc - DELTA_BRANCH) {
	  child_tmp[kids] = assign;
	  //child_tmp[kids] = new QAPAssignment(assign);
	  child_tmp[kids].fix(row_in_assign,col_in_assign);
	  child_tmp[kids].predicted_bound = bound + 
	    U->me[row_in_U][col_in_U];
	  if (bp->reorder_children) {
	    d_vals[kids] = U->me[row_in_U][col_in_U];
	  }
	  kids++; loc_kept++;
	} /*
	else
	  printf("%f + %f >= %f,  (%d %d)  (%d %d)\n",
		 U->me[row_in_U][col_in_U] , 
		 bound , qpbp->inc - DELTA_BRANCH,
		 row_in_assign, col_in_assign,
		 row_in_U, col_in_U);  */
      }
  }
  else {
    row_in_U = -1;
    loc_in_row = 0;
    loc_kept = 0;
    for(row_in_assign=0;row_in_assign < assign.n;row_in_assign++) {
      if(assign.p[row_in_assign] > NONE) continue;
	  
      row_in_U++; 
      
      loc_in_row++;
      // if(ceil(U->me[row_in_U][col_in_U] + bound) < qpbp->inc) 
      if(U->me[row_in_U][col_in_U] + bound < qpbp->inc - DELTA_BRANCH) 
	{
	  child_tmp[kids] = assign;
	  //child_tmp[kids] = new QAPAssignment(assign);
	  child_tmp[kids].fix(row_in_assign,col_in_assign);
	  child_tmp[kids].predicted_bound = bound + 
	    U->me[row_in_U][col_in_U];
	  if (bp->reorder_children) {
	    d_vals[kids] = U->me[row_in_U][col_in_U];
	  }
	  kids++; loc_kept++;
	}/*
	   else
	   printf("%f + %f >= %f,  (%d %d)  (%d %d)\n",
	   U->me[row_in_U][col_in_U] , 
	   bound , qpbp->inc - DELTA_BRANCH,
	   row_in_assign, col_in_assign,
	   row_in_U, col_in_U); */
    }
  }

  if (!bp->reorder_children) {
    for (int i = 0 ; i < kids ; i++){
      child[i] = child_tmp[i];
    }
  }
  else {
    reorder(child_tmp,child,d_vals,kids);
    delete [] d_vals;
  }
}

/*
void binary_branch(QAP *qap,QAPAssignment &assign,MAT *X,MAT *U,
		   QAPAssignment *child, int &kids,
		   double bound, QPBParameters *qpbp, int use_best,
		   int how,bool use_J,bool reorder_children,
		   int &loc_in_row, int &loc_kept)
{
  int n = assign.n;
  int i,j;

  for(i=0;i<n;i++) {
    for(j=0;j<n;j++)
      if(!assign.disallowed(i,j)) break;
    if(!assign.disallowed(i,j)) break;
  }

  printf("I want to disallow %d,%d\n",i,j);
  
  kids = 2;
  loc_in_row = n;
  loc_kept = 2;

  child[0] = assign;
  child[0].fix(i,j);

  child[1] = assign;
  child[1].disallow(i,j);

}
*/

//====================================================================



