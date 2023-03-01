/*
============================================================


============================================================
*/

#ifndef UTIL_H
#define UTIL_H

#include <stdarg.h>
#include <iostream>
#include <time.h>
#include <fstream>

/*
============================================================
  Constants and options:
============================================================
*/
static const char *qap_version_info = "[Jan 1, 2023]";


// If we can only solve LAPs with integer input, define this
//#define INTEGER_LAP

// If the EVB3 code is available (in dir EVB3), define this
//#define HAVE_EVB3

// if solution is even, can fathom more aggressively
#define SOLUTION_IS_EVEN          

#ifdef SOLUTION_IS_EVEN
const double DELTA        = 1.95;
const double DELTA_BRANCH = 1.95;
#else
const double DELTA        = 0.95;
const double DELTA_BRANCH = 0.95;
#endif

const double ROLLOVER = 4295.0;

// for branching rules 1, 3 and 4
const bool use_weighted_branching = true; 
const double branching_weight = 0.50;

// print cpu time every PRINT_MESSAGE_EVERY nodes
const int PRINT_MESSAGE_EVERY = 10000000;  
//const int PRINT_MESSAGE_EVERY = 100000;  
//const int PRINT_MESSAGE_EVERY = 1;  

// if algorithm is run only to a given depth,
// are nodes greater than that depth saved?
const bool save_too_deep = true;  // est. doesn't work if false

// define this if you want to be able to disallow assignments from
// being made -- downside: O(n^2) storage at each node, doesn't seem
// to help.
//#define STORE_J_MATRIX

// bound types
typedef int bound_t;
const bound_t QP_BND       = 0;
const bound_t QP_BND_PARAM = 1;
const bound_t QP_GLB_BND   = 2;
const bound_t QP_BND_IMP   = 3;
const bound_t GLB_BND      = 4;
const bound_t PB_BND       = 5;
const bound_t EVB3_BND     = 6;
const bound_t BUFF_BND     = 7;

// branch types
const int BRANCH_U_SUM = 1;
const int BRANCH_U_FATHOM = 2;
const int BRANCH_STRONG = 3;
const int BRANCH_LOOKAHEAD = 4;
const int BRANCH_BUFF = 5;
const int BRANCH_FIXED = -1;


/*
  changed 12/15/99.  Added a relative_gap parameter.  Use
  this branching rule if
        (opt-bnd)/opt >= relative_gap.
 */
struct BNBParameters 
{
  double relative_gap;
  int depth;

  // bound-related:
  int which_bound;         // use parametric or nonparametric QPB?
  bool compute_lapU;       // whether or not to solve LAP(U)
  int maxier_bnd;          // maximum # Frank-Wolfe iterations
  int maxier_nofathom;     // maximum # Frank-Wolfe iterations, when
                           //  we know we can't fathom current node
  int param_step;          // update S and T every step iterations

  // branching-related:
  int maxier_br;           // max FW iers, when computing B
  int rule;                // row selection rule
  int use_best;            // num of rows of U to look at

  // estimator-related
  double est_factor;
  double score_exponent;

  char *log_file;            // console log file (optional)
  std::string log_prefix;          // prefix for log entries
};

// for testing purposes...
class stopwatch
{
public:

  stopwatch() {t_total = 0;};
  
  void go() {t_start = clock();};
  double stop() {t_end = clock(); double d=diff(t_start,t_end); 
                 t_total +=d; return d;};
  double elapsed() {return t_total;};

private:
  double diff(double start_time, double stop_time) {
    double d = (stop_time - start_time)/CLOCKS_PER_SEC;
    if(d < 0) d += ROLLOVER;
    return d;
  }

  double t_start;
  double t_end;

  double t_total;
};


inline bool is_qpbnd(bound_t bnd)
{
  return ( 
          (bnd==QP_BND) ||
          (bnd==QP_BND_PARAM) ||
          (bnd==QP_GLB_BND) ||
          (bnd==QP_BND_IMP) ||
          (bnd==EVB3_BND) 
          );
}

void qap_read_params(char *filename, BNBParameters *&params, int &nlevel,
		     char *&sym_file,int &max_depth);

void print_bnb_params(BNBParameters *&params,int nlevel);
void print_time(std::ostream &outs=std::cout);

void ins_sort_ascending(double *v, int *order, int n);
void ins_sort_descending(double *v, int *order, int n);

void nice_array_print(double *d, int n, char delim= ' ');
void nice_int_array_print(int *d,int n);
void nice_matrix_print(double **d, int m, int n);
void write_matrix_flat(std::ofstream &f, double **a, int m, int n, char delim);
void matrix_to_file(std::ofstream &f, double **a, int m, int n);
void matrix_to_ampl_file(std::ofstream &f, const char *name, double **a, int m, int n);
bool find(int *list, int n, int item);

inline double computeRelativeGap(double opt,double root,double bnd)
{ return (opt-bnd)/(opt-root); }

double drand();

double score(double opt,double root_bound,double predicted_bound,
	     double score_exponent);

void shuffle(int *p, int n);

/*
============================================================
  T *new_zero(int m)
  allocate an array of objects, and set them to 0.
  T must support assignment to 0.

  T *copy_array(T *d,int m)
  copy the first m items of d into a new array.

  since they are templated functions, they need to go here.
============================================================
 */

template <class T>
inline T *new_zero(int m)
{
  /*
  T *x;
  
  x = new T[m];
  for(int i=0;i<m;i++)
    x[i] = 0;
  */
  T *x, *y;
  
  x = new T[m]; y = x;
  for(int i=0;i<m;i++)
    (*y++) = 0;

  return x;
}

template <class T>
inline T *new_perm(int m)
{
  /*
  T *x;
  
  x = new T[m];
  for(int i=0;i<m;i++)
    x[i] = i;
  */
  T *x, *y;
  
  x = new T[m]; y = x;
  for(int i=0;i<m;i++)
    (*y++) = i;

  return x;
}

template <class T>
inline T *copy_array(T *d,int m)
{
  T *x;
  T *y, *z;
  
  x = new T[m]; y = x; z = d;
  for(int i=0;i<m;i++)
    (*y++) = (*z++);

  return x;
}

template <class T>
void copy_array(T *d,T *x,int m)
{
  T *y, *z;
  
  y = x; z = d;
  for(int i=0;i<m;i++)
    (*y++) = (*z++);
}

// return true if first m elements of a and b are equal.
// assumes != operator defined for class T.
template <class T>
bool equal_array(T *a,T *b,int m)
{
  for(int i=0;i<m;i++)
    if(*a++ != *b++)
      return false;
  return true;
}

#endif
