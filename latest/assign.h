#ifndef ASSIGN_H
#define ASSIGN_H

#include <iostream>
#include "meschach.h"
#undef catch
#include "qap.h"
#include "util.h"

/*
============================================================

QAPAssignment:

Allows for storage and operation on 'partial assignments'
or 'partial permutations'.

p:    An array of integers.  If p[i] >=0, then the assignment
      i-->p[i] has been made.  If p[i]<0, then no assignment 
      has been made yet, and row i of J has some unassigned
      columns.

pinv: "p inverse".  If pinv[i] >= 0, then the assignment
      pinv[i]-->i has been made.

n:    The size of arrays p, pinv, and the number of rows and
      columns in J.

nfix: The number of assignments made to this point.  Should
      be equal to the number of nonnegative elements in p.

depth: depth of this node in the branch and bound tree.  This
       is not necessarily the same as nfix.

bound: lower bound on this problem.

predicted_bound: an estimate of the lower bound on this problem;
                 lower bound on the lower bound!  We use this 
                 value to help decide which branching strategy 
		 to use.

history:  order in which assignments are made...the first
          facility to be fixed is in history[0].

OPTIONAL:  if STORE_J_MATRIX is defined:
J:    A matrix showing which assignments have been made, and
      which are disallowed.  
         J[i][j] == ASSIGNED means i-->j.
         J[i][j] == DISALLOWED means i cannot be assigned to j.
         J[i][j] == ALLOWED means i may be assigned to j.


============================================================
*/

typedef int Index;
typedef char Jtype;

#define ALLOWED     '-'
#define DISALLOWED  '0'
#define ASSIGNED    '1'
#define REALLY_BIG  (1.0e5)

// constants for p and pinv.  If p[i] == NONE, then all assignments
// in row i are permitted.  If p[i] == PARTIAL, some of them have
// been disallowed -- scan J[i] to find out which ones.
const int NONE = -1;
const int PARTIAL = -2;

/*
============================================================
WarmStart: the Frank-Wolfe point X of a node that was branched on,
shared by its children. A child starts its own bound there
(init_X_from_parent, with its fixed row and column removed) instead of
at the uniform matrix. It holds the parent's assignment p so the child
can tell which row and column it fixed, and is reference counted, so
it goes away with the last child. It lives in the node (and in state
files) rather than in a cache, so a node's bound does not depend on
which process computed its parent.
============================================================
*/
class QAPAssignment;

struct WarmStart
{
  int refs;
  int n;          // instance size
  Index *p;       // the parent's assignment
  MAT *X;         // the parent's X: its free rows x free columns, in order
  // The parent's final DNN state (a DNNWarm, see dnn.h) when its DNN bound
  // ran and missed, else NULL; owned, freed by dnn_free. Not kept in state
  // files: it only changes where a child's ADMM starts.
  void *dnn;
  void (*dnn_free)(void *);
};

/*
 LowerSum: a sum of products and terms, rounded to a double that is a
 lower bound on the exact sum. Accumulates in long double. The result is
 exact when every term is a multiple of 1/4 (operands multiples of 1/2:
 integer data, or integer data symmetrized by halves) and the sum of
 magnitudes stays below 2^(LDBL_MANT_DIG - 2): then every product and
 partial sum is representable. Otherwise it subtracts an explicit bound
 on the long double rounding, (k + 2) u_ld sum|terms| for k terms, and
 rounds down. Bounds on subproblems built from lower-bounded linear terms
 and shift are valid for the exact subproblem, since assignment variables
 are nonnegative.
*/
struct LowerSum
{
  long double sum, abs;
  int terms;
  bool dyadic;
  LowerSum() : sum(0), abs(0), terms(0), dyadic(true) {}
  void add(double t);
  void add_product(double a, double b);
  double lower() const;
};

// position of the row and column node fixed among its WarmStart parent's
// free ones, if node differs from that parent by exactly one assignment
bool warm_child_offsets(const QAPAssignment &node, int &di, int &dj);

WarmStart *new_warm_start(const Index *p, int n, MAT *X);
void warm_start_retain(WarmStart *w);
void warm_start_release(WarmStart *&w);

class QAPAssignment
{
  

public:
  QAPAssignment();
  QAPAssignment(int n);
  QAPAssignment(const QAPAssignment &s);
  ~QAPAssignment();

  void clear();
  void operator =(const QAPAssignment &source);

  void fix(int i,int j);      // fix X_ij to 1
  void unfix(int i,int j);      // fix X_ij to 1

  bool isPermutation();
  bool isInvalid();
  int adjust();               // find additional vars. to be fixed to 1.

  void print(std::ostream& outs = std::cout, bool printJ = false);
  friend std::ostream& operator <<(std::ostream& outs, const QAPAssignment &a);
  friend bool operator ==(const QAPAssignment &a,const QAPAssignment &b);

  Index *p;                   // permutation p
  Index *pinv;                // inv(p)

  Index *history; // first assignment was history[0]-->history[p[0]]

  int n;
  int nfix;

  int depth;
  double bound;
  double predicted_bound;

  WarmStart *warm;            // parent's X, or NULL (see WarmStart)
  void set_warm(WarmStart *w);

  int *get_unassigned(bool rows) const;
  
#ifdef STORE_J_MATRIX
  Jtype **J;
  // fix X_ij to 0
  inline void disallow(int i,int j) {J[i][j] = DISALLOWED;}; 

  // we don't know what X_ij should be.
  inline void allow(int i,int j) {J[i][j] = ALLOWED;}; // fix X_ij to 0

  inline bool disallowed(int i,int j) const {
    return (pinv[j]!=NONE)||(p[i]!=NONE)||(J[i][j]==DISALLOWED);}; 
            // is X_ij = 0?
  //  inline bool disallowed(int i,int j) const {return false;}; // is X_ij = 0?
#else
  inline bool disallowed(int i,int j) const {return false;}; // is X_ij = 0?
#endif

private:
  void pack(int i,int j);
};

void write_permutation(std::ostream &outs,Index *p,int n);
void read_permutation(std::istream &outs,Index *p,int n);

QAP *reduceQAP(QAP *q,const QAPAssignment &a);
void reduceQAP(QAP *q,QAP *q1,const QAPAssignment &a);
void reduceA(MAT *A,MAT *&A1,int *p,int m,int n);
void reduceB(MAT *B,MAT *&B1,int *pinv,int m,int n);
void reduceC(QAP *q, QAP *q1,const QAPAssignment &a);

int numberUnassigned(int *p,int i);
MAT *submatrix(MAT *A,int *row,int *col,int m,int n);
void submatrix(MAT *A,MAT *&Asub,int *row,int *col,int m,int n);

double QAPobjective(QAP *q,const QAPAssignment &a);
double QAPobjective(QAP *q,Index *p,int m);
double QAPobjective(MAT *A, MAT *B, MAT *C, double shift, Index *p, int m);

QAPAssignment randomly_fix(const QAPAssignment &node, int count);

#endif
