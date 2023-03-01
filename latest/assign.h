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
