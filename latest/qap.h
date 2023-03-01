/*
============================================================
QAP:

Standard form for QAP:
        n   n                          n
       --- ---                        ---
   min \   \    a(i,j) b(p(i),p(j)) + \   c(i,p(i))
       /   /                          /
       --- ---                        ---
       i=1 j=1                        i=1

The QAP structure simply holds matrices A, B, C, which 
are of type MAT *.  MAT is a matrix data structure defined
in the Meschach library -- see the Meschach library for
details.  If A is of type MAT *, then A->m and A->n give
the dimensions of the matrix, and A->me[i][j] is the (i,j)
component of A.

============================================================
*/

#ifndef QAP_H
#define QAP_H

#include "meschach.h"
#undef catch
#include <string>

class QAP
{
public:
  QAP();
  QAP(int n);
  QAP(const QAP &qap);
  ~QAP();
  void clear();

  MAT *A;
  MAT *B;
  MAT *C;
  double shift;
};

void printQAP(QAP *qap);
void negateQAP(QAP *qap);
QAP *readQAPLIB(char *filename);
bool writeQAPLIB(char *filename, QAP *q);
void applyPerturbations(QAP *qap, double *r, double *s, 
		       double *g = NULL, double *h = NULL);

#endif 
