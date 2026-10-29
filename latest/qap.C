#include "qap.h"
#include <iostream>
#include <fstream>
#include "linalg.h"
#include "util.h"
#include <math.h>

QAP::QAP()
{
  shift = 0.0;
  A = NULL;
  B = NULL;
  C = NULL;
}

QAP::QAP(int n)
{
  A = m_get(n,n);
  B = m_get(n,n);
  C = m_get(n,n);
  shift = 0.0;
}

QAP::QAP(const QAP &qap)
{
  A = m_copy(qap.A,MNULL);
  B = m_copy(qap.B,MNULL);
  C = m_copy(qap.C,MNULL);
  shift = qap.shift;
}

QAP::~QAP()
{
  clear();
}

void QAP::clear()
{
  M_FREE(A); A = MNULL;
  M_FREE(B); B = MNULL;
  M_FREE(C); C = MNULL;
}

void printQAP(QAP *qap)
{
  std::cout << std::endl;

  //  m_output(qap->A);
  //  m_output(qap->B);
  //  m_output(qap->C);

  my_nice_m_output(qap->A,"A = ",qap->A->m,qap->A->n);
  my_nice_m_output(qap->B,"B = ",qap->B->m,qap->B->n);
  my_nice_m_output(qap->C,"C = ",qap->C->m,qap->C->n);
  
  std::cout << std::endl;

}

/*
============================================================
 q = readQAPLIB(filename);

 Read a file with name 'filename' in QAPLIB format into
 the QAP object q.

============================================================
 */
QAP *readQAPLIB(char *filename)
{
  std::ifstream f;
  int i,j,n,in;
  double val;
  QAP *qap;

  //  cout << "Input File: [" << filename << "]" << std::endl;

  f.open(filename);

  if(f)
    {
      f >> n;

      qap = new QAP(n);

      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	  {
	    f >> val;
	    qap->A->me[i][j] = val;
	  }

      if(!f) {
	std::cout << "readQAPLIB: error reading A" << std::endl; return qap;
      }

      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	  {
	    f >> val;
	    qap->B->me[i][j] = val;
	  }
      if(!f) { std::cout << "readQAPLIB: error reading B" << std::endl; return qap; }

      // read C, if it has been supplied.  Otherwise C = 0.
      f >> val;
      if(f)
	{
	  for(i=0;i<n;i++)
	    for(j=0;j<n;j++)
	      {
		qap->C->me[i][j] = val;
		f >> val;
	      }
	}
    }
  else
    {
      std::cout << "readQAPLIB: could not open file [" 
	   << filename << "]" << std::endl;
      return NULL;
    }
  f.close();
  
  return qap;
}

static bool is_symmetric(MAT *M)
{
  for (unsigned int i = 0; i < M->m; i++)
    for (unsigned int j = i + 1; j < M->n; j++)
      if (M->me[i][j] != M->me[j][i]) return false;
  return true;
}

static bool is_integer(MAT *M)
{
  for (unsigned int i = 0; i < M->m; i++)
    for (unsigned int j = 0; j < M->n; j++)
      if (M->me[i][j] != floor(M->me[i][j])) return false;
  return true;
}

static bool is_even(double v) { return fmod(fabs(v), 2.0) == 0.0; }

/*
============================================================
  set_fathom_tolerance(qap);

  Set DELTA / DELTA_BRANCH (see util.h) from the data. Every
  objective value is an integer when A, B, C and the shift are
  integers. It is also even when A and B are symmetric (the
  off-diagonal terms come in equal pairs), every a_ii * b_kk is
  even, and C and the shift are even. Call before symmetrize_QAP,
  which can introduce halves.
============================================================
*/
void set_fathom_tolerance(QAP *qap)
{
  const int n = qap->A->m;
  // Integer data only gives integer computed values if every product and
  // partial sum is exact in double: |objective| <= sum|A| max|B| + sum|C| +
  // |shift| must stay below 2^53. Otherwise rounding can exceed DELTA.
  double sumA = 0, maxB = 0, sumC = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      sumA += fabs(qap->A->me[i][j]);
      maxB = std::max(maxB, fabs(qap->B->me[i][j]));
      sumC += fabs(qap->C->me[i][j]);
    }
  bool exact = sumA * maxB + sumC + fabs(qap->shift) < 9007199254740992.0;  // 2^53
  bool integer = exact && is_integer(qap->A) && is_integer(qap->B) &&
                 is_integer(qap->C) && qap->shift == floor(qap->shift);
  bool even = integer && is_symmetric(qap->A) && is_symmetric(qap->B) &&
              is_even(qap->shift);
  for (int i = 0; even && i < n; i++)
    for (int k = 0; even && k < n; k++)
      even = is_even(qap->A->me[i][i] * qap->B->me[k][k]) &&
             is_even(qap->C->me[i][k]);

  DELTA = DELTA_BRANCH = even ? 1.95 : (integer ? 0.95 : 0.0);
  std::cout << "objective values are "
            << (even ? "even integers" : (integer ? "integers" :
                (exact ? "not integers" : "too large to compute exactly")))
            << ": fathoming tolerance " << DELTA << std::endl;
}

/*
============================================================
  symmetrize_QAP(qap);

  The QP bound, the eigenvalue computations and the gradient
  (2(AXB - SX - XT) + C) all assume A and B are symmetric. When
  one of them is, replacing the other by (M + M')/2 leaves every
  objective value unchanged: sum a_ij b_p(i)p(j) pairs (i,j) with
  (j,i). Returns false when neither is symmetric, in which case
  the bounds are not valid.
============================================================
*/
bool symmetrize_QAP(QAP *qap)
{
  MAT *M = MNULL;
  if (is_symmetric(qap->A) && !is_symmetric(qap->B))
    M = qap->B;
  else if (is_symmetric(qap->B) && !is_symmetric(qap->A))
    M = qap->A;
  else
    return is_symmetric(qap->A);

  for (unsigned int i = 0; i < M->m; i++)
    for (unsigned int j = i + 1; j < M->n; j++)
      M->me[i][j] = M->me[j][i] = 0.5 * (M->me[i][j] + M->me[j][i]);
  std::cout << "symmetrized " << (M == qap->A ? "A" : "B")
            << " (objective values are unchanged)" << std::endl;
  return true;
}

bool writeQAPLIB(char *filename, QAP *qap)
{
  std::ofstream f;
  int i,j,n,in;
  double val;

  //  cout << "Output File: [" << filename << "]" << std::endl;
  n = qap->A->n;

  if(n<=0) return false; // failure

  f.open(filename);

  if(f)
    {
      f << n << std::endl << std::endl;

      for(i=0;i<n;i++) {
	for(j=0;j<n;j++)
	  {
	    f << qap->A->me[i][j] << " ";
	  }
	f << std::endl;
      }
      f << std::endl;
      if(!f) { std::cout << "writeQAPLIB: error writing A" << std::endl; return false; }

      for(i=0;i<n;i++) {
	for(j=0;j<n;j++)
	  {
	    f << qap->B->me[i][j] << " ";
	  }
	f << std::endl;
      }
      f << std::endl;
      if(!f) { std::cout << "writeQAPLIB: error writing B" << std::endl; return false; }

      for(i=0;i<n;i++) {
	for(j=0;j<n;j++)
	  {
	    f << qap->C->me[i][j] << " ";
	  }
	f << std::endl;
      }
      f << std::endl;
      if(!f) { std::cout << "writeQAPLIB: error writing C" << std::endl; return false; }

    }
  else
    {
      std::cout << "writeQAPLIB: could not open file [" 
	   << filename << "]" << std::endl;
      return false;
    }
  f.close();
  return true;
}

/*
A' = A + eg' + ge' + Diag(r)
B' = B + eh' + he' + Diag(s)
C' = C - 2[Aeh' + ge'B + gs' + rh' + ngh' + (e'g)eh'] - [as' + rb' + rs']

where a = diag(A), b = diag(B).
*/
void applyPerturbations(QAP *qap, double *r, double *s, 
		       double *g, double *h)
{
  int i,j,n;

  n = qap->A->m;

  for(i=0; i<n; i++) {
    qap->A->me[i][i] += r[i];
    qap->B->me[i][i] += s[i];
  }

  if(g&&h) {
    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
	qap->A->me[i][j] += g[i] + g[j];
	qap->B->me[i][j] += h[i] + h[j];
      }
    }
  }

  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      qap->C->me[i][j] -= qap->A->me[i][i]*s[j] + r[i]*qap->B->me[j][j] 
	                  + r[i]*s[j];
      // plus some other stuff from g and h...
    }
  }
  
}


// unused.
void negateQAP(QAP *qap)
{
  int i,j;
  for(i=0;i<qap->A->m;i++)
    for(j=0;j<qap->A->n;j++)
      {
	qap->A->me[i][j] = -qap->A->me[i][j];
	qap->B->me[i][j] = -qap->B->me[i][j];
	qap->C->me[i][j] = -qap->C->me[i][j];
      }
}
