#include "qap.h"
#include <iostream>
#include <fstream>
#include "linalg.h"

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
