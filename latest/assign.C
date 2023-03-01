#include "assign.h"
#include "util.h"

QAPAssignment::QAPAssignment()
{
  n = 0;
  nfix = 0;
  depth = 0;
  bound = 0.0;
  predicted_bound = 0.0;

  p = NULL;
  pinv = NULL;
  history = NULL;

#ifdef STORE_J_MATRIX
  J  = NULL;
#endif
}

/*
============================================================

 QAPAssignment a(n);

 A new QAP assignment of size n is created.  All assignments
 are permitted initially.
============================================================
 */
QAPAssignment::QAPAssignment(int m)
{
  int i;
  n = m;
  nfix = 0;
  depth = 0;
  bound = 0.0;
  predicted_bound = 0.0;

  p = new Index[n];
  pinv = new Index[n];
  history = new Index[n];

  for(i=0;i<n;i++)
    {
      p[i] = NONE;
      pinv[i] = NONE;
      history[i] = NONE;
    }

#ifdef STORE_J_MATRIX
  J  = new Jtype*[n];
  for(i=0;i<n;i++)
    {
      J[i] = new Jtype[n];
      for(int j=0;j<n;j++)
	J[i][j] = ALLOWED;
    }
#endif
}

/*
============================================================
 QAPAssignment a1(a);
============================================================
 */
QAPAssignment::QAPAssignment(const QAPAssignment &s)
{
  n = s.n;
  nfix = s.nfix;
  depth = s.depth;
  bound = s.bound;
  predicted_bound = s.predicted_bound;

  p = copy_array<Index>(s.p,n);
  pinv = copy_array<Index>(s.pinv,n);
  history = copy_array<Index>(s.history,n);

#ifdef STORE_J_MATRIX
  J  = new Jtype*[n];
  for(int i=0;i<n;i++)
    {
      J[i] = copy_array<Jtype>(s.J[i],n);
    }
#endif
}

/*
============================================================
  a = a1;
============================================================
 */
void QAPAssignment::operator =(const QAPAssignment &s)
{
  if(s.p != p)
  {
    clear();

    n = s.n;
    nfix = s.nfix;
    depth = s.depth;
    bound = s.bound;
    predicted_bound = s.predicted_bound;

    p = copy_array<Index>(s.p,n);
    pinv = copy_array<Index>(s.pinv,n);
    history = copy_array<Index>(s.history,n);

#ifdef STORE_J_MATRIX
    J  = new Jtype*[n];
    for(int i=0;i<n;i++)
    {
      J[i] = copy_array<Jtype>(s.J[i],n);
    }
#endif 
  }
}

/*
============================================================
  A.clear();

  resets A to a size zero QAPAssignment.
============================================================
 */
void QAPAssignment::clear()
{
  delete [] p;   p = NULL;
  delete [] pinv;   pinv = NULL;
  delete [] history;   history = NULL;

  n = 0; nfix = 0; depth = 0;

#ifdef STORE_J_MATRIX
  for(int i=0;i<n;i++)
    {
      delete [] J[i];
    }
  delete [] J;
  J = NULL;
#endif
}

/*
============================================================
(destructor)
============================================================
 */
QAPAssignment::~QAPAssignment()
{
  clear();
}

/*
============================================================
 a.fix(i,j);

 fix X_ij to 1, i.e. make the assignment i-->j.
============================================================
 */
void QAPAssignment::fix(int i,int j) 
{
  p[i] = j;
  pinv[j] = i;

  history[nfix] = i;
  depth++;
  nfix++;

  // we could loop thru row i and col j of J, filling in with
  // DISALLOWED OR FIXED, but it won't make any difference.

}

void QAPAssignment::unfix(int i,int j) 
{
  p[i] = NONE;
  pinv[j] = NONE;

  depth--;
  nfix--;
  history[nfix] = NONE;

  // we could loop thru row i and col j of J, filling in with
  // DISALLOWED OR FIXED, but it won't make any difference.

}

/*
============================================================
 ok = a.isPermutation();

 Returns true when all facilities have been assigned 
 locations.

============================================================
 */
bool QAPAssignment::isPermutation()
{
  for(int i=0;i<n;i++)
    if(p[i] < 0)
      return false;
  return true;
}

/*
============================================================
 ok = a.isInvalid();

 Returns false if the assignment does not correspond to
 a valid partial solution to a QAP.  That is, if all 
 assignments in a given row or column are disallowed.

 (currently unused)
============================================================
 */
bool QAPAssignment::isInvalid()
{

  return false;
}

/*
============================================================
 a.adjust();
 (unused)
============================================================
 */
int QAPAssignment::adjust()
{
  return 1;
}

int *get_unassigned_core(int *p, int n, int m) {
  int *q = new int[m];
  int k = 0;
  for (int i = 0; i < n; i++) {
    if (p[i] == NONE) {
      q[k++] = i;
    }
  }
  return q;
}

int *QAPAssignment::get_unassigned(bool rows) const
{
  return get_unassigned_core(rows ? p : pinv, n, n - nfix);
}

void write_permutation(std::ostream &outs,Index *p,int n)
{
  int i;
  outs << "[";
  
  for(i=0;i<n;i++)
    {
      if(p[i]==NONE)
	outs << "- ";
      else if(p[i]==PARTIAL)
	outs << "? ";
      else
	outs << (int)p[i] << " ";
    }
  outs << "]";
}

// assumes p is allocated to appropriate size
void read_permutation(std::istream &ins,Index *p,int n)
{
  char c = '\0';
  char str[10];
  int i;

  while(ins &&(c != '[')) ins >> c;
  
  if(!ins) return;

  for(i=0;i<n;i++)
    {
      ins >> str;

      if(str[0] == '-')
	p[i] = NONE;
      else if(str[0] == '?')
	p[i] = PARTIAL;
      else
	p[i] = atoi(str);
    }
  while(ins &&(c != ']')) ins >> c;
}

void QAPAssignment::print(std::ostream& outs, bool printJ)
{
  outs << "depth = " << depth 
       << "  bnd = " << bound
       << "  est.bnd = " << predicted_bound
       << std::endl;

  outs << "  ";  write_permutation(outs,p,n);  outs << std::endl;
  outs << "  ";  write_permutation(outs,pinv,n);  outs << " (inv)" <<std::endl;
  outs << "  ";  write_permutation(outs,history,n);  outs << " (history)" <<std::endl;

#ifdef STORE_J_MATRIX
  int i,j;
  if(printJ)
    for(i=0;i<n;i++)
      {
	for(j=0;j<n;j++)
	  outs << J[i][j] << " ";
	outs << std::endl;
      }
#endif
}

std::ostream& operator <<(std::ostream& outs, const QAPAssignment &a)
{
  outs << "  ";  write_permutation(outs,a.p,a.n);  outs << std::endl;
  outs << "  ";  write_permutation(outs,a.pinv,a.n);  outs << " (inv)" <<std::endl;
  outs << "  ";  write_permutation(outs,a.history,a.n);  outs << " (history)" <<std::endl;

#ifdef STORE_J_MATRIX
  int i,j;
  for(i=0;i<a.n;i++)
    {
      for(j=0;j<a.n;j++)
	outs << a.J[i][j] << " ";
      outs << std::endl;
    }
#endif
  
  return outs;
}

/*
============================================================
  a==b
  
  Tests equality of two assignments.  Unused at this point.
============================================================
 */
bool operator ==(const QAPAssignment &a,const QAPAssignment &b)
{
  return true;
}


/*
  find the size of the set
   { j | p[j] == NONE, j = 1..i}
  that is, the number of unassigned entries of i, including i itself.
 */
int numberUnassigned(int *p, int i)
{
  int ctr,j;
  for(ctr=0,j=0;j<=i;j++)
    if(p[j]==NONE) ctr++;
  return ctr;
}



/*
============================================================

   Return a matrix B = A(r,c) where
    r = {i | row[i]<0 }
    c = {j | col[j]<0 }
  
   m: (size of matrix to create) = (# true vals in row and col)
   n: size of A = length of row and col.
============================================================
 */
MAT *submatrix(MAT *A,int *row,int *col,int m,int n)
{
  MAT *B;
  int i,ii,j,jj;

  B = m_get(m,m);
  for(ii=-1,i=0;i<n;i++) {
    if(row[i]<=NONE) {
      ii++;
      for(jj=-1,j=0;j<n;j++) {
	if(col[j]<=NONE) {
	  jj++;
	  //printf("ii jj %d %d   i j %d %d  size %d %d   sizeA %d %d\n", ii, jj, i, j, B->m, B->n, A->m, A->n);
	  B->me[ii][jj] = A->me[i][j];
	}
      }
    }
  }
  return B;
}

void submatrix(MAT *A,MAT *&Asub,int *row,int *col,int m,int n)
{
  int i,ii,j,jj;

  Asub = m_resize(Asub,m,m);
  for (ii=-1,i=0; i<n; i++) {
    if (row[i] <= NONE) {
      ii++;
      for (jj=-1,j=0; j<n; j++) {
	if (col[j] <= NONE) {
	  jj++;
	  Asub->me[ii][jj] = A->me[i][j];
	}
      }
    }
  }
}

/*
============================================================
addPenalties(C,assign);

for all disallowed assignments in assign, add a large
constant to the corresponding location in the matrix C,
ensuring that these assignments will not be considered
by the lower bounding procedure.
============================================================
*/
#ifdef STORE_J_MATRIX
inline void addPenalties(MAT *C, const QAPAssignment &a)
{
  int i,ii,j,jj,n=a.n;
  // add on penalties for disallowed assignments
  for(ii=-1,i=0;i<n;i++)
    {
      if(a.p[i] > NONE) continue;
      ii++; // the ii_th free row.
      for(jj=-1,j=0;j<n;j++)
	{
	  if(a.pinv[j] > NONE) continue;
	  jj++;// the jj_th free column
	  if(a.disallowed(i,j))
	    {
	      C->me[ii][jj] = REALLY_BIG;
	    }
	}
    }
}
#endif

/*
============================================================
 add linear term due to fixing variables
============================================================
 */
inline void addLinearTerm(QAP *q, QAP *q1,const QAPAssignment &a)
{
  int i,ii,j,jj,k;

  for(k=0;k<a.n;k++)
    {
      if(a.p[k]>NONE) {
	for(ii=-1,i=0;i<a.n;i++)
	  {
	    if(a.p[i]<=NONE) {
	      ii++; // the ii_th free row.
	      for(jj=-1,j=0;j<a.n;j++)
		{
		  if(a.pinv[j]<=NONE) {
		    jj++;// the ii_th free column
		    q1->C->me[ii][jj] += 
		      q->A->me[i][k]*q->B->me[j][a.p[k]] +
		      q->A->me[k][i]*q->B->me[a.p[k]][j];
		  }
		}
	    }
	  }
      }
    }
}

/*
============================================================
 compute constant term due to fixing variables
============================================================
 */
inline void computeShift(QAP *q, QAP *q1,const QAPAssignment &a)
{
  int i,j;
  // compute shift
      
  for(i=0;i<a.n;i++)
    {
      if(a.p[i]>NONE) {
	
	q1->shift += q->C->me[i][a.p[i]];
	for(j=0;j<a.n;j++)
	  {
	    if(a.p[j]>NONE) {
	      q1->shift += q->A->me[i][j]*q->B->me[a.p[i]][a.p[j]];
	    }
	  }
      }
    }
}

/*
============================================================
 qr = reduceQAP(q,a);

 Given a QAP q = (A,B,C), produce a (possibly smaller) 
 QAP qr = (A',B',C') so that ...

============================================================
 */
QAP *reduceQAP(QAP *q,const QAPAssignment &a)
{
  QAP *qr;
  qr = new QAP;

  reduceQAP(q,qr,a);

  return qr;
}


/*
============================================================
 reduceQAP(q,qr,a);

 Given a QAP q = (A,B,C), produce a (possibly smaller) 
 QAP qr = (A',B',C') that takes into account the assignments
 made in the given QAPAssignment a.  For each assignment
 i-->j, row and column i are removed from A, row and column
 j are removed from B, and a linear term is added on to C.
 For each disallowed assignment i-\->j, a large positive
 term is added to C_ij.

============================================================
 */
void reduceQAP(QAP *q,QAP *q1,const QAPAssignment &a)
{
  int ii,jj,i,j,k,n,m;

  n = a.n; 
  m = n - a.nfix; 
  q1->shift = 0.0;

  if(a.nfix>0)
    {
      // compute shift

      for(i=0;i<n;i++)
	{
	  if(a.p[i]>NONE) {

	    q1->shift += q->C->me[i][a.p[i]];
	    for(j=0;j<n;j++)
	      {
		if(a.p[j]>NONE) {
		  q1->shift += q->A->me[i][j]*q->B->me[a.p[i]][a.p[j]];
		}
	      }
	  }
	}
      
      // get q1 = (A1,B1,C1)
      M_FREE(q1->A);M_FREE(q1->B);M_FREE(q1->C);

      q1->A = submatrix(q->A,a.p,a.p,m,n);
      q1->B = submatrix(q->B,a.pinv,a.pinv,m,n);
      q1->C = submatrix(q->C,a.p,a.pinv,m,n);
      // add linear term due to fixing variables
      for(k=0;k<n;k++)
	{
	  if(a.p[k]>NONE) {
	    for(ii=-1,i=0;i<n;i++)
	      {
		if(a.p[i]<=NONE) {
		  ii++; // the ii_th free row.
		  for(jj=-1,j=0;j<n;j++)
		    {
		      if(a.pinv[j]<=NONE) {
			jj++;// the ii_th free column
			q1->C->me[ii][jj] += 
			  q->A->me[i][k]*q->B->me[j][a.p[k]] +
			  q->A->me[k][i]*q->B->me[a.p[k]][j];
		      }
		    }
		}
	      }
	  }
	}
    }
  else
    {
      M_FREE(q1->A);M_FREE(q1->B);M_FREE(q1->C);
      q1->A = m_copy(q->A,MNULL);
      q1->B = m_copy(q->B,MNULL);
      q1->C = m_copy(q->C,MNULL);
    }

#ifdef STORE_J_MATRIX
  // add on penalties for disallowed assignments
  for(ii=-1,i=0;i<n;i++)
    {
      if(a.p[i] > NONE) continue;
      ii++; // the ii_th free row.
      for(jj=-1,j=0;j<n;j++)
	{
	  if(a.pinv[j] > NONE) continue;
	  jj++;// the jj_th free column
	  if(a.disallowed(i,j))
	    {
	      q1->C->me[ii][jj] = REALLY_BIG;
	    }
	}
    }
#endif

}

// m = size of A, n = size of A1
void reduceA(MAT *A, MAT *&A1, int *p, int m, int n)
{
  M_FREE(A1);
  if (m != n) {
    A1 = submatrix(A, p, p, m, n);
  }
  else {
    A1 = m_copy(A, MNULL);
  }
}

void reduceB(MAT *B, MAT *&B1, int *pinv, int m, int n)
{
  M_FREE(B1);
  if (m != n) {
    B1 = submatrix(B,pinv,pinv,m,n);
  }
  else {
    B1 = m_copy(B,MNULL);
  }
}

// assumes A',B' already set.
void reduceC(QAP *q, QAP *q1,const QAPAssignment &a)
{
  int n,m;

  n = a.n; 
  m = n - a.nfix; 
  q1->shift = 0.0;

  M_FREE(q1->C);
  if(a.nfix>0)
    {
      computeShift(q,q1,a);
      q1->C = submatrix(q->C,a.p,a.pinv,m,n);
      addLinearTerm(q,q1,a);
    }
  else
    {
      M_FREE(q1->C);
      q1->C = m_copy(q->C,MNULL);
    }
#ifdef STORE_J_MATRIX
  addPenalties(q1->C,a);
#endif
}

/*
============================================================
  obj = QAPobjective(q,a);

  Compute tr AXBX' + CX', where the permutation matrix X
  is determined by a.
============================================================
 */
double QAPobjective(QAP *q,const QAPAssignment &a)
{
  double obj = q->shift;
  int i,j,m=a.n;

  for(i=0;i<m;i++)
    {
      obj += q->C->me[i][a.p[i]];
      for(j=0;j<m;j++)
	{
	  obj += q->A->me[i][j]*q->B->me[a.p[i]][a.p[j]];
	}
    }
  return obj;
}

/*
============================================================
  obj = QAPobjective(q,a);

  Compute tr AXBX' + CX', where the permutation matrix X
  is determined by the array p: X(i,p[i]) = 1.
============================================================
 */
double QAPobjective(QAP *q,Index *p,int m)
{
  return QAPobjective(q->A, q->B, q->C, q->shift, p, m);
}


double QAPobjective(MAT *A, MAT *B, MAT *C, double shift, Index *p, int m)
{
  double obj = shift;
  for (int i = 0;i < m; i++)
    {
      obj += C->me[i][p[i]];
      for (int j = 0; j < m; j++)
	{
	  obj += A->me[i][j] * B->me[p[i]][p[j]];
	}
    }
  return obj;
}

QAPAssignment randomly_fix(const QAPAssignment &node, int count)
{
  int m = node.n - node.nfix;
  int *u_r = node.get_unassigned(true);
  shuffle(u_r, m);
  int *u_c = node.get_unassigned(false);
  shuffle(u_c, m);
  
  QAPAssignment child = node;
  for (int k = 0; k < count; k++) {
    child.fix(u_r[k], u_c[k]);
  }
  delete [] u_r;
  delete [] u_c;
  return child;
}

