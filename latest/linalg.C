#include <iostream>
#include "linalg.h"

/*
============================================================
copies all rows but di, all columns but dj from X into X1.

[currently unused]
============================================================
*/
void warmstart(MAT *X, MAT *X1, int di, int dj)
{
  int i,j;

  for(i=0;i<di;i++)
    for(j=0;j<dj;j++)
      X1->me[i][j] = X->me[i][j];
  for(i=di+1;i<X->m;i++)
    for(j=0;j<dj;j++)
      X1->me[i-1][j] = X->me[i][j];
  for(j=dj+1;j<X->n;j++)
    for(i=0;i<di;i++)
      X1->me[i][j-1] = X->me[i][j];
  for(i=di+1;i<X->m;i++)
    for(j=dj+1;j<X->n;j++)
      X1->me[i-1][j-1] = X->me[i][j];

}

/*
========================================================================

 V = projectionMatrix(n);

 A matrix V of size n x n-1 is returned.  The columns of V form
 an orthonormal basis for the nullspace of e', that is,
 if V*a = b then e'b = 0, and V'*V = I.

 [currently unused, no need to explicitly compute V,
  use vtav, vatv instead.]
========================================================================
*/
MAT *projectionMatrix(int n)
{
  int i,j;
  MAT *V;

  V = m_get(n,n-1);

  for(i=0; i < n-1; i++)
    V->me[0][i] = -1.0/sqrt(n);
  for(i=1; i < n; i++)
    for(j=0;j<n-1;j++)
      V->me[i][j] = -1.0/(n+sqrt(n));
  for(i=0; i < n-1; i++)
    V->me[i+1][i] += 1.0;

  return V;
}

// compute V'AV
void vtav(MAT *A, MAT *OUT)
{
  int i,j;
  double k1,k2;
  int n = A->m, m = n-1;
  VEC *w;
  MAT *R;

  // ********

  w = v_get(n);
  R = m_get(n,m);

  k1 = -1.0/sqrt(n);
  k2 = -1.0/(n+sqrt(n));

  for(i=0;i<n;i++)
    w->ve[i] = k1*A->me[i][0];

  for(i=0;i<n;i++)
    for(j=1;j<n;j++)
      w->ve[i] += k2*A->me[i][j];

  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      R->me[j][i] = w->ve[j] + A->me[j][i+1];

  // ********

  for(i=0;i<m;i++)
    w->ve[i] = k1*R->me[0][i];

  for(i=0;i<m;i++)
    for(j=1;j<n;j++)
      w->ve[i] += k2*R->me[j][i];

  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      OUT->me[i][j] = w->ve[j] + R->me[i+1][j];

  V_FREE(w);
  M_FREE(R);
}

// compute VAV'
void vavt(MAT *A, MAT *OUT)
{
  double c_j,k1,k2;
  int i,j,n = A->m, m = n+1;
  MAT *R;

  R = m_get(m,n);

  k1 = -1.0/sqrt(m);
  k2 = -1.0/(m+sqrt(m));

  for(j=0;j<n;j++) {

    for(c_j=0.0,i=0;i<n;i++)
      c_j += A->me[i][j];

    R->me[0][j] = k1*c_j;
    c_j = k2*c_j;
    for(i=0;i<n;i++)
      R->me[i+1][j] = A->me[i][j] + c_j;
    
  }

  for(j=0;j<m;j++) {

    for(c_j=0.0,i=0;i<n;i++)
      c_j += R->me[j][i];
    
    OUT->me[0][j] = k1*c_j;
    c_j = k2*c_j;
    for(i=0;i<n;i++)
      OUT->me[i+1][j] = R->me[j][i] + c_j;
  }
    
  M_FREE(R);
}

void nice_m_output(MAT *A)
{
  for(int i=0; i<A->m; i++) {
    for(int j=0; j<A->n; j++) {
      printf("%.2f ",A->me[i][j]);
    }
    std::cout << std::endl;
  }
}

/*
========================================================================

 my_nice_m_output(X,"Here is X",m,n)

 Print the first m rows and first n columns of X in a nicer format
 than Meschach.

========================================================================
*/
void my_nice_m_output(MAT *X,char *title,int m,int n)
{
  int i,j;
  printf("   %s  = \n",title);
  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
	{
	  printf("%6.3f ",X->me[i][j]);
	}
      printf("\n");
    }
  printf("\n");
}

/* my_m_mlt -- matrix-matrix multiplication */
/*
========================================================================

(basically the meschach m_mlt routine, with the call to __mltadd__
replaced with its definition, and some of the error-check code
removed)

========================================================================
*/
MAT	*my_m_mlt(MAT *A,MAT *B,MAT *OUT)
{
	u_int	i, k, m, n, p;
	Real	**A_v, **B_v /*, *B_row, *OUT_row, sum, tmp */;

	m = A->m;	n = A->n;	p = B->n;
	A_v = A->me;		B_v = B->me;

	if ( OUT==(MAT *)NULL || OUT->m != A->m || OUT->n != B->n )
		OUT = m_resize(OUT,A->m,B->n);

	m_zero(OUT);
	for ( i=0; i<m; i++ )
		for ( k=0; k<n; k++ )
		{
		  if ( A_v[i][k] != 0.0 ) {
		    // replaced call to __mltadd__ with its code below...
		    //__mltadd__(OUT->me[i],B_v[k],A_v[i][k],(int)p);

		    register double *dp1, *dp2, s;
		    register int	l;
		    dp1 = OUT->me[i];
		    dp2 = B_v[k];
		    s = A_v[i][k];
#ifdef VUNROLL
		    register int        p4;
		    
		    p4 = p / 4;
		    p  = p % 4;
		    for ( l = 0; l < p4; l++ )
		      {
			dp1[4*l]   += s*dp2[4*l];
			dp1[4*l+1] += s*dp2[4*l+1];
			dp1[4*l+2] += s*dp2[4*l+2];
			dp1[4*l+3] += s*dp2[4*l+3];
		      }
		    dp1 += 4*p4;	dp2 += 4*p4;
#endif
		    
		    for ( l = 0; l < p; l++ )
		      dp1[l] += s*dp2[l];
		  }

		}

	return OUT;
}

/* mtrm_mlt -- matrix transposed-matrix multiplication
	-- A^T.B is returned, result stored in OUT */
MAT *my_mtrm_mlt(MAT *A,MAT *B,MAT *OUT)
{
	int	i, k, l, limit;

	if ( ! OUT || OUT->m != A->n || OUT->n != B->n )
		OUT = m_resize(OUT,A->n,B->n);

	limit = B->n;
	m_zero(OUT);
	for ( k = 0; k < A->m; k++ )
		for ( i = 0; i < A->n; i++ )
		{
		  if ( A->me[k][i] != 0.0 ) {
		    //__mltadd__(OUT->me[i],B->me[k],A->me[k][i],(int)limit);


		    register double *dp1, *dp2, s;
		    register int	l;
		    dp1 = OUT->me[i];
		    dp2 = B->me[k];
		    s = A->me[k][i];
#ifdef VUNROLL
		    register int        p4;
		    
		    p4 = limit / 4;
		    limit  = limit % 4;
		    for ( l = 0; l < p4; l++ )
		      {
			dp1[4*l]   += s*dp2[4*l];
			dp1[4*l+1] += s*dp2[4*l+1];
			dp1[4*l+2] += s*dp2[4*l+2];
			dp1[4*l+3] += s*dp2[4*l+3];
		      }
		    dp1 += 4*p4;	dp2 += 4*p4;
#endif
		    
		    for ( l = 0; l < limit; l++ )
		      dp1[l] += s*dp2[l];
		  }
		}

	return OUT;
}

/*
========================================================================

 R = mtvm_mlt(A,b,R);

 Computes A'*diag(b)*A, where b is a vector.

 If R is not allocated to be of an appropriate size, it is resized
 automatically.

========================================================================
*/
MAT *mtvm_mlt0(MAT *A, VEC *b, MAT *R)
{
  int i,j;
  MAT *temp=MNULL;

  temp = m_get(A->m,A->n);
  for(i=0;i<temp->m;i++)
    for(j=0;j<temp->n;j++)
      temp->me[i][j] = b->ve[i]*A->me[i][j];
  R = my_mtrm_mlt(A,temp,R);
  M_FREE(temp);
  return R;
}

// Computes A'*diag(b)*A, where b is a vector.
MAT *mtvm_mlt1(MAT *A, VEC *b, MAT *R)
{
  double v;
  int i,j,k;
  
  if(R==NULL)
    R = m_get(A->n,A->n);

  for(i=0;i<A->n;i++) {
    for(j=0;j<A->n;j++) {
      v = 0.0;
      for(k=0;k<A->m;k++) {
	v += b->ve[k]*A->me[k][i]*A->me[k][j];
      }
      R->me[i][j] = v;
    }
  }
  return R;
}

// compute R = A'*diag(R)*A
// result is a symmetric matrix.
MAT *mtvm_mlt(MAT *A, VEC *b, MAT *R)
{
  double v;
  int i,j,k;
  register double *Ak,*Aki,bk;

  if(R==NULL)
    R = m_get(A->n,A->n);
  else
    m_zero(R);

  for(k=0;k<A->m;k++) {
    Ak = A->me[k];
    Aki = Ak;
    bk = b->ve[k];
    for(i=0;i<A->n;i++) {
      /*
      for(j=0;j<A->n;j++) {
	R->me[i][j] += bk*(*Aki)*Ak[j];
      }
      */
      for(j=0;j<i;j++) {
	R->me[i][j] += bk*(*Aki)*Ak[j];
	R->me[j][i] += bk*(*Aki)*Ak[j];
      }
      R->me[i][i] += bk*(*Aki)*(*Aki);
      Aki++;
    }
  }
  return R;
}


// permute rows of A, storing them in OUT.
void permrows0(MAT *A,int *p,MAT *OUT)
{
  u_int i,j;

  for(i=0;i<A->m;i++)
    for(j=0;j<A->n;j++)
      OUT->me[i][j] = A->me[p[i]][j];
}

// permute rows of A, storing them in OUT.
void permrows(MAT *A,int *p,MAT *OUT)
{
  register double *op, *op2;
  register int i,j;

  for(i=0;i<A->m;i++) {
    op = OUT->me[i];
    op2 = A->me[p[i]];
    for(j=0;j<A->n;j++)
      *op++ = *op2++;
  }
}

// permute cols of A, storing them in OUT.
void permcols_old(MAT *A,int *p,MAT *OUT)
{
  u_int i,j;

  for(i=0;i<A->m;i++)
    for(j=0;j<A->n;j++)
      OUT->me[i][p[j]] = A->me[i][j];
}

// permute cols of A, storing them in OUT.
void permcols(MAT *A,int *p,MAT *OUT)
{
  register double *op,*op2;
  register u_int i,j;
  
  op2 = A->me[0];
  for(i=0;i<A->m;i++) {
    op = OUT->me[i];
    for(j=0;j<A->n;j++) {
      op[p[j]] = *(op2++);
    }
  }
}

/***********************************************************************/
/*
========================================================================

 prod = rowdot(A,ra,B,rb);

 Takes dot product of A(ra,:) with B(rb,:).

========================================================================
*/
inline double rowdot(MAT *A,int ra,MAT *B,int rb)
{
  int i;
  register double *a = A->me[ra];
  register double *b = B->me[rb];
  register double sum = 0.0;
  for(i=0;i<A->n;i++)
    sum += (*a++) * (*b++);
  return sum;
}

/*
======================================================================

 my_svd(A,d);

 Given a square, symmetric matrix A, calculate W and d so that

 A = W*diag(d)*W', 

 and all elements of d are nonnegative.

 On return from this function, the first argument will hold W,
 i.e. the contents of A are overwritten.

======================================================================
*/
void my_svd(MAT *A, VEC *d)
{
  MAT *U,*A0;
  int i,j;
  U = m_get(A->m,A->n);

  /*
  A0 = m_get(1,A->n);

  for(i=0;i<A->n;i++)
    A0->me[0][i] = A->me[0][i];
  */

  A0 = m_copy(A,MNULL);

  /* 
     call to svd in svd.C will destroy A. 
     normal Meschach svd will not destroy A.
  */

  svd(A,MNULL,U,d); 

  /*
  for(i=0;i<U->m;i++)
  {
    if(fabs(rowdot(A0,0,U,i) + d->ve[i]*U->me[i][0]) < 1e-6)
    {
      d->ve[i] = -d->ve[i];
      for(j=0;j<U->n;j++)
      {
	A->me[i][j] = -U->me[i][j]; 
      }
    }
    else {
      for(j=0;j<U->n;j++)
	A->me[i][j] = U->me[i][j];
    }
  }
  */

  for(i=0;i<U->m;i++)
  {
    if(fabs(rowdot(A0,i,U,i) + d->ve[i]*U->me[i][i]) < 1e-6)
    {
      d->ve[i] = -d->ve[i];
      for(j=0;j<U->n;j++)
      {
	A->me[i][j] = -U->me[i][j]; 
      }
    }
    else {
      for(j=0;j<U->n;j++)
	A->me[i][j] = U->me[i][j];
    }
  }
  M_FREE(A0);
  M_FREE(U);
}
