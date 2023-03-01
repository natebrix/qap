#include <algorithm>
#include	<stdio.h>
#include    "meschach.h"
#include	<math.h>

#define	sgn(x)	((x) >= 0 ? 1 : -1)
#define	MAX_STACK	30

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
inline void __mltadd__(Real *dp1,Real *dp2,double s,int len)
{
    register int	i;
#ifdef VUNROLL
    register int        len4;
    
    len4 = len / 4;
    len  = len % 4;
    for ( i = 0; i < len4; i++ )
    {
	dp1[4*i]   += s*dp2[4*i];
	dp1[4*i+1] += s*dp2[4*i+1];
	dp1[4*i+2] += s*dp2[4*i+2];
	dp1[4*i+3] += s*dp2[4*i+3];
    }
    dp1 += 4*len4;	dp2 += 4*len4;
#endif
    
    for ( i = 0; i < len; i++ )
	dp1[i] += s*dp2[i];
}


/* hhvec -- calulates Householder vector to eliminate all entries after the
	i0 entry of the vector vec. It is returned as out. May be in-situ */
void hhvec(VEC *out,u_int i0,Real *beta,Real *newval)
{
	Real	norm;
	Real *outi;
	norm = 0.0;
	outi = out->ve;
	//	for(int i=0;i<i0;i++) {
	for(int i=i0;i<out->dim;i++) {
	  norm += outi[i]*outi[i];
	}
	norm = sqrt(norm);
	/*
	*/
	//	norm = sqrt(_in_prod(out,out,i0));
	//	if ( norm <= 0.0 )
	if ( norm <= MACHEPS*MACHEPS ) // changed, 7/17/2000
	{
		*beta = 0.0;
		return;
	}
	*beta = 1.0/(norm * (norm+fabs(out->ve[i0])));
	if ( out->ve[i0] > 0.0 )
		*newval = -norm;
	else
		*newval = norm;
	out->ve[i0] -= *newval;
}

/* hhtrcols2 -- transform a matrix by a Householder vector by columns
	starting at row i0 from column j0 -- in-situ */
MAT	*hhtrcols2(MAT *M,u_int i0,u_int j0,VEC *hh,double beta)
{
  /* Real	ip, scale; */
  int	i /*, k */;
  static	VEC	*w = VNULL;
  
  if ( M==(MAT *)NULL || hh==(VEC *)NULL )
    error(E_NULL,"hhtrcols2");
  if ( M->m != hh->dim )
    error(E_SIZES,"hhtrcols2");
  if ( i0 > M->m || j0 > M->n )
    error(E_BOUNDS,"hhtrcols2");
  
  if ( beta == 0.0 )	return (M);
  
  w = v_resize(w,M->n);
  MEM_STAT_REG(w,TYPE_VEC);
  v_zero(w);
  
  for ( i = i0; i < M->m; i++ )
    if ( hh->ve[i] != 0.0 )
      __mltadd__(&(w->ve[j0]),&(M->me[i][j0]),hh->ve[i],
		 (int)(M->n-j0));
  for ( i = i0; i < M->m; i++ )
    if ( hh->ve[i] != 0.0 )
      __mltadd__(&(M->me[i][j0]),&(w->ve[j0]),-beta*hh->ve[i],
		 (int)(M->n-j0));
  return (M);
}

/* hhtrrows2 -- transform a matrix by a Householder vector by rows
	starting at row i0 from column j0 -- in-situ */
MAT	*hhtrrows2(MAT *M,u_int i0,u_int j0,VEC *hh,double beta)
{
  Real	ip, scale;
  int	i , j ;
  
  if ( M==(MAT *)NULL || hh==(VEC *)NULL )
    error(E_NULL,"hhtrrows2");
  if ( M->n != hh->dim )
    error(E_RANGE,"hhtrrows2");
  if ( i0 > M->m || j0 > M->n )
    error(E_BOUNDS,"hhtrrows2");
  
  if ( beta == 0.0 )	return (M);
  
  /* for each row ... */
  for ( i = i0; i < M->m; i++ )
    {	/* compute inner product */
      //ip = __ip__(&(M->me[i][j0]),&(hh->ve[j0]),(int)(M->n-j0));
      ip = 0.0;
      for ( j = j0; j < M->n; j++ )
	ip += M->me[i][j]*hh->ve[j];
      /**************************************************
       **************************************************/
      scale = beta*ip;
      if ( scale == 0.0 )
	continue;
      
      /* do operation */
      __mltadd__(&(M->me[i][j0]),&(hh->ve[j0]),-scale,
		 (int)(M->n-j0));
      /**************************************************
	for ( j = j0; j < M->n; j++ )
	  M->me[i][j] -= scale*hh->ve[j];
      **************************************************/
    }
  
  return (M);
}


void rot_rows(MAT *out,u_int i,u_int k,double c,double s)
{
  u_int	j;
  Real	temp;
  double outij, outkj;

  /*
  if ( mat==(MAT *)NULL )
    error(E_NULL,"rot_rows");
  if ( i >= mat->m || k >= mat->m )
    error(E_RANGE,"rot_rows");
  */

  
  for ( j=0; j<out->n; j++ )
    {
      outij = out->me[i][j];
      outkj = out->me[k][j];

      // temp = c*out->me[i][j] + s*out->me[k][j]; 
      // out->me[k][j] = -s*out->me[i][j] + c*out->me[k][j];

      temp = c*outij + s*outkj; 
      out->me[k][j] = -s*outij + c*outkj;

      out->me[i][j] = temp;
    }
}

inline void givens(double x,double y,Real &c,Real &s)
{
	Real	norm;

	norm = sqrt(x*x+y*y);
	if ( norm != 0.0 )
	{  c = x/norm;	s = y/norm;	}
	else
	{  c = 1.0;	s = 0.0;	}	/* identity */
}

/* fixsvd -- fix minor details about SVD
	-- make singular values non-negative
	-- sort singular values in decreasing order
	-- variables as for bisvd()
	-- no argument checking */
static void	fixsvd(VEC *d,MAT *U,MAT *V)
{
  int i,j;

  /* make singular values non-negative */
  for ( i = 0; i < d->dim; i++ )
    if ( d->ve[i] < 0.0 )
      {
	d->ve[i] = - d->ve[i];
	if ( U != MNULL )
	  for ( j = 0; j < U->m; j++ )
	    U->me[i][j] = - U->me[i][j];
      }

#ifdef SORT_SING_VALS
    int sp,stack[MAX_STACK];
    int		k, l, r;
    Real	tmp, v;
    
    /* sort singular values */
    /* nonrecursive implementation of quicksort due to R.Sedgewick,
       "Algorithms in C", p. 122 (1990) */
    sp = -1;
    l = 0;	r = d->dim - 1;
    for ( ; ; )
    {
	while ( r > l )
	{
	    /* i = partition(d->ve,l,r) */
	    v = d->ve[r];

	    i = l - 1;	    j = r;
	    for ( ; ; )
	    {	/* inequalities are "backwards" for **decreasing** order */
		while ( d->ve[++i] > v )
		    ;
		while ( d->ve[--j] < v )
		    ;
		if ( i >= j )
		    break;
		/* swap entries in d->ve */
		tmp = d->ve[i];	  d->ve[i] = d->ve[j];	d->ve[j] = tmp;
		/* swap rows of U & V as well */
		if ( U != MNULL )
		    for ( k = 0; k < U->n; k++ )
		    {
			tmp = U->me[i][k];
			U->me[i][k] = U->me[j][k];
			U->me[j][k] = tmp;
		    }
		if ( V != MNULL )
		    for ( k = 0; k < V->n; k++ )
		    {
			tmp = V->me[i][k];
			V->me[i][k] = V->me[j][k];
			V->me[j][k] = tmp;
		    }
	    }
	    tmp = d->ve[i];    d->ve[i] = d->ve[r];    d->ve[r] = tmp;
	    if ( U != MNULL )
		for ( k = 0; k < U->n; k++ )
		{
		    tmp = U->me[i][k];
		    U->me[i][k] = U->me[r][k];
		    U->me[r][k] = tmp;
		}
	    if ( V != MNULL )
		for ( k = 0; k < V->n; k++ )
		{
		    tmp = V->me[i][k];
		    V->me[i][k] = V->me[r][k];
		    V->me[r][k] = tmp;
		}
	    /* end i = partition(...) */
	    if ( i - l > r - i )
	    {	stack[++sp] = l;    stack[++sp] = i-1;	l = i+1;    }
	    else
	    {	stack[++sp] = i+1;  stack[++sp] = r;	r = i-1;    }
	}
	if ( sp < 0 )
	    break;
	r = stack[sp--];	l = stack[sp--];
    }
#endif
}


/* bisvd -- svd of a bidiagonal m x n matrix represented by d (diagonal) and
			f (super-diagonals)
	-- returns with d set to the singular values, f zeroed
	-- if U, V non-NULL, the orthogonal operations are accumulated
		in U, V; if U, V == I on entry, then SVD == U^T.A.V
		where A is initial matrix
	-- returns d on exit */
VEC	*bisvd(VEC *d,VEC *f,MAT *U,MAT *V)
{
	int	i, j, n;
	int	i_min, i_max, split;
	Real	c, s, shift, size, z;
	Real	d_tmp, diff, t11, t12, t22, *d_ve, *f_ve;

	if ( ! d || ! f )
		error(E_NULL,"bisvd");
	if ( d->dim != f->dim + 1 )
		error(E_SIZES,"bisvd");
	n = d->dim;
	if ( ( U && U->n < n ) || ( V && V->m < n ) )
		error(E_SIZES,"bisvd");
	if ( ( U && U->m != U->n ) || ( V && V->m != V->n ) )
		error(E_SQUARE,"bisvd");


	if ( n == 1 )
	  return d;
	d_ve = d->ve;	f_ve = f->ve;

	size = v_norm_inf(d) + v_norm_inf(f);

	i_min = 0;
	while ( i_min < n )	/* outer while loop */
	{
	    /* find i_max to suit;
		submatrix i_min..i_max should be irreducible */
	    i_max = n - 1;
	    for ( i = i_min; i < n - 1; i++ )
		if ( d_ve[i] == 0.0 || f_ve[i] == 0.0 )
		{   i_max = i;
		    if ( f_ve[i] != 0.0 )
		    {
			/* have to ``chase'' f[i] element out of matrix */
			z = f_ve[i];	f_ve[i] = 0.0;
			for ( j = i; j < n-1 && z != 0.0; j++ )
			{
			  // givens(d_ve[j+1],z, &c, &s);
			  givens(d_ve[j+1],z, c, s);
			    s = -s;
			    d_ve[j+1] =  c*d_ve[j+1] - s*z;
			    if ( j+1 < n-1 )
			    {
				z         = s*f_ve[j+1];
				f_ve[j+1] = c*f_ve[j+1];
			    }
			    if ( U )
			      //rot_rows(U,i,j+1,c,s,U);
			      rot_rows(U,i,j+1,c,s);
			}
		    }
		    break;
		}
	    if ( i_max <= i_min )
	    {
		i_min = i_max + 1;
		continue;
	    }
	    /* printf("bisvd: i_min = %d, i_max = %d\n",i_min,i_max); */

	    split = FALSE;
	    while ( ! split )
	    {
		/* compute shift */
		t11 = d_ve[i_max-1]*d_ve[i_max-1] +
			(i_max > i_min+1 ? f_ve[i_max-2]*f_ve[i_max-2] : 0.0);
		t12 = d_ve[i_max-1]*f_ve[i_max-1];
		t22 = d_ve[i_max]*d_ve[i_max] + f_ve[i_max-1]*f_ve[i_max-1];
		/* use e-val of [[t11,t12],[t12,t22]] matrix
				closest to t22 */
		diff = (t11-t22)/2;
		shift = t22 - t12*t12/(diff +
			sgn(diff)*sqrt(diff*diff+t12*t12));

		/* initial Givens' rotation */
		//givens(d_ve[i_min]*d_ve[i_min]-shift,
		//d_ve[i_min]*f_ve[i_min], &c, &s);
		givens(d_ve[i_min]*d_ve[i_min]-shift,
			d_ve[i_min]*f_ve[i_min], c, s);

		/* do initial Givens' rotations */
		d_tmp         = c*d_ve[i_min] + s*f_ve[i_min];
		f_ve[i_min]   = c*f_ve[i_min] - s*d_ve[i_min];
		d_ve[i_min]   = d_tmp;
		z             = s*d_ve[i_min+1];
		d_ve[i_min+1] = c*d_ve[i_min+1];
		if ( V )
		  //rot_rows(V,i_min,i_min+1,c,s,V);
		    rot_rows(V,i_min,i_min+1,c,s);
		/* 2nd Givens' rotation */
		//givens(d_ve[i_min],z, &c, &s);
		givens(d_ve[i_min],z, c, s);
		d_ve[i_min]   = c*d_ve[i_min] + s*z;
		d_tmp         = c*d_ve[i_min+1] - s*f_ve[i_min];
		f_ve[i_min]   = s*d_ve[i_min+1] + c*f_ve[i_min];
		d_ve[i_min+1] = d_tmp;
		if ( i_min+1 < i_max )
		{
		    z             = s*f_ve[i_min+1];
		    f_ve[i_min+1] = c*f_ve[i_min+1];
		}
		if ( U )
		  //rot_rows(U,i_min,i_min+1,c,s,U);
		    rot_rows(U,i_min,i_min+1,c,s);

		for ( i = i_min+1; i < i_max; i++ )
		{
		    /* get Givens' rotation for zeroing z */
		  //givens(f_ve[i-1],z, &c, &s);
		    givens(f_ve[i-1],z, c, s);
		    f_ve[i-1] = c*f_ve[i-1] + s*z;
		    d_tmp     = c*d_ve[i] + s*f_ve[i];
		    f_ve[i]   = c*f_ve[i] - s*d_ve[i];
		    d_ve[i]   = d_tmp;
		    z         = s*d_ve[i+1];
		    d_ve[i+1] = c*d_ve[i+1];
		    if ( V )
		      //rot_rows(V,i,i+1,c,s,V);
			rot_rows(V,i,i+1,c,s);
		    /* get 2nd Givens' rotation */
		    //givens(d_ve[i],z, &c, &s);
		    givens(d_ve[i],z, c, s);
		    d_ve[i]   = c*d_ve[i] + s*z;
		    d_tmp     = c*d_ve[i+1] - s*f_ve[i];
		    f_ve[i]   = c*f_ve[i] + s*d_ve[i+1];
		    d_ve[i+1] = d_tmp;
		    if ( i+1 < i_max )
		    {
			z         = s*f_ve[i+1];
			f_ve[i+1] = c*f_ve[i+1];
		    }
		    if ( U )
		      //rot_rows(U,i,i+1,c,s,U);
			rot_rows(U,i,i+1,c,s);
		}
		/* should matrix be split? */
		for ( i = i_min; i < i_max; i++ )
		    if ( fabs(f_ve[i]) <
				MACHEPS*(fabs(d_ve[i])+fabs(d_ve[i+1])) )
		    {
			split = TRUE;
			f_ve[i] = 0.0;
		    }
		    else if ( fabs(d_ve[i]) < MACHEPS*size )
		    {
			split = TRUE;
			d_ve[i] = 0.0;
		    }
		    /* printf("bisvd: d =\n");	v_output(d); */
		    /* printf("bisvd: f = \n");	v_output(f); */
		}
	}
	fixsvd(d,U,V);

	return d;
}

/* bifactor -- perform preliminary factorisation for bisvd
	-- updates U and/or V, which ever is not NULL */
MAT	*bifactor(MAT *A,MAT *U,MAT *V)
{
	int	k;
	static VEC	*tmp1=VNULL, *tmp2=VNULL;
	Real	beta;

	if ( ! A )
		error(E_NULL,"bifactor");
	if ( ( U && ( U->m != U->n ) ) || ( V && ( V->m != V->n ) ) )
		error(E_SQUARE,"bifactor");
	if ( ( U && U->m != A->m ) || ( V && V->m != A->n ) )
		error(E_SIZES,"bifactor");
	tmp1 = v_resize(tmp1,A->m);
	tmp2 = v_resize(tmp2,A->n);
	MEM_STAT_REG(tmp1,TYPE_VEC);
	MEM_STAT_REG(tmp2,TYPE_VEC);

	if ( A->m >= A->n )
	    for ( k = 0; k < A->n; k++ )
	    {
		get_col(A,k,tmp1);
		//hhvec(tmp1,k,&beta,tmp1,&(A->me[k][k]));
		hhvec(tmp1,k,&beta,&(A->me[k][k]));
		hhtrcols2(A,k,k+1,tmp1,beta);
		if ( U )
		    hhtrcols2(U,k,0,tmp1,beta);
		if ( k+1 >= A->n )
		    continue;
		get_row(A,k,tmp2);
		//hhvec(tmp2,k+1,&beta,tmp2,&(A->me[k][k+1]));
		hhvec(tmp2,k+1,&beta,&(A->me[k][k+1]));
		hhtrrows2(A,k+1,k+1,tmp2,beta);
		if ( V )
		    hhtrcols2(V,k+1,0,tmp2,beta);
	    }
	else
	    for ( k = 0; k < A->m; k++ )
	    {
		get_row(A,k,tmp2);
		//hhvec(tmp2,k,&beta,tmp2,&(A->me[k][k]));
		hhvec(tmp2,k,&beta,&(A->me[k][k]));
		hhtrrows2(A,k+1,k,tmp2,beta);
		if ( V )
		    hhtrcols2(V,k,0,tmp2,beta);
		if ( k+1 >= A->m )
		    continue;
		get_col(A,k,tmp1);
		//hhvec(tmp1,k+1,&beta,tmp1,&(A->me[k+1][k]));
		hhvec(tmp1,k+1,&beta,&(A->me[k+1][k]));
		hhtrcols2(A,k+1,k+1,tmp1,beta);
		if ( U )
		    hhtrcols2(U,k+1,0,tmp1,beta);
	    }

	return A;
}

/* svd -- returns vector of singular values in d
	-- also updates U and/or V, if one or the other is non-NULL
	-- destroys A */
VEC	*svd(MAT *A_tmp,MAT *U,MAT *V,VEC *d)
{
  static VEC	*f=VNULL;
  int	i, limit;
  
  if ( U != MNULL )
    m_ident(U);
  if ( V != MNULL )
    m_ident(V);
  limit = std::min(A_tmp->m,A_tmp->n);
  d = v_resize(d,limit);
  f = v_resize(f,limit-1);
  MEM_STAT_REG(f,TYPE_VEC);
  
  bifactor(A_tmp,U,V);
  if ( A_tmp->m >= A_tmp->n )
    for ( i = 0; i < limit; i++ )
      {
	d->ve[i] = A_tmp->me[i][i];
	if ( i+1 < limit )
	  f->ve[i] = A_tmp->me[i][i+1];
      }
  else
    for ( i = 0; i < limit; i++ )
      {
	d->ve[i] = A_tmp->me[i][i];
	if ( i+1 < limit )
	  f->ve[i] = A_tmp->me[i+1][i];
      }
  
  
  if ( A_tmp->m >= A_tmp->n )
    bisvd(d,f,U,V);
  else
    bisvd(d,f,V,U);
  
  return d;
}


VEC *svd0(MAT *A,MAT *U,MAT *V,VEC *d)
{
	static VEC	*f=VNULL;
	int	i, limit;
	MAT	*A_tmp;

	if ( ! A )
		error(E_NULL,"svd");
	if ( ( U && ( U->m != U->n ) ) || ( V && ( V->m != V->n ) ) )
		error(E_SQUARE,"svd");
	if ( ( U && U->m != A->m ) || ( V && V->m != A->n ) )
		error(E_SIZES,"svd");

	A_tmp = m_copy(A,MNULL);
	if ( U != MNULL )
	    m_ident(U);
	if ( V != MNULL )
	    m_ident(V);
	limit = std::min(A_tmp->m,A_tmp->n);
	d = v_resize(d,limit);
	f = v_resize(f,limit-1);
	MEM_STAT_REG(f,TYPE_VEC);

	bifactor(A_tmp,U,V);
	if ( A_tmp->m >= A_tmp->n )
	    for ( i = 0; i < limit; i++ )
	    {
		d->ve[i] = A_tmp->me[i][i];
		if ( i+1 < limit )
		    f->ve[i] = A_tmp->me[i][i+1];
	    }
	else
	    for ( i = 0; i < limit; i++ )
	    {
		d->ve[i] = A_tmp->me[i][i];
		if ( i+1 < limit )
		    f->ve[i] = A_tmp->me[i+1][i];
	    }


	if ( A_tmp->m >= A_tmp->n )
	    bisvd(d,f,U,V);
	else
	    bisvd(d,f,V,U);

	M_FREE(A_tmp);

	return d;
}

