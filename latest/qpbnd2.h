#ifndef QPBND2_H
#define QPBND2_H

#include "qap.h"
#include "fw.h"
#include "qpbnd.h"

/*
  computes the same lower bound as in qpbnd().  The only difference
  is that some of the computation required to calculate S and T
  has been done already.  This information is passed into qpbnd2()
  to save time.

  When we compute B_ij, we reuse the same A' and B' matrices for
  several calls to qpbnd.  Basically, we are trying to avoid 
  duplicate work.

  [fill this out later]
 */
void init_w0(MAT *A,VEC *&w,MAT *&W,bool ascending);
void init_w(MAT *A,VEC *&w,MAT *&W,bool ascending);

double qpbnd2(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
	      VEC *w,VEC *u,MAT *W,MAT *U,QPBParameters *qpbp);

double qpbnd_parametric2(MAT *A,MAT *B,MAT *C,MAT *X,MAT *U1,
			 VEC *w,VEC *u,MAT *W,MAT *U,
			 QPBParameters *qpbp);

#endif


