/*
============================================================

 Linear assignment solver of Jonker and Volgenant.

 --I have written a driver for the code, in lap_jv().
 --A floating-point version is in dlap.[hC]

============================================================
*/

#ifndef LAP2_H
#define LAP2_H

/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/*************** CONSTANTS  *******************/

//#define LAP_USES_INT

//#define BIG 100000
// for tai problems, big needs to be:
//#define BIG 1000000
//temp
#define BIG 100000000

/*************** TYPES      *******************/

typedef int row;
typedef int col;
typedef int cost;

/*************** FUNCTIONS  *******************/

extern cost lap(int dim, cost **assigncost,
               int *rowsol, int *colsol, cost *u, cost *v);

extern cost lap_noinit(int dim, 
		       cost **assigncost,
		       col *rowsol, 
		       row *colsol, 
		       cost *u, 
		       cost *v);

extern void checklap(int dim, int **assigncost,
                     int *rowsol, int *colsol, int *u, int *v);

/* driver by NB */

double lap_jv(double **D,int n,double mult);
double lap_jv(double **D,int *perm,int n,double mult);
double lap_jv(double **D,int **D_int,int *perm,int *invperm,
	      double *ys,int *ys_int,double *yt,int *yt_int,int n,double mult);
double lap_jv(double **D,int *perm,
	      double *ys,double *yt,int n,double mult);

//void reducedCosts(double **D, double **U, double *ys,double *yt,int n)
inline void reducedCosts(double **D, double **U, double *ys,double *yt,int n)
{
  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++) {
      U[i][j] = D[i][j] - (yt[j] + ys[i]);
    }
  }
}
double dual_obj(double *ys,double *yt,int n);

#endif

