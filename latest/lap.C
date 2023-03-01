/************************************************************************
*
*  lap.cpp
   version 1.0 - 4 September 1996
   author: Roy Jonker @ MagicLogic Optimization Inc.
   e-mail: roy_jonker@magiclogic.com

   Code for Linear Assignment Problem, according to 
   
   "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear   
    Assignment Problems," Computing 38, 325-340, 1987
   
   by
   
   R. Jonker and A. Volgenant, University of Amsterdam.
*
*************************************************************************/
/*************** CONSTANTS *******************/

#if !defined TRUE
#define	 TRUE		1
#endif
#if !defined FALSE
#define  FALSE		0
#endif

/*************** DATA TYPES *******************/

typedef int boolean;
#include "lap.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define UNROLL_4
cost lap(int dim, 
        cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v)

// input:
// dim        - problem size
// assigncost - cost matrix

// output:
// rowsol     - column assigned to row in solution
// colsol     - row assigned to column in solution
// u          - dual variables, row reduction numbers
// v          - dual variables, column reduction numbers

{
  int *ap;
  //  col *mp;
  row *cp;
  boolean unassignedfound;
  register row i;
  register col j;
  row imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
  col j1, j2, endofpath, last, low, up, *collist, *matches;
  cost min, h, umin, usubmin, v2, *d;

  free = new row[dim];       // list of unassigned rows.
  collist = new col[dim];    // list of columns to be scanned in various ways.
  matches = new col[dim];    // counts how many times a row could be assigned.
  d = new cost[dim];         // 'cost-distance' in augmenting path calculation.
  pred = new row[dim];       // row-predecessor of column in augmenting/alternating path.

  // init how many times a row will be assigned in the column reduction.
  for (i = 0; i < dim; i++)  
    matches[i] = 0;

#ifdef UNROLL_1
  // COLUMN REDUCTION 
  //cp = &colsol[dim-1];
  cp = colsol + dim - 1;
  for (j = dim-1; j >= 0; j--)    // reverse order gives better results.
  {
    // find minimum cost over rows.
    ap = &assigncost[0][j];
    min = *ap;
    //    min = assigncost[0][j]; 
    imin = 0;
    for (i = 1; i < dim; i++) {
      //      ap += dim;
      ap = &assigncost[i][j];
      if(*ap < min) {
	min = *ap;
	imin = i;
	} 
    }
    v[j] = min; 

    if (++matches[imin] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin] = j; 
      //      colsol[j] = imin; 
      *cp = imin; 
    }
    else {
      *cp = -1;   // row already assigned, column not assigned.
    }
    cp--;
  }
#endif
#ifdef UNROLL_2
  cost min2;  row imin2; int *ap2;
  // COLUMN REDUCTION 
  cp = colsol + dim - 1;
  for (j = dim-1; j >= 1; j-=2)
  {
    // find minimum cost over rows.
    ap = &assigncost[0][j];    ap2 = &assigncost[0][j-1];
    min = *ap;                 min2 = *ap2;
    imin = 0;                  imin2 = 0;
    for (i = 1; i < dim; i++) {
      ap += dim;
      if(*ap < min) {
	min = *ap;
	imin = i;
	} 
      ap2 += dim;
      if(*ap2 < min2) {
	min2 = *ap2;
	imin2 = i;
	} 
    }
    v[j  ] = min; 
    v[j-1] = min2; 

    if (++matches[imin] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin] = j; 
      *cp = imin; 
    }
    else {
      *cp = -1;   // row already assigned, column not assigned.
    }
    cp--;
    if (++matches[imin2] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin2] = j-1; 
      *cp = imin2; 
    }
    else {
      *cp = -1;   // row already assigned, column not assigned.
    }
    cp--;
  }

  //  for (; j >= 0; j--,cp--)    // reverse order gives better results.
  if(j==0)
  {
    // find minimum cost over rows.
    ap = &assigncost[0][j];
    min = *ap;
    //    min = assigncost[0][j]; 
    imin = 0;
    for (i = 1; i < dim; i++) {
      ap += dim;
      if(*ap < min) {
	min = *ap;
	imin = i;
	} 
    }
    v[j] = min; 

    if (++matches[imin] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin] = j; 
      *cp = imin; 
    }
    else
      *cp = -1;   // row already assigned, column not assigned.
  }
#endif
#ifdef UNROLL_4
  cost min2;  row imin2; int *ap2;
  cost min3;  row imin3; int *ap3;
  cost min4;  row imin4; int *ap4;
  // COLUMN REDUCTION 
  cp = colsol + dim - 1;
  for (j = dim-1; j >= 3; j-=4)
  {
    // find minimum cost over rows.
    ap  = &assigncost[0][j];   ap2 = &assigncost[0][j-1];
    ap3 = &assigncost[0][j-2]; ap4 = &assigncost[0][j-3];
    min  = *ap;                min2 = *ap2;
    min3 = *ap3;               min4 = *ap4;
    imin  = 0;                 imin2 = 0;
    imin3 = 0;                 imin4 = 0;
    for (i = 1; i < dim; i++) {
      ap = &assigncost[i][j];
      //ap2 = ap - 1;
      //ap3 = ap2 - 1;
      //ap4 = ap3 - 1;
      //      ap += dim;
      if(*ap < min) {
	min = *ap;
	imin = i;
	} 
      //      ap2 = &assigncost[i][j-1];
      ap2 = ap - 1;
      if(*ap2 < min2) {
	min2 = *ap2;
	imin2 = i;
	} 
      //      ap3 = &assigncost[i][j-2];
      ap3 = ap - 2;
      if(*ap3 < min3) {
	min3 = *ap3;
	imin3 = i;
	} 
      //      ap4 = &assigncost[i][j-3];
      ap4 = ap - 3;
      if(*ap4 < min4) {
	min4 = *ap4;
	imin4 = i;
	} 
    }
    v[j  ] = min; 
    v[j-1] = min2; 
    v[j-2] = min3; 
    v[j-3] = min4; 

    if (++matches[imin] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin] = j; 
      *cp = imin; 
    }
    else
      *cp = -1;   // row already assigned, column not assigned.
    cp--;
    if (++matches[imin2] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin2] = j-1; 
      *cp = imin2; 
    }
    else
      *cp = -1;   // row already assigned, column not assigned.
    cp--;
    if (++matches[imin3] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin3] = j-2; 
      *cp = imin3; 
    }
    else
      *cp = -1;   // row already assigned, column not assigned.
    cp--;
    if (++matches[imin4] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin4] = j-3; 
      *cp = imin4; 
    }
    else
      *cp = -1;   // row already assigned, column not assigned.
    cp--;
  }

  for (; j >= 0; j--)    // reverse order gives better results.
  {
    // find minimum cost over rows.
    ap = &assigncost[0][j];
    min = *ap;
    //    min = assigncost[0][j]; 
    imin = 0;
    for (i = 1; i < dim; i++) {
      ap = &assigncost[i][j];
      //ap += dim;
      if(*ap < min) {
	min = *ap;
	imin = i;
	} 
    }
    v[j] = min; 

    if (++matches[imin] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin] = j; 
      *cp = imin; 
    }
    else
      *cp = -1;   // row already assigned, column not assigned.
    cp--;
  }
#endif

  // REDUCTION TRANSFER
  for (i = 0; i < dim; i++) 
    if (matches[i] == 0)     // fill list of unassigned 'free' rows.
      free[numfree++] = i;
    else
      if (matches[i] == 1)   // transfer reduction from rows that are assigned once.
      {
        j1 = rowsol[i]; 
        min = BIG;
        for (j = 0; j < dim; j++)  
          if (j != j1)
            if (assigncost[i][j] - v[j] < min) 
              min = assigncost[i][j] - v[j];
        v[j1] = v[j1] - min;
      }

  // AUGMENTING ROW REDUCTION 
  //  int lapeps = 0;
  int loopcnt = 0;           // do-loop to be done twice.
  do
  {
    loopcnt++;

    // scan all free rows.
    // in some cases, a free row may be replaced with another one to be scanned next.
    k = 0; 
    prvnumfree = numfree; 
    numfree = 0;             // start list of rows still free after augmenting row reduction.
    while (k < prvnumfree)
    {
      i = free[k]; 
      k++;

      // find minimum and second minimum reduced cost over columns.
      umin = assigncost[i][0] - v[0]; 
      j1 = 0; 
      usubmin = BIG;

      for (j = 1; j < dim; j++) 
      {
	h = assigncost[i][j] - v[j];
	if (h < usubmin)
          if (h >= umin) 
          { 
            usubmin = h; 
            j2 = j;
          }
          else 
          { 
            usubmin = umin; 
            umin = h; 
            j2 = j1; 
            j1 = j;
          }
      }

      i0 = colsol[j1];
      if (umin < usubmin) 
      //      if (umin < usubmin - lapeps) 
        // change the reduction of the minimum column to increase the minimum
        // reduced cost in the row to the subminimum.
        v[j1] = v[j1] - (usubmin - umin);
      else                   // minimum and subminimum equal.
	//if (umin >= usubmin) // minimum and subminimum equal.
        if (i0 >= 0)         // minimum column j1 is assigned.
        { 
          // swap columns j1 and j2, as j2 may be unassigned.
          j1 = j2; 
          i0 = colsol[j2];
        }

      // (re-)assign i to j1, possibly de-assigning an i0.
      rowsol[i] = j1; 
      colsol[j1] = i;

      if (i0 >= 0)           // minimum column j1 assigned earlier.
        if (umin < usubmin) 
	  //if (umin < usubmin - lapeps) 
          // put in current k, and go back to that k.
          // continue augmenting path i - j1 with i0.
          free[--k] = i0; 
        else 
          // no further augmenting reduction possible.
          // store i0 in list of free rows for next phase.
          free[numfree++] = i0; 
    }
  }
  while (loopcnt < 2);       // repeat once.

  // AUGMENT SOLUTION for each free row.
  for (f = 0; f < numfree; f++) 
  {
    freerow = free[f];       // start row of augmenting path.

    // Dijkstra shortest path algorithm.
    // runs until unassigned column added to shortest path tree.
    for (j = 0; j < dim; j++)  
    { 
      d[j] = assigncost[freerow][j] - v[j]; 
      pred[j] = freerow;
      collist[j] = j;        // init column list.
    }

    low = 0; // columns in 0..low-1 are ready, now none.
    up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
             // columns in up..dim-1 are to be considered later to find new minimum, 
             // at this stage the list simply contains all columns 
    unassignedfound = FALSE;
    do
    {
      if (up == low)  // no more columns to be scanned for current minimum.
      {
        last = low - 1; 

        // scan columns for up..dim-1 to find all indices for which new minimum occurs.
        // store these indices between low..up-1 (increasing up). 
        min = d[collist[up++]]; 
        for (k = up; k < dim; k++) 
        {
          j = collist[k]; 
          h = d[j];
          if (h <= min)
          {
            if (h < min)     // new minimum.
            { 
              up = low;      // restart list at index low.
              min = h;
            }
            // new index with same minimum, put on undex up, and extend list.
            collist[k] = collist[up]; 
            collist[up++] = j; 
          }
        }

        // check if any of the minimum columns happens to be unassigned.
        // if so, we have an augmenting path right away.
        for (k = low; k < up; k++) 
          if (colsol[collist[k]] < 0) 
          {
            endofpath = collist[k];
            unassignedfound = TRUE;
            break;
          }
      }

      if (!unassignedfound) 
      {
        // update 'distances' between freerow and all unscanned columns, via next scanned column.
        j1 = collist[low]; 
        low++; 
        i = colsol[j1]; 
        h = assigncost[i][j1] - v[j1] - min;

        for (k = up; k < dim; k++) 
        {
          j = collist[k]; 
          v2 = assigncost[i][j] - v[j] - h;
          if (v2 < d[j])
          {
            pred[j] = i;
            if (v2 == min)   // new column found at same minimum value
              if (colsol[j] < 0) 
              {
                // if unassigned, shortest augmenting path is complete.
                endofpath = j;
                unassignedfound = TRUE;
                break;
              }
              // else add to list to be scanned right away.
              else 
              { 
                collist[k] = collist[up]; 
                collist[up++] = j; 
              }
            d[j] = v2;
          }
        }
      } 
    }
    while (!unassignedfound);

    // update column prices.
    for (k = 0; k <= last; k++)  
    { 
      j1 = collist[k]; 
      v[j1] = v[j1] + d[j1] - min;
    }

    // reset row and column assignments along the alternating path.
    do
    {
      i = pred[endofpath]; 
      colsol[endofpath] = i; 
      j1 = endofpath; 
      endofpath = rowsol[i]; 
      rowsol[i] = j1;
    }
    while (i != freerow);
  }

  // calculate optimal cost.
  cost lapcost = 0;
  for (i = 0; i < dim; i++)  
  {
    j = rowsol[i];
    u[i] = assigncost[i][j] - v[j];
    lapcost = lapcost + assigncost[i][j]; 
  }

  // free reserved memory.
  delete[] pred;
  delete[] free;
  delete[] collist;
  delete[] matches;
  delete[] d;

  return lapcost;
}

cost lap0(int dim, 
        cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v)

// input:
// dim        - problem size
// assigncost - cost matrix

// output:
// rowsol     - column assigned to row in solution
// colsol     - row assigned to column in solution
// u          - dual variables, row reduction numbers
// v          - dual variables, column reduction numbers

{
  boolean unassignedfound;
  row  i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
  col  j, j1, j2, endofpath, last, low, up, *collist, *matches;
  cost min, h, umin, usubmin, v2, *d;

  free = new row[dim];       // list of unassigned rows.
  collist = new col[dim];    // list of columns to be scanned in various ways.
  matches = new col[dim];    // counts how many times a row could be assigned.
  d = new cost[dim];         // 'cost-distance' in augmenting path calculation.
  pred = new row[dim];       // row-predecessor of column in augmenting/alternating path.

  // init how many times a row will be assigned in the column reduction.
  for (i = 0; i < dim; i++)  
    matches[i] = 0;

  // COLUMN REDUCTION 
  for (j = dim-1; j >= 0; j--)    // reverse order gives better results.
  {
    // find minimum cost over rows.
    min = assigncost[0][j]; 
    imin = 0;
    for (i = 1; i < dim; i++)  
      if (assigncost[i][j] < min) 
      { 
        min = assigncost[i][j]; 
        imin = i;
      }
    v[j] = min; 

    if (++matches[imin] == 1) 
    { 
      // init assignment if minimum row assigned for first time.
      rowsol[imin] = j; 
      colsol[j] = imin; 
    }
    else
      colsol[j] = -1;        // row already assigned, column not assigned.
  }

  // REDUCTION TRANSFER
  for (i = 0; i < dim; i++) 
    if (matches[i] == 0)     // fill list of unassigned 'free' rows.
      free[numfree++] = i;
    else
      if (matches[i] == 1)   // transfer reduction from rows that are assigned once.
      {
        j1 = rowsol[i]; 
        min = BIG;
        for (j = 0; j < dim; j++)  
          if (j != j1)
            if (assigncost[i][j] - v[j] < min) 
              min = assigncost[i][j] - v[j];
        v[j1] = v[j1] - min;
      }

  // AUGMENTING ROW REDUCTION 
  int loopcnt = 0;           // do-loop to be done twice.
  do
  {
    loopcnt++;

    // scan all free rows.
    // in some cases, a free row may be replaced with another one to be scanned next.
    k = 0; 
    prvnumfree = numfree; 
    numfree = 0;             // start list of rows still free after augmenting row reduction.
    while (k < prvnumfree)
    {
      i = free[k]; 
      k++;

      // find minimum and second minimum reduced cost over columns.
      umin = assigncost[i][0] - v[0]; 
      j1 = 0; 
      usubmin = BIG;
      for (j = 1; j < dim; j++) 
      {
        h = assigncost[i][j] - v[j];
        if (h < usubmin)
          if (h >= umin) 
          { 
            usubmin = h; 
            j2 = j;
          }
          else 
          { 
            usubmin = umin; 
            umin = h; 
            j2 = j1; 
            j1 = j;
          }
      }

      i0 = colsol[j1];
      if (umin < usubmin) 
        // change the reduction of the minimum column to increase the minimum
        // reduced cost in the row to the subminimum.
        v[j1] = v[j1] - (usubmin - umin);
      else                   // minimum and subminimum equal.
        if (i0 >= 0)         // minimum column j1 is assigned.
        { 
          // swap columns j1 and j2, as j2 may be unassigned.
          j1 = j2; 
          i0 = colsol[j2];
        }

      // (re-)assign i to j1, possibly de-assigning an i0.
      rowsol[i] = j1; 
      colsol[j1] = i;

      if (i0 >= 0)           // minimum column j1 assigned earlier.
        if (umin < usubmin) 
          // put in current k, and go back to that k.
          // continue augmenting path i - j1 with i0.
          free[--k] = i0; 
        else 
          // no further augmenting reduction possible.
          // store i0 in list of free rows for next phase.
          free[numfree++] = i0; 
    }
  }
  while (loopcnt < 2);       // repeat once.

  // AUGMENT SOLUTION for each free row.
  for (f = 0; f < numfree; f++) 
  {
    freerow = free[f];       // start row of augmenting path.

    // Dijkstra shortest path algorithm.
    // runs until unassigned column added to shortest path tree.
    for (j = 0; j < dim; j++)  
    { 
      d[j] = assigncost[freerow][j] - v[j]; 
      pred[j] = freerow;
      collist[j] = j;        // init column list.
    }

    low = 0; // columns in 0..low-1 are ready, now none.
    up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
             // columns in up..dim-1 are to be considered later to find new minimum, 
             // at this stage the list simply contains all columns 
    unassignedfound = FALSE;
    do
    {
      if (up == low)         // no more columns to be scanned for current minimum.
      {
        last = low - 1; 

        // scan columns for up..dim-1 to find all indices for which new minimum occurs.
        // store these indices between low..up-1 (increasing up). 
        min = d[collist[up++]]; 
        for (k = up; k < dim; k++) 
        {
          j = collist[k]; 
          h = d[j];
          if (h <= min)
          {
            if (h < min)     // new minimum.
            { 
              up = low;      // restart list at index low.
              min = h;
            }
            // new index with same minimum, put on undex up, and extend list.
            collist[k] = collist[up]; 
            collist[up++] = j; 
          }
        }

        // check if any of the minimum columns happens to be unassigned.
        // if so, we have an augmenting path right away.
        for (k = low; k < up; k++) 
          if (colsol[collist[k]] < 0) 
          {
            endofpath = collist[k];
            unassignedfound = TRUE;
            break;
          }
      }

      if (!unassignedfound) 
      {
        // update 'distances' between freerow and all unscanned columns, via next scanned column.
        j1 = collist[low]; 
        low++; 
        i = colsol[j1]; 
        h = assigncost[i][j1] - v[j1] - min;

        for (k = up; k < dim; k++) 
        {
          j = collist[k]; 
          v2 = assigncost[i][j] - v[j] - h;
          if (v2 < d[j])
          {
            pred[j] = i;
            if (v2 == min)   // new column found at same minimum value
              if (colsol[j] < 0) 
              {
                // if unassigned, shortest augmenting path is complete.
                endofpath = j;
                unassignedfound = TRUE;
                break;
              }
              // else add to list to be scanned right away.
              else 
              { 
                collist[k] = collist[up]; 
                collist[up++] = j; 
              }
            d[j] = v2;
          }
        }
      } 
    }
    while (!unassignedfound);

    // update column prices.
    for (k = 0; k <= last; k++)  
    { 
      j1 = collist[k]; 
      v[j1] = v[j1] + d[j1] - min;
    }

    // reset row and column assignments along the alternating path.
    do
    {
      i = pred[endofpath]; 
      colsol[endofpath] = i; 
      j1 = endofpath; 
      endofpath = rowsol[i]; 
      rowsol[i] = j1;
    }
    while (i != freerow);
  }

  // calculate optimal cost.
  cost lapcost = 0;
  for (i = 0; i < dim; i++)  
  {
    j = rowsol[i];
    u[i] = assigncost[i][j] - v[j];
    lapcost = lapcost + assigncost[i][j]; 
  }

  // free reserved memory.
  delete[] pred;
  delete[] free;
  delete[] collist;
  delete[] matches;
  delete[] d;

  return lapcost;
}

void checklap(int dim, cost **assigncost,
              col *rowsol, row *colsol, cost *u, cost *v)
{
  row  i;
  col  j;
  cost lapcost = 0, redcost = 0;
  boolean *matched;
  char wait;
  
  matched = new boolean[dim];
  
  for (i = 0; i < dim; i++)  
    for (j = 0; j < dim; j++)  
      if ((redcost = assigncost[i][j] - u[i] - v[j]) < 0)
      {
        printf("negative reduced cost i %d j %d redcost %d\n", i, j, redcost);
        break; 
      }

  for (i = 0; i < dim; i++)  
    if ((redcost = assigncost[i][rowsol[i]] - u[i] - v[rowsol[i]]) != 0)
    {
      printf("non-null reduced cost i %d soli %d redcost %d\n", 
           i, rowsol[i], redcost);
      break; 
    }
  
  for (j = 0; j < dim; j++)  
    matched[j] = FALSE;
    
  for (i = 0; i < dim; i++)  
    if (matched[rowsol[i]])
    {
      printf("column matched more than once - i %d soli %d\n", i, rowsol[i]);
      break; 
    }
    else
      matched[rowsol[i]] = TRUE;
      
    
  for (i = 0; i < dim; i++)  
    if (colsol[rowsol[i]] != i)
    {
      printf("error in row solution i %d soli %d solsoli %d\n", i, 
	     rowsol[i], colsol[rowsol[i]]);
      break; 
    }

  for (j = 0; j < dim; j++)  
    if (rowsol[colsol[j]] != j)
    {
      printf("error in col solution j %d solj %d solsolj %d\n", j, 
      	     colsol[j], rowsol[colsol[j]]);
      break; 
    }

  delete[] matched;
  return;
}

/************************** LAP DRIVER ***************************/

double lap_jv(double **D,int n,double mult)
{
  int i;
  int **D_int;
  int *ys_int, *yt_int, *perm, *invperm;
  double *ys, *yt;
  double obj;

  ys_int = new int[n];
  yt_int = new int[n];
  ys = new double[n];
  yt = new double[n];
  perm = new int[n];
  invperm = new int[n];

  D_int = new int*[n];
  for(i=0;i<n;i++)
    D_int[i] = new int[n];

  obj = lap_jv(D,D_int,perm,invperm,ys,ys_int,yt,yt_int,n,mult);    

  for(i=0;i<n;i++)
    delete [] D_int[i];
  delete [] D_int;
  delete [] ys;  delete [] yt;
  delete [] ys_int;  delete [] yt_int;
  delete [] perm;
  delete [] invperm;
  return obj;
}

double lap_jv(double **D,int *perm,int n,double mult)
{
  int i;
  int **D_int;
  int *ys_int, *yt_int, *invperm;
  double *ys, *yt;
  double obj;

  ys_int = new int[n];
  yt_int = new int[n];
  ys = new double[n];
  yt = new double[n];
  invperm = new int[n];

  D_int = new int*[n];
  for(i=0;i<n;i++)
    D_int[i] = new int[n];

  obj = lap_jv(D,D_int,perm,invperm,ys,ys_int,yt,yt_int,n,mult);    

  for(i=0;i<n;i++)
    delete [] D_int[i];
  delete [] D_int;
  delete [] ys;  delete [] yt;
  delete [] ys_int;  delete [] yt_int;
  delete [] invperm;
  return obj;
}

double lap_jv(double **D,int **D_int,int *perm,int *invperm,
	      double *ys,int *ys_int,double *yt,int *yt_int,int n,double mult)
{
  int i,j;
  int obj;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	D_int[i][j] = (int)(D[i][j]*mult);
	//D_int[i][j] = int(floor(D[i][j]*mult));
      }
  if(n> 4)
    obj = lap(n,D_int,perm,invperm,ys_int,yt_int);
  else
    obj = lap0(n,D_int,perm,invperm,ys_int,yt_int);
  //obj = lap0(n,D_int,perm,invperm,ys_int,yt_int);

  /*
  for(i=0;i<n;i++)
    for(int j=0;j<n;j++)
      if(D_int[i][j] - ys_int[i] - yt_int[j] == 0)
	if(perm[i]==j)
	  printf("(%d,%d) basic arc, flow 1\n",i,j);
	else
	  printf("(%d,%d) basic arc, flow 0\n",i,j);
  printf("\n");
  */
  for(i=0;i<n;i++) {
    ys[i] = ys_int[i] / mult;
    yt[i] = yt_int[i] / mult;
  }

  return ((double)obj)/mult;
}

double lap_jv(double **D,int *perm,
	      double *ys,double *yt,int n,double mult)
{
  int i;
  int **D_int;
  int *ys_int, *yt_int, *ip;
  double bnd;

  D_int = new int*[n];
  for(i=0;i<n;i++)
    D_int[i] = new int[n];

  ys_int = new int[n];
  yt_int = new int[n];
  ip = new int[n];

  bnd =  lap_jv(D,D_int,perm,ip,ys,ys_int,yt,yt_int,n,mult);

  for(i=0;i<n;i++)
    delete [] D_int[i];
  delete [] D_int;
  delete [] ys_int;  delete [] yt_int;
  delete [] ip;

  return bnd;
}

double dual_obj(double *ys,double *yt,int n)
{
  double v = 0.0;
  int i,j;

  for(i=0;i<n;i++)
    v += ys[i] + yt[i];
  return v;
}

