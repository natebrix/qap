#ifndef BRANCH_H
#define BRANCH_H

/*
============================================================

void branch(QAP *qap,QAPAssignment &assign,MAT *X,MAT *U,
	    QAPAssignment *child, int &kids,
	    double bound, double opt,
	    int maxier_br,int use_best,int how,bool use_J3,
	    int &loc_in_row, int &loc_kept);

Notes:
   The basic procedure is:
       1) choose a row i (or column j) on which to branch.
       2) for each X_ij, j=1..n,
          create a child problem where X_ij=1.

   To do step 1, we either 
    a) look at U, and choose the row i with the maximum row sum
    b) look at U, and choose the row i that has the largest
       number of U_ij greater than the current upper bound.
    c) compute a few rows of B.  B_ij is obtained by making the
       assignment i->j (equivalently, X_ij=1) and computing
       the QP bound.  Choose the row with maximum row sum of B.
    d) a lookahead strategy where we compute B, then look at the
       resulting U matrices that are obtained when we compute
       each U_ij.  In some sense this is like looking two levels
       down on this tree.  
    e) choose the row in some really clever way not described here :)

Arguments to branch:
--------------------

X: solution to the QP solved to obtain the bound for this node.

U: obtained from the dual solution to QP.  If the QP bound for
   the current problem is QPB, and we fix X_ij = 1,
   then the resulting QP bound will be at least QPB + U_ij.

child: (output) an array of QAPAssignments to be added to the
       todo list.

bound: the current QP bound

opt: the current upper bound for the QAP.

maxFWier_br: if any QP bounds are evaluated to make branching
             decisions, no more than maxFWier_br iterations
             of the bound will be made.

use_best: if we decide to compute QP bounds to make our 
          branching decision, we will choose the best 

how: the branching strategy used.  For now there are 4 strategies,
     corresponding to strategies a) -- d) outlined above.

loc_in_row:  (output) once we have made the choice of row i,
             the number of X_ij in the parent that were undetermined.

loc_kept:    (output) the number of X_ij that were either a) assigned
             to 1, or b) left alone over all the children.  
             loc_kept /loc_in_row gives us a measure of how much 
             we pruned this subtree during the branching phase.  These
             last two arguments are needed for bookkeeping only...they
             aren't used to make any decisions in the code.

----------------------------------------------------------

IMPORTANT NOTES:
----------------

branch_init(n) must be called before the first call to bnb(),
branch_shutdown() must be called after the last call to bnb().  To 
handle symmetry, branch() uses some global variables.

Symmetry breaking:
------------------

Most of the nugent problems possess symmetry which can 
be exploited to discard redundant problems.  For example,
here is the distance matrix for nug06:

  0 1 2 1 2 3
  1 0 1 2 1 2
  2 1 0 3 2 1  = D
  1 2 3 0 1 2
  2 1 2 1 0 1
  3 2 1 2 1 0

The 6 locations are arranged in a grid structure:

  1 2 3
  4 5 6

[The distance from location 4 to location 3 is 3 for example.]

Suppose I have decided to assign facility 2 to location 6.
Because of the grid structure, this is equivalent to assigning
2 --> 1.  In fact, assignments x-->3, x-->4, x-->6 are all
equivalent to assigning x-->1.  So at the very top level of
my branch and bound tree, I need only consider branches for
locations 1 and 2.  This reasoning is only valid at the
root of the tree.  We (Kurt and I) call this set of locations

J1.  So for nug06, J1 = {1,2}, and instead of creating children
for each location 1..n, I need only create children for the elements
of J1.

There is additional symmetry farther down the tree for some
problems.  Suppose I have already made the assignments
1-->2 and 2-->5.

  ? x ? 
  ? x ?

Then it is clear (at least to me) that if I make the assignment
3-->3, that is the same as assigning 3-->1.  Also 3-->6 is equiv.
to 3-->4. Note that this is valid only because the only two
locations that have been assigned are 2 and 5.


Summarizing:
  J1: the root problem has only |J1| children.
  J2,J3:  if all fixed locations are in the set J2, then the current
          node will only have |J3| children.

At the beginning of the program, we read in a "symmetry file" to
initialize J1,J2,J3.  Sometimes there are several J2,J3 sets.

The format is 
---------------------
|J1|
J1
number of J2,J3 sets

first |J2|
first J2
second |J3|
second J3
...
---------------------
The symmetry files end in .sym.

============================================================
*/

#include "assign.h"
#include "qap.h"
#include "qpbnd.h"

class BranchParameters
{
public:
  int use_best;
  int how;
  bool use_J;
  bool reorder_children;
  bound_t which_bound;
  
  // if you want to forgo branching and just pick a fixed row, col,
  // set how = BRANCH_FIXED, and set thsse fields
  int fixed_index_in_assign;
  int fixed_index_in_U;
  bool choose_row;

  
  BranchParameters() {
    use_best = 9999;
    how = 3;
    use_J = false;
    reorder_children = false;
    which_bound = 2;
    fixed_index_in_assign = -1;
    fixed_index_in_U = -1;
    choose_row = false;
  }
};

void set_branch_parameters(BranchParameters *bp, BNBParameters *param,
			   bool use_J, bool reorder_children);


// utility functions
void find_best_match(MAT *X, int *p);
void merge_U(MAT *U1, double bound1, MAT *U2, double bound2);

// branching functions
void branch(QAP *qap,QAPAssignment &assign,MAT *X,MAT *U,
            QAPAssignment *child, int &kids,
            double bound, QPBParameters *qpbp,
	    BranchParameters *bp,
            int &loc_in_row, int &loc_kept);

void binary_branch(QAP *qap,QAPAssignment &assign,MAT *X,MAT *U,
		   QAPAssignment *child, int &kids,
		   double bound, QPBParameters *qpbp, int use_best,
		   int how,bool use_J,bool reorder_children,
		   int &loc_in_row, int &loc_kept);

bool read_symmetry(char *filename);
void branch_init(int n);
void branch_shutdown();

#endif
