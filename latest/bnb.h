/*
============================================================

bnb()          : the branch-and-bound code.
bnb_console()  : testing program.
benchmark()    : used to compare cpu times over different
                 processors.
find_suboptimal(): heuristic technique that uses branch-and-
                   bound code.

============================================================
*/

#ifndef QAPBNB_H
#define QAPBNB_H

#include "qap.h"
#include "assign.h"
#include "link1.h"
#include "solution.h"
#include "stats.h"
#include "branch.h"
#include "util.h"
#include "console.h"
#include <fstream>

/*
  qs_int is the type used by QAPStatistics to record the number
  of nodes at each level, and total number of nodes.  Anything
  that counts nodes in bnb() should be defined as a qs_int.
*/
typedef long long qs_int;   // 64-bit (review B6, Nov 9): int wrapped at 2^31 in the
                            // dive estimates' l_level / e_level sums on nug28-30

void list_clean_up(Node *&todo,double val,QAPStatistics *s);
void pop(Node *&stack,QAPAssignment &a,QAPStatistics *stats);
void push(Node *&stack, const QAPAssignment &a,QAPStatistics *stats);

void log_matrix(BNBParameters *p, const char *data, MAT *U);

BNBParameters *pick_param(const QAPAssignment &node,double root,double opt,
			  BNBParameters *params,int &i,double &rel_gap,
			  int nlevel);


// The root's QPB bound with the first level's parameters, as the search sets
// stats->root_bound at depth 0. A -F run (driver.C) starts below the root and
// needs it: the relative gaps that pick the parameter level and gate the DNN
// bound are measured against it.
double root_qpb_bound(QAP *q, BNBParameters *param);

double QAPobjective(QAP *q,const QAPAssignment &a);
double QAPobjective(QAP *q,Index *p,int m);
void solveSize3QAP(QAP *q,const QAPAssignment &a,double &obj_save,
			  int *p_save);

void bnb(QAP *q,Node *&stack,QAPSolution *s,QAPStatistics *stats,
         int max_CPU, int max_computed_bounds,
	 BNBParameters *params,int nlevel,int max_depth,
	 bool dive,int max_returned_depth);

void benchmark(QAP *q,Node *&stack,QAPSolution *s,
	       QAPStatistics *stats,
	       int max_CPU, int max_computed_bounds,
	       BNBParameters *params,int nlevel,int max_depth,
	       bool dive,int max_returned_depth);

void stack_setup(Node *&stack,QAPStatistics *stats,int &id,int n);



void find_suboptimal(QAP *q,Node *&stack,QAPSolution *s);

#endif
