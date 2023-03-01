/*
============================================================
QAPStatistics

Records information about the branch and bound solution 
process.  In particular:

total_time      : Total time spent running the branch and bound alg.
bound_time      : Time spent computing bounds (*)
ier_start_time  : when current iteration started.

stack_size      : current size of todo list.
update_number   : (used by MW)

node      : Total number of nodes in the branch and bound tree,
            including those fathomed.
fathom    : Total number of nodes fathomed by the B&B algorithm.
enumerate : Total number of subproblems ...

b_level     : An array containing the number of bounds computed
              at each level of the tree
f_level     : An array containing the number of times bound exceeded
              the current upper bound. (fathoms)
t_level     : amount of CPU time spent at each level of the tree.

(These last two are described in branch.h: loc_in_row and loc_kept.
Basically these help measure how the tree is pruned during the branching
phase...they are used for bookkeeping only.)
l_level     : total number of locations at this level
e_level     : number of locations eliminated at this level

*** This class is now templated (1/3/99).  QAPStatistics<double>
stores all quantities related to nodes in doubles, QAPStatistics<int>
uses ints, and so on.  Any class may be used so long as it
supports the usual arithmetic operations and assignment to 0.

(*) unimplemented
============================================================
*/

#ifndef QAPSTAT_H
#define QAPSTAT_H

#include <iostream>
#include <time.h>
#define CLOCK_ROLLOVER 4295

#define RECORD_STDDEV

//#define ctr_type int 

//template <class ctr_type>
class QAPStatistics
{
public:

  typedef unsigned int ctr_type;
  typedef unsigned int ctr_type1;

  QAPStatistics();
  QAPStatistics(int n,int m);
  QAPStatistics(const QAPStatistics &s);
  ~QAPStatistics();

  void clear();

  void anotherNode(int depth, int strat,double factor = 1.0);
  void anotherFathom(int depth,double factor = 1.0);
  void anotherEnumeration() {enumerate++;};
  void print(std::ostream& outs=std::cout,int print_level=1) const;

  friend std::ostream& operator <<(std::ostream& outs, const QAPStatistics &s);

  void start();
  void stop();
  void update();
  void updateRelativeGap(int depth,double rel_gap,double factor = 1.0);
  void set_root_bound(double value){root_bound = value;};

  void startIteration(int depth);
  void stopIteration(int depth,int strat,int fw,double factor=1.0);

  void merge_stats(QAPStatistics *new_stats,double factor=1.0);
  void scaleTimes(double factor);
  void scale(double factor);

  /************************************************************/

  double total_time;
  double bound_time;
  double ier_start_time;

  double root_bound;

  int stack_size;
  int update_number; // used only in distributed implementation.
  
  ctr_type node;
  ctr_type fathom;
  ctr_type enumerate;

  ctr_type fw_iter; // number of FW iers

  ctr_type *n_level;  // number of nodes at each level of the tree
  ctr_type *f_level;  // number of times nodes fathomed
  double *t_level;    // cpu time at each level
  ctr_type *l_level;  // total number of locations at this level
  ctr_type *e_level;  // number of locations eliminated at this level

  ctr_type **n_strat; // number of nodes used at each level, broken
                      // down by strategy used.
  double *t_strat;

#ifdef RECORD_STDDEV
  double *rg_sum;
  double *rg_sumsq;
#endif

  int n; // number of levels
  int m; // number of branching stategies

  double start_time;
  double stop_time;

  static const int DEFAULT=0;
  static const int LATEX=1;

  void setPrintMode(int mode) { print_mode = mode; };
private:
  int print_mode;
  void print_default(std::ostream& outs=std::cout,int print_level=1) const;
  void print_latex(std::ostream& outs=std::cout,int print_level=1) const;
  
};


#endif




