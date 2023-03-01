/*
============================================================
QAPSolution:

Stores the objective value of a QAP, and the permutation
that evaluates to this obj. value.

============================================================
*/

#ifndef QAPSOL_H
#define QAPSOL_H

#include <iostream>

class QAPSolution
{
public:

  QAPSolution();
  QAPSolution(int dim);
  QAPSolution(const QAPSolution &s);
  ~QAPSolution();

  friend std::ostream& operator <<(std::ostream& outs, const QAPSolution &s);
  void print(std::ostream& outs=std::cout);

  double getObjective() { return obj;};
  void setObjective(double val);
  void setPermutation(int *p);
  void update(double val,int *p);
  void update(double val,char *p);

  int ifs_found; // integer feasible solution found
  int *perm;     // array containing solution
  double obj;    // objective value of this solution
  int n;         // problem size
};


#endif
