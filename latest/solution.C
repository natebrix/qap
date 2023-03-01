#include "solution.h"

QAPSolution::QAPSolution()
{
  obj = 9.9e10;
  ifs_found = 0;
}
QAPSolution::QAPSolution(int dim)
{
  obj = 9.9e10;
  ifs_found = 0;
  perm = new int[dim];
  n = dim;
  for(int i=0; i<n; i++)
    perm[i] = -1;
}

QAPSolution::QAPSolution(const QAPSolution &s)
{

  this->obj = s.obj;
  this->ifs_found = s.ifs_found;
  this->n = s.n;
  this->perm = new int[this->n];
  for(int i=0; i<this->n; i++)
    this->perm[i] = s.perm[i];
}

QAPSolution::~QAPSolution()
{
  delete [] perm;
}

void QAPSolution::setObjective(double val)
{
  obj = val;
}

void QAPSolution::setPermutation(int *p)
{
  for(int i=0;i<n;i++)
    perm[i] = p[i];
}

void QAPSolution::update(double val,int *p)
{
  if(val < obj)
    {
      std::cout << "New objective value: " << val << std::endl;
      obj = val;
      setPermutation(p);
      ifs_found = 1;
    }
}

void QAPSolution::update(double val,char *p)
{
  if(val < obj)
    {
      obj = val;
      for(int i=0;i<n;i++)
	perm[i] = p[i];
      ifs_found = 1;
    }
}

std::ostream& operator <<(std::ostream& outs, const QAPSolution &s)
{
  outs << "Objective: " << s.obj << std::endl;
  if(s.ifs_found)
    outs << "incumbent solution: ";
  else
    outs << "incumbent solution: ";
  for(int i=0;i<s.n;i++)
    outs << s.perm[i] << " ";
  outs << std::endl;

}

void QAPSolution::print(std::ostream& outs)
{
  outs << "Objective: " << obj << std::endl;
  if(ifs_found)
    outs << "incumbent solution: ";
  else
    outs << "incumbent solution: ";
  for(int i=0;i<n;i++)
    outs << perm[i] << " ";
  outs << std::endl;

}
