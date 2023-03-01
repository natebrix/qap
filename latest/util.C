#include "util.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double drand()
{
  return ((double)rand())/RAND_MAX;
}

double score(double opt,double root_bound,double predicted_bound,
	     double scoring_exponent)
{
  double gap = computeRelativeGap(opt,root_bound,predicted_bound);
  //    return 1.0;
  //  return gap;
  return pow(gap,scoring_exponent);
}

void print_time(std::ostream &outs)
{
  char *t;
  time_t ct = time(0);
  t = ctime( &ct );
  t[19] = '\0';
  outs << t << " ";
}

void print_bnb_params(BNBParameters *&params,int nlevel)
{
  
  std::cout << "gap\tdepth\trule \t    FW   \tbest\tupdate" << std::endl;
  for(int i=0;i<nlevel;i++) {

    std::cout << params[i].relative_gap    << "\t";
    std::cout << params[i].depth  << "\t";
    
    std::cout << params[i].rule << "\t";

    if(params[i].which_bound==GLB_BND)
      printf("---/---/--\t");
    else if((params[i].rule==1)||(params[i].rule==2))
      printf("%3d/%3d/--\t",
		params[i].maxier_bnd,
		params[i].maxier_nofathom);
    else
      printf("%3d/%3d/%2d\t",
		params[i].maxier_bnd,
		params[i].maxier_nofathom,
		params[i].maxier_br);

    if((params[i].rule==1)||(params[i].rule==2))
      std::cout << "--" << "\t";
    else
      std::cout << params[i].use_best << "\t";

    if((params[i].which_bound == QP_BND) || 
       (params[i].which_bound == GLB_BND) || 
       (params[i].which_bound == PB_BND))      
      std::cout << "--" << "\t";
    else
      std::cout << params[i].param_step << "\t";

    // identify bound used (default: QPB_PARAM)
    if(params[i].which_bound == GLB_BND)
      std::cout << "(GLB)";
    if(params[i].which_bound == PB_BND)
      std::cout << "(PB)";
    if(params[i].which_bound == EVB3_BND)
      std::cout << "(EVB3)";
    if(params[i].which_bound == QP_GLB_BND)
      std::cout << "(QP/GLB)";
    if(params[i].which_bound == QP_BND_IMP)
      std::cout << "(improved QPB)";

    std::cout << std::endl;
  }

  
}

/*
============================================================
qap_read_params(filename,params,nlevel,sym_file,max_depth);

read parameters from file 'filename', storing them in params.
nlevel contains the number of branching strategies read.
If the parameter file specifies the name of a symmetry file
to read from, it is stored in sym_file.  The additional
parameter max_depth is stored in a separate variable.

============================================================
 */
void qap_read_params(char *filename, BNBParameters *&params, int &nlevel,
		     char *&sym_file,int &max_depth)
{
  int i;
  const int len = 80;
  const int filename_len = 80;
  std::ifstream f;
  char junk[len];
  int sym_file_given;

  f.open(filename);

  f >> sym_file_given;
  f.getline(junk, len);

  if(sym_file_given)
    {
      sym_file = new char[filename_len];
      f.getline(sym_file,filename_len);
    }

  f >> max_depth;
  f.getline(junk, len);

  f >> nlevel;
  f.getline(junk, len);

  params = new BNBParameters[nlevel];

  for(i=0;i<nlevel;i++) {

    if(i < nlevel-1) {
      f >> params[i].relative_gap;
      f.getline(junk, len);

      f >> params[i].depth;
      f.getline(junk, len);
    }

    f >> params[i].rule;
    f.getline(junk, len);

    f >> params[i].maxier_bnd;
    f.getline(junk, len);
    
    f >> params[i].maxier_nofathom;
    f.getline(junk, len);

    f >> params[i].maxier_br;
    f.getline(junk, len);

    f >> params[i].which_bound;
    f.getline(junk, len);

    if(params[i].which_bound == PB_BND) {
      params[i].maxier_bnd = 1;
      params[i].maxier_nofathom = 1;
    }

    f >> params[i].param_step;
    f.getline(junk, len);

    f >> params[i].use_best;
    f.getline(junk, len);

  }
  f.close();
  
  params[nlevel-1].relative_gap = -1.0;
  params[nlevel-1].depth = 1000;


  // --------------------------------------------------------
  std::cout << std::endl;

  print_bnb_params(params,nlevel);
  if(sym_file_given)
    std::cout << "Symmetry file = " << sym_file << std::endl;
  if(max_depth >= 0)
    std::cout << "max_depth = " << max_depth << std::endl;

  std::cout << std::endl;

}

void ins_sort_ascending(double *v, int *order,int n)
{
  double save;
  int j,k,save_p;
  for(k=n-1;k>0;k--)
    {
      j = k;
      save = v[k-1];
      if ( order != NULL )
	save_p = order[k-1];
      while((j<n)&&(save > v[j]))
	{
	  v[j-1] = v[j];
	  if ( order != NULL )
	    order[j-1] = order[j];
	  j = j + 1;
	}
      v[j-1] = save;
      if ( order != NULL )
	order[j-1] = save_p;
    }
}

void ins_sort_descending(double *v, int *order, int n)
{
  double save;
  int j,k,save_p;
  for(k=n-1;k>0;k--)
    {
      j = k;
      save = v[k-1];
      if ( order != NULL )
	save_p = order[k-1];
      while((j<n)&&(save < v[j]))
	{
	  v[j-1] = v[j];
	  if ( order != NULL )
	    order[j-1] = order[j];
	  j = j + 1;
	}
      v[j-1] = save;
      if ( order != NULL )
	order[j-1] = save_p;
    }
}

void nice_array_print(double *d, int n, char delim)
{
  for (int j = 0; j < n; j++) {
    printf("%.3f%c",d[j], delim);
  }
  std::cout << std::endl;
}

void matrix_to_file(std::ofstream &f, double **a, int m, int n) {
  char *s = new char[50]; // todo
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      snprintf(s, 50, "%.8f", a[i][j]);
      f << s << " ";
    }
    f << std::endl;
  }
  f << std::endl;
}

void write_matrix_flat(std::ofstream &f, double **a, int m, int n, char delim) {
  char *s = new char[50]; // todo
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      snprintf(s, 50, "%.8f", a[i][j]);
      f << s << delim;
    }
  }
}



void matrix_to_ampl_file(std::ofstream &f, const char *name, double **a, int m, int n) {
  char *s = new char[50]; // todo
  f << "param " << name << std::endl;
  f << ": ";
  for (int i = 0; i < n; i++) {
    f << (i+1) << " ";
  }
  f << " := " << std::endl;
  for (int i = 0; i < m; i++) {
    f << (i+1) << " ";
    for (int j = 0; j < n; j++) {
      snprintf(s, 50, "%.8f", a[i][j]);
      f << s << " ";
    }
    f << std::endl;
  }
    f << ";" << std::endl;
}

void nice_matrix_print(double **d, int m, int n)
{
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      printf("%.2f ",d[i][j]);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void nice_int_array_print(int *d,int n)
{
  for(int j=0; j<n; j++) {
    printf("%2d ",d[j]);
  }
  std::cout << std::endl;
}

// return true if item occurs in the array list, false otherwise.
bool find(int *list, int n, int item)
{
  for(int i=0;i<n;i++)
    if(list[i] == item)
      return true;
  return false;
}

void shuffle(int *p, int n) {
  for (int k = 0; k < n; k++) {
    int i = (rand() % n);
    int j = (rand() % n);
    int t = p[i];
    p[i] = p[j];
    p[j] = t;
  }
}

