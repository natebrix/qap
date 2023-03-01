#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

/****************************************************************/
/*
    This programme implement a simulated annealing for the
    quadratic assignment problem along the lines describes in
    the article D. T. Connoly, "An improved annealing scheme for 
    the QAP", European Journal of Operational Research 46, 1990,
    93-100.

    Compiler : g++ should work. 

    Author : E. Taillard, 
             IDSIA, Corso Elvezia 36, 6900 Lugano, Switzerland

    Date : 16. 3. 98

    Format of data file : Example for problem nug5 :

5

0 1 1 2 3
1 0 2 1 2
1 2 0 1 2
2 1 1 0 1
3 2 2 1 0

0 5 2 4 1
5 0 3 0 2
2 3 0 0 0
4 0 0 0 5
1 2 0 5 0

   Additionnal parameters : Number of iterations, number of runs

*/

/********************************************************************/
typedef double qap_type ;
    
const int n_max = 151;
const qap_type infini = 1399999999;
const int nb_iter_initialisation = 1000; // Connolly propose nb_iterations/100

qap_type optimum;

//typedef qap_type type_matrice[n_max][n_max];
typedef qap_type **type_matrice;

/*--------------- choses manquantes -----------------*/

int max(int a, int b) {if (a > b) return(a); else return(b);};
double max(double a, double b) {if (a > b) return(a); else return(b);}
int min(int a, int b) {if (a < b) return(a); else return(b);}
double min(double a, double b) {if (a < b) return(a); else return(b);}
/*

*/
template <class Item>
void swap(Item &a,  Item &b) {Item temp = a; a = b; b = temp;}

void a_la_ligne(std::ifstream & fichier_donnees)
{char poubelle[1000]; fichier_donnees.getline(poubelle, sizeof(poubelle));}
/*-------------------------------------------------*/

/************* random number generators ****************/

const int m = 2147483647; const int m2 = 2145483479; 
const int a12 = 63308; const int q12 = 33921; const int r12 = 12979; 
const int a13 = -183326; const int q13 = 11714; const int r13 = 2883; 
const int a21 = 86098; const int q21 = 24919; const int r21 = 7417; 
const int a23 = -539608; const int q23 = 3976; const int r23 = 2071;
const double invm = 4.656612873077393e-10;
int x10 = 12345, x11 = 67890, x12 = 13579, 
     x20 = 24680, x21 = 98765, x22 = 43210;

double mon_rand()
{
  int h, p12, p13, p21, p23;
  h = x10/q13; p13 = -a13*(x10-h*q13)-h*r13;
  h = x11/q12; p12 = a12*(x11-h*q12)-h*r12;
  if (p13 < 0) p13 = p13 + m; if (p12 < 0) p12 = p12 + m;
  x10 = x11; x11 = x12; x12 = p12-p13; if (x12 < 0) x12 = x12 + m;
  h = x20/q23; p23 = -a23*(x20-h*q23)-h*r23;
  h = x22/q21; p21 = a21*(x22-h*q21)-h*r21;
  if (p23 < 0) p23 = p23 + m2; if (p21 < 0) p21 = p21 + m2;
  x20 = x21; x21 = x22; x22 = p21-p23; if(x22 < 0) x22 = x22 + m2;
  if (x12 < x22) h = x12 - x22 + m; else h = x12 - x22;
  if (h == 0) return(1.0); else return(h*invm);
}

int unif(int low, int high)
{
  return(low + long(double(high - low + 1) *  mon_rand() - 0.5));
}

/************************** sa for qap ********************************/

void lire(char *nom_fichier,int &n, qap_type **&a, qap_type **&b)
{
  std::ifstream fichier_donnees;
  // char nom_fichier[30];
  int i, j;
  
  //  cout << "nom du fichier de donnees : \n";
  //  cin >> nom_fichier;
  //  cout << nom_fichier << '\n';

  fichier_donnees.open(nom_fichier);
  fichier_donnees >> n; a_la_ligne(fichier_donnees);

  a = new qap_type*[n];
  b = new qap_type*[n];
  for(i = 0; i < n; i = i + 1) 
    {
      a[i] = new qap_type[n];
      b[i] = new qap_type[n];
    }

  for (i = 0; i < n; i = i+1) for (j = 0; j < n; j = j+1)
    fichier_donnees >> a[i][j];
  for (i = 0; i < n; i = i+1) for (j = 0; j < n; j = j+1)
    fichier_donnees >> b[i][j];
  fichier_donnees.close();
}

qap_type calc_delta_complet2(int n, qap_type **& a, qap_type **& b,
			     int * & p, int r, int s)
{
  qap_type d;
  d = (a[r][r]-a[s][s])*(b[p[s]][p[s]]-b[p[r]][p[r]]) +
    (a[r][s]-a[s][r])*(b[p[s]][p[r]]-b[p[r]][p[s]]);
  for (int k = 0; k < n; k = k + 1) 
    if (k!=r && k!=s)
      d = d + (a[k][r]-a[k][s])*(b[p[k]][p[s]]-b[p[k]][p[r]]) +
	(a[r][k]-a[s][k])*(b[p[s]][p[k]]-b[p[r]][p[k]]);
  return d;
}

inline qap_type calcule_cout(qap_type n, qap_type **a, qap_type **b, 
			     int * p)
{
  qap_type c = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      c += a[i][j] * b[p[i]][p[j]];
  return c;
}

inline void tire_solution_aleatoire(int n, int * p)
{
  int i;
  for (i = 0; i < n; i = i+1) p[i] = i;
  for (i = 1; i < n; i = i+1) swap(p[i-1], p[unif(i-1, n-1)]);
}

void recuit(int n, qap_type **a, qap_type **b,
            int * meilleure_sol, qap_type & meilleur_cout,
            int nb_iterations)

 {
   int * p;
   int i, r, s;
   qap_type delta;
   qap_type k = n*(n-1)/2, mxfail = k, nb_fail, no_iteration;
   qap_type dmin = infini, dmax = 0;
   qap_type Cout;
   double t0, tf, beta, tfound, temperature;

   p = new int[n];
   for (i = 0; i < n; i = i + 1) 
     p[i] = meilleure_sol[i];

   Cout = calcule_cout(n, a, b, p);
   meilleur_cout = Cout;

  for(no_iteration = 1; no_iteration <= nb_iter_initialisation;
      no_iteration++)
    {
      r = unif(0, n-1);
      s = unif(0, n-2);
      if (s >= r) s = s+1;

      delta = calc_delta_complet2(n,a,b,p,r,s);
      if (delta > 0)
	{dmin = min(dmin, delta); dmax = max(dmax, delta);}; 
      Cout = Cout + delta;
      swap(p[r], p[s]);
    }
  t0 = dmin + (dmax - dmin)/10.0;
  tf = dmin;
  beta = (t0 - tf)/(nb_iterations*t0*tf);
  
  nb_fail = 0;
  tfound = t0;
  temperature = t0;
  r = 0; s = 1;
  //  r = 1; s = 2;
  for (no_iteration = 1; 
       no_iteration <= nb_iterations - nb_iter_initialisation; 
       no_iteration = no_iteration + 1)
    {
      if (meilleur_cout > optimum)
	{temperature = temperature / (1.0 + beta*temperature);
	
	s = s + 1;
	if (s >= n)
	  {r = r + 1; 
	  if (r > n - 2) r = 0;
	  s = r + 1;
	  };
	
	delta = calc_delta_complet2(n,a,b,p,r,s);
	if ((delta < 0) || (mon_rand() < exp(-double(delta)/temperature)) ||
	    mxfail == nb_fail)
	  {Cout = Cout + delta; swap(p[r], p[s]); nb_fail = 0;}
	else nb_fail = nb_fail + 1;
	
	if (mxfail == nb_fail) {beta = 0; temperature = tfound;};
	if (Cout < meilleur_cout)
	  {
	    meilleur_cout = Cout;
	    for (i = 0; i < n; i++) meilleure_sol[i] = p[i];
	    tfound = temperature;
	  }
	}
    }
  delete [] p;
 }

void sa_qap(double **a, double **b, int n, 
	    int *p, double &Cost,int nb_iterations)
{
  tire_solution_aleatoire(n, p);
  recuit(n,a,b,p,Cost, nb_iterations);
}


#if 0

int main(int argc, char *argv[])
{
  int  i,n, nb_iterations, nb_res, no_res;
  qap_type Cost;
  qap_type **a, b;
  int * p;

  nb_iterations = 100000;
  nb_res = 1;
  if(argc < 2) {
    cout << "First argument should be QAPLIB input file." << std::endl;
    return 0;
  }
  if(argc > 2) 
    nb_iterations = atoi(argv[2]);
  if(argc > 3) 
    nb_res = atoi(argv[3]);

  lire(argv[1],n, a, b);
  p = new int[n];

  //  cout << "# iterations, # resolutions : \n";
  //  cin >> nb_iterations >> nb_res;
  for (no_res = 1; no_res <= nb_res; no_res = no_res + 1)
    {
      sa_qap(a,b,n,p,Cost,nb_iterations);

      cout << "Meilleure solution trouvee : " << Cost << std::endl;
      for(i = 0; i < n; i = i + 1) 
	cout << p[i] << ' ';
      cout << '\n';
    };


  for(i = 0; i < n; i = i + 1) 
    {
      delete [] a[i];
      delete [] b[i];
    }
  delete [] a;
  delete [] b;
}
#endif

