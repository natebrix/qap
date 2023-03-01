#include "util.h"
#include "linalg.h"
#include "lap.h"
#include "bnb.h"
#include "console.h"
#include "stats.h"
#include "solution.h"
#include "qpbnd.h"
#include "branch.h"

#ifdef HAVE_EVB3
#include "evb3.h"
#endif

#include <fstream>
#include <math.h>
//#include <time.h>

// compute root bound, then fix every variable to 1 and compute resulting
// bound.
double level1_bounds(QAP *q,QPBParameters *qpbp,int maxier_bnd,int maxier_br)
{
  int i,j,ii,jj;
  double bound;
  MAT *X,*X1,*U,*B;
  VEC *row_sum;
  QAP *qr;

  int n = q->A->m;
  QAPAssignment a(n),a0(n);
  BranchParameters *bp = new BranchParameters();
  
  B = m_get(n,n);
  X = m_get(n,n);
  X1 = m_get(n-1,n-1);
  U = m_get(n,n);
  row_sum = v_get(n);

  init_X(X,n);

  qpbp->maxier = maxier_bnd;
  qpbp->maxier_nofathom = maxier_bnd;

  bound = qpbnd(q->A,q->B,q->C,X,U,qpbp);
  std::cout << "Root: " << ceil(bound) << std::endl;
  n--;

  qpbp->maxier = maxier_br;
  qpbp->maxier_nofathom = maxier_br;
  for(ii=0;ii<=n;ii++)
    {
      std::cout << ".";
      std::cout.flush();
      for(jj=0;jj<=n;jj++)
	{
	  //init_X(X1,n);
	  warmstart(X,X1,ii,jj);
  
	  a = a0;
	  a.fix(ii,jj);
	  qr = reduceQAP(q,a);

	  B->me[ii][jj] =   glbnd(qr->A,qr->B,qr->C,U,qr->shift);
	  //	  B->me[ii][jj] = qpbnd(qr->A,qr->B,qr->C,X1,U,qpbp);
	  /*
	  B->me[ii][jj] = qpbnd(qr->A,qr->B,qr->C,X1,U,maxier_kid,maxier_kid,
				qr->shift,10000.0,false);
	  */
	  row_sum->ve[ii] += B->me[ii][jj];
	}
    }
  std::cout << std::endl;
  for(ii=0;ii<=n;ii++)
    {
      for(jj=0;jj<=n;jj++)
	{
	  printf("%4.0f ",ceil(B->me[ii][jj]));
	}
      printf("\n");
    }

  std::cout << "Row sum: " << std::endl;
  for(jj=0;jj<=n;jj++)
    printf("%4.0f ",ceil(row_sum->ve[jj]));
  printf("\n");

  int *order = new int[n+1];
  for(jj=0;jj<=n;jj++)
    order[jj] = jj;
  ins_sort_descending(row_sum->ve,order,n);

  std::cout << "Rank: " << std::endl;
  for(jj=0;jj<=n;jj++)
    printf("%4d ",order[jj]);
  printf("\n");


  delete [] order;
  M_FREE(B);
  M_FREE(X);
  M_FREE(U);
  V_FREE(row_sum);

  return bound;
}

VEC *my_variance(MAT *X)
{
  VEC *mean, *var;
  int i,j;
  int m = X->m, n = X->n;

  mean = v_get(m);
  var  = v_get(m);

  for(i=0;i<m;i++) {
    for(j=0;j<n;j++) {
      mean->ve[i] += X->me[i][j];
    }
    mean->ve[i] = mean->ve[i] / n;
  }

  for(i=0;i<m;i++) {
    for(j=0;j<n;j++) {
      var->ve[i] += (X->me[i][j] - mean->ve[i])*(X->me[i][j] - mean->ve[i]);
    }
  }

  V_FREE(mean);
  return var;
}

void test_vector(VEC *v, const char* name)
{
  printf("%s\n", name);
  nice_array_print(v->ve, v->dim);
  double norm = sqrt(in_prod(v, v));
  printf("%f\n", norm);
}

void qpb_test(QAP *q)
{
  int status;
  int iter,total_iter;
  double eig;
  double lo_bound, up_bound;
  MAT *S,*T;
  MAT *G;

  MAT *W,*U;
  VEC *s,*t,*w,*u;

  MAT *A = q->A;
  MAT *B = q->B;
  MAT *C = q->C;
  
  int *xstar;

  double C_X, G_X;

  up_bound =  1.0e99;
  lo_bound = -1.0e99;

  S = MNULL;  T = MNULL;
  W = MNULL;  U = MNULL;
  w = VNULL;  u = VNULL;  s = VNULL;  t = VNULL;

  eig = init_ST(A, B, C, w, u, W, U, s, t, S, T, true);

  test_vector(w, "w");
  test_vector(u, "u");
  printf("W\n");
  nice_matrix_print(W->me, 3, 3);
  printf("U\n");
  nice_matrix_print(U->me, 3, 3);

  test_vector(s, "s");
  test_vector(t, "t");

  printf("S\n");
  nice_matrix_print(S->me, 3, 3);
  printf("T\n");
  nice_matrix_print(T->me, 3, 3);

  MAT *V = projectionMatrix(A->n);
  printf("V\n");
  nice_matrix_print(V->me, 3, 3);
  
  MAT *A1 = m_get(A->n-1, A->n-1);
  vtav(A, A1);
  printf("VTAV(A)\n");
  nice_matrix_print(A1->me, 3, 3);

  A1 = m_get(A->n+1, A->n+1);
  vavt(A, A1);
  printf("VAVT(A)\n");
  nice_matrix_print(A1->me, 3, 3);
}


Node *create_stack_from_node(const QAPAssignment &node, QAPStatistics *stats)
{
  stats->clear();
  Node *stack = NULL;
  push(stack, node, stats);
  return stack;
}

Node *create_stack_from_nodes(QAPAssignment *nodes, int k, QAPStatistics *stats)
{
  stats->clear();
  Node *stack = NULL;
  for (int i = 0; i < k; i++) {
    push(stack, nodes[i], stats);
  }
  return stack;
}

// need a variant that takes a list.
void run_at_node(QAP *q, QAPAssignment &node, QAPSolution *sln,
		 QAPStatistics *stats, BNBParameters *params,
		 ConsoleParameters *c) {
  Node *stack = create_stack_from_node(node, stats);
  bnb(q,stack,sln,stats,1000,1000,params, c->nlevel, -1, false, -1);
  
  std::cout << "  ---- RESULTS ----" << std::endl;
  stats->print(std::cout,1);
  sln->print();

  c->set_log_prefix("g");
  log_bnb_run(c, node, sln, stats);
}

// run one node and set its bound. Returns true if children.
bool run_one_node(QAP *q, QAPAssignment &node, 
		  QAPSolution *sln, QAPStatistics *stats,
		  BNBParameters *params, ConsoleParameters *c) {
  double inc = sln->getObjective(); // reset this every time.
  Node *stack = create_stack_from_node(node, stats);
  params->log_file = c->log_file;
  printf("log prefix = %s\n", c->log_prefix.c_str());
  params->log_prefix = c->log_prefix;

  bnb(q, stack, sln, stats, 1000, 1000, params, c->nlevel, node.depth, false, -1);
  params->log_file = NULL;
  sln->setObjective(inc);

  if (stack) {
    node.bound = stack->data.bound; // todo may not have anything on stack.
  }
  return (stack != NULL);
}

void run_for_unassigned(QAP *q, QAPAssignment &node,
			QAPSolution *sln, QAPStatistics *stats,
			BNBParameters *params, ConsoleParameters *c,
			QPBParameters *qpbp, BranchParameters *bp, 
			int *unassigned, const char *prefix) {
{
  double inc = sln->getObjective(); // reset this every time.
  
  int m = node.n - node.nfix;
  QAPAssignment child[node.n];
  int kids;
  int loc_in_row, loc_kept;
  MAT *U = m_get(m, m);
  char *run_prefix = new char[10];
  printf("Trying %d different branching decisions for %s\n", m, c->log_prefix.c_str());
  for (int i = 0; i < m; i++) {
    m_zero(U);
    bp->fixed_index_in_assign = unassigned[i];
    bp->fixed_index_in_U = i;
    branch(q, node, U, U, child, kids, node.bound, qpbp, bp,
	   loc_in_row, loc_kept);


    for (int j = 0; j < kids; j++) {
      child[i].bound = node.bound;
    }
       
    // now run BNB on these kids.
    Node *stack = create_stack_from_nodes(child, kids, stats);
    bnb(q, stack, sln, stats, 1000, 1000, params, c->nlevel, -1, false, -1);
    snprintf(run_prefix, 10, "%s|%d", prefix, i);
    c->set_log_prefix(run_prefix);
    log_bnb_run(c, child[0], sln, stats);
    sln->setObjective(inc);
  }
  M_FREE(U);
 }
}

void run_for_nodes(QAP *q, QAPAssignment &node, 
		   QAPSolution *sln, QAPStatistics *stats,
		   BNBParameters *params, ConsoleParameters *c,
		   QPBParameters *qpbp, BranchParameters *bp,
		   bool rows, const char *prefix) {
  int *unassigned = node.get_unassigned(rows);
  bp->how = BRANCH_FIXED;
  bp->choose_row = rows;
  run_for_unassigned(q, node, sln, stats, params, c, qpbp, bp, unassigned, prefix);
  delete[] unassigned;
}

void run_all_branches(QAP *q, QAPAssignment &node, 
		      QAPSolution *sln, QAPStatistics *stats,
		      BNBParameters *params, ConsoleParameters *c,
		      const char* prefix) {
  const int prefix_len = 35;
  char run_prefix[prefix_len];
  snprintf(run_prefix, prefix_len, "%s|_", prefix);
  c->set_log_prefix(run_prefix);
  // compute the bound first.
  if (run_one_node(q, node, sln, stats, params, c)) {
    int m = node.n - node.nfix;
    BranchParameters *bp = new BranchParameters();
    QPBParameters *qpbp = new QPBParameters();
    qpbp->compute_lapU = false;  qpbp->lapU = 0.0; qpbp->feasible = true;
    qpbp->lap_scale = 1.0;
    qpbp->inc = sln->getObjective();
    qpbp->shift = q->shift;     
    qpbp->debug = false;

    // These are all the row assignments.
    snprintf(run_prefix, prefix_len, "t|r|%s", prefix);
    run_for_nodes(q, node, sln, stats, params, c, qpbp, bp, true, run_prefix);

    // These are all the column assignments.
    snprintf(run_prefix, prefix_len, "t|c|%s", prefix);
    run_for_nodes(q, node, sln, stats, params, c, qpbp, bp, false, run_prefix);
  }
}

void the_time() {
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  printf("%d-%02d-%02d-%02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

}

void start_log(ConsoleParameters *c) {
  the_time();
  if (c->log_file && c->log_header) {
    std::ofstream outfile;
    outfile.open(c->log_file, std::ios_base::app);
    outfile << "*"
	    << "|" << c->input_file 
	    << "|" << c->param_file
	    << "|" << c->sym_file;
    outfile << std::endl;
  }
}

void log_fix(ConsoleParameters *c, int i, int j) {
  if (c->log_file) {
    std::ofstream outfile;
    outfile.open(c->log_file, std::ios_base::app);
    outfile << "f"
	    << "|" << i
	    << "|" << j;
    outfile << std::endl;
  }
}

void log_int(ConsoleParameters *c, const char *str, int i) {
  if (c->log_file) {
    std::ofstream outfile;
    outfile.open(c->log_file, std::ios_base::app);
    outfile << str
	    << "|" << i;
    outfile << str << std::endl;
  }
}

void log_string(ConsoleParameters *c, const char *str) {
  if (c->log_file) {
    std::ofstream outfile;
    outfile.open(c->log_file, std::ios_base::app);
    outfile << str << std::endl;
  }
}

void log_bnb_run(ConsoleParameters *c,
		 QAPAssignment &node,
		 QAPSolution *sln,
		 QAPStatistics *stats)
{
  if (c->log_file) {
    std::ofstream outfile;
    outfile.open(c->log_file, std::ios_base::app);
    outfile << c->log_prefix << "|";
    write_permutation(outfile, node.p, node.n);
    outfile 
      << "|" << sln->getObjective()
      << "|" << node.bound
      << "|" << stats->total_time 
      << "|" << stats->node
      << "|" << stats->fathom
      << "|" << stats->enumerate
      << "|" << stats->root_bound;
    outfile << std::endl;
  }
}

// todo: I want to also tag which "choice" I am making
// the choice is an in U

void one_branch_trial(QAP *q, QAPAssignment &node,
		      QAPSolution *sln,
		      QAPStatistics *stats,
		      BNBParameters *params,
		      ConsoleParameters *c,
		      QPBParameters *qpbp,
		      MAT *X,
		      MAT *U, int m,
		      const char *prefix) {
  // started with n, want m. But nfix have already been assigned, and we
  // will try all children.
  int count = node.n - m - node.nfix - 1;
  QAPAssignment child = randomly_fix(node, count);
  std::cout << child << std::endl;

  run_all_branches(q, child, sln, stats, params, c, prefix);
  c->set_log_prefix("g");
  log_bnb_run(c, node, sln, stats);
}

void branch_trial(QAP *q, QAPAssignment &node,
	      QAPSolution *sln,
	      QAPStatistics *stats,
	      BNBParameters *params,
	      ConsoleParameters *c,
	      QPBParameters *qpbp,
	      int m, int trials) {
  MAT *X = NULL;
  MAT *U = NULL;
  const int prefix_len = 5;
  char prefix[prefix_len];
  for (int trial = 0; trial < trials; trial++) {
    printf("Trial %d/%d starting\n", trial+1, trials);
    // this is the trial. Can I uniquify it?
    snprintf(prefix, prefix_len, "%d", trial+1);
    one_branch_trial(q, node, sln, stats, params, c, qpbp, X, U, m, prefix);
  }
}


//------------------

char next_choice(ConsoleParameters *c, QAPAssignment &node) {
  char choice;
  std::cout << std::endl << std::endl;
  std::cout << "a) compute B matrix" << std::endl;
  std::cout << "b) compute bound" << std::endl;
  std::cout << "c) clear assignments" << std::endl;
  std::cout << "d) disallow assignment" << std::endl;
  std::cout << "f) fix variable " << std::endl;
  std::cout << "g) go run bnb on this node" << std::endl;
  std::cout << "h) heuristic" << std::endl;
  std::cout << "l) change depth of this node" << std::endl;
  std::cout << "m) create branching training data" << std::endl;
  std::cout << "o) output current subproblem to QAPLIB file" << std::endl;
  std::cout << "p) change parameters" << std::endl;
  std::cout << "r) branching decision" << std::endl;
  std::cout << "s) evaluate a permutation" << std::endl;
  std::cout << "t) try bnb on all possible branching decisions" << std::endl;
  //    std::cout << "t) trace path to solution" << std::endl;
  std::cout << "u) show U matrix " << std::endl;
  std::cout << "v) unit tests " << std::endl;
  std::cout << "/) change debug mode " << std::endl;
  std::cout << "q) quit" << std::endl;
  std::cout << "current assignments:" << std::endl;
  node.print();
  
  std::cout << "bnb> ";
  std::cin >> choice;
  return choice;
}

int get_int_value(const char *prompt) {
  int i;
  std::cout << prompt << "> ";
  std::cin >> i;
  return i;
}

void fix(ConsoleParameters *c, QAPAssignment &node) {
  int i, j;
  std::cout << "row> ";
  std::cin >> i;
  std::cout << "col> ";
  std::cin >> j;
  node.fix(i,j);
  log_fix(c, i, j);
}

void branch_training_data(QAP *q, QAPAssignment &node,
			  BNBParameters *params,
			  QAPSolution *sln,
			  QAPStatistics *stats,
			  ConsoleParameters *c,
			  QPBParameters *qpbp) {
  int trials = get_int_value("trials");
  int m = get_int_value("size");
  branch_trial(q, node, sln, stats, params, c, qpbp, m, trials);
}


void change_param(const QAPAssignment &node, BNBParameters *params,
		  QAPSolution *sln, QAPStatistics *stats) {
  char choice2;
  int strat;
  double rel_gap;
  int nlevel;
  BNBParameters *param;
  double lap_scale;


  std::cout << "Current Strategy: #" << param - params << std::endl;
  std::cout << "1) maxier (bound) = " << param->maxier_bnd << std::endl;
  std::cout << "2) maxier (fathom) = " << param->maxier_nofathom << std::endl;
  std::cout << "3) maxier (branch) = " << param->maxier_br << std::endl;
  std::cout << "4) step = " << param->param_step << std::endl;
  std::cout << "5) branching rule = " << param->rule << std::endl;
  std::cout << "6) incumbent = " << sln->obj << std::endl;
  std::cout << "7) LAP scaling = " << lap_scale << std::endl;
  std::cout << "0) never mind..." << std::endl;
  std::cout << "bnb (param)> ";
  
  std::cin >> choice2;
  switch(choice2) {
  case '1':
    std::cout << "maxier (bound): ";	std::cin >> param->maxier_bnd;
    break;
  case '2':
    std::cout << "maxier (fathom): ";	std::cin >> param->maxier_nofathom;
    break;
  case '3':
    std::cout << "maxier (branch)> ";	std::cin >> param->maxier_br;
    break;
  case '4':
    std::cout << "step> ";	std::cin >> param->param_step;
    break;
  case '5':
    std::cout << "rule> ";	std::cin >> param->rule;
    break;
  case '6':
    std::cout << "incumbent> ";	std::cin >> sln->obj;
    break;
  case '7':
    std::cout << "scale> ";	std::cin >> lap_scale;
    break;
  default:
    break;
  }
}


void print_B(QAP *q, QAPAssignment &node, QAPSolution *sln,
	     QAPStatistics *stats, BNBParameters *params,
	     ConsoleParameters *c, BranchParameters *bp,
	     BNBParameters *param, QPBParameters *qpbp) {
  int n = q->A->m;
  MAT *B, *X, *U, *X1, *prevX;
  int i, j;
  double bound;
  int kids;
  QAPAssignment kid[n];

  B = m_get(n,n);
  X = m_get(n,n);
  X1 = m_get(n,n);
  prevX = m_get(n,n);
  U = m_get(n,n);

  QAP *qr = reduceQAP(q,node);
  
  printf("qr->n = %d\n",qr->A->n);
  
  bound = level1_bounds(qr,qpbp,param->maxier_bnd,param->maxier_br);
  
  qpbp->inc = 1e99; // sln->getObjective();
  qpbp->maxier = param->maxier_br;
  qpbp->maxier_nofathom = param->maxier_br;
  qpbp->step = param->param_step;
  
  set_branch_parameters(bp, param, false, false);
  bp->how = 3;
  branch(q,node,X,U,kid,kids,bound,qpbp,bp,i,j);
  
  nice_m_output(U);
  log_matrix(params, "B", U);
  delete qr;
}

double console_qpb(QAP *q, QAPAssignment &node,
		   QAPSolution *sln,
		   QAPStatistics *stats,
		   QPBParameters *qpbp, MAT *X, MAT *U) {
  QAP *qr = reduceQAP(q, node);
  int m = qr->A->m;

  X = m_resize(X, m, m);
  U = m_resize(U, m, m);
  
  init_X(X, m);
  
  qpbp->inc = sln->getObjective();
  qpbp->shift = qr->shift;     
  double bound = qpbnd(qr->A,qr->B,qr->C,X,U,qpbp);
  delete qr;
  return bound;
}

void bnb_console(QAP *q, Node *&stack,
		 QAPSolution *sln,
		 QAPStatistics *stats,
		 BNBParameters *params, ConsoleParameters *c)
{
#ifdef HAVE_EVB3
  EVB3 evb3;
#endif
  QAP *q_evb3=NULL;

  stopwatch s_bnd;
  int *p;
  int n = q->A->m,m;
  int i,j;
  bool done = false;
  char input[80];
  double bound, heuristic_obj;
  int strat;
  double rel_gap;
  QPBParameters *qpbp;
  QAPAssignment node(n),init;
  BNBParameters *param;

  int kids;
  int loc_in_row, loc_kept;
  QAPAssignment kid[n];

  QAP *qr=NULL;
  bool have_prevX = false;
  MAT *B, *X, *U, *X1, *prevX;
  int *fr,*fc;
  VEC *row_sum;

  B = m_get(n,n);
  X = m_get(n,n);
  X1 = m_get(n,n);
  prevX = m_get(n,n);
  U = m_get(n,n);
  row_sum = v_get(n);
  p = new int[n];
  fr = new int[n];
  fc = new int[n];

  qpbp = new QPBParameters();
  qpbp->compute_lapU = false;  qpbp->lapU = 0.0; qpbp->feasible = true;
  qpbp->lap_scale = 1.0;
  BranchParameters *bp = new BranchParameters();
  
  start_log(c);
  
  pop(stack,node,stats);
  init = node;
  param = &(params[0]);
  params->log_file = c->log_file;
  params->log_prefix = c->log_prefix;

  do {
    char choice = next_choice(c, node);
    
    qpbp->inc = 1e99;  // don't want to fathom early
    qpbp->maxier = param->maxier_bnd;
    qpbp->maxier_nofathom = param->maxier_nofathom;
    qpbp->step = param->param_step;

    switch (choice) {

    case '/':
      qpbp->debug = !qpbp->debug;
      std::cout << "debug mode = " << qpbp->debug << std::endl;
      break;

    case 'a':     case 'A':
      print_B(q, node, sln, stats, params, c, bp, param, qpbp);
      break;

    case 'b':     case 'B':     
      delete qr;
      //      qr = new QAP(n);
      qr = reduceQAP(q, node);
      m = qr->A->m;

      if(m <= 5) {
	printQAP(qr);
      }
      
      X = m_resize(X,m,m);
      U = m_resize(U,m,m);

      init_X(X,m);

      qpbp->inc = sln->getObjective();
      qpbp->shift = qr->shift;     
      bound = qpbnd(qr->A,qr->B,qr->C,X,U,qpbp);


      std::cout << "QPB \t\t= " << bound << std::endl;

      /*      
      for(double lam = 0.0;lam<=1.0;lam+=0.1) {
	delete qr;
	qr = reduceQAP(q,node);
	qpbp->shift = qr->shift;     
	init_X(X,m);

	std::cout << lam << std::endl;
	bound = qpbnd0(qr->A,qr->B,qr->C,X,U,qpbp,lam);
	std::cout << "QPB0\t\t= " << bound << std::endl;
	//	std::cout << mm_trace(qr->C,X) + qpbp->shift << std::endl;
      }
      */

      if(have_prevX) {

	if(prevX->m > m) {
	  // we fixed some variables since the last bound computation.
	  // need to compensate by extracting some rows and cols from X.
	  int just_fixed = prevX->m - m; 
	  int r,r1,c,c1;

	  for(i=0;i<prevX->m;i++) {
	    fr[i] = fc[i] = -1;
	  }

	  for(i=node.nfix-just_fixed;i<node.nfix;i++) {
	    r = node.history[i]; r1 = r;
	    for(j=0;j<node.nfix-just_fixed;j++)
	      if(node.history[j] < r)
		r1--;
	    c = node.p[node.history[i]]; c1 = c;
	    for(j=0;j<node.nfix-just_fixed;j++)
	      if(node.p[node.history[j]] < c)
		c1--;
	    fr[r1] = 0;
	    fc[c1] = 0;
	  }

	  for(i=0;i<prevX->m;i++)
	    printf("%2d ",fr[i]);
	  printf("\n");
	  for(i=0;i<prevX->m;i++)
	    printf("%2d ",fc[i]);
	  printf("\n");

	  submatrix(prevX,X,fr,fc,m,prevX->m);
	  find_best_match(X,p);
	  m_zero(X);
	  for(i=0;i<X->m;i++)
	    X->me[i][p[i]] = 1.0;
	}

	bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);

	std::cout << "QPBpf\t\t= " << bound << std::endl;
      }

      s_bnd.go();
      //for(int ii=0;ii<100;ii++) {
       init_X(X,m);
	bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);
	//}
      s_bnd.stop();
      
      m_copy(X,prevX);  have_prevX = true;

      if(node.nfix==0) {
	stats->root_bound = bound;
	std::cout << "[set root bound]" << std::endl;
      }

      std::cout << "QPBp \t\t= " << bound << std::endl;
      //std::cout << "time = " << s_bnd.elapsed()/100.0 << std::endl;

#ifdef HAVE_EVB3
      delete q_evb3; q_evb3 = new QAP(*qr);
      evb3.perturb(qr,q_evb3);
      init_X(X,m);
      bound = qpbnd_parametric(q_evb3->A,q_evb3->B,q_evb3->C,X,U,qpbp);
      std::cout << "QPEVB3 \t\t= " << bound << std::endl;
#endif

      /*
      m_output(X);
      find_best_match(X,p);
      nice_int_array_print(p,X->m);
      std::cout <<	QAPobjective(qr,p,X->m) << std::endl;
      */
      bound =  pbnd(qr->A,qr->B,qr->C,qr->shift);
      std::cout << "PB  \t\t= " << bound << std::endl;

      bound =  glbnd(qr->A,qr->B,qr->C,U,qr->shift);
      std::cout << "GLB \t\t= " << bound << std::endl;

      break;
      
    case 'd':     case 'D': 
#ifdef USE_J_MATRIX
      std::cout << "row> ";
      std::cin >> i;
      std::cout << "col> ";
      std::cin >> j;
      node.disallow(i,j);
      //std::cout.form("%d->%d has been disallowed\n",i,j);
      printf("%d-->%d has been disallowed\n",i,j);
#else
      std::cout << "this version does not support that operation" << std::endl;
#endif
      break;
      
    case 'f':     case 'F':
      fix(c, node);
      break;
      
    case 'g':    case 'G':
      run_at_node(q, node, sln, stats, params, c);
      break;

    case 'h':    case 'H':

      // uses: q, node, sln, qpbp
      
      qr = reduceQAP(q, node);
      m = qr->A->m;
      X = m_resize(X,m,m);
      U = m_resize(U,m,m);
      init_X(X, m);

      qpbp->inc = sln->getObjective();
      qpbp->shift = qr->shift;     
      bound = qpbnd(qr->A, qr->B, qr->C, X, U, qpbp);

      p = qpb_heuristic(X);
      heuristic_obj = QAPobjective(qr, p, m);
      printf("HEURISTIC = %f\n", heuristic_obj);
      break;
      
    case 'l':    case 'L':
      std::cout << "current depth is " << node.depth << std::endl;
      i = get_int_value("depth");
      if (i>=0) node.depth = i;
      log_int(c, "l", i);
      break;

    case 'm': case 'M':
      printf("%s\n", c->input_file);
      qr = reduceQAP(q, node);
      //      branch_training_data(qr, node, params, sln, stats, c, qpbp);
      branch_training_data(q, node, params, sln, stats, c, qpbp);
      break;
      
    case ' ':     case 'c':    case 'C':
      node = init;
      log_string(c, "c");
      break;
      
    case 'o':    case 'O':
      std::cout << "filename> ";
      std::cin >> input;
      
      if (!strcmp(input,"")) {
	std::cout << "(no file given)" << std::endl;
	break;
      }

      delete qr;
      qr = reduceQAP(q,node);

      writeQAPLIB(input,qr);
      std::cout << "Wrote QAP to [" << input << "]" << std::endl;

      break;

    case 'p':    case 'P':
      change_param(node, param, sln, stats);
      break;
      
    case 'r': case 'R':
      delete qr;
      qr = reduceQAP(q,node);
      m = qr->A->m;

      X = m_resize(X,m,m);
      U = m_resize(U,m,m);

      param = pick_param(node,stats->root_bound,sln->getObjective(),
			 params,strat,rel_gap, c->nlevel);

      /* bound computation */

      qpbp->shift = qr->shift;
      qpbp->inc = sln->getObjective();

      qpbp->maxier = param->maxier_bnd;
      qpbp->maxier_nofathom = param->maxier_nofathom;
      qpbp->step = param->param_step;
      qpbp->fw_iter = 0;
      qpbp->improve = false;

      init_X(X,m);

      bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);

      node.bound = std::max(std::max(bound, node.bound), node.predicted_bound);

      if(node.nfix==0) {
        stats->root_bound = bound;
        std::cout << "[set root bound]" << std::endl;
      }
      std::cout << "bound = " << bound << ", upper bnd = " << sln->obj << std::endl;

      /* branching decision */

      qpbp->maxier = param->maxier_br;
      qpbp->maxier_nofathom = param->maxier_br;
      qpbp->step = param->param_step;
      qpbp->save_best = true;
      qpbp->improve   = false;

      set_branch_parameters(bp, param, false, false);
      branch(q,node,X,U,kid,kids,bound,qpbp,bp,loc_in_row,loc_kept);

      //      m_output(U);
      nice_matrix_print(U->me, U->m, U->n);
      //std::cout.form("strat %d, rule = %d, %d kids created, (%d/%d) locations kept.\n",
      printf("strat %d, rule = %d, %d kids created, (%d/%d) locations kept.\n",
		strat,param->rule, kids, loc_kept, loc_in_row);

      for(i=0; i<kids; i++) {
	std::cout << "Child #" << i+1 << std::endl;
	kid[i].print();
      }

      i = get_int_value("child");
      if (i && (i>0) && (i<=kids)) {
	node = kid[i-1];
      }
      else {
	std::cout << "(leaving problem as it was)" << std::endl;
      }

      break;

    case 's': case 'S':
      bool zero_based;
      std::cout << "Enter permutation (1-based): ";

      for(i=0;i<n;i++) {
	std::cin >> p[i];
      }

      zero_based = false;
      for(i=0;i<n;i++) {
	if(p[i]	== 0) zero_based = true;
      }

      if (!zero_based) {
	for(i=0;i<n;i++) p[i]--;
      }
      
      std::cout << std::endl << "value = " << QAPobjective(q,p,n) << std::endl;
      break;
      
    case 't': case 'T':
      run_all_branches(q, node, sln, stats, params, c, "0");
      break;

    case 'u':    case 'U':
      std::cout << "(using parametric bound)" << std::endl;

      delete qr;
      qr = reduceQAP(q,node);
      m = qr->A->m;

      X = m_resize(X,m,m);
      U = m_resize(U,m,m);

      init_X(X,m);

      qpbp->shift = qr->shift;

      bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);

      std::cout << "bound: " << bound << std::endl;

      std::cout << "U = " << std::endl;
      for(i=0; i<U->m; i++) {
        for(j=0; j< U->n; j++)
            //std::cout.form("%5.2f ",U->me[i][j]);
            printf("%5.2f ",U->me[i][j]);
        std::cout << std::endl;
      }

      std::cout << "X = " << std::endl;
      for(i=0; i<X->m; i++) {
	for(j=0; j< X->n; j++)
	  //std::cout.form("%5.2f ",X->me[i][j]);
	  printf("%5.2f ",X->me[i][j]);
	std::cout << std::endl;
      }
      VEC *var;
      var = my_variance(X);

      nice_array_print(var->ve,var->dim);
      log_matrix(params, "U", U);
      log_matrix(params, "X", X);
      break;
      
    case 'v': case 'V':
      qpb_test(q);
      break;
      
    case 'w': case 'W':

      std::cout << "row> ";
      std::cin >> i;
      std::cout << "col> ";
      std::cin >> j;


      delete qr;
      qr = reduceQAP(q,node);
      m = qr->A->m;

      X = m_resize(X,m,m);
      X1 = m_resize(X1,m-1,m-1);
      U = m_resize(U,m,m);

      init_X(X,m);
      bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);
      std::cout << "QPBp \t\t= " << bound << std::endl;

      node.fix(i,j);

      delete qr;
      qr = reduceQAP(q,node);
      qpbp->shift = qr->shift;


      U = m_resize(U,m-1,m-1);
      init_X(X1,m-1);
      bound = qpbnd_parametric(qr->A,qr->B,qr->C,X1,U,qpbp);
      std::cout << "QPBp \t\t= " << bound << std::endl;

      warmstart(X,X1,i,j);
      qpbp->feasible = false;
      bound = qpbnd_parametric(qr->A,qr->B,qr->C,X1,U,qpbp);
      std::cout << "QPBpw \t\t= " << bound << std::endl;

      break;

    case 'q':    case 'Q':
      done = true;
      break;
    default:
      break;
    }
  } while(!done);

  delete [] p;
  M_FREE(X);
  M_FREE(X1);
  M_FREE(U);

}


