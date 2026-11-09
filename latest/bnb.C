#include <cstdlib>
#include "util.h"
#include "linalg.h"
#include "lap.h"
#include "bnb.h"
#include "stats.h"
#include "solution.h"
#include "qpbnd.h"
#include "branch.h"

#ifdef HAVE_EVB3
#include "evb3.h"
#endif

#include <fstream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include "state.h"
#include <csignal>

// Saving state (driver flags -o, -k, -t, -r; see state.h).
const char *bnb_state_file = NULL;   // where to write state (-o)
double bnb_checkpoint_every = 0;     // CPU seconds between checkpoints (-k); 0 = never
double bnb_time_limit = 0;           // stop after this many CPU seconds (-t); 0 = never
bool bnb_resumed = false;            // stats were loaded from a state file (-r)
bool bnb_stopped = false;            // set when bnb() stops on the time limit
// Set by SIGUSR1 (see driver.C): save state and stop at the next node, like -t.
// bench/parallel.py uses it to take work back from a busy worker.
volatile sig_atomic_t bnb_stop_requested = 0;
// -U FILE (docs/distributed/plan.md 5.2): an incumbent published elsewhere,
// "obj [p0 .. p(n-1)]" in the state file's convention, re-read at node
// boundaries every bnb_incumbent_every nodes or at once after SIGUSR2, and
// adopted when lower. Any attained objective is a safe cutoff at any time.
const char *bnb_incumbent_file = NULL;
int bnb_incumbent_every = 50;
volatile sig_atomic_t bnb_incumbent_poke = 0;

static void read_incumbent_file(QAPSolution *s, int n)
{
  FILE *f = fopen(bnb_incumbent_file, "r");
  if (!f) return;
  double obj;
  if (fscanf(f, "%lf", &obj) == 1 && obj < s->getObjective()) {
    std::vector<int> p(n), seen(n, 0);
    bool perm = true;
    for (int i = 0; i < n && perm; i++)
      perm = fscanf(f, "%d", &p[i]) == 1 && p[i] >= 0 && p[i] < n && !seen[p[i]]++;
    s->setObjective(obj);
    if (perm) s->setPermutation(p.data());
    std::cout << "incumbent from file: " << std::setprecision(17) << obj
              << std::setprecision(6) << (perm ? "" : " (no permutation)") << std::endl;
  }
  fclose(f);
}

// DNN bound (-x dmin:dmax[:maxgap[:hopeless[:warm]]]): applied to nodes QPB does not
// fathom at depths dmin..dmax whose relative gap (inc - QPB)/(inc - root QPB)
// is at most maxgap. Off while dnn_max_depth < 0.
int dnn_min_depth = 0, dnn_max_depth = -1;
double dnn_max_gap = 1.0;
int dnn_hopeless_after = 0;          // DNNParams::hopeless_after; 0 = never give up
int dnn_warm = 0;                    // start a child's DNN from its parent's (1) or not (0)
#ifdef HAVE_DNN
#include "dnn.h"
#include "dnncut.h"
#include <unistd.h>
static void free_dnn_warm(void *p) { delete (DNNWarm *)p; }
#endif
#include <math.h>
//#include <time.h>

void list_clean_up(Node *&todo, double val, QAPStatistics *s);

//extern stopwatch s_init,s_fw,s_update,s_lap,s_dir;

#ifdef HAVE_EVB3
void use_EVB3(QAP *&q)
{
  QAP *q_evb3;
  EVB3 evb3;

  q_evb3 = new QAP(*q);
  evb3.set_fm(-750);
  evb3.set_maxit(100);
  evb3.perturb(q,q_evb3);
  delete q;
  q = q_evb3;
}
#endif

/*
============================================================

 solveSmallQAP(q,a,obj,p);

 Like solveSize3QAP, for nodes with fewer than three unassigned
 facilities. The search never creates such nodes itself (size-3
 nodes are enumerated), but the root of an n <= 2 problem is
 one, and the bound code cannot handle it.

============================================================
 */
void solveSmallQAP(QAP *q,const QAPAssignment &a,double &obj_save,
		   int *p_save)
{
  int rows[3], cols[3], m = 0, k = 0;
  Index *p = copy_array<Index>(a.p,a.n);

  for(int i=0;i<a.n;i++) if(a.p[i] < 0)    rows[m++] = i;
  for(int j=0;j<a.n;j++) if(a.pinv[j] < 0) cols[k++] = j;

  obj_save = 1e100;
  std::sort(cols, cols + m);
  do {
    for(int i=0;i<m;i++) p[rows[i]] = cols[i];
    double obj = QAPobjective(q,p,a.n);
    if(obj < obj_save) {
      obj_save = obj;
      for(int i=0;i<a.n;i++) p_save[i] = p[i];
    }
  } while(std::next_permutation(cols, cols + m));
  delete [] p;
}

/*
============================================================

 solveSize3QAP(q,a,obj,p);

 Assume that a corresponds to a QAPAssignment where all 
 facilities have been assigned to locations, except three.
 Then there are six different assignments that could be
 made to the remaining variables.

 solveSize2QAP makes the six remaining assignments, compares
 objective values, and returns the permutation corresponding
 to the smallest objective value in p.

============================================================
 */
void solveSize3QAP(QAP *q,const QAPAssignment &a,double &obj_save,
			  int *p_save)
{
  Index *p;
  double obj[6];
  int i;
  int r1,r2,r3,c1,c2,c3;

  p = copy_array<Index>(a.p,a.n);

  for(i=0;(i<a.n)&&(a.p[i]>=0);i++);
  r1 = i;
  for(i++;(i<a.n)&&(a.p[i]>=0);i++);
  r2 = i;
  for(i++;(i<a.n)&&(a.p[i]>=0);i++);
  r3 = i;

  for(i=0;(i<a.n)&&(a.pinv[i]>=0);i++);
  c1 = i;
  for(i++;(i<a.n)&&(a.pinv[i]>=0);i++);
  c2 = i;
  for(i++;(i<a.n)&&(a.pinv[i]>=0);i++);
  c3 = i;

  p[r1] = c1;  p[r2] = c2;  p[r3] = c3;
  obj[0] = QAPobjective(q,p,a.n);
  p[r1] = c1;  p[r2] = c3;  p[r3] = c2;
  obj[1] = QAPobjective(q,p,a.n);

  p[r1] = c2;  p[r2] = c1;  p[r3] = c3;
  obj[2] = QAPobjective(q,p,a.n);
  p[r1] = c2;  p[r2] = c3;  p[r3] = c1;
  obj[3] = QAPobjective(q,p,a.n);

  p[r1] = c3;  p[r2] = c1;  p[r3] = c2;
  obj[4] = QAPobjective(q,p,a.n);
  p[r1] = c3;  p[r2] = c2;  p[r3] = c1;
  obj[5] = QAPobjective(q,p,a.n);

  int best = 0;
  for(i=1;i<6;i++) 
    {
      if(obj[i] < obj[best])
	best = i;
    }

  obj_save = obj[best];
  for(i=0;i<a.n;i++)
    p_save[i] = p[i];

  // yuck
  switch(best) {
  case 0:
    p_save[r1] = c1;  p_save[r2] = c2;  p_save[r3] = c3;  
    break;
  case 1:
    p_save[r1] = c1;  p_save[r2] = c3;  p_save[r3] = c2;
    break;
  case 2:
    p_save[r1] = c2;  p_save[r2] = c1;  p_save[r3] = c3;
    break;
  case 3:
    p_save[r1] = c2;  p_save[r2] = c3;  p_save[r3] = c1;
    break;
  case 4:
    p_save[r1] = c3;  p_save[r2] = c1;  p_save[r3] = c2;
    break;
  case 5:
    p_save[r1] = c3;  p_save[r2] = c2;  p_save[r3] = c1;
    break;
  }


  delete [] p;
  return;

}


/*
============================================================

  param = pick_param(...);

  select the appropriate branching strategy for this node and return
  a pointer to the corresponding BNBParameters object.  No memory
  is allocated.

============================================================
 */
BNBParameters *pick_param(const QAPAssignment &node,double root,double opt,
			  BNBParameters *params,int &i,double &rel_gap,
			  int nlevel)
{
  double bnd;
  BNBParameters *param;

  bnd = std::max(node.bound, node.predicted_bound);

  rel_gap = computeRelativeGap(opt,root,bnd);

  for(i=0;(i<nlevel-1)&&(rel_gap <= params[i].relative_gap);i++);
  for(;(i<nlevel-1)&&(node.depth > params[i].depth);i++);

  /*
  i=0;
  while((i<nlevel)&&(rel_gap > params[i].relative_gap)||
	(node.depth < params[i].depth))
    i++;
  i--;
  */

  if(node.depth==0) i = 0;
  param = &(params[i]);

  /*

  if(1)
    std::cout.form("depth = %d  (%f-%f)/(%f-%f) = %f, rule %d\n",node.depth,opt,
	      node.predicted_bound,opt,root,
	      (opt-node.predicted_bound)/(opt-root),i);
  */
  return param;
}


/*
============================================================
 init_X_warm(node, X);

 Starting point for node's Frank-Wolfe bound: its parent's X with the
 row and column the node fixed removed (init_X_from_parent), when the
 node carries one (see WarmStart in assign.h) and differs from that
 parent by exactly one assignment; otherwise the uniform matrix.
 Fathomed nodes need far fewer iterations from there.
============================================================
 */
static void init_X_warm(const QAPAssignment &node, MAT *X)
{
  int di, dj;
  if (node.warm && (int)node.warm->X->m == (int)X->m + 1 &&
      warm_child_offsets(node, di, dj))
    init_X_from_parent(X, node.warm->X, di, dj);
  else
    init_X(X, X->m);
}

void log_matrix(BNBParameters *p, const char *data, MAT *U)
{
  if (p->log_file) {
    std::ofstream outfile;
    outfile.open(p->log_file, std::ios_base::app);
    outfile << data << "|";
    if (!p->log_prefix.empty()) {
      outfile << p->log_prefix << "|";
    }
    write_matrix_flat(outfile, U->me, U->m, U->n, '|');
    outfile << std::endl;
  }
}


/*
============================================================

 bnb(q,stack,soln,stats,max_CPU,max_bnds,params,nlevel,
     max_depth, dive);

 Solve a quadratic assignment problem q using a branch-and-
 bound algorithm.  See qap.h for a description of the 
 format of the QAP data structure q.

 soln may contain an initial feasible point to the QAP,
 along with an associated objective value.  Upon return
 from bnb(), it will contain the optimal solution to the
 QAP.

 stats contains information about the number of nodes 
 processed, see stats.h for more info.

 max_CPU: if using MW, code will return after max_CPU seconds
          of CPU time.

 max_bnds: if using MW, code will return after max_bnds
          bounds have been computed.

 params: An array of BNBParameters structures.  params[i]
         contains branching rules and parameters to be used
         for a certain section of the tree.

 nlevel: size of params

 max_depth: if nonnegative, do not process any node deeper
            than max_depth.  A list of nodes of depth max_depth
            will be returned in stack.

 dive: if true, perform a single 'dive' down the tree, i.e.
       rather than adding all children to the todo list, 
       randomly select one.  Used to estimate size of branch
       and bound tree.

============================================================
 */
void bnb(QAP *q,Node *&stack,QAPSolution *s,QAPStatistics *stats,
         int max_CPU, int max_computed_bounds,
	 BNBParameters *params, int nlevel, int max_depth,
	 bool dive, int max_returned_depth)
{
  int i,j;

  int kids, strat;
  int loc_in_row, loc_kept;
  int *p;
  int n = q->A->m, m;
  double obj, bound, rel_gap;
  bool finishing_up, finishing_up_second_recourse;
  bool use_J,reorder;

  BNBParameters *param;
  QAPAssignment a(n),kid[n],node;
  QPBParameters *qpbp;    
  BranchParameters *bp = new BranchParameters();
  
  QAP *qr;
  MAT *X,*U;

  Node *too_deep = NULL; // list of nodes that were too deep to process.

  double factor = params->est_factor;

  double val, sum_of_scores, sum_so_far, importance;

  // ************ initialization **********************
  if (max_depth < 0)
    max_depth = 100000;

  {
    QAPStatistics::ctr_type fw0 = stats->fw_iter;
    stats->start();                    // (which zeroes fw_iter)
    if (bnb_resumed) stats->fw_iter = fw0;
  }
  const double session_start = stats->total_time;
  double next_checkpoint = stats->total_time + bnb_checkpoint_every;
  bnb_stopped = false;

  qpbp = new QPBParameters();
  qpbp->compute_lapU = false;  
  qpbp->lapU = 0.0;   qpbp->lap_scale = 1.0;  
  qpbp->save_best = true;
  qpbp->feasible = true;

#ifdef BNB_STANDALONE
  reorder = false;
#else
  reorder = true;
  finishing_up = false; 
  finishing_up_second_recourse = false;  
#endif

  X = m_get(n,n);
  U = m_get(n,n);
  p = new int[n];
  qr = new QAP(n);

  // *********** branch and bound *********************
#ifdef BNB_STANDALONE
  while (stack != NULL)
#else
  while (         
	(stack != NULL) &&      
	(s->ifs_found == 0)
	)
#endif
    {
      // Checkpoints and the time limit are taken here, before the next node
      // is popped, so a saved stack never loses a node.
      if (bnb_state_file) {
	if (bnb_stop_requested ||
	    (bnb_time_limit > 0 &&
	     stats->total_time - session_start >= bnb_time_limit)) {
	  while (too_deep != NULL) {         // keep nodes set aside by -d
	    list_head_insert(stack, too_deep->data);
	    list_head_remove(too_deep);
	  }
	  write_state(bnb_state_file, stack, s, stats);
	  bnb_stopped = true;
	  break;
	}
	if (bnb_checkpoint_every > 0 && stats->total_time >= next_checkpoint) {
	  // the nodes set aside by -d belong to the saved state too (the stop
	  // path above moves them for good; here they go back afterwards so
	  // the run continues unchanged, in the same order)
	  int moved = 0;
	  for (Node *t = too_deep; t != NULL; t = t->link) {
	    list_head_insert(stack, t->data);
	    moved++;
	  }
	  write_state(bnb_state_file, stack, s, stats);
	  for (; moved > 0; moved--) list_head_remove(stack);
	  next_checkpoint = stats->total_time + bnb_checkpoint_every;
	}
      }
      if (bnb_incumbent_file) {
	static long inc_tick = 0;
	if (bnb_incumbent_poke || inc_tick++ % std::max(1, bnb_incumbent_every) == 0) {
	  bnb_incumbent_poke = 0;
	  read_incumbent_file(s, n);
	}
      }

      pop(stack, node, stats);
      stats->startIteration(node.depth);

#ifndef BNB_STANDALONE
      if(( stats->total_time >= max_CPU ) || 
         ( stats->node >= max_computed_bounds))
        finishing_up = true;
      
      if(finishing_up && (node.depth <= max_returned_depth) ) {
	push(stack,node,stats);
        break;
      }

      if (finishing_up && (stats->total_time >= 2.5*max_CPU) ){
        finishing_up_second_recourse = true;
      }
      
      if(finishing_up_second_recourse && 
	 ((node.depth <=  (max_returned_depth + 1) )
	  || (stats->total_time >= 3.5*max_CPU))) {
	push(stack,node,stats);
        break;
      }
#endif

      if (node.depth > max_depth) {
	if (save_too_deep) {
	    list_head_insert(too_deep,node);
	}
	continue;
      }

      // It may happen that the incumbent value was updated
      // since this node was put on the stack.  Check to see if the
      // new incumbent value is low enough to fathom this node.
      if (node.bound > s->getObjective() - DELTA) {
	//std::cout << "Fathomed due to updated incumbent (" 
	//     << s->getObjective() << ")" << std::endl;
	stats->anotherFathom(node.depth);
	continue;
      }

      // determine how far down the tree we are, and choose
      // the appropriate parameters object.
      param = pick_param(node,stats->root_bound,s->getObjective(),
			 params,strat,rel_gap, nlevel);

      if (!dive) {
	stats->updateRelativeGap(node.depth,rel_gap);
	stats->anotherNode(node.depth,strat);
      }
      else {
	/*	std::cout << "depth = " << node.depth
	     << ", rel gap = " << rel_gap 
	     << ", factor = " << factor << std::endl;
	*/
	stats->updateRelativeGap(node.depth,rel_gap,factor);
	stats->anotherNode(node.depth,strat,factor);
      }

      if (stats->node % PRINT_MESSAGE_EVERY == 0) {
	print_time();
        std::cout << "nodes = " << stats->node << "  "
             << "stack size = " << stats->stack_size << "  "
             << "cpu time = " << stats->total_time << "  "
             << std::endl;
      }

      /*
      if(node.depth >= 0) {
	std::cout << "*** Node #" << stats->node 
		  << ", Level = " << node.n - node.nfix 
		  << ", Depth = " << node.depth
		  << "  stack size = " << stats->stack_size 
		  << std::endl;
	node.print(std::cout,false);
      }
      */

      // If the assignment corresponds to a permutation, there is no
      // need to compute a bound: just compute the objective value
      // corresponding to the permutation, update, and continue.
      if (node.isPermutation()) {
	stats->anotherEnumeration();
	obj = QAPobjective(q,node);
	s->update(obj,node.p);
	
	continue;
      }

      if (node.n - node.nfix < 3) {
	stats->anotherEnumeration();
	solveSmallQAP(q,node,obj,p);
	s->update(obj,p);
	continue;
      }

      // If the assignment has only three remaining unassigned variables,
      // then we can enumerate all six possible solutions and pick the one
      // with the best objective value.
      if (node.n - node.nfix == 3) {
	stats->anotherEnumeration();
	solveSize3QAP(q,node,obj,p); 
	s->update(obj,p);
	/*
	  node.print(std::cout,false);
	  std::cout << "  ";	  write_permutation(std::cout,p,n);	  std::cout << std::endl;
	*/
	continue;
      }

      // If k assignments have already been made, then we can
      // reduce the problem of computing a bound on (A,B,C) to
      // the problem of computing a bound on (A',B',C') where
      // (A',B',C') are (n-k) x (n-k) matrices.

      reduceQAP(q, qr, node);
#ifdef HAVE_EVB3
      if(param->which_bound == EVB3_BND) {
	use_EVB3(qr);
      }
#endif

      m  = qr->A->m;             // m = size of reduced QAP = n - k
      X = m_resize(X,m,m);
      U = m_resize(U,m,m);

      // ************* bound computation *************

      qpbp->shift = qr->shift;
      qpbp->inc = s->getObjective();

      qpbp->maxier = param->maxier_bnd;
      qpbp->maxier_nofathom = param->maxier_nofathom;
      // Patience only where the rule scores children with their own bounds
      // (strong, lookahead). With rule 1 or 2, U is the only thing that fixes
      // children out, so those nodes keep their full iteration budget. Never at
      // the root: its bound sets the relative gaps that pick_param uses for
      // every node, and it is a single node.
      qpbp->nofathom_patience =
	(node.depth > 0 &&
	 (param->rule == BRANCH_LOOKAHEAD || param->rule == BRANCH_LOOKAHEAD_ECC ||
	  param->rule == BRANCH_STRONG))
	? fw_nofathom_patience : 0;
      qpbp->step = param->param_step;
      qpbp->fw_iter = 0;
      qpbp->improve = false;

      use_J = false;
      switch(param->which_bound) {

      case PB_BND:
	init_X(X,m);
	bound = qpbnd(qr->A,qr->B,qr->C,X,U,qpbp);
	break;

      case QP_BND_IMP:
	qpbp->improve = true;    // drop thru to QP_BND_PARAM
      case QP_BND_PARAM:
      case EVB3_BND:
	init_X_warm(node, X);
	bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);

	break;
      case QP_BND:

	init_X(X,m);
	bound = qpbnd(qr->A,qr->B,qr->C,X,U,qpbp);
	break;
      case QP_GLB_BND:

	init_X(X,m);
	bound = qp_glb_bnd(qr->A,qr->B,qr->C,X,U,qpbp);
	break;

      case GLB_BND:

	use_J = true;
	bound = glbnd(qr->A,qr->B,qr->C,U,qr->shift);
	break;
      }


      if (params->log_file) {
	log_matrix(params, "U|_", U);
	log_matrix(params, "X|_", X);
      }

      /*
      if((node.depth >= 8)&&(node.p[0]==11)) {
	printf("%d Bnd = %.2f, Parent = %.2f, Predicted = %.2f, Depth = %d\n",
	       strat,bound,node.bound,node.predicted_bound,node.depth);
      }
      */
      node.bound = std::max(std::max(bound, node.bound), node.predicted_bound);
      if (node.depth == 0) stats->root_bound = bound;

      // ************** DNN bound (-x) **************
      // A much stronger, much more expensive bound for nodes near the top
      // that QPB cannot fathom (see dnn.h, experiments/dnn_bound). It is used
      // only to fathom and in node.bound: branch() below still gets the QPB
      // bound, since bound + U_ij is only valid with the QPB value and U.
      // (node.bound also feeds pick_param, the children's inherited bounds and
      // some branching scores, so the search differs from a QPB-only run even
      // where the DNN bound does not fathom.) Skipped until there is a real
      // incumbent: against 1e100 nothing fathoms and every call runs to its limit.
      bool dnn_fathomed = false, dnn_ran = false;
      void *dnn_keep = NULL;   // this node's DNN state, for its children
      std::vector<double> dnn_xbar;   // its primal marginals, for rule 8
      std::string dnn_pool_out;       // its cut pool for its children (QAP_DNN_POOL)
      // QAP_DNN_EXT=k: below dnn_max_depth a node whose parent's DNN bound ran
      // and did not fathom keeps the DNN bound, up to k levels down; under a
      // fathomed or a QPB-only node nothing changes. The flag travels with the
      // parent's warm start, so not across a checkpoint restore. Without it a
      // last-level miss hands its whole subtree to QPB alone (580,000 nodes
      // under two depth-7 misses in the laptop nug30 run; Goux's D8).
      static const int dnn_ext = getenv("QAP_DNN_EXT") ? atoi(getenv("QAP_DNN_EXT")) : 0;
      const int dnn_ext_in = node.warm ? node.warm->dnn_ext : 0;
#ifdef HAVE_DNN
      if (bound <= s->getObjective() - DELTA && m >= 2 && s->getObjective() < 1e99 &&
	  node.depth >= dnn_min_depth && (node.depth <= dnn_max_depth || dnn_ext_in > 0)) {
	const double inc = s->getObjective();
	const double gap = (inc - bound) / std::max(inc - stats->root_bound, 1e-9);
	if (gap <= dnn_max_gap) {
	  DNNParams dp;
	  // giving up early on a node that would fathom costs its whole subtree
	  dp.hopeless_after = dnn_hopeless_after;
	  // prototype switches (experiments/eigensolver): Peaceman-Rachford first
	  // dual step and Chebyshev filter degree of the subspace eigensolver
	  if (const char *e = getenv("QAP_DNN_GAMMA1")) dp.gamma1 = atof(e);
	  if (const char *e = getenv("QAP_DNN_GAMMA")) dp.gamma = atof(e);
	  if (const char *e = getenv("QAP_DNN_CHEB")) dp.cheb_degree = atoi(e);
	  // tight lam_hi for the filter (experiments/wil50_feasibility/solvers/accel):
	  // with QAP_DNN_CHEB=3 the strength of degree 4 at 85% of the cost
	  if (const char *e = getenv("QAP_DNN_TIGHT_LAM")) dp.tight_lam = atoi(e);
	  // gap-ratio give-up rule (experiments/wil50_feasibility/tuning), opt-in
	  if (const char *e = getenv("QAP_DNN_GAP_RATIO")) dp.hopeless_gap_ratio = atof(e);
	  // a give-up within this fraction of the cutoff must be confirmed 100
	  // iterations later (Goux's D6), opt-in
	  if (const char *e = getenv("QAP_DNN_GU_ZONE")) dp.hopeless_zone = atof(e);
	  // fast ADMM momentum with Goux's safeguards (restart factor), opt-in
	  if (const char *e = getenv("QAP_DNN_MOMENTUM")) dp.momentum = atof(e);
	  // ADMM penalty factor (0.3 = default; Goux: x 0.7 with momentum)
	  if (const char *e = getenv("QAP_DNN_BETA")) dp.beta = atof(e);
	  // full eigendecomposition every this many iterations (20; Goux: 120)
	  if (const char *e = getenv("QAP_DNN_REFRESH")) dp.refresh = atoi(e);
	  // cut phase (QAP_DNN_CUT_ZONE, experiments/wil50_feasibility/stronger_bounds):
	  // when the uncut bound ends within this fraction of |incumbent| of the
	  // cutoff, regional triangle cuts from the node's own iterate, rounds of
	  // QAP_DNN_CUT_ITERS iterations with QAP_DNN_CUT_K new cuts each, Polyak
	  // step toward best + QAP_DNN_CUT_DELTA, stopping at the cutoff. Regions:
	  // the prefixes of the QAP_DNN_CUT_GRID=R:C grid over the original indices
	  // of side QAP_DNN_CUT_SIDE (0: the first matrix; wil50 5:10, nug30 5:6
	  // on side 1). Every value is certified; the children keep the uncut
	  // multiplier.
	  static const double cut_zone = getenv("QAP_DNN_CUT_ZONE") ? atof(getenv("QAP_DNN_CUT_ZONE")) : 0;
	  // QAP_DNN_CUT_MAXDEPTH: run the cut phase only at nodes of at most this
	  // depth (default: every DNN depth); on nug30 the depth-2 phases pay and
	  // the depth-3 ones do not (experiments/wil50_feasibility/tuning, refine/)
	  static const int cut_max_depth = getenv("QAP_DNN_CUT_MAXDEPTH") ? atoi(getenv("QAP_DNN_CUT_MAXDEPTH")) : 1 << 30;
	  // warm start from the parent's final multiplier, when the parent's
	  // DNN bound ran and missed and this node fixes one more facility
	  DNNWarm start, *save = (dnn_warm || param->rule == BRANCH_DNN_XBAR || param->rule == BRANCH_XBAR_LOOKAHEAD || cut_zone > 0) ? new DNNWarm : NULL;
	  const DNNWarm *sp = NULL;
	  int di, dj;
	  // QAP_DNN_WARM_Y=2 (experiments/wil50_feasibility/warm_y): the child also
	  // starts from the parent's final Y conditioned on its new assignment (1:
	  // restricted); 0 (default) the parent's Z only. Falls back to Z only when
	  // the parent's Y is missing.
	  static const int warm_y = getenv("QAP_DNN_WARM_Y") ? atoi(getenv("QAP_DNN_WARM_Y")) : 0;
	  if (dnn_warm && node.warm && node.warm->dnn &&
	      ((DNNWarm *)node.warm->dnn)->m == m + 1 && warm_child_offsets(node, di, dj)) {
	    const DNNWarm &pw = *(DNNWarm *)node.warm->dnn;
	    if (warm_y > 0 && dnn_warm_child(pw, di, dj, start, warm_y) && !start.Y.empty()) {
	      sp = &start;
	      dp.warm_primal = 1;
	    } else if (dnn_warm_child(pw, di, dj, start)) {
	      sp = &start;
	    }
	  }
	  // QAP_DNN_POOL=1 (experiments/pool_inheritance, wil50_feasibility/
	  // pool_search): a node whose cut phase ran hands its cut pool (the
	  // cuts with a positive multiplier in the solve its children start from,
	  // at those multipliers) to its children with that solve's state. A child
	  // with a pool solves with it frozen in place of the uncut solve (same
	  // iteration limit and give-up rule); its cut phase, if it runs, starts
	  // from the pool and that solve, and hands on its own pool; without a cut
	  // phase the child hands on the inherited pool. The pool travels in state
	  // files (the DNN state does not). QAP_DNN_POOL_CHECK=K checks each
	  // assembled term on K random permutations; QAP_DNN_TRACE=file appends a
	  // line per DNN node.
	  static const int dnn_pool = getenv("QAP_DNN_POOL") ? atoi(getenv("QAP_DNN_POOL")) : 0;
	  static const int pool_check = getenv("QAP_DNN_POOL_CHECK") ? atoi(getenv("QAP_DNN_POOL_CHECK")) : 0;
	  static const char *dnn_trace = getenv("QAP_DNN_TRACE");
	  static const int cut_side = getenv("QAP_DNN_CUT_SIDE") ? atoi(getenv("QAP_DNN_CUT_SIDE")) : 0;
	  std::vector<int> free_orig, other_orig, fixed_at(n, -1);
	  for (int a = 0; a < n; a++) {
	    if ((cut_side == 0 ? node.p[a] : node.pinv[a]) == NONE) free_orig.push_back(a);
	    if ((cut_side == 0 ? node.pinv[a] : node.p[a]) == NONE) other_orig.push_back(a);
	    else fixed_at[a] = cut_side == 0 ? node.pinv[a] : node.p[a];
	  }
	  DNNPool pool_in, pool_mid;
	  bool have_pool = dnn_pool && node.warm && !node.warm->dnn_pool.empty() &&
	    dnn_pool_decode(node.warm->dnn_pool, pool_in) && pool_in.side == cut_side && !pool_in.cuts.empty();
	  int pool_kept = 0, pool_dropped = 0;
	  // QAP_DNN_POOL_CP=1 (roadmap/admm/cpadmm, "In the search"): a pooled child's
	  // first solve in CP form (dnn_cp_bound with DNNCPParams::child): the inherited
	  // pool inside the PRSM's Y-step from nu = mu / (scale beta) instead of frozen
	  // at mu in the objective; it hands on (and a cut phase starts from) the cuts
	  // with nu > 0 at its final mu_D. Either side; off by default
	  static const int pool_cp = getenv("QAP_DNN_POOL_CP") ? atoi(getenv("QAP_DNN_POOL_CP")) : 0;
	  bool pool_cp_ran = false;
	  DNNResult r;
	  if (have_pool) {
	    DNNCutParams pc;
	    pc.side = cut_side; pc.pool_in = &pool_in; pc.pool_out = &pool_mid;
	    pc.other_orig = &other_orig; pc.fixed_at = &fixed_at; pc.pool_check = pool_check;
	    if (pool_cp) {
	      DNNCPParams cpp;
	      cpp.child = 1; cpp.nu_from_pool = 1;
	      r = dnn_cp_bound(qr->A, qr->B, qr->C, qr->shift, inc - DELTA, dp, free_orig, n, pc, cpp, sp, save,
	                       &pool_kept, &pool_dropped);
	      pool_cp_ran = pool_kept > 0 && cpp.pool_final > 0;   // not when it ran dnn_pool_bound instead
	    } else
	      r = dnn_pool_bound(qr->A, qr->B, qr->C, qr->shift, inc - DELTA, dp, free_orig, n, pc, sp, save,
	                         &pool_kept, &pool_dropped);
	    if (pool_kept == 0) have_pool = false;
	  } else
	    r = dnn_bound(qr->A, qr->B, qr->C, qr->shift, inc - DELTA, dp, sp, save);
	  dp.warm_primal = 0;   // the cut phase below sets its own
	  DNNWarm *cut_state = NULL;   // the cut phase's state for the children (QAP_DNN_POOL)
	  DNNCutResult cr_trace;
	  bool cut_ran = false, cut_cp = false;
	  dnn_ran = true;
	  stats->dnn_nodes++;
	  stats->dnn_time += r.cpu;
	  node.bound = std::max(node.bound, r.lb);
	  if (r.status == DNN_FATHOMED) {
	    dnn_fathomed = true;
	    stats->dnn_fathoms++;
	  }
	  if (!dnn_fathomed && cut_zone > 0 && node.depth <= cut_max_depth && save && save->m == m && r.lb > -HUGE_VAL &&
	      (inc - DELTA) - r.lb <= cut_zone * std::max(1.0, fabs(inc))) {
	    DNNCutParams cp;
	    if (const char *e = getenv("QAP_DNN_CUT_GRID")) sscanf(e, "%d:%d", &cp.grid_rows, &cp.grid_cols);
	    if (const char *e = getenv("QAP_DNN_CUT_SIDE")) cp.side = atoi(e);
	    if (const char *e = getenv("QAP_DNN_CUT_ROUNDS")) cp.rounds = atoi(e);
	    if (const char *e = getenv("QAP_DNN_CUT_ITERS")) cp.iters = atoi(e);
	    if (const char *e = getenv("QAP_DNN_CUT_K")) cp.kcuts = atoi(e);
	    if (const char *e = getenv("QAP_DNN_CUT_K5")) cp.kcuts5 = atoi(e);   // 5-point count inequalities per round (0)
	    if (const char *e = getenv("QAP_DNN_CUT_K7")) cp.kcuts7 = atoi(e);   // local-7 inequalities per round (0)
	    if (const char *e = getenv("QAP_DNN_CUT_SUB7")) cp.sub7 = atoi(e);   // random 7-subsets per region and round (1000)
	    // box (count) inequalities per round, QAP_DNN_CUT_KBOX=K[:QMAX] (0)
	    if (const char *e = getenv("QAP_DNN_CUT_KBOX")) {
	      cp.kbox = atoi(e);
	      if (const char *c = strchr(e, ':')) cp.box_qmax = atoi(c + 1);
	    }
	    if (const char *e = getenv("QAP_DNN_CUT_DELTA")) cp.delta = atof(e);
	    // stop rule: after QAP_DNN_CUT_STALL rounds (3), stop when the gap
	    // exceeds 3x the gain of the last that many rounds (0 = off)
	    cp.stall_rounds = 3;
	    if (const char *e = getenv("QAP_DNN_CUT_STALL")) cp.stall_rounds = atoi(e);
	    // primal gate (wil50_feasibility/trims): QAP_DNN_CUT_PGATE=R:M stops the
	    // phase after round R when its best primal is more than M below the
	    // cutoff (0 = off; the replay on 88 logged phases supports 4:20)
	    if (const char *e = getenv("QAP_DNN_CUT_PGATE")) {
	      cp.pgate_rounds = atoi(e);
	      if (const char *c = strchr(e, ':')) cp.pgate_margin = atof(c + 1);
	    }
	    // the dual method's two fixes (research_bound_quality.md, measured in
	    // stronger_bounds/dual): the Polyak target best + delta, uncapped, so
	    // the step does not vanish at the cutoff (QAP_DNN_CUT_NOCAP=0 restores
	    // the cap), and a quiet round of QAP_DNN_CUT_FREEZE iterations (600;
	    // 0 = off) with the multipliers frozen when the phase ends unfathomed
	    // with its primal above the cutoff, so the certificate catches up
	    cp.cap_at_cutoff = 0;
	    cp.freeze_iters = 600; cp.freeze_if_primal_above = 1;
	    if (const char *e = getenv("QAP_DNN_CUT_NOCAP")) cp.cap_at_cutoff = atoi(e) ? 0 : 1;
	    if (const char *e = getenv("QAP_DNN_CUT_FREEZE")) cp.freeze_iters = atoi(e);
	    // QAP_DNN_CUT_WARMY=1 (roadmap/convergence2): every cut round and the
	    // quiet round restart from the previous round's Y as well as its Z
	    // (dnncuts -warmY; same_node_continuation). Off by default: today each
	    // round restarts from Y = e0 e0', and the certificate spends the first
	    // part of every round rebuilding it
	    if (const char *e = getenv("QAP_DNN_CUT_WARMY")) cp.warm_primal = atoi(e) ? 1 : 0;
	    // QAP_DNN_CUT_DAMP1=F[:G]: the first Polyak step of a phase that
	    // starts from an inherited pool is F times the usual one, its later
	    // steps G times (1 = off; pipeline_ab/README.md)
	    if (const char *e = getenv("QAP_DNN_CUT_DAMP1")) {
	      cp.damp1 = atof(e);
	      if (const char *c = strchr(e, ':')) cp.damp_rest = atof(c + 1);
	    }
	    // QAP_DNN_CUT_EQUIP=G:R (tail_probe/README.md, "equip phases"): a phase
	    // whose uncut value is more than G below the cutoff cannot fathom (no
	    // phase has gained more than 143) and runs only to equip the children
	    // with a pool: at most R Polyak rounds (default 3). Off by default.
	    if (const char *e = getenv("QAP_DNN_CUT_EQUIP")) {
	      double gap_min = atof(e); int eq_rounds = 3;
	      if (const char *c = strchr(e, ':')) eq_rounds = atoi(c + 1);
	      if (gap_min > 0 && eq_rounds > 0 && (inc - DELTA) - r.lb > gap_min && eq_rounds < cp.rounds) cp.rounds = eq_rounds;
	    }
	    if (cp.grid_rows == 0) cp.prefix = true;
	    DNNPool pool_out;
	    if (dnn_pool && cp.side == cut_side) {
	      cp.pool_in = have_pool ? (pool_cp_ran ? &pool_mid : &pool_in) : NULL;   // QAP_DNN_POOL_CP: at mu_D
	      cp.pool_out = &pool_out;
	      cp.other_orig = &other_orig; cp.fixed_at = &fixed_at; cp.pool_check = pool_check;
	    }
	    // QAP_DNN_CUT_CP=k:K (roadmap/admm/cpadmm): the phase as one CP-ADMM solve
	    // (dnn_cut_cp) instead of the rounds: the pool inside the PRSM's Y-step,
	    // K new triangles separated every k iterations on its own iterate, from the
	    // inherited pool (its multipliers) plus round 0's QAP_DNN_CUT_K cuts on the
	    // first solve's Y, warm from its Z (and Y under QAP_DNN_CUT_WARMY), at most
	    // QAP_DNN_CUT_CP_ITERS iterations (default the rounds' ROUNDS x ITERS +
	    // FREEZE), stopped by the stall rule per check over STALL x ITERS
	    // iterations. QAP_DNN_CUT_CP_SWEEPS=S[:P]: Dykstra sweeps and passes per
	    // iteration (1:1). Off by default; regional triangles on either side
	    // (atomic cuts would run the rounds; bnbdrv never sets them)
	    if (const char *e = getenv("QAP_DNN_CUT_CP"))
	      if (sscanf(e, "%d:%d", &cp.cp_every, &cp.cp_k) == 2 && cp.cp_every > 0 && cp.cp_k > 0) {
	        cut_cp = true;
	        if (const char *f = getenv("QAP_DNN_CUT_CP_ITERS")) cp.cp_iters = atoi(f);
	        if (const char *f = getenv("QAP_DNN_CUT_CP_SWEEPS")) {
	          cp.cp_sweeps = atoi(f);
	          if (const char *c = strchr(f, ':')) cp.cp_passes = atoi(c + 1);
	        }
	      }
	    DNNWarm *wc = new DNNWarm(*save);
	    DNNCutResult cr = cut_cp ? dnn_cut_cp(qr->A, qr->B, qr->C, qr->shift, inc - DELTA, dp,
	                                          free_orig, n, cp, wc, r.lb, r.primal)
	                             : dnn_cut_rounds(qr->A, qr->B, qr->C, qr->shift, inc - DELTA, dp,
	                                              free_orig, n, cp, wc, r.lb, r.primal);
	    cut_ran = true; cr_trace = cr;
	    if (dnn_pool && cr.status != DNN_FATHOMED && !pool_out.cuts.empty() && wc->m == m) {
	      dnn_pool_out = dnn_pool_encode(pool_out);
	      cut_state = wc;
	    } else
	      delete wc;
	    stats->dnn_cut_nodes++;
	    stats->dnn_cut_time += cr.cpu;
	    node.bound = std::max(node.bound, cr.best_lb);
	    if (cr.status == DNN_FATHOMED) {
	      dnn_fathomed = true;
	      stats->dnn_cut_fathoms++;
	    }
	  }
	  if (save && !dnn_fathomed && save->m == m &&
	      (param->rule == BRANCH_DNN_XBAR || param->rule == BRANCH_XBAR_LOOKAHEAD)) {
	    // x_ij = Y_0,(ij) of the final primal iterate (row 0 of the lifted
	    // matrix, after the leading 1), m x m row-major in local indices
	    const int n1 = m * m + 1;
	    dnn_xbar.assign(save->Y.begin() + 1, save->Y.begin() + n1);
	  }
	  if (have_pool && dnn_pool_out.empty() && !dnn_fathomed && !pool_mid.cuts.empty())
	    dnn_pool_out = dnn_pool_encode(pool_mid);   // no cut phase: hand on the inherited pool
	  if (cut_state && save && !dnn_fathomed) { *save = std::move(*cut_state); }
	  delete cut_state;
	  if (dnn_trace) {
	    if (FILE *tf = fopen(dnn_trace, "a")) {
	      fprintf(tf, "pid %d depth %d fixes", (int)getpid(), node.depth);
	      for (int a = 0; a < n; a++) if (node.p[a] != NONE) fprintf(tf, " %d:%d", a, (int)node.p[a]);
	      fprintf(tf, " | first lb %.4f primal %.4f status %s iters %d cpu %.2f warm %d pool %d kept %d dropped %d",
	              r.lb, r.primal, dnn_status_name(r.status), r.iterations, r.cpu, sp ? 1 : 0,
	              have_pool ? 1 : 0, pool_kept, pool_dropped);
	      if (pool_cp_ran) fprintf(tf, " pcp 1 handed %zu", pool_mid.cuts.size());
	      fprintf(tf, " | cut %d", cut_ran ? 1 : 0);
	      if (cut_ran) fprintf(tf, " best %.4f rounds %d iters %d cpu %.2f status %s", cr_trace.best_lb,
	                           cr_trace.rounds, cr_trace.iterations, cr_trace.cpu, dnn_status_name(cr_trace.status));
	      if (cut_ran && cut_cp) fprintf(tf, " cp %d active %d stalled %d", cr_trace.cp_form, cr_trace.active, cr_trace.stalled);
	      fprintf(tf, " | cutoff %.4f fathomed %d pool_out %zu\n", inc - DELTA, dnn_fathomed ? 1 : 0, dnn_pool_out.size());
	      fclose(tf);
	    }
	  }
	  if (save && !dnn_fathomed && save->m == m && dnn_warm) dnn_keep = save;
	  else delete save;
	}
      }
#endif
      node.set_warm(NULL);    // done with the parent's X and DNN state
    
      // ************** fathom / branch **************

      kids = 0;

      /*
      */
      if (bound  <= s->getObjective() - DELTA && !dnn_fathomed) {
	
	// kids = the number of child nodes formed,
	// kid[0]..kid[kids-1] hold the child nodes.
	
	qpbp->inc = s->getObjective();
	qpbp->maxier = param->maxier_br;
	qpbp->maxier_nofathom = param->maxier_br;
	qpbp->step = param->param_step;
	qpbp->save_best = true;
	qpbp->improve   = false;
	qpbp->nofathom_patience = 0;   // only node bounds use patience
	set_branch_parameters(bp, param, use_J, reorder);
	// every bound but GLB leaves its Frank-Wolfe point in X; the
	// children start from it (copied first: branch() may use X)
	WarmStart *warm = (param->which_bound != GLB_BND)
	  ? new_warm_start(node.p, n, X) : NULL;
#ifdef HAVE_DNN
	if (warm && dnn_keep) {
	  warm->dnn = dnn_keep;
	  warm->dnn_free = free_dnn_warm;
	  dnn_keep = NULL;
	}
	// the DNN extension below the last DNN depth: k levels from a miss at
	// that depth, one fewer per level below it
	if (warm && dnn_ran && !dnn_fathomed)
	  warm->dnn_ext = node.depth < dnn_max_depth ? 0
	    : (node.depth == dnn_max_depth ? dnn_ext : dnn_ext_in - 1);
	if (warm && !dnn_pool_out.empty()) warm->dnn_pool.swap(dnn_pool_out);
#endif
#ifdef HAVE_DNN
	branch_set_dnn_marginals(dnn_xbar.empty() ? NULL : dnn_xbar.data(), m);
#endif
	branch(q, node, X, U, kid, kids, bound + qpbp->lapU,
	       qpbp, bp, loc_in_row, loc_kept);
#ifdef HAVE_DNN
	branch_set_dnn_marginals(NULL, 0);
#endif
	for (i = 0; i < kids; i++)
	  kid[i].set_warm(warm);
	warm_start_release(warm);

	if (params->log_file && (bp->how == BRANCH_LOOKAHEAD || bp->how == BRANCH_LOOKAHEAD_ECC)) {
	  log_matrix(params, "B|_", U);
	}
	
	if (dive) {
	  stats->l_level[node.depth] += qs_int(factor*(loc_in_row));
	  stats->e_level[node.depth] += 
	    qs_int(factor*(loc_in_row - loc_kept));
	  
	  if(kids > 0) {
	    // Sampling probability of each child. Default: the QPB score of its
	    // predicted bound. QAP_EST_XBAR=a: when the node has DNN marginals
	    // (rules 8 and 9), child (i,j) gets (1-mix) x_ij^a / sum + mix/kids
	    // (QAP_EST_XBAR_MIX, default 0.1; the uniform part keeps every
	    // child's probability positive, so the estimate stays unbiased).
	    static const double est_xbar_pow =
	      getenv("QAP_EST_XBAR") ? atof(getenv("QAP_EST_XBAR")) : 0.0;
	    static const double est_xbar_mix =
	      getenv("QAP_EST_XBAR_MIX") ? atof(getenv("QAP_EST_XBAR_MIX")) : 0.1;
	    std::vector<double> prob(kids, 0.0);
	    bool by_xbar = false;
#ifdef HAVE_DNN
	    if (est_xbar_pow > 0.0 && (int)dnn_xbar.size() == m * m) {
	      // local index of a free row/column = its rank among the free ones
	      std::vector<int> lrow(n, -1), lcol(n, -1);
	      for (int r = 0, k = 0; r < n; r++) if (node.p[r] <= NONE) lrow[r] = k++;
	      for (int c = 0, k = 0; c < n; c++) if (node.pinv[c] <= NONE) lcol[c] = k++;
	      double wsum = 0.0;
	      for (int k = 0; k < kids; k++) {
		double x = 0.0;
		for (int r = 0; r < n; r++)   // the pair this child fixes
		  if (node.p[r] <= NONE && kid[k].p[r] > NONE) {
		    x = dnn_xbar[lrow[r] * m + lcol[kid[k].p[r]]];
		    break;
		  }
		prob[k] = pow(std::max(x, 0.0), est_xbar_pow);
		wsum += prob[k];
	      }
	      if (wsum > 0.0) {
		by_xbar = true;
		for (int k = 0; k < kids; k++)
		  prob[k] = (1.0 - est_xbar_mix) * prob[k] / wsum + est_xbar_mix / kids;
	      }
	    }
#endif
	    if (!by_xbar) {
	      for(sum_of_scores=0.0,i=0;i<kids;i++) {
		prob[i] = score(s->getObjective(),stats->root_bound,
				kid[i].predicted_bound,
				params->score_exponent);
		sum_of_scores += prob[i];
	      }
	      for(i=0;i<kids;i++) prob[i] /= sum_of_scores;
	    }
	    val = drand();

	    sum_so_far = 0.0; i=0;
	    do {
	      sum_so_far += prob[i];
	      if(sum_so_far >= val) break;
	      i++;
	    } while(i<kids-1);

	    importance = prob[i];

	    push(stack,kid[i],stats);
	  }
	}
	else {
	  stats->l_level[node.depth] += loc_in_row;
	  stats->e_level[node.depth] += loc_in_row - loc_kept;
	  
	  for(i=0;i<kids;i++) {
	    push(stack,kid[i],stats);
	  }
	}
	for (i = 0; i < kids; i++)
	  kid[i].set_warm(NULL);    // the stack holds the references now
      }
      else {
	if(dive)
	  stats->anotherFathom(node.depth,factor);
	else
	  stats->anotherFathom(node.depth);
      }
#ifdef HAVE_DNN
      if (dnn_keep) free_dnn_warm(dnn_keep);   // not handed to children
#endif
      
      if (dive) {
	stats->stopIteration(node.depth,strat,qpbp->fw_iter,factor);
	/*
	  std::cout << "factor: " << factor 
	  << "\tnew_factor: " << factor/importance
	  << "\tkids: " << kids
	  << "\timp: " << importance
	  << std::endl;
	*/
	
	factor = factor / importance;
      }
      else {
	stats->stopIteration(node.depth,strat,qpbp->fw_iter);
      }
    }

  // *************** shutdown *************************
  /*
  printf("init up fw lap dir [%.4f %.4f %.4f %.4f %.4f]\n",
	 s_init.elapsed(),s_update.elapsed(),
	 s_fw.elapsed()-(s_lap.elapsed()+s_dir.elapsed()),
	 s_lap.elapsed(),s_dir.elapsed());
  */

  if (dive) {
    params->est_factor = factor;
  }

  list_clean_up(stack, s->getObjective(), stats);

  if (stack == NULL) {
    stack = too_deep;
    stats->stack_size = list_length(stack);
  }

  M_FREE(X);  M_FREE(U);
  delete [] p;
  delete qr;
  delete qpbp;
  delete bp;

  stats->stop();  
}

void benchmark(QAP *q,Node *&stack,QAPSolution *s,
	       QAPStatistics *stats,
	       int max_CPU, int max_computed_bounds,
	       BNBParameters *params,int nlevel,int max_depth,
	       bool dive,int max_returned_depth)
{
  int benchmark_fixed = 7;

  for(int i=0;i<benchmark_fixed;i++) {
    stack->data.fix(i,i);
  }  
 
  s->setObjective(1e10);
  
  bnb(q,stack,s,stats,max_CPU,max_computed_bounds,params,nlevel,
      benchmark_fixed +2,false,50);

  stats->print();
  
}


/*
============================================================
 remove all nodes whose bound is >= val.
============================================================
 */
void list_clean_up(Node *&todo, double val, QAPStatistics *stats)
{
  if (todo == NULL) {
    return;
  }
  
  while ((todo != NULL) && (todo->data.bound >= val)) {
    list_head_remove(todo);
    stats->stack_size--;
  }

  if ((todo == NULL) || (todo->link == NULL)) {
    return;
  }
  
  Node *pred = todo;
  while (pred->link != NULL) {
    if(pred->link->data.bound >= val) {
      list_remove(pred);
      stats->stack_size--;
    }
    else {
      pred = pred->link;
    }
  }
}

/*
============================================================
  stack_setup(stack,id,n);

  Initializes the stack.  Normally, we initialize the stack
  to a single problem, in which all possible assignments
  are permitted.

============================================================
 */
void stack_setup(Node *&stack,QAPStatistics *stats,int &id,int n)
{
  QAPAssignment a(n);

  list_head_insert(stack,a);
  stats->stack_size = 1;
}

/*
============================================================
Below are some functions that can be used to test various
parts of the system.  These functions are not used in the
branch and bound algorithm.
============================================================

 */


/*
  return true if the set of permutations determined by p1
  is a subset of the set of permutations determined by p2.

  e.g.  [1 - -] subset of [- - -]
        [- - -] not subset of [1 - -]
        ([1 - -] might be subset of [? - -])
        [? 2 -] not a subset of [- 2 -]
        [0 - -] not a subset of [1 - -] or [? - -]

  *** FOR NOW WE WILL NOT GUARANTEE CORRECTNESS IF P1 OR P2
      CONTAINS 'PARTIAL' ASSIGNMENTS, I.E. DISALLOWED ASSIGNMENTS ***

      this means that we return false only when p1[i] has been 
      assigned a value, and p1[i] != p2[i].

 */
bool subset_permutation(int *p1, int *p2, int n)
{
  for(int i=0; i<n; i++) {

    if(p1[i]==p2[i]) continue;

    if(p2[i]>NONE)
      return false;

  }

  return true;
}

void bnb_trace(QAP *q,QAPSolution *sln,
	       BNBParameters *params,int nlevel,char *inputfile)
{
      /*

	1) read history file.-- given permutation p.
	   -- read thru file, line by line.
	   1 2 - - 4, Bound = 123.456
	   -- when p matches perm. in file, print out line.
	   -- requires *all* nodes be saved.

	2) read history file.-- given permutation p.
	   -- read thru file, line by line.
	   1 2 3 - 4    4 3 1 2 -
	   -- when p matches perm. in file, read history 
	   -- requires *leaf* nodes be saved.
       */
  int p1[sln->n];
  int p2[sln->n];

  std::cout << "Enter permutation #1: ";
  read_permutation(std::cin,p1,sln->n);

  std::cout << "Enter permutation #2: ";
  read_permutation(std::cin,p2,sln->n);

  write_permutation(std::cout,p1,sln->n);
  std::cout << "\t";
  write_permutation(std::cout,p2,sln->n);
  std::cout << std::endl;

  if(subset_permutation(p1,p2,sln->n))
    std::cout << "Match" << std::endl;
  else
    std::cout << "No match" << std::endl;
}


void print_stack(QAP *q,Node *stack,BNBParameters *params, int nlevel,
		 bool compute_bounds,double opt)
{
  Node *current;
  int len, div = 20;
  int i,n = q->A->m,m;
  double bound,root_bound,*relgap;

  QPBParameters *qpbp;    
  BNBParameters *param;
  
  QAP *qr;
  MAT *X,*U;

  // wow, this is a waste of O(n) time...
  len = list_length(stack);
  relgap = new double[len];

  // ************ initialization **********************
  qpbp = new QPBParameters;  qpbp->feasible = true;

  X = m_get(n,n);
  U = m_get(n,n);
  qr = new QAP(n);

  param = params;

  init_X(X,n);
  
  qpbp->shift = q->shift;
  qpbp->inc = 1e10;
  
  qpbp->maxier = param->maxier_bnd;
  qpbp->maxier_nofathom = param->maxier_nofathom;
  qpbp->step = param->param_step;
  
  if(param->which_bound==0)
    root_bound = qpbnd(q->A,q->B,q->C,X,U,qpbp);
  else
    root_bound = qpbnd_parametric(q->A,q->B,q->C,X,U,qpbp);
  

  for(i=0,current=stack; current != NULL; current = current->link,i++) {

    relgap[i] = (opt - current->data.predicted_bound)/(opt - root_bound);
    /*
    std::cout.form("%f %f %f\n",current->data.bound,
	      current->data.predicted_bound,relgap[i]);
    */

  }
  /*  
  ins_sort_descending(relgap, NULL, len);
  for(i=0;i<div;i++) {
    std::cout.form("%d %f\n",i,relgap[i*len/div]);
  }
  std::cout.form("%d %f\n",i,relgap[len-1]);
  */
  delete [] relgap;
}

// called by find_suboptimal()
void pop_alt(Node *&stack,QAPAssignment &a)
{
  if(stack!=NULL)
    {
      a = stack->data;
      list_head_remove(stack);
    }
}
void find_suboptimal(QAP *q,Node *&stack,QAPSolution *sln)
{
  int *p, m, n = sln->n;
  int kids, loc_in_row, loc_kept;
  QAPAssignment node;//, kid[n];
  double obj, bound;
  QPBParameters *qpbp;
  MAT *X, *U;
  QAP *qr;

  //  QAPAssignment **kid = new (QAPAssignment *)[n];
  QAPAssignment kid[n];

  qpbp = new QPBParameters;
  qpbp->compute_lapU = false;  qpbp->lapU = 0.0; qpbp->feasible = true;
  BranchParameters *bp = new BranchParameters();
  
  X = m_get(n,n);
  U = m_get(n,n);
  p = new int[n];
  qr = new QAP(n);

  qpbp->inc = sln->getObjective();
  
  qpbp->maxier = 200;
  qpbp->maxier_nofathom = 200;
  qpbp->step = 30;

  sln->obj = 1e10;

  while(1) {

  pop_alt(stack,node);

  if(node.n - node.nfix == 3) {
    solveSize3QAP(q,node,obj,p); 
    sln->update(obj,p);
    break;   // (was return, which skipped the clean-up below)
  }

  reduceQAP(q,qr,node);
  m  = qr->A->m;
  
  X = m_resize(X,m,m);
  U = m_resize(U,m,m);
  init_X(X,m);

  qpbp->shift = qr->shift;

  bound = qpbnd_parametric(qr->A,qr->B,qr->C,X,U,qpbp);
  std::cout << "Level " << node.depth << ", bound = " << bound << std::endl;
  
  qpbp->inc = sln->getObjective();
  qpbp->maxier = 200;
  qpbp->maxier_nofathom = 100;
  qpbp->step = 30;
  bp->use_best = n;
  bp->how = 3;
  bp->use_J = false;
  bp->reorder_children = true;
  bp->which_bound = QP_BND_PARAM;
  branch(q,node,X,U,kid,kids,bound,qpbp,bp,loc_in_row,loc_kept);

  list_head_insert(stack,kid[0]); 

  }

  M_FREE(X);
  M_FREE(U);
  delete [] p;
  delete qr;
  delete qpbp;
  delete bp;
}

/*
============================================================
  pop(stack,node,stats);

  Remove the first item in the stack and store it in the
  given QAPAssignment.
============================================================
 */
void pop(Node *&stack,QAPAssignment &a, QAPStatistics *stats)
{
  if(stack!=NULL)
    {
      a = stack->data;
      list_head_remove(stack);
      stats->stack_size--;
    }
  else
    stats->stack_size=0;
}

/*
============================================================
  push(stack,node,stats);

============================================================
 */
void push(Node *&stack,const QAPAssignment &node, QAPStatistics *stats)
{
  list_head_insert(stack,node);
  stats->stack_size++;
}
