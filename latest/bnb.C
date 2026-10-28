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

// DNN bound (-x dmin:dmax[:maxgap[:hopeless]]): applied to nodes QPB does not
// fathom at depths dmin..dmax whose relative gap (inc - QPB)/(inc - root QPB)
// is at most maxgap. Off while dnn_max_depth < 0.
int dnn_min_depth = 0, dnn_max_depth = -1;
double dnn_max_gap = 1.0;
int dnn_hopeless_after = 0;          // DNNParams::hopeless_after; 0 = never give up
long dnn_calls = 0, dnn_fathoms = 0;
double dnn_cpu = 0;
#ifdef HAVE_DNN
#include "dnn.h"
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
  const WarmStart *w = node.warm;
  const int n = node.n;
  int extra = -1;

  if (w && w->n == n && (int)w->X->m == (int)X->m + 1) {
    for (int i = 0; i < n; i++) {
      if (w->p[i] == node.p[i]) continue;
      if (w->p[i] == NONE && node.p[i] >= 0 && extra < 0) extra = i;
      else { extra = -1; break; }
    }
  }
  if (extra < 0) {
    init_X(X, X->m);
    return;
  }

  // position of the fixed row and column among the parent's free ones
  const int col = node.p[extra];
  int di = 0, dj = 0;
  for (int i = 0; i < extra; i++)
    if (node.p[i] == NONE) di++;
  for (int j = 0; j < col; j++)
    if (node.pinv[j] == NONE) dj++;
  init_X_from_parent(X, w->X, di, dj);
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
	  write_state(bnb_state_file, stack, s, stats);
	  next_checkpoint = stats->total_time + bnb_checkpoint_every;
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

      node.set_warm(NULL);    // done with the parent's X

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
      bool dnn_fathomed = false;
#ifdef HAVE_DNN
      if (bound <= s->getObjective() - DELTA && m >= 2 &&
	  node.depth >= dnn_min_depth && node.depth <= dnn_max_depth) {
	const double inc = s->getObjective();
	const double gap = (inc - bound) / std::max(inc - stats->root_bound, 1e-9);
	if (gap <= dnn_max_gap) {
	  DNNParams dp;
	  // giving up early on a node that would fathom costs its whole subtree
	  dp.hopeless_after = dnn_hopeless_after;
	  DNNResult r = dnn_bound(qr->A, qr->B, qr->C, qr->shift, inc - DELTA, dp);
	  dnn_calls++;
	  dnn_cpu += r.cpu;
	  node.bound = std::max(node.bound, r.lb);
	  if (r.status == DNN_FATHOMED) {
	    dnn_fathomed = true;
	    dnn_fathoms++;
	  }
	}
      }
#endif
    
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
	branch(q, node, X, U, kid, kids, bound + qpbp->lapU,
	       qpbp, bp, loc_in_row, loc_kept);
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
	    for(sum_of_scores=0.0,i=0;i<kids;i++)
	      sum_of_scores+=score(s->getObjective(),stats->root_bound,
				   kid[i].predicted_bound,
				   params->score_exponent);
	    val = sum_of_scores*drand();
	    
	    sum_so_far = 0.0; i=0;
	    do {
	      sum_so_far += score(s->getObjective(),stats->root_bound,
				  kid[i].predicted_bound,
				  params->score_exponent);
	      if(sum_so_far >= val) break;
	      i++;
	    } while(i<kids-1);
	    
	    importance = score(s->getObjective(),stats->root_bound,
			       kid[i].predicted_bound,
			       params->score_exponent)/sum_of_scores;
	    
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

  if (dnn_calls > 0)
    printf("DNN bound: %ld nodes, %ld fathomed, %.1f s\n", dnn_calls, dnn_fathoms, dnn_cpu);

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
