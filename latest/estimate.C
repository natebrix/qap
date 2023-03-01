#include "estimate.h"
#include "util.h"

#include <fstream>

const char *dive_hist_file = "dive_history_nug25_2.txt";

void standard_errs(double n_mean,double n_sumsq,
		   double t_mean,double t_sumsq,
		   double trials);


void print_dive(std::ostream &outs,
		QAPStatistics *top_stats,
		QAPStatistics *path_stats,
		double time_scale)
{
  // terminal depth, node estimate, and time estimate for each dive

  int depth=0;
  qs_int node;
  double time;

  node = top_stats->node + path_stats->node;
  time = top_stats->total_time + path_stats->total_time;

  // calculate depth of dive
  while(path_stats->n_level[depth++]==0); 
  while(path_stats->n_level[depth++]>0);
  depth -= 2;

  //  path_stats->print();

  outs << depth << " ";
  outs << node << " ";
  outs << time_scale*time;

  outs << std::endl;

}

QAPStatistics *estimate(QAP *q,Node *&stack,QAPSolution *s,
			BNBParameters *params,int nlevel,
			EstimateParameters *ep)
{
  std::ofstream dive_hist;

  Node *p, *path;
  QAPStatistics *path_stats, *top_stats, *total_stats;
  QAPStatistics *estimated_stats;
  int too_deep = 0;
  double sum_of_bounds, sum_so_far, importance, p_score;
  double val, root_bound, opt, size;
  double n_predicted=0.0, t_predicted=0.0; // node, time sum squared
  double n_sumsq=0.0, t_sumsq=0.0; // node, time sum squared

  top_stats  = new QAPStatistics(s->n,nlevel);
  path_stats = new QAPStatistics(s->n,nlevel);
  total_stats = new QAPStatistics(s->n,nlevel);
  estimated_stats = NULL;

  opt = s->getObjective();

  params->score_exponent = ep->score_exponent;
  std::cout << "Scoring exponent = " << params->score_exponent << std::endl;

  // step 1: run branch and bound code on the given problem,
  // returning all nodes of level start_dive_level.

  bnb(q,stack,s,top_stats,0,0,params,nlevel,ep->start_dive_level-1,false,-1);

  root_bound = top_stats->root_bound;

  top_stats->print();

  // yuck.
  int sdl = ep->start_dive_level;
  top_stats->node = top_stats->node - top_stats->n_level[sdl+1];
  top_stats->n_level[sdl+1] = 0;
  top_stats->f_level[sdl+1] = 0;
  top_stats->t_level[sdl+1] = 0;
  top_stats->l_level[sdl+1] = 0;
  top_stats->e_level[sdl+1] = 0;
  for(int k=0;k<top_stats->m;k++)
    top_stats->n_strat[sdl+1][k] = 0;

  // step 2: compute sum of bounds
  sum_of_bounds = 0.0;
  for(p=stack;p!=NULL;p=p->link)
    sum_of_bounds += score(opt,root_bound,p->data.predicted_bound,
			   params->score_exponent);
  //  sum_of_bounds += 1.0;

  std::cout << "sum = " << sum_of_bounds << std::endl;
  std::cout << "stack = " << top_stats->stack_size << std::endl;

  size = double(top_stats->stack_size);
  if(top_stats->stack_size == 0) {
    total_stats->merge_stats(top_stats);
    return total_stats;
  } 

  for(int i=1;i<=ep->trials;i++) {
    path_stats->clear();

    // step 3: get a random number between 0 and sum_of_bounds.
    //         select a node for estimation based on this value.

    val = sum_of_bounds*drand();
    
    sum_so_far = 0.0; p = stack;
    do {
      sum_so_far += score(opt,root_bound,p->data.predicted_bound,
			  params->score_exponent);
      if(sum_so_far >= val) break;
      p = p->link;
    } while(p->link!=NULL);

    // p now points to a node we want to estimate.

    path = new Node;
    path->data = p->data;
    path->link = NULL;

    p_score = score(opt,root_bound,path->data.predicted_bound,
		    params->score_exponent);
    importance  = p_score / sum_of_bounds;
    
    // step 4: random walk down path of tree.
    path_stats->root_bound = root_bound;
    params->est_factor = 1.0/importance;
    bnb(q,path,s,path_stats,0,0,params,nlevel,-1,true,-1);

    if(dive_depth(path_stats) >= ep->ignore_depth) {
      //      cout << "too deep: " << i << std::endl;
      //      path_stats->print();
      too_deep++;
      i--; continue;
    }
    
    n_predicted = top_stats->node + path_stats->node;
    t_predicted = ep->time_scale*(top_stats->total_time + 
				  path_stats->total_time);
    n_sumsq += n_predicted*n_predicted;
    t_sumsq += t_predicted*t_predicted;

    total_stats->merge_stats(path_stats);
    total_stats->merge_stats(top_stats);

    if(i % ep->print_every == 0) {

      dive_hist.open(dive_hist_file,std::ofstream::app);
      print_dive(dive_hist,top_stats,path_stats,ep->time_scale);
      dive_hist.close();
    }
    if((i % ep->print_everything_every == 0)&&(i<ep->trials)) {
      print_time();
      printf("[ %5d dives  %12d nodes  %12.2f seconds %d ignored]\n",
		i,int(double(total_stats->node)/i),
		ep->time_scale*total_stats->total_time/i,too_deep);

      estimated_stats = new QAPStatistics(*total_stats);
      estimated_stats->scale(1.0/i);
      estimated_stats->scaleTimes(ep->time_scale);
      estimated_stats->print();
      delete estimated_stats;
    }
  }

  total_stats->scale(1.0/ep->trials);
  total_stats->scaleTimes(ep->time_scale);
  total_stats->print();

  double t_mean,n_mean;

  t_mean = total_stats->total_time;
  n_mean = total_stats->node;
  standard_errs(n_mean,n_sumsq,t_mean,t_sumsq,ep->trials);

  std::cout << too_deep << " dives with depth >= " << ep->ignore_depth 
       << " were thrown out." << std::endl;

  // free up memory...
  delete path_stats;
  delete top_stats;
  // estimated_stats already freed.

  return total_stats;
}

void standard_errs(double n_mean,double n_sumsq,
		   double t_mean,double t_sumsq,
		   double trials)
{
  double t_var, n_var;
  double t_stddev, n_stddev;

  std::cout << "tmean = " << t_mean <<  "\t nmean = " << n_mean << std::endl;
  std::cout << "tssq  = " << t_sumsq << "\t nssq  = " << n_sumsq << std::endl;

  t_var  = t_sumsq/(trials) - t_mean*t_mean;
  n_var  = n_sumsq/(trials) - n_mean*n_mean;

  t_stddev = sqrt(t_var);
  n_stddev = sqrt(n_var);

  std::cout << "tvar  = " << t_var << "\t nvar  = " << n_var << std::endl;
  std::cout << "time stddev = " << t_stddev << ", err = " << t_stddev/sqrt(trials) << std::endl;
  std::cout << "node stddev = " << n_stddev << ", err = " << n_stddev/sqrt(trials) << std::endl;
}

void readEstimateParameters(const char *filename,EstimateParameters *&ep)
{
  if(!ep) ep = new EstimateParameters();

  std::ifstream f;

  f.open(filename);

  if(!f)
    {
      std::cout << "sorry chief, couldn't find [" << filename << "].\n";
      f.close();
      return;
    }

  if(f)
    f >> ep->start_dive_level;
  if(f)
    f >> ep->trials;
  if(f)
    f >> ep->print_every;
  if(f)
    f >> ep->print_everything_every;
  if(f)
    f >> ep->score_exponent;
  if(f)
    f >> ep->ignore_depth;

  if(!f)
    {
      std::cout << "warning: couldn't read estimator parameters, using defaults" 
	   << std::endl;
    }
  f.close();

  printf("start dives at %d, %d trials, print every %d/%d, score exp = %f\n",
	    ep->start_dive_level,ep->trials,ep->print_every,
	    ep->print_everything_every,ep->score_exponent);
  printf("throw out dives with depth >= %d\n",ep->ignore_depth);
}

/*
  return depth of deepest node in s.
 */
int dive_depth(QAPStatistics *s)
{
  int depth = 0;

  // find first level with positive number of nodes.
  while((depth <= s->n)&&(s->n_level[depth] == 0)) depth++;

  if(depth > s->n) return depth;

  // loop until we find a level with no nodes.
  while((depth <= s->n)&&(s->n_level[depth] != 0)) depth++;
  
  if(depth > s->n) return depth;

  // depth-1 is the last level reached.
  return depth - 1;
}



