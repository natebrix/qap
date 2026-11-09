#include "stats.h"
#include "util.h"
#include <iomanip>
#include <math.h>
#include <stdio.h>

//template <class ctr_type>
QAPStatistics::QAPStatistics()
{
  print_mode = DEFAULT;

  total_time = 0;
  bound_time = 0;
  stack_size = 0;
  update_number = 0;

  root_bound = 0.0;

  n = 0;
  node = 0;
  fathom = 0;
  enumerate = 0;
  fw_iter = 0;
  dnn_nodes = 0;  dnn_fathoms = 0;  dnn_time = 0.0;
  dnn_cut_nodes = 0;  dnn_cut_fathoms = 0;  dnn_cut_time = 0.0;

  n_level = NULL;
  f_level = NULL;
  t_level = NULL;
  l_level = NULL;
  e_level = NULL;

  n_strat = NULL;
  t_strat = NULL;
  pop_f_level = NULL;
  dnn_n_level = dnn_f_level = dnn_g_level = dnn_i_level = NULL;
  dnn_t_level = NULL;
  cut_n_level = cut_f_level = NULL;
  cut_t_level = NULL;

#ifdef RECORD_STDDEV
  rg_sum  = NULL;
  rg_sumsq= NULL;
#endif
}

void QAPStatistics::alloc_dnn_levels()
{
  pop_f_level = new_zero<ctr_type>(n+1);
  dnn_n_level = new_zero<ctr_type>(n+1);
  dnn_f_level = new_zero<ctr_type>(n+1);
  dnn_g_level = new_zero<ctr_type>(n+1);
  dnn_i_level = new_zero<ctr_type>(n+1);
  dnn_t_level = new_zero<double>(n+1);
  cut_n_level = new_zero<ctr_type>(n+1);
  cut_f_level = new_zero<ctr_type>(n+1);
  cut_t_level = new_zero<double>(n+1);
}

void QAPStatistics::free_dnn_levels()
{
  delete [] pop_f_level;
  delete [] dnn_n_level;  delete [] dnn_f_level;  delete [] dnn_g_level;
  delete [] dnn_i_level;  delete [] dnn_t_level;
  delete [] cut_n_level;  delete [] cut_f_level;  delete [] cut_t_level;
}

QAPStatistics::QAPStatistics(int n_in,int m_in)
{
  print_mode = DEFAULT;

  total_time = 0;
  bound_time = 0;
  stack_size = 0;
  update_number = 0;

  root_bound = 0.0;

  n = n_in;
  m = m_in;
  node = 0;
  fathom = 0;
  enumerate = 0;
  fw_iter = 0;
  dnn_nodes = 0;  dnn_fathoms = 0;  dnn_time = 0.0;
  dnn_cut_nodes = 0;  dnn_cut_fathoms = 0;  dnn_cut_time = 0.0;

  n_level = new_zero<ctr_type>(n+1);
  f_level = new_zero<ctr_type>(n+1);
  t_level = new_zero<double>(n+1);
  l_level = new_zero<ctr_type>(n+1);
  e_level = new_zero<ctr_type>(n+1);

  n_strat = new ctr_type *[n+1];
  for(int i=0;i<=n;i++)
    n_strat[i] = new_zero<ctr_type>(m);
  t_strat = new_zero<double>(m);
  alloc_dnn_levels();

#ifdef RECORD_STDDEV
  rg_sum = new_zero<double>(n+1);
  rg_sumsq = new_zero<double>(n+1);
#endif

}

QAPStatistics::QAPStatistics(const QAPStatistics &s)
{
  this->print_mode = s.print_mode;

  this->total_time = s.total_time;
  this->bound_time = s.bound_time;
  this->stack_size = s.stack_size;
  this->update_number = s.update_number;
  this->root_bound = s.root_bound;

  this->n = s.n;
  this->m = s.m;
  this->node = s.node;
  this->fathom = s.fathom;
  this->enumerate = s.enumerate;
  this->dnn_nodes = s.dnn_nodes;
  this->dnn_fathoms = s.dnn_fathoms;
  this->dnn_time = s.dnn_time;
  this->dnn_cut_nodes = s.dnn_cut_nodes;
  this->dnn_cut_fathoms = s.dnn_cut_fathoms;
  this->dnn_cut_time = s.dnn_cut_time;
  this->fw_iter = s.fw_iter;

  this->n_level = copy_array<ctr_type>(s.n_level,s.n+1);
  this->f_level = copy_array<ctr_type>(s.f_level,s.n+1);
  this->t_level = copy_array<double>(s.t_level,s.n+1);
  this->l_level = copy_array<ctr_type>(s.l_level,s.n+1);
  this->e_level = copy_array<ctr_type>(s.e_level,s.n+1);

  this->n_strat = new ctr_type *[s.n + 1];
  for(int i=0;i<=n;i++)
    this->n_strat[i] = copy_array<ctr_type>(s.n_strat[i],s.m);
  this->t_strat = copy_array<double>(s.t_strat,s.m);
  this->pop_f_level = copy_array<ctr_type>(s.pop_f_level,s.n+1);
  this->dnn_n_level = copy_array<ctr_type>(s.dnn_n_level,s.n+1);
  this->dnn_f_level = copy_array<ctr_type>(s.dnn_f_level,s.n+1);
  this->dnn_g_level = copy_array<ctr_type>(s.dnn_g_level,s.n+1);
  this->dnn_i_level = copy_array<ctr_type>(s.dnn_i_level,s.n+1);
  this->dnn_t_level = copy_array<double>(s.dnn_t_level,s.n+1);
  this->cut_n_level = copy_array<ctr_type>(s.cut_n_level,s.n+1);
  this->cut_f_level = copy_array<ctr_type>(s.cut_f_level,s.n+1);
  this->cut_t_level = copy_array<double>(s.cut_t_level,s.n+1);

#ifdef RECORD_STDDEV
  rg_sum = copy_array<double>(s.rg_sum,s.n+1);
  rg_sumsq = copy_array<double>(s.rg_sumsq,s.n+1);
#endif
}

QAPStatistics::~QAPStatistics()
{
  delete [] n_level;
  delete [] f_level;
  delete [] t_level;
  delete [] e_level;
  delete [] l_level;

  for(int i=0;i<=n;i++)
    delete [] n_strat[i];
  delete [] n_strat;
  delete [] t_strat;
  free_dnn_levels();
#ifdef RECORD_STDDEV
  delete [] rg_sum;
  delete [] rg_sumsq;
#endif
}

void QAPStatistics::clear()
{
  total_time = 0;
  bound_time = 0;
  stack_size = 0;
  update_number = 0;

  node = 0;
  fathom = 0;
  enumerate = 0;
  fw_iter = 0;
  dnn_nodes = 0;  dnn_fathoms = 0;  dnn_time = 0.0;
  dnn_cut_nodes = 0;  dnn_cut_fathoms = 0;  dnn_cut_time = 0.0;

  for(int i=0; i<=n; i++)
    {
      n_level[i] = 0;
      f_level[i] = 0;
      t_level[i] = 0.0;
      e_level[i] = 0;
      l_level[i] = 0;
      for(int j=0;j<m;j++)
	n_strat[i][j] = 0;
      pop_f_level[i] = 0;
      dnn_n_level[i] = dnn_f_level[i] = dnn_g_level[i] = dnn_i_level[i] = 0;
      dnn_t_level[i] = 0.0;
      cut_n_level[i] = cut_f_level[i] = 0;
      cut_t_level[i] = 0.0;
    }
  for(int i=0; i<m; i++)
    t_strat[i] = 0.0;
#ifdef RECORD_STDDEV
  for(int i=0; i<=n; i++)
    {
      rg_sum[i] = 0.0;
      rg_sumsq[i] = 0.0;
    }
#endif
}

//template <class ctr_type>
std::ostream& operator <<(std::ostream& outs, const QAPStatistics &s)
{
  s.print(outs,1);
  return outs;
}

//template <class ctr_type>
void QAPStatistics::print(std::ostream& outs, int print_level) const
{
  switch(print_mode) {
  case DEFAULT:
    print_default(outs,print_level);
    break;
  case LATEX:
    print_latex(outs,print_level);
    break;
  }
}

//template <class ctr_type>
void QAPStatistics::print_default(std::ostream& outs,int print_level) const
{
  int i;
  ctr_type sum;

  outs << std::endl;
  outs << "    Total Time: " << total_time << " seconds " << std::endl;
  outs << "         Nodes: " << node << " / " 
       << fathom << " / " << enumerate << std::endl;
  outs << "    root bound: " << root_bound << std::endl;
  outs << " fw iterations: " << fw_iter << std::endl;
  if (dnn_nodes > 0)
    outs << "     DNN bound: " << dnn_nodes << " nodes, " << dnn_fathoms
	 << " fathomed, " << dnn_time << " s" << std::endl;
  if (dnn_cut_nodes > 0)
    outs << "     cut phase: " << dnn_cut_nodes << " nodes, " << dnn_cut_fathoms
	 << " fathomed, " << dnn_cut_time << " s" << std::endl;

  if (print_level>0) {
    printf("\n");
    printf("level    nodes");
    printf("\t%%fathom");
    printf("\t%%c.elim");
    printf("\t    time\n");
      
    for(i=0;(i<=n) && (n_level[i]==0);i++);
    for(;(i<=n) && (n_level[i]>0);i++)
      {
	printf("%2d",i);
	printf("  %10llu",n_level[i]);
	//	outs.width(12); outs << n_level[i];
	
	if(n_level[i] > 0) {
	  //outs << "\t";
	  //outs.setf(ios_base::fixed,ios_base::floatfield);
	  //outs.precision(4);
	  //outs << 1.0*f_level[i] / n_level[i];
	  //printf("\t%.4f",1.0*f_level[i] / n_level[i]);
	  printf("\t%.4f",1.0*f_level[i] / n_level[i]);
	}
	else
	  printf("\tN/A");
	
	if(l_level[i] > 0) {
	  //outs << "\t";
	  //outs.setf(ios_base::fixed,ios_base::floatfield);
	  //outs.precision(4);
	  //outs << 1.0*e_level[i]/l_level[i];
	  printf("\t%.4f",1.0*e_level[i]/l_level[i]);
	}
	else
	  printf("\tN/A");
	
	printf("\t%8.2f\n",t_level[i]);
      }
  }
  
  if(1)
    {
      printf("\nStrategy breakdown:\n");
      for(int i=0;(i<=n) && (n_level[i]>0);i++) {
	printf("%2d | ",i);
	for(int j=0;j<m;j++)
	  printf("%8.4f ",1.0*n_strat[i][j]/n_level[i]);
	//printf(" | %d ",n_level[i]);
	printf(" | %llu \n", n_level[i]);
      }

      
      printf("nodes");
      for(int j=0;j<m;j++) {
	sum = 0;
	for(int i=0; i<=n; i++) {
	  sum += n_strat[i][j];
	}
	printf("%8llu ",sum);
      }
      printf("\n");
      
      printf("time ");
      for(int j=0;j<m;j++) {
	printf("%8.2f ",t_strat[j]);
      }
      printf("\n");
      printf("\n");
  }

#ifdef RECORD_STDDEV
  double mean,var;
  printf("Relative gap information:\n");
  printf("      mean  \tstd. dev\n");
  for(int i=0;(i<n) && (n_level[i]>0); i++) {
    mean = rg_sum[i] / n_level[i];
    var  = rg_sumsq[i] / n_level[i] - mean*mean;
    printf("%2d| %f\t%f\r\n",i,mean,sqrt(var));
  }

  /*
  outs << std::endl;
  for(int i=0;(i<n);i++) {
    printf("%2d| %f\t%f\n",i,rg_sum[i],rg_sumsq[i]);
  }
  */
#endif

  if (dnn_nodes > 0) print_dnn_levels();

  outs << std::endl;
}

// The DNN bound (and its cut phase) by level, and every level's fathoms by
// cause, after the tables of a run without -x, which stay as they were.
void QAPStatistics::print_dnn_levels() const
{
  int i;
  fflush(stdout);
  printf("\nDNN bound by level:\n");
  printf("level    nodes\t%%fathom\t%%gaveup\t      iters\t it/node\t    time\n");
  for (i = 0; i <= n; i++) {
    if (dnn_n_level[i] == 0) continue;
    printf("%2d  %10llu\t%.4f\t%.4f\t%11llu\t%8.1f\t%8.2f\n", i, dnn_n_level[i],
	   1.0*dnn_f_level[i]/dnn_n_level[i], 1.0*dnn_g_level[i]/dnn_n_level[i],
	   dnn_i_level[i], 1.0*dnn_i_level[i]/dnn_n_level[i], dnn_t_level[i]);
  }
  if (dnn_cut_nodes > 0) {
    printf("\nCut phase by level:\n");
    printf("level    nodes\t%%fathom\t    time\n");
    for (i = 0; i <= n; i++) {
      if (cut_n_level[i] == 0) continue;
      printf("%2d  %10llu\t%.4f\t%8.2f\n", i, cut_n_level[i],
	     1.0*cut_f_level[i]/cut_n_level[i], cut_t_level[i]);
    }
  }
  // f_level less the counted causes is QPB's share (a run resumed from a
  // state file before version 7 has no causes for its first part)
  printf("\nFathoms by cause:\n");
  printf("level  fathoms\t     qpb\t     dnn\t    cuts\t  popped\n");   // as wide as the level table
  for (i = 0; (i <= n) && (n_level[i] == 0); i++);
  for (; (i <= n) && (n_level[i] > 0); i++) {
    long long qpb = (long long)f_level[i] - (long long)pop_f_level[i]
      - (long long)dnn_f_level[i] - (long long)cut_f_level[i];
    printf("%2d  %10llu\t%8lld\t%8llu\t%8llu\t%8llu\n", i, f_level[i], qpb,
	   dnn_f_level[i], cut_f_level[i], pop_f_level[i]);
  }
  fflush(stdout);
}

//template <class ctr_type>
void QAPStatistics::print_latex(std::ostream& outs,int print_level) const
{
  int i;
  ctr_type sum;

  outs << std::endl;
  outs << "\\begin{tabular}{rl}" << std::endl;
  outs << "    Total Time: &" << total_time << " seconds \\\\" << std::endl;
  outs << "         Nodes: &" << node << " / " 
       << fathom << " / " << enumerate << "\\\\" << std::endl;
  outs << "    root bound: &" << root_bound << "\\\\" << std::endl;
  outs << " fw iterations: &" << fw_iter << "\\\\" << std::endl;
  outs << "\\end{tabular}" << std::endl;

  if(print_level>0)
    {
      outs << std::endl;
      outs << "\\begin{tabular}{rrlll}" << std::endl;
      outs << "  & nodes";
      outs << "&\\%fathom";
      outs << "&\\%c.elim";
      outs << "&    time";
      outs << "\\\\" << std::endl;
      
      for(i=0;(i<=n) && (n_level[i]==0);i++);
      for(;(i<=n) && (n_level[i]>0);i++)
	{
	  outs.width(2); outs << i;
	  outs << "&";
	  outs.width(10); outs << n_level[i];

	  if(n_level[i] > 0) {
	    printf("\t& %.4f",1.0*f_level[i] / n_level[i]);
	  }
	  else
	    outs << "\t& N/A";

	  if(l_level[i] > 0) {
	    printf("\t& %.4f",1.0*e_level[i]/l_level[i]);
	  }
	  else
	   outs << "\t& N/A";

	  printf("\t& %8.2f",t_level[i]);
	  outs << "\\\\" << std::endl;
	}
      outs << "\\end{tabular}" << std::endl;
    }

  if(1)
    {
      outs << std::endl << "Breakdown of strategies used:" << std::endl;

      outs << "\\begin{tabular}{r";
      for(int kk=0;kk<m;kk++)
	outs << "r";
      outs << "l}" << std::endl;

      for(int i=0;(i<=n) && (n_level[i]>0);i++) {
	printf("%2d ",i);
	for(int j=0;j<m;j++)
	  printf("& %8.4f ",1.0*n_strat[i][j]/n_level[i]);
	outs << " & " << n_level[i] << " ";
	outs << "\\\\" << std::endl;
      }

      outs << "nodes ";
      for(int j=0;j<m;j++)
	{
	  sum = 0;
	  for(int i=0; i<=n; i++)
	    sum += n_strat[i][j];
	  outs.width(8); outs << "&" << sum << " ";
	}
      outs << "&\\\\" << std::endl;

      outs << "time ";
      for(int j=0;j<m;j++) {
	printf("& %8.2f ",t_strat[j]);
      }
      //      outs << "&\\\\" << std::endl;
      outs << std::endl;
      outs << "\\end{tabular}" << std::endl;
  }

#ifdef RECORD_STDDEV
  double mean,var;

  outs << std::endl << "Relative gap information:" << std::endl;
  outs << "\\begin{tabular}{rrr}";
  outs << "  &    mean  & \tstd. dev";
  outs << "\\\\" << std::endl;

  for(int i=0;(i<n) && (n_level[i]>0);i++) {
    mean = rg_sum[i]/n_level[i];
    var  = rg_sumsq[i]/n_level[i]-mean*mean;
    printf("%2d & %f\t & %f \\\\ \n",i,mean,sqrt(var));
  }
#endif

  outs << "\\end{tabular}" << std::endl;
  outs << std::endl;
}

//template <class ctr_type>
void QAPStatistics::anotherNode(int depth,int strat,double factor) 
{
  node += ctr_type(factor);

  if(depth>n) depth = n;

  n_level[depth] += ctr_type(factor);
  n_strat[depth][strat] += ctr_type(factor);

}

//template <class ctr_type>
void QAPStatistics::anotherFathom(int depth,double factor) 
{
  fathom += ctr_type(factor);  
  if(depth>n) depth = n;

  f_level[depth] += ctr_type(factor);
}

void QAPStatistics::anotherDiscard(int depth)
{
  anotherFathom(depth);
  if (depth > n) depth = n;
  pop_f_level[depth]++;
}

void QAPStatistics::anotherDNN(int depth, bool fathomed, bool gave_up, int iterations, double cpu)
{
  dnn_nodes++;
  dnn_time += cpu;
  if (fathomed) dnn_fathoms++;
  if (depth > n) depth = n;
  dnn_n_level[depth]++;
  if (fathomed) dnn_f_level[depth]++;
  if (gave_up) dnn_g_level[depth]++;
  dnn_i_level[depth] += iterations > 0 ? ctr_type(iterations) : 0;
  dnn_t_level[depth] += cpu;
}

void QAPStatistics::anotherCut(int depth, bool fathomed, double cpu)
{
  dnn_cut_nodes++;
  dnn_cut_time += cpu;
  if (fathomed) dnn_cut_fathoms++;
  if (depth > n) depth = n;
  cut_n_level[depth]++;
  if (fathomed) cut_f_level[depth]++;
  cut_t_level[depth] += cpu;
}

//template <class ctr_type>
void QAPStatistics::start()
{
 start_time = clock();
 fw_iter = 0;
}

//template <class ctr_type>
void QAPStatistics::stop()
{
  double diff;
  stop_time = clock();
  diff = (stop_time - start_time)/CLOCKS_PER_SEC;
  if(diff < 0) diff += CLOCK_ROLLOVER;
  total_time += diff;
  start_time = stop_time;
}

//template <class ctr_type>
void QAPStatistics::startIteration(int depth)
{
  start_time = clock();
  ier_start_time = start_time;
}

//template <class ctr_type>
void QAPStatistics::stopIteration(int depth,int strat,
					    int fw,double factor)
{
  double ier_time,diff;
  stop_time = clock();

  diff = (stop_time - start_time)/CLOCKS_PER_SEC;
  if(diff < 0) diff += CLOCK_ROLLOVER;

  ier_time = (stop_time - ier_start_time)/CLOCKS_PER_SEC;
  if(ier_time < 0) ier_time += CLOCK_ROLLOVER;

  if(depth > n) depth = n;

  total_time += ier_time*factor;
  t_level[depth] += ier_time*factor;
  t_strat[strat] += ier_time*factor;

  fw_iter += ctr_type(factor*fw);

  start_time = clock();
}

//template <class ctr_type>
void QAPStatistics::update()
{
  double diff;
  stop_time = clock();
  diff = (stop_time - start_time)/CLOCKS_PER_SEC;
  if(diff < 0) diff += CLOCK_ROLLOVER;
  total_time += diff;
  start_time = stop_time;
}

void QAPStatistics::merge_stats(QAPStatistics *new_stats,
				double factor)
{
 
  total_time += factor*new_stats->total_time;
  bound_time += factor*new_stats->bound_time;

  stack_size = ( ( ( stack_size * update_number) + new_stats->stack_size ) 
		 / (update_number + 1) );
  update_number += 1;

  node += new_stats->node;
  fathom += new_stats->fathom;
  enumerate += new_stats->enumerate;
  fw_iter += new_stats->fw_iter;
  dnn_nodes += new_stats->dnn_nodes;
  dnn_fathoms += new_stats->dnn_fathoms;
  dnn_time += new_stats->dnn_time;
  dnn_cut_nodes += new_stats->dnn_cut_nodes;
  dnn_cut_fathoms += new_stats->dnn_cut_fathoms;
  dnn_cut_time += new_stats->dnn_cut_time;
  
  root_bound = new_stats->root_bound;

  for(int i=0;i<=n;i++)
    n_level[i] += new_stats->n_level[i];

  for(int i=0;i<=n;i++)
    f_level[i] += new_stats->f_level[i];

  for(int i=0;i<=n;i++)
    t_level[i] += factor*new_stats->t_level[i];

  for(int i=0;i<=n;i++)
    e_level[i] += new_stats->e_level[i];

  for(int i=0;i<=n;i++)
    l_level[i] += new_stats->l_level[i];

  for(int i=0;i<=n;i++)
    for(int j=0;j<m;j++)
      n_strat[i][j] += new_stats->n_strat[i][j];

  for(int i=0;i<m;i++)
    t_strat[i] += factor*new_stats->t_strat[i];

  // like the DNN totals above: summed, times not scaled by factor
  for(int i=0;i<=n;i++) {
    pop_f_level[i] += new_stats->pop_f_level[i];
    dnn_n_level[i] += new_stats->dnn_n_level[i];
    dnn_f_level[i] += new_stats->dnn_f_level[i];
    dnn_g_level[i] += new_stats->dnn_g_level[i];
    dnn_i_level[i] += new_stats->dnn_i_level[i];
    dnn_t_level[i] += new_stats->dnn_t_level[i];
    cut_n_level[i] += new_stats->cut_n_level[i];
    cut_f_level[i] += new_stats->cut_f_level[i];
    cut_t_level[i] += new_stats->cut_t_level[i];
  }

#ifdef RECORD_STDDEV
  for(int i=0;i<=n;i++) {
    rg_sum[i] += new_stats->rg_sum[i];
    rg_sumsq[i] += new_stats->rg_sumsq[i];
  }
#endif

}

//template <class ctr_type>
void QAPStatistics::updateRelativeGap(int depth, double rel_gap,
						double factor)
{
  if(depth > n) depth = n;
  rg_sum[depth] += factor*rel_gap;
  rg_sumsq[depth] += factor*rel_gap*rel_gap;
}

//template <class ctr_type>
void QAPStatistics::scaleTimes(double factor)
{
  total_time = factor*total_time;
  bound_time = factor*bound_time;
  for(int i=0;i<=n;i++)
    t_level[i] = factor*t_level[i];
  for(int i=0;i<m;i++)
    t_strat[i] = factor*t_strat[i];

}


void QAPStatistics::scale(double factor)
{
 
  total_time = factor*total_time;
  bound_time = factor*bound_time;

  node      = ctr_type1(factor*node);
  fathom    = ctr_type1(factor*fathom);
  enumerate = ctr_type1(factor*enumerate);
  fw_iter   = ctr_type1(factor*fw_iter);  
  
  for(int i=0;i<=n;i++)
    n_level[i] = ctr_type1(factor*n_level[i]);

  for(int i=0;i<=n;i++)
    f_level[i] = ctr_type1(factor*f_level[i]);

  for(int i=0;i<=n;i++)     // part of f_level; the DNN arrays, like the DNN totals, are not scaled
    pop_f_level[i] = ctr_type1(factor*pop_f_level[i]);

  for(int i=0;i<=n;i++)
    t_level[i] = factor*t_level[i];

  for(int i=0;i<=n;i++)
    e_level[i] = ctr_type1(factor*e_level[i]);

  for(int i=0;i<=n;i++)
    l_level[i] = ctr_type1(factor*l_level[i]);

  for(int i=0;i<=n;i++)
    for(int j=0;j<m;j++)
      n_strat[i][j] = ctr_type1(factor*n_strat[i][j]);

  for(int i=0;i<m;i++)
    t_strat[i] = factor*t_strat[i];

#ifdef RECORD_STDDEV
  for(int i=0;i<=n;i++) {
    rg_sum[i]   = factor*rg_sum[i];
    rg_sumsq[i] = factor*rg_sumsq[i];
  }
#endif

}
