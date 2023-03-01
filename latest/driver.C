/*A
  usage:

  -b      benchmark
  -c      enter console (do not run branch and bound)
  -d X    run to depth X (overrides parameter file)
  -e      estimate (do not run branch and bound)
  -f X    scale all times by X
  -h      use heuristic (do not run branch and bound)
[  -o      change output mode]
  -p F    use parameters file P
  -s F    use symmetry file F
  -u X    set initial upper bound X
  (The path to a QAPLIB input file must also be provided.)

 */


#include "bnb.h"
#include "estimate.h"
#ifdef HAVE_EVB3
#include "evb3.h"
#endif

#include <fstream>
#include <math.h>
#include <time.h>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

//char *sym_file = NULL, *input_file = NULL, *param_file = NULL,
//  *log_file=NULL, *command_file=NULL;
bool run_benchmark, scale, heuristic, console, run_estimate;
bool set_max_depth = false;
bool print_deep_nodes = false;
int max_depth = -1;
double scale_factor = 1.0;
double qap_upper_bound = 100000;
EstimateParameters *ep = NULL;

/*******************************************************************/

void driver_init()
{
  run_benchmark = false;
  run_estimate = false;
  scale = false;
  heuristic = false;
  console = false;

  scale_factor = 1.0;
  qap_upper_bound = 100000;
}

void print_help()
{
  cerr << std::endl;
  cerr << "QPB Branch and Bound Code" << std::endl;
  cerr << "Author: nathan.brixius@gmail.com / @natebrix" << std::endl;
  cerr << "--------------------------------------------" << std::endl;
  cerr << std::endl;
  cerr << "Options:" << std::endl;
  cerr << "--------" << std::endl;
  cerr << "-a F\trun console commands from script file F" << std::endl
       << "-b\tbenchmark" << std::endl
       << "-c\tenter console (do not run branch and bound)"  << std::endl
       << "-d X\trun to depth X (overrides parameter file)" << std::endl
       << "-e\testimate (do not run branch and bound)"  << std::endl
       << "-f X\tscale all times by X"  << std::endl
       << "-h\tuse heuristic (do not run branch and bound)"  << std::endl
       << "-l F\tuse log file F" << std::endl
       << "-p F\tuse parameters file P"  << std::endl
       << "-s F\tuse symmetry file F"  << std::endl
       << "-u X\tset initial upper bound X"  << std::endl;
  cerr << "(The path to a QAPLIB input file must also be provided.)" << std::endl;
}

void print_results(std::ostream &outs,QAPStatistics *stats,
		   QAPSolution *solution,double scale_factor)
{
    
    outs << "  ---- RESULTS ----" << std::endl;
    if (scale_factor != 1.0) {
      stats->scaleTimes(scale_factor);
      outs << "[times have been scaled by " << scale_factor << "]" << std::endl;
    }
    stats->print(outs, 1);
    //solution->print(outs);
}

void process_command_line(int argc, char *argv[], ConsoleParameters *c)
{
  int i = 1;

  while(i < argc) {

    if((argv[i][0] == '-')&&(strlen(argv[1])>1)) {
      
      switch(argv[i][1]) {
      case 'a': // auto console
	i++;
	if (!c->command_file) {
	  c->command_file = strdup(argv[i]);
	}
	cout << "command file: " << c->command_file << std::endl;
	break;
	
      case 'b': // benchmark
	run_benchmark = true;
	cout << "benchmark ON, branch and bound OFF" << std::endl;
	break;

      case 'c': // console
	console = true;
	cout << "console ON, branch and bound OFF" << std::endl;
	break;

      case 'd': // depth
	i++;
	max_depth = atoi(argv[i]);
	set_max_depth = true;
	cout << "max_depth has been set to " << max_depth << std::endl;
	break;

      case 'e': // estimator
	run_estimate = true;
	cout << "estimator ON, branch and bound OFF" << std::endl;
	break;

      case 'f': // statistics scaling factor
	i++;
	scale_factor = atof(argv[i]);
	cout << "scale times by: " << scale_factor << std::endl;
	break;

      case 'h': // heuristic
	cout << "upper bound ON, branch and bound OFF" << std::endl;
	heuristic = true;
	break;

      case 'l': // log file
	i++;
	if (!c->log_file) {
	  c->log_file = strdup(argv[i]);
	}
	cout << "log file: " << c->log_file << std::endl;
	break;
	
      case 'p': // parameters file
	i++;
	if(!c->param_file)
	  c->param_file = strdup(argv[i]);
	cout << "parameters file: " << c->param_file << std::endl;
	break;

      case 's': // symmetry file
	i++;
	if(!c->sym_file)
	  c->sym_file = strdup(argv[i]);
	cout << "symmetry file: " << c->sym_file << std::endl;
	break;

      case 'u': // upper bound
	i++;
	qap_upper_bound = atof(argv[i]);
	cout << "incumbent: " << qap_upper_bound << std::endl;
	break;

      default:
	// unrecognized option
	cout << "ignoring unrecognized option " << argv[i] << std::endl;
	break;
      }

    }
    else   // must be the input file
      {
	if (!c->input_file) {
	  c->input_file = strdup(argv[i]);
	  cout << "input file: " << c->input_file << std::endl;
	}
      }

    i++;
  }

}

// try to get name of file, stripping away path to the file, e.g.
// "../QAPLIB/nug12.dat" --> "nug12.dat"
char *get_suffix(char *input_file)
{
  char *ptr;
  ptr = input_file+strlen(input_file)-1;
  while((ptr!=input_file)&&((*ptr != '/')&&(*ptr != '\\')))
    ptr--;
  if((*ptr == '/')||(*ptr == '\\'))
    ptr++;
  return ptr;
}

// assumes input_file ends in ".dat"
void set_sym_file(char *input_file,char *&sym_file)
{
  char *suffix;

  sym_file = new char[80];
  strcpy(sym_file,"./SYM/");

  suffix = new char[strlen(input_file)];
  strcpy(suffix, get_suffix(input_file));

  int len = strlen(suffix);
  suffix[len-3] = 's';
  suffix[len-2] = 'y';
  suffix[len-1] = 'm';

  sym_file = strcat(sym_file,suffix);
  delete [] suffix;
}

// assumes input_file ends in ".dat"
void set_param_file(char *input_file,char *&param_file)
{
  char *suffix;

  param_file = new char[80];
  strcpy(param_file,"./PARAM/");

  suffix = new char[strlen(input_file)];
  strcpy(suffix, get_suffix(input_file));

  int len = strlen(suffix);
  suffix[len-3] = 'p';
  suffix[len-2] = 'a';
  suffix[len-1] = 'r';

  param_file = strcat(param_file,suffix);
  delete [] suffix;
}

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

int main(int argc, char *argv[])
{
  signal(SIGSEGV, handler);   // install our handler
    
  int n;
  int param_max_depth;  // max depth param read in from param file
  QAP *q;
  QAPSolution *solution;
  QAPStatistics *stats;
  QAPStatistics *estimated_stats; // natbr todo; 6/9/2017 this should be double
  BNBParameters *params;
  ConsoleParameters *c = new ConsoleParameters();
  Node *stack = NULL;

  cout << qap_version_info << std::endl;
  system("hostname");
  driver_init();
  process_command_line(argc, argv, c);

  if (!c->input_file) {
    print_help();
    return -1;
  }
  else {
    q = readQAPLIB(c->input_file);
  }
  
  if (q==NULL) {
    printf("qpb: couldn't read file, bailing out.\n");
    return -1;
  } 

  if (!c->sym_file) {
    set_sym_file(c->input_file, c->sym_file);
  }
  
  if (!c->param_file) {
    set_param_file(c->input_file, c->param_file);
  }
  
  n = q->A->m;
  solution = new QAPSolution(n);

  solution->setObjective(qap_upper_bound);
  
  qap_read_params(c->param_file,params, c->nlevel, c->sym_file,param_max_depth);

  if (!set_max_depth) {
    max_depth = param_max_depth;
  }
  
  stats = new QAPStatistics(n, c->nlevel);
  //  stats->setPrintMode(QAPStatistics<qs_int>::LATEX);

  // Setup the initial state of the stack. 

  int id;
  stack_setup(stack,stats,id,n);
  branch_init(n);

  if (read_symmetry(c->sym_file)) {
    cout << "read symmetry info from " << c->sym_file  << std::endl;
  }
  else {
    cout << "couldn't find [" << c->sym_file << "], assuming no symmetry"
	 << std::endl;
  }
  
  if(run_benchmark) {
    benchmark(q,stack,solution,stats,0,0,params, c->nlevel,
	      max_depth,false,-1);
  }
  else if(console) {
    bnb_console(q, stack, solution, stats, params, c);
  }
  else if(heuristic) {
      find_suboptimal(q,stack,solution);
      solution->print();
  }
  else if(run_estimate) {

    ep =  new EstimateParameters();
    ep->time_scale = scale_factor;
    readEstimateParameters("estimate.par",ep);

    cout << "!! ESTIMATE !! " << std::endl;
    //    estimated_stats = 
    estimate(q,stack,solution,params, c->nlevel,ep);
  }
  else {
    //    use_EVB3(q);
    bnb(q,stack,solution,stats,0,0,params, c->nlevel,
	max_depth,false,-1);

    print_results(cout,stats,solution,scale_factor);

    // nb 1/11/2016: this no longer works because output goes to stdout only.
    //ofstream outf;
    //outf.open("results.out");
    //print_results(outf,stats,solution,scale_factor);
    //outf.close();


    if((max_depth >= 0) && print_deep_nodes) {
      print_stack(q,stack,params, c->nlevel,true,solution->getObjective());
    }
  }
  
  branch_shutdown();
  list_clear(stack);
  delete q;
  delete stats;
  delete solution;
  delete [] params;

  return 0;
}

