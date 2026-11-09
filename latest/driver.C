#include <cstdio>
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
  -o F    write state to F (checkpoints, -t, end of run); -k X checkpoint every
          X CPU s; -t X stop after X CPU s (exit 3); -r F resume from F
  -w W    Frank-Wolfe patience at lookahead/strong levels once a node cannot
          be fathomed (default 0 = off)
  (The path to a QAPLIB input file must also be provided.)

 */


#include "bnb.h"
#include "estimate.h"
#include "state.h"
#include "settings.h"
#include <cerrno>
#ifdef HAVE_EVB3
#include "evb3.h"
#endif

#include <fstream>
#include <math.h>
#include <cmath>
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
int max_depth = -1;
double scale_factor = 1.0;
// No incumbent by default. This used to be 100000, which is below every
// objective value of instances such as tai12b (~39M), so without -u the
// search pruned everything and found nothing.
double qap_upper_bound = 1e100;
EstimateParameters *ep = NULL;

// state files (state.h); bnb() reads the first four
extern const char *bnb_state_file;
extern double bnb_checkpoint_every;
extern double bnb_time_limit;
extern int dnn_min_depth, dnn_max_depth, dnn_hopeless_after, dnn_warm;
extern double dnn_max_gap;
#ifdef HAVE_DNN
#include "dnn.h"
#endif
extern bool bnb_resumed;
extern bool bnb_stopped;
extern volatile sig_atomic_t bnb_stop_requested;
static void request_stop(int) { bnb_stop_requested = 1; }
extern const char *bnb_incumbent_file;
extern volatile sig_atomic_t bnb_incumbent_poke;
static void poke_incumbent(int) { bnb_incumbent_poke = 1; }
char *resume_file = NULL;
char *root_fixes = NULL;   // -F i:j,...: the search starts at this node (its depth = the number of fixes)

/*******************************************************************/

void driver_init()
{
  run_benchmark = false;
  run_estimate = false;
  scale = false;
  heuristic = false;
  console = false;

  scale_factor = 1.0;
  qap_upper_bound = 1e100;
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
       << "-u X\tset initial upper bound X"  << std::endl
       << "-o F\twrite branch-and-bound state to F (at checkpoints, on -t, at the end)" << std::endl
       << "-k X\twith -o, checkpoint every X CPU seconds" << std::endl
       << "-U F\tincumbent file 'obj [perm]', re-read every 50 nodes and on SIGUSR2;\n\ta lower value is adopted" << std::endl
       << "-x D1:D2[:G[:H[:W]]]\tDNN bound on nodes QPB does not fathom at depths D1-D2\n"
       << "\twith relative gap at most G (default 1); H: give up after H\n"
       << "\titerations when hopeless (default 0, never); W=1: start a node's\n"
       << "\tDNN from its parent's when that missed (default 0)" << std::endl
       << "-t X\twith -o, stop after X CPU seconds, saving state (exit code 3)" << std::endl
       << "-r F\tresume from state file F instead of the root" << std::endl
       << "-F i:j,...\tstart at the node with these assignments fixed (depth = number of fixes)" << std::endl
       << "-w W\tat lookahead/strong levels, stop a node's FW W iterations after it\n"
       << "\tlast fixed out more children (default 0 = off)" << std::endl
       << "--print-settings\twith the other options, print the resolved settings and exit" << std::endl
       << "--settings-table\tthe settings reference table (Markdown) and exit" << std::endl
       << "--settings-map\tthe retired environment names and their keys and exit" << std::endl
       << "--allow-settings-change\tresume (-r) although the settings differ from the state file's" << std::endl
       << "(-x and -w are the keys dnn.depths/maxgap/giveup/warm and fw.patience of the\n"
       << "params file's [settings] section; a flag overrides the file. RUNNING.md, Settings.)" << std::endl;
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
    {
      long hits = 0, misses = 0;
      init_ST_cache_stats(hits, misses);
      if (hits + misses > 0)
        outs << "eigendecompositions for W and U: " << misses << " computed, " << hits
             << " reused from the per-size cache" << std::endl;
    }
    //solution->print(outs);
}

// ConsoleParameters frees its file names with delete[], so allocate them
// with new[] rather than strdup() (malloc).
static char *new_string(const char *s)
{
  char *t = new char[strlen(s) + 1];
  strcpy(t, s);
  return t;
}

static bool print_settings_only = false;
extern bool state_allow_settings_change;   // state.C: --allow-settings-change

// Option arguments (review C5): a flag without its argument, a number that
// does not parse, and an unknown option are errors, so that a run never
// silently does something other than what was asked.
static const char *need_arg(int &i, int argc, char *argv[])
{
  if (i + 1 >= argc) {
    cerr << argv[i] << " needs an argument (./bnbdrv with no arguments prints help)" << std::endl;
    exit(1);
  }
  return argv[++i];
}

static long arg_int(const char *flag, const char *s)
{
  errno = 0;
  char *end = NULL;
  const long v = strtol(s, &end, 10);
  if (!*s || *end || errno) {
    cerr << flag << ": " << s << " is not an integer" << std::endl;
    exit(1);
  }
  return v;
}

static double arg_num(const char *flag, const char *s)
{
  errno = 0;
  char *end = NULL;
  const double v = strtod(s, &end);
  if (!*s || *end || errno || !std::isfinite(v)) {
    cerr << flag << ": " << s << " is not a number" << std::endl;
    exit(1);
  }
  return v;
}

void process_command_line(int argc, char *argv[], ConsoleParameters *c)
{
  int i = 1;

  while(i < argc) {
    const char *a = argv[i];

    if (strcmp(a, "--print-settings") == 0) {        // dump the settings and exit (needs no instance)
      print_settings_only = true;
      i++;
      continue;
    }
    if (strcmp(a, "--allow-settings-change") == 0) {  // -r with settings other than the file's
      state_allow_settings_change = true;
      i++;
      continue;
    }

    if (a[0] == '-' && a[1] != '\0' && a[2] == '\0') {

      switch(a[1]) {
      case 'a': { // auto console
	const char *v = need_arg(i, argc, argv);
	if (!c->command_file) {
	  c->command_file = new_string(v);
	}
	cout << "command file: " << c->command_file << std::endl;
	break;
      }

      case 'b': // benchmark
	run_benchmark = true;
	cout << "benchmark ON, branch and bound OFF" << std::endl;
	break;

      case 'c': // console
	console = true;
	cout << "console ON, branch and bound OFF" << std::endl;
	break;

      case 'd': // depth
	max_depth = (int)arg_int("-d", need_arg(i, argc, argv));
	set_max_depth = true;
	cout << "max_depth has been set to " << max_depth << std::endl;
	break;

      case 'e': // estimator
	run_estimate = true;
	cout << "estimator ON, branch and bound OFF" << std::endl;
	break;

      case 'f': // statistics scaling factor
	scale_factor = arg_num("-f", need_arg(i, argc, argv));
	cout << "scale times by: " << scale_factor << std::endl;
	break;

      case 'h': // heuristic
	cout << "upper bound ON, branch and bound OFF" << std::endl;
	heuristic = true;
	break;

      case 'l': { // log file
	const char *v = need_arg(i, argc, argv);
	if (!c->log_file) {
	  c->log_file = new_string(v);
	}
	cout << "log file: " << c->log_file << std::endl;
	break;
      }

      case 'p': { // parameters file
	const char *v = need_arg(i, argc, argv);
	if(!c->param_file)
	  c->param_file = new_string(v);
	cout << "parameters file: " << c->param_file << std::endl;
	break;
      }

      case 's': { // symmetry file
	const char *v = need_arg(i, argc, argv);
	if(!c->sym_file)
	  c->sym_file = new_string(v);
	cout << "symmetry file: " << c->sym_file << std::endl;
	break;
      }

      case 'u': // upper bound
	qap_upper_bound = arg_num("-u", need_arg(i, argc, argv));
	cout << "incumbent: " << qap_upper_bound << std::endl;
	break;

      case 'o': // write state here (checkpoints, time limit, end of run)
	bnb_state_file = new_string(need_arg(i, argc, argv));
	cout << "state file: " << bnb_state_file << std::endl;
	break;

      case 'U': // incumbent file, re-read at node boundaries (SIGUSR2: at once)
	bnb_incumbent_file = new_string(need_arg(i, argc, argv));
	signal(SIGUSR2, poke_incumbent);
	cout << "incumbent file: " << bnb_incumbent_file << std::endl;
	break;

      case 'k': // checkpoint every X CPU seconds
	bnb_checkpoint_every = arg_num("-k", need_arg(i, argc, argv));
	cout << "checkpoint every " << bnb_checkpoint_every << " s" << std::endl;
	break;

      case 't': // stop after X CPU seconds, saving state
	bnb_time_limit = arg_num("-t", need_arg(i, argc, argv));
	cout << "time limit " << bnb_time_limit << " s" << std::endl;
	break;

      case 'r': // resume from a state file
	resume_file = new_string(need_arg(i, argc, argv));
	cout << "resume from: " << resume_file << std::endl;
	break;

      case 'F': // start at the node with these assignments fixed (i:j = first index i to second index j)
	root_fixes = new_string(need_arg(i, argc, argv));
	cout << "root fixes: " << root_fixes << std::endl;
	break;

      case 'x': { // DNN bound at depths dmin..dmax: -x dmin:dmax[:maxgap[:hopeless[:warm]]]
	const char *v = need_arg(i, argc, argv);
#ifdef HAVE_DNN
	{
	  // dmin:dmax[:maxgap[:hopeless[:warm]]]; an empty field leaves that
	  // setting to the params file or its default. The flag is recorded in the
	  // registry (dnn.depths, dnn.maxgap, dnn.giveup, dnn.warm); the globals
	  // are set from the registry once the params file has been read
	  // (apply_policy_settings), so a flag overrides the file.
	  char *f[5] = {NULL, NULL, NULL, NULL, NULL};
	  char *arg = new_string(v);
	  int n = 0;
	  for (char *p = arg; n < 5; n++) {
	    f[n] = p;
	    char *colon = strchr(p, ':');
	    if (!colon) { n++; break; }
	    *colon = '\0';
	    p = colon + 1;
	  }
	  int dmin = -1, dmax = -1;
	  if (n >= 1 && *f[0]) dmin = (int)arg_int("-x", f[0]);
	  if (n >= 2 && *f[1]) dmax = (int)arg_int("-x", f[1]);
	  if (n < 2 || n > 5 || dmin < 0 || dmax < 0 || dmin > dmax ||
	      (n == 5 && strchr(f[4], ':'))) {
	    cerr << "-x needs dmin:dmax[:maxgap[:hopeless[:warm]]] with dmin <= dmax, e.g. -x 2:3" << std::endl;
	    exit(1);
	  }
	  char buf[64];
	  snprintf(buf, sizeof buf, "%d:%d", dmin, dmax);
	  settings_set("dnn.depths", buf, SRC_FLAG);
	  if (n >= 3 && *f[2]) settings_set("dnn.maxgap", f[2], SRC_FLAG);
	  if (n >= 4 && *f[3]) settings_set("dnn.giveup", f[3], SRC_FLAG);
	  if (n >= 5 && *f[4]) settings_set("dnn.warm", f[4], SRC_FLAG);
	  delete [] arg;
	}
#else
	cerr << "-x: this build has no DNN bound (make with HAVE_DNN)" << std::endl;
	exit(1);
#endif
	break;
      }

      case 'w': // Frank-Wolfe patience once a node cannot be fathomed (fw.patience)
	settings_set("fw.patience", need_arg(i, argc, argv), SRC_FLAG);
	break;

      default:
	cerr << "unknown option " << a << " (./bnbdrv with no arguments prints help)" << std::endl;
	exit(1);
      }

    }
    else if (a[0] == '-' && a[1] != '\0') {
      cerr << "unknown option " << a << " (./bnbdrv with no arguments prints help)" << std::endl;
      exit(1);
    }
    else   // must be the input file
      {
	if (!c->input_file) {
	  c->input_file = new_string(a);
	  cout << "input file: " << c->input_file << std::endl;
	}
	else {
	  cerr << "two input files given: " << c->input_file << " and " << a << std::endl;
	  exit(1);
	}
      }

    i++;
  }

  // combinations that cannot mean what was asked
  if (root_fixes && resume_file) {
    cerr << "-F and -r together: a state file already fixes its nodes' assignments" << std::endl;
    exit(1);
  }
  if (bnb_time_limit > 0 && !bnb_state_file) {
    cerr << "-t without -o: the time limit would stop the run with nowhere to save its state" << std::endl;
    exit(1);
  }
  if (bnb_checkpoint_every > 0 && !bnb_state_file) {
    cerr << "-k without -o: checkpoints need a state file" << std::endl;
    exit(1);
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

  // (both buffers used to be one byte short: suffix had no room for the
  // NUL when input_file has no directory, and sym_file was a fixed 80)
  suffix = new char[strlen(input_file) + 1];
  strcpy(suffix, get_suffix(input_file));

  int len = strlen(suffix);
  suffix[len-3] = 's';
  suffix[len-2] = 'y';
  suffix[len-1] = 'm';

  sym_file = new char[strlen("./SYM/") + len + 1];
  strcpy(sym_file,"./SYM/");
  strcat(sym_file,suffix);
  delete [] suffix;
}

// assumes input_file ends in ".dat"
// The policy settings the solver keeps in globals, from the registry once
// the environment, the flags and the params file have all been read
// (settings.h: a flag overrides the file, the file its default).
static void apply_policy_settings()
{
  const char *depths = setting_s("dnn.depths");
  if (*depths) {
#ifdef HAVE_DNN
    int dmin, dmax;
    char tail;
    if (sscanf(depths, "%d:%d%c", &dmin, &dmax, &tail) != 2 || dmin < 0 || dmin > dmax) {
      cerr << "dnn.depths needs dmin:dmax with 0 <= dmin <= dmax, not '" << depths << "'" << std::endl;
      exit(2);
    }
    dnn_min_depth = dmin;
    dnn_max_depth = dmax;
    dnn_max_gap = setting_d("dnn.maxgap");
    dnn_hopeless_after = (int)setting_i("dnn.giveup");
    dnn_warm = (int)setting_i("dnn.warm");
    // one core per solver process. Accelerate (macOS) has no such call;
    // it reads VECLIB_MAXIMUM_THREADS when it is loaded, so that must be
    // set in the environment (bench/parallel.py sets it for its workers).
    if (!dnn_set_blas_threads(1) && !getenv("VECLIB_MAXIMUM_THREADS"))
      cerr << "warning: the BLAS may use several threads per process; "
	   << "export VECLIB_MAXIMUM_THREADS=1 when running several solvers"
	   << std::endl;
    cout << "DNN bound at depths " << dnn_min_depth << "-" << dnn_max_depth
	 << ", relative gap <= " << dnn_max_gap << ", give up after "
	 << dnn_hopeless_after << " iterations (0 = never), "
	 << (dnn_warm ? "warm-started from the parent, " : "")
	 << "extended " << setting_s("dnn.ext") << " levels below a last-level miss, "
	 << "give-ups within " << setting_s("dnn.gu_zone") << " of the cutoff confirmed, "
	 << dnn_lapack_name() << std::endl;
#else
    cerr << "dnn.depths: this build has no DNN bound (make with HAVE_DNN)" << std::endl;
    exit(2);
#endif
  }
  fw_nofathom_patience = (int)setting_i("fw.patience");
  if (fw_nofathom_patience)
    cout << "FW patience after cannot-fathom: " << fw_nofathom_patience << std::endl;
}

void set_param_file(char *input_file,char *&param_file)
{
  char *suffix;

  // (both buffers used to be one byte short: suffix had no room for the
  // NUL when input_file has no directory, and param_file was a fixed 80)
  suffix = new char[strlen(input_file) + 1];
  strcpy(suffix, get_suffix(input_file));

  int len = strlen(suffix);
  suffix[len-3] = 'p';
  suffix[len-2] = 'a';
  suffix[len-1] = 'r';

  param_file = new char[strlen("./PARAM/") + len + 1];
  strcpy(param_file,"./PARAM/");
  strcat(param_file,suffix);
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
  signal(SIGUSR1, request_stop);  // with -o: save state and stop at the next node
    
  int n;
  int param_max_depth;  // max depth param read in from param file
  QAP *q;
  QAPSolution *solution;
  QAPStatistics *stats;
  QAPStatistics *estimated_stats; // natbr todo; 6/9/2017 this should be double
  BNBParameters *params;
  ConsoleParameters *c = new ConsoleParameters();
  Node *stack = NULL;

  // the two listings are documentation (tests/settings_doc.py, bench/envpar.py
  // read them): nothing else on stdout
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--settings-table") == 0) { settings_table(cout); return 0; }
    if (strcmp(argv[i], "--settings-map") == 0) { settings_env_map(cout); return 0; }
  }
  cout << qap_version_info << std::endl;
  system("hostname");
  driver_init();
  process_command_line(argc, argv, c);

  if (print_settings_only) {
    // the resolved settings (environment, flags, the par file if given), then exit
    if (c->param_file) {
      settings_set("par.file", c->param_file, SRC_FLAG);
      settings_set("par.hash", settings_hash_file(c->param_file), SRC_PAR);
      settings_read_par(c->param_file);
    }
    settings_dump(cout);
    return 0;
  }

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

  set_fathom_tolerance(q);   // from the original data, before symmetrizing
  if (!symmetrize_QAP(q)) {
    cout << "WARNING: neither A nor B is symmetric; the QP bounds assume one is,"
	 << " so results for this instance are not reliable." << std::endl;
  }

  if (!c->sym_file) {
    set_sym_file(c->input_file, c->sym_file);
  }
  
  const bool param_file_given = c->param_file != NULL;
  if (!c->param_file) {
    set_param_file(c->input_file, c->param_file);
  }
  settings_set("par.file", c->param_file, param_file_given ? SRC_FLAG : SRC_DEFAULT);
  settings_set("par.hash", settings_hash_file(c->param_file), SRC_PAR);
  
  n = q->A->m;
  solution = new QAPSolution(n);

  // The fathoming tolerance DELTA assumes the incumbent is an attained value:
  // "bound > inc - DELTA" means no attainable value below inc fits. A -u off
  // the objective grid (6125 on nug30, whose values are even) would fathom
  // nodes holding strictly better solutions. Round it up to the grid.
  if (std::isfinite(qap_upper_bound)) {
    double g = DELTA >= 1.9 ? 2 * ceil(qap_upper_bound / 2) :
               DELTA >= 0.9 ? ceil(qap_upper_bound) : qap_upper_bound;
    if (g != qap_upper_bound) {
      cout << "incumbent " << qap_upper_bound << " is not an attainable objective value; using "
	   << g << std::endl;
      qap_upper_bound = g;
    }
  }
  solution->setObjective(qap_upper_bound);
  
  qap_read_params(c->param_file,params, c->nlevel, c->sym_file,param_max_depth);
  settings_read_par(c->param_file);   // its [settings] section (settings.h)
  apply_policy_settings();            // -x/-w or the file into the solver's globals
  settings_dump(cout);   // every setting with its value and source (settings.h)

  if (!set_max_depth) {
    max_depth = param_max_depth;
  }
  
  stats = new QAPStatistics(n, c->nlevel);
  //  stats->setPrintMode(QAPStatistics<qs_int>::LATEX);

  // Setup the initial state of the stack. 

  int id;
  if (resume_file) {
    if (!read_state(resume_file, n, stack, solution, stats)) return 1;
    bnb_resumed = true;
    // an explicit -u below the saved incumbent wins (the default is +inf)
    if (qap_upper_bound < solution->getObjective())
      solution->setObjective(qap_upper_bound);
    cout << "resumed: " << stats->stack_size << " open nodes, " << stats->node
	 << " nodes and " << stats->total_time << " s so far, incumbent "
	 << solution->getObjective() << std::endl;
  }
  else {
    stack_setup(stack,stats,id,n);
    if (root_fixes) {
      // the subtree below a node: its depth is the number of fixes (fix()
      // counts them), so -x and -d count from the root, not from this node
      // (dnntest -x uses the same i:j spelling)
      for (char *tok = strtok(root_fixes, ","); tok; tok = strtok(NULL, ",")) {
	int fi, fj;
	if (sscanf(tok, "%d:%d", &fi, &fj) != 2 || fi < 0 || fi >= n || fj < 0 || fj >= n ||
	    stack->data.p[fi] != NONE || stack->data.pinv[fj] != NONE) {
	  cerr << "-F: bad or repeated assignment " << tok << std::endl;
	  exit(1);
	}
	stack->data.fix(fi, fj);
      }
      // The relative gaps that pick the parameter level and gate the DNN
      // bound are measured against the root's bound, which this run never
      // computes in its search: compute it here, as a full run would at
      // depth 0 (it was left at 0 before Nov 9, so a -F run picked its
      // parameters as if every node were far from the root).
      stats->set_root_bound(root_qpb_bound(q, params));
      cout << "root bound " << stats->root_bound << " (for the relative gaps of a -F run)" << std::endl;
    }
  }
  branch_init(n);

  if (read_symmetry(c->sym_file, n)) {
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
    estimated_stats = estimate(q,stack,solution,params, c->nlevel,ep);
    delete estimated_stats;   // estimate() reports as it goes; we own the result
    delete ep;  ep = NULL;
  }
  else {
    //    use_EVB3(q);
    bnb(q,stack,solution,stats,0,0,params, c->nlevel,
	max_depth,false,-1);

    if (bnb_state_file) {
      // with -d this is the frontier; after -t the unfinished nodes;
      // otherwise an empty stack with the final statistics
      write_state(bnb_state_file, stack, solution, stats);
    }
    if (bnb_stopped)
      cout << "stopped at the time limit with " << list_length(stack)
	   << " open nodes; resume with -r " << bnb_state_file << std::endl;

    print_results(cout,stats,solution,scale_factor);

    // nb 1/11/2016: this no longer works because output goes to stdout only.
    //ofstream outf;
    //outf.open("results.out");
    //print_results(outf,stats,solution,scale_factor);
    //outf.close();


  }
  
  branch_shutdown();
  list_clear(stack);
  delete q;
  delete stats;
  delete solution;
  delete [] params;

  delete c;
  delete [] resume_file;
  return bnb_stopped ? 3 : 0;
}

