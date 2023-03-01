/*
============================================================
bnb_console()  : testing program.
                   bound code.

============================================================
*/

#ifndef CONSOLE_H
#define CONSOLE_H

#include "qap.h"
#include "assign.h"
#include "link1.h"
#include "solution.h"
#include "stats.h"
#include "branch.h"
#include "util.h"
#include "bnb.h"
#include <fstream>

class ConsoleParameters
{
public:
  int nlevel;
  char *log_file;
  char *command_file;
  char *sym_file;
  char *input_file;
  char *param_file;
  std::string log_prefix;
  bool log_header;

  ConsoleParameters() {
    param_file = NULL;
    log_file = NULL;
    command_file = NULL;
    sym_file = NULL;
    input_file = NULL;
    log_header = false;
  }

  ~ConsoleParameters() {
    if (param_file) delete[] param_file;
    if (log_file) delete[] log_file;
    if (command_file) delete[] command_file;
    if (sym_file) delete[] sym_file;
    if (input_file) delete[] input_file;
  }

  void set_log_prefix(const char *prefix) {
    log_prefix = prefix;
  }
};

void start_log(ConsoleParameters *c);
void log_bnb_run(ConsoleParameters *c, QAPAssignment &node, QAPSolution *sln, QAPStatistics *stats);
void run_at_node(QAP *q, QAPAssignment &node, QAPSolution *sln,
		 QAPStatistics *stats, BNBParameters *params,
		 ConsoleParameters *c);
void run_all_branches(QAP *q, QAPAssignment &node, Node *&stack,
		      QAPSolution *sln, QAPStatistics *stats,
		      BNBParameters *params, ConsoleParameters *c, const char *prefix);

void level1_bounds(QAP *q,int maxFWier_bnd);
void bnb_console(QAP *q,Node *&stack,QAPSolution *solution,
		 QAPStatistics *stats,BNBParameters *params,
		 ConsoleParameters *c);

#endif
