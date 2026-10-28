/*
============================================================

 Saving and restoring branch-and-bound state.

 A state file holds everything bnb() needs to carry on: the open
 nodes in stack order, the incumbent, the root bound (which
 pick_param uses for every node) and all statistics. Doubles are
 written with 17 significant digits, so a resumed search makes the
 same decisions bit for bit.

 It is plain text, so other tools (bench/parallel.py) can split a
 state into tasks and merge the statistics of finished ones.

 Format (version 3; versions 1 and 2, which lack the warm starts and
 the DNN counts, are still read):
   QAP-STATE 3
   n <n> levels <stats.n> strategies <stats.m>
   incumbent <obj> <found> <perm[0..n-1]>
   root_bound <value>
   counters <node> <fathom> <enumerate> <fw_iter>
   times <total_time> <bound_time>
   n_level ... / f_level ... / l_level ... / e_level ...  (stats.n+1 each)
   t_level ...                                            (stats.n+1)
   n_strat ...                                            ((stats.n+1)*stats.m)
   t_strat ...                                            (stats.m)
   rg_sum ... / rg_sumsq ...                              (stats.n+1 each)
   dnn_counts <nodes> <fathomed>                          (DNN bound, bnb -x)
   dnn_time <CPU seconds>
   warm_starts <count>
   warm <m> <p[0..n-1]> <X[0][0..m-1]> ... <X[m-1][0..m-1]>
   ...                         (a branched node's X, see WarmStart in assign.h)
   nodes <count>
   node <depth> <bound> <predicted_bound> <nfix> <p[0..n-1]> <history[0..n-1]> <warm>
                               (warm: index into the warm starts, or -1)
   ...
   end

============================================================
*/

#ifndef STATE_H
#define STATE_H

#include "link1.h"
#include "solution.h"
#include "stats.h"

// Write the state atomically (to file.tmp, then rename). Returns false on error.
bool write_state(const char *file, Node *stack, QAPSolution *s,
		 QAPStatistics *stats);

// Read a state written by write_state. n must match the instance; the
// statistics must have the same number of levels and strategies.
// stack must be empty; on return it holds the nodes in the saved order.
bool read_state(const char *file, int n, Node *&stack, QAPSolution *s,
		QAPStatistics *stats);

#endif
