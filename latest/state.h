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

 Format (version 9; versions 1-8, which lack the warm starts, the DNN
 counts, dnn_ext, the cut pools, the DNN state, the per-level DNN counts,
 the settings record or the @path DNN references, are still read):
   QAP-STATE 9
   n <n> levels <stats.n> strategies <stats.m>
   incumbent <obj> <found> <perm[0..n-1]>
   root_bound <value>
   counters <node> <fathom> <enumerate> <fw_iter>
   times <total_time> <bound_time>
   settings <hash> [NAME=VALUE ...]   (version 8: settings_record() in settings.h;
                                       read_state refuses a file whose hash differs
                                       from this run's unless --allow-settings-change)
   n_level ... / f_level ... / l_level ... / e_level ...  (stats.n+1 each)
   t_level ...                                            (stats.n+1)
   n_strat ...                                            ((stats.n+1)*stats.m)
   t_strat ...                                            (stats.m)
   rg_sum ... / rg_sumsq ...                              (stats.n+1 each)
   dnn_counts <nodes> <fathomed>                          (DNN bound, bnb -x)
   dnn_time <CPU seconds>
   pop_f_level ... / dnn_n_level ... / dnn_f_level ... / dnn_g_level ...
   dnn_i_level ... / dnn_t_level ... / cut_n_level ... / cut_f_level ...
   cut_t_level ...                    (stats.n+1 each, version 7; see stats.h:
                                       popped, DNN computed/fathomed/given
                                       up/iterations/CPU, cut phase run/
                                       fathomed/CPU, by level)
   warm_starts <count>
   warm <m> <p[0..n-1]> <X[0][0..m-1]> ... <X[m-1][0..m-1]> <dnn_ext> <pool> <dnn>
   ...                         (a branched node's X, see WarmStart in assign.h;
                                pool: dnn_pool_encode's token, "-" for none;
                                dnn: the node's final DNN state for its
                                children, written with ship.bits = 16|32|64,
                                "-" otherwise; format in state.C. Version 9:
                                with QAP_DNN_WARM_DIR set the field is
                                "@<path>", the token being the first line of
                                that file, path relative to this file's
                                directory unless absolute; a checkpoint and
                                its store directory travel together)
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
