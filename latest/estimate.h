/*
============================================================

est_stats = estimate(q,stack,s,params,nlevel,ep);

Estimates a full branch-and-bound run on QAP q using the
parameters params.  ep is a set of "Estimate Parameters"
which determine how accurate the estimate will be, and
how long it will take.

EstimateParameters:
 start_dive_level:   all nodes in the B&B tree at this level
                     will be collected before starting
                     the estimation process.  Therefore,
                     the information at the first start_dive_
                     level levels of the tree will be 'exact'.

 print_every:        print a simple status message every x nodes.

 print_everything_every:  print the current estimated statistics
                     every x nodes.

 trials:             number of dives to make to compute estimate.

 time_scale:         scale estimated times by this factor to account
                     for platform differences.

 score_exponent:     the higher this is, the more likely that deeper
                     paths will be sampled.  1.5 is usually a good
                     value.

 ignore_depth:       throw out all dives greater than or equal
                     to this depth.  prevents overestimates.
============================================================
*/

#ifndef ESTIMATE_H
#define ESTIMATE_H

#include "bnb.h"
#include "qpbnd.h"
#include "branch.h"

// maybe bnbparameters should have a pointer to one of these.
class EstimateParameters 
{
public:
  EstimateParameters() {
    start_dive_level = 2;
    print_every = 100;
    print_everything_every = 1000;
    trials = 10000;
    time_scale = 1.0;
    score_exponent = 1.0;
    ignore_depth = 1000;
  }

  int ignore_depth;     
  int start_dive_level; 
  int print_every;
  int print_everything_every;
  int trials;
  double time_scale;
  double score_exponent;
};

QAPStatistics *estimate(QAP *q,Node *&stack,QAPSolution *s,
				BNBParameters *params,int nlevel,
				EstimateParameters *ep);

void readEstimateParameters(const char *filename,EstimateParameters *&ep);

int dive_depth(QAPStatistics *s);

#endif
