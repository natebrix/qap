# Common ground for the branching_lab agents (Nov 19 2026)

You are one of several agents working in parallel for Nathan Brixius on the QAP
branch-and-bound solver. Read `/home/user/qap/CLAUDE.md` first. The question for all of
us: how should the solver choose where to branch, and which nodes can it drop without
bounding, now that its bound is the certified DNN bound with a cut phase? Nathan's own
framing (Nov 19): in 2000 "rule 2" read the QP bound's dual matrix U and told us for free
how each child would behave; the DNN bound has no such thing, and strong branching is
expensive.

## Where things are

- `/home/user/wt-cognac`: a read-only checkout of the main branch `cognac` (commit
  0e6c2ecf). **Read there, never write there.** It has the solver (`latest/`), every run's
  node trace, and the main session's branching studies. Read these before you start, so
  you do not repeat them:
  - `experiments/line_signals/README.md` (Nov 18): cheap signals do not rank lines (the
    marginal score, rule 2, U sum, QPB lookahead); the children's first DNN solve does
    (Spearman 0.8); 30 warm pooled iterations per child are enough; a lookahead rule
    `QAP_LA` exists in an experiment build (`experiment_build.patch`); the DNN certificate's
    child bounds are flat and useless.
  - `experiments/tai60b_postmortem/README.md`: tai60b took 331,149 nodes and 3,597
    core-hours against an estimate of 29,357 and 279. The tree is two chains, the best
    known solution's path `44:45` (155,549 nodes) and its mirror `44:58` (175,311), and the
    subtrees hanging off them; two branchings (depth 3 and depth 7) decide two thirds of
    the cost.
  - `experiments/tai60b_lines/`, `experiments/rule9_bottom/` (rule 9 at every depth:
    -19% nodes on tai50b), `experiments/tho30_solve/line_truth/`,
    `experiments/tho40_line_truth/` (the "truth" sets: nodes where 8 candidate lines each
    had every child's subtree solved to the end), `experiments/xbar_lines_sweep/`,
    `experiments/wil50_feasibility/tuning/lines/README.md`.
- `/home/user/qap`: this lab's checkout (branch `claude/cool-babbage-ilnogj`). Write only
  under `experiments/branching_lab/<your directory>/`. `experiments/bound_lab/` here has
  the Nov 15-17 bound work (`README.md` is its index); `bound_lab/tai60b/nearopt.py` and
  `bound_lab/swapcuts/` are relevant to some briefs.
- Node traces (one line per node), in `/home/user/wt-cognac/infra/phase1/runs/<run>/nodes.txt.gz`:
  `tai60b-261117` (331,149 nodes), `tai50b-261116` (8,216, rule 2 below depth 6),
  `tai50b-r9` (6,629, rule 9 everywhere), `tho40ab-32` (a tho40 depth-2 subtree, 22,650
  nodes), `nug30-trace-a`, and others. Format:
  `node <depth> <fixes> <qpb> <dnn> <cut> <cutoff> <outcome> <rule> <line> <children> <cpu s>`,
  with `-` where a field does not apply. Outcomes: B branched, D fathomed by the DNN first
  solve, C by the cut phase, Q by QPB, P discarded when popped, E enumerated (a leaf).
  A fix is `row:column`; a line is `rN` or `cN`. **Work out for yourself, against
  `QAPSOLN/<instance>.sln` and the instance file, which matrix the rows index and what a
  fix means; do not take it from this page.** The recipes are in
  `infra/phase1/recipes/` (`tai60b.par`, `tai50b.par`, ...), region files beside them or in
  `experiments/bound_lab/`.
- Primary sources: `docs/sources/` (the thesis `nb_thesis.txt`: branching rules in
  chapter 3, Knuth's estimator in chapter 4). Read the source before reasoning from
  memory about the 2000 design.

## Ground rules

- This machine has 4 cores and 15 GB, shared by three agents. Python analysis of traces is
  cheap; keep it to one or two processes. If your brief involves solver runs, one
  `bnbdrv` at a time, single-threaded, and build only in your own worktree
  (`git -C /home/user/wt-cognac worktree add --detach /home/user/wt-<name> 0e6c2ecf`), never
  in `/home/user/qap/latest` or `/home/user/wt-cognac`.
- Mark every number **[m]** measured (name the script and its output file) or **[e]**
  estimated. Keep scripts and their outputs in your directory so a result can be rerun.
- A result that would let the solver drop nodes is a claim about a proof. Validate it the
  hard way (brute force on small instances, against every optimal solution) before you
  measure what it would save, and say exactly what you validated.
- Report what does not work as plainly as what does. A clean negative with the reason is
  a result. Give a kill rule before each experiment and keep to it.
- Do not commit: state files, anything over 2 MB, anything containing an AWS account
  number (an S3 bucket's name in this project holds one). Do not push; the coordinating
  session commits. No pull requests, no AWS, nothing outside this container.
- Deliver `experiments/branching_lab/<your directory>/README.md`: the answer first (two or
  three sentences), then the evidence in tables, then what you did not test, then the
  cheapest next step if it is worth one. Write for a reader who knows the solver but has
  not seen your work. Your final message to the coordinator should be that README's
  answer and caveats, not a log of what you did.
