# Branching lines at wil50 clear misses (set up Oct 31 2026; to run on the laptop)

Question: at the depth-3 nodes the DNN bound misses by a lot (about 95% of the
estimated wil50 tree is below them, ../../NOTES.md section 4), which line
(facility = row, location = column) leaves the smallest subtree, and do the
rules we have pick it: the most eccentric location (rule 6), QPB lookahead
(rule 4), or a new candidate, the DNN primal's own marginals x_ij?

`lines.py` (see its docstring) bounds, for each parent of `parents.txt`, every
child of about 9 candidate lines (eccentric; the top 2 rows and 2 columns by U
sum, rule 4's candidates, with lookahead's own pick among them marked; the row
and column whose x_ij are most concentrated; 2 random rows and 2 random
columns), with a fixed DNN budget per child, and scores each line by
sum over children of 47 ** (gap / 130), the size of the subtree a child
leaves if the bound gains about 130 per level (wil50, ../../NOTES.md). It also
reports the rank correlation of x_ij (known before branching) and of the
child's QPB with the child's gap.

The parent run uses `dnntest -X`, which prints the marginals x_ij of the last
primal iterate (new, Oct 31).

## Running it

    VECLIB_MAXIMUM_THREADS=1 python3 experiments/wil50_feasibility/tuning/lines/lines.py \
        qaplib/wil50.dat --ub 48816 \
        --parents experiments/wil50_feasibility/tuning/lines/parents.txt \
        --work experiments/wil50_feasibility/tuning/lines/work --jobs 8

`--dry-run` lists the candidate lines (after the parent runs) without bounding
children. Cost: about 9 lines x 47 children x 3 parents = 1,270 child bounds
at 300 iterations each (about 150 s each on the 2.1 GHz cloud Xeon), plus the
child QPB bounds through the console: about 50 CPU-hours there, so roughly
2-3 hours on the M2 Pro with 8 workers. It caches every child in `--work`
and resumes after an interruption.

## Check on nug30 (cloud)

%%NUG30%%
