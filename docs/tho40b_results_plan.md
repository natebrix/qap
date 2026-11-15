# tho40-b: plan for using the results

Nov 15 2026. The shared copy of this plan is a Claude doc that the group comments on; this
file is its record in the repository.

tho40-b will leave the first complete n = 40 tree under the production settings, with every
node traced. This plan turns it into three things: corrected estimates and settings for the
other open instances, a measured target for bound improvements, and a tested path to a larger
run on Spot.

At 02:15 UTC on Nov 15 the run was 9.4 h in, at 559,283 nodes and 3.75 CPU s per node on c7a.
It projects to about 3e6 nodes (2e6 to 5e6) and about 2 days of wall.

## What the run leaves us

Two traces carry almost everything below, and the collect saves neither of them whole.

- **Node trace:** a line per node with its fixes, QPB and DNN bounds, cutoff, outcome,
  branching line, child count and CPU. About 3e6 lines; format in `latest/bnb.C`.
- **DNN trace:** a line per DNN solve with its iterations, whether it gave up or fathomed, and
  whether it started warm.
- **Run records:** per-unit timing, disk and store figures, the coordinator logs, and the
  published per-depth tables.

Before cleanup:

- [ ] Keep the run's traces in S3. `qapctl purge-states` leaves them; nothing else should
  delete the run's prefix.
- [ ] Build `qapctl export` to pull both traces to a durable copy. The collect keeps the node
  trace only under 25 MB gzipped (tho40-b's will be about 60 to 80 MB), and it gathers no DNN
  traces.

## (a) Tuning the other instances

The complete tree is ground truth for the methods that sized wil50, sko42, sko49 and tai50b,
so the first job is to score those methods on it.

1. **Score the estimators.** Cut the finished trace off at depths 2, 3 and 4 and rerun each
   method: the Nov 8 calibrations (`experiments/tho40_estimate/estimate.py`), the wil50
   sampling stages, and the Knuth dives used for tai50b. Each method's error at n = 40
   restates the other estimates with honest error bars, the main input to the quota request.
2. **Fix the cost model.** A node costs about 3.9 CPU s on c7a against 1.39 s on the M2 Pro:
   2.8x, not the 2x assumed. Re-price each instance from a short calibration subtree on c7a;
   nodes at m = 50 cost more again.
3. **Settings tho40 can inform.**
    - Give-up (200 iterations): when a node gave up, did its subtree turn out tiny? That gives
      a per-depth break-even.
    - Branching levels: whether rule 9's levels and their gap thresholds paid, from child
      counts and child fathom rates by level.
    - Deep DNN: whether the DNN bound below depth 27, new in tho40-b, paid for itself.
    - Warm starts: the gain from shipping parent states by depth, at full scale (2x in CPU per
      node on the tho40ab subtree, `infra/phase1/runs/tho40ab-*`).
4. **A reusable tuning set.** 6 to 10 subtrees at depth 3 or 4 across the range of gaps, each
   with its exact size known from the trace and each sized to run in about 30 min on one c7a. A
   settings change then becomes an A/B against known numbers, as tho40ab was. tho40 runs no
   cut phase, so cut settings for wil50, sko42 and sko49 need their own subtrees.

## (b) Bound improvements

The trace measures what a better bound would save, and where in the tree, before we build one.

1. **Value of a better bound.** For each branched node we will know how far its bound fell
   short of the cutoff and how big its subtree turned out. The curve: nodes saved if the bound
   rose by d, at each depth. It says how much a new cut must gain to matter.
2. **Top of the tree.** tho40's root gap is 5.8% (nug30's is 2.8%), and the median gap at
   depth 2 is 2.9%. Subtree size grows about 0.7 orders of magnitude per 1% of gap at depth 2,
   so a bound 0.45% tighter there (about 1,100 units) would roughly halve the tree. Costly
   bounds may pay when used only near the top, for example the cut phase at depth 4 or less.
3. **A bound lab.** Freeze a sample of near-miss nodes as state files: nodes that branched with
   a small shortfall and grew big subtrees. Test candidates offline, scored by fathoms gained
   per CPU second:
    - the DNN bound run to convergence, without the give-up
    - the grid cut phase (tho40 is a 5 x 8 grid and has never run it)
    - CP-ADMM and more cut rounds
    - the newer ideas on our research list
4. **Branching at the top.** An earlier study put choosing the root line by its children's
   bounds at 2 to 5x (`experiments/wil50_feasibility/tuning/lines/`). The finished tree
   measures what each line actually cost.

## (c) Before a larger run, on Spot

Workers run on Spot from here on, decided Nov 15. The work before a large run is making the
launcher and coordinator fit a fleet that grows and shrinks, then proving it on a subtree.

1. **Spot readiness.**
    - A worker already drains on a Spot notice: running units return their open nodes, or
      their last checkpoint within 90 s. That path is tested by SIGTERM in the coordinator
      tests, never by a real Spot notice.
    - A run's slot count is fixed when the run is created, and the coordinator sizes its queue
      reserve and its work-stealing by it. A Spot fleet needs both to follow the workers
      actually present.
    - Capacity: spread the fleet over several instance types and Availability Zones, so one
      Spot pool running dry does not stall the run.
    - Test: a mid-size subtree on 2 to 4 Spot instances with interruptions injected (AWS Fault
      Injection Service can send Spot interruption notices), its totals checked against a
      clean run.
2. **Several instances at once.** Nothing at scale has run on more than one instance yet. Watch
   the SQS, DynamoDB and S3 rates, the coordinator's merge rate, the parent-state store's
   cleanup under more churn, publish time with millions of nodes, and the tree viewer.
3. **Coordinator scaling.** The known gaps (`coord/README.md`): counting units without a
   DynamoDB scan, ledger snapshots, and an auto-scaling group in place of `launch.sh`. With
   Spot, the auto-scaling group moves up the list.
4. **Resume at the time limit.** A run that reaches `--limit` shuts down, and launching it
   again resumes it from the store. tho40-b may test this; if not, a small run will.
5. **Quota and budget.** Request Spot and on-demand vCPU quota in the chosen region, sized from
   the restated estimates. Budget alerts, and per-run cost from instance and disk tags plus a
   recorded hourly price.
6. **Data volume.** Full traces at 1e7 nodes run to GBs. We keep them whole: both the node
   and the DNN trace (decided Nov 15), so the export and the S3 budget must allow for that.
7. **Standard of proof.** Already recorded: the commit, the settings and a verified incumbent.
   A possible addition is an independent recheck of a random sample of fathoms and give-ups,
   agreed with the group before the run.

## Order of work

The Spot test and the quota request come last, so both are sized from estimates restated
against the finished tree.

1. **Now, while tho40-b runs:** build the trace export (node and DNN traces); write the
   analysis scripts against the tho40ab traces; design Spot support (a live slot count, a
   changing fleet).
   - Gate: tho40-b finished and its traces exported.
2. **When tho40-b ends:** collect and export the run, keeping its traces in S3; score the
   estimators against the complete tree; draw the bound-value curve and pick the tuning set.
   - Gate: estimates restated for wil50, sko42, sko49 and tai50b.
3. **Then, toward a larger run:** the bound lab on frozen near-miss nodes; the Spot test on
   several instances with injected interruptions; the quota request and budget, and the choice
   of the next instance.

## Questions for the group

- Which instance comes next: wil50, sko42, sko49 or tai50b?
- What standard of proof do we want for a published result?
- Should the coordinator's own instance stay on-demand while the workers run on Spot? Losing it
  pauses the run until a relaunch resumes it.
- What did the 2000 nug30 run teach that we should plan for, beyond disk?
