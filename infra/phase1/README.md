# Phase 1: the coordinator on one instance (Nov 8)

Phase 1 of `docs/distributed/plan.md`: the coordinator and worker of `coord/`
on S3, SQS and DynamoDB, on one EC2 instance. Its acceptance test is nug30
with `-u 6124` on one 16-core instance: the optimum proven, and node, fathom,
FW and DNN totals equal to `bench/parallel.py` on the same machine and
settings. Then nug30 without `-u`, which must finish on 6124 with the
incumbents passed through S3.

`launch.sh` (on the laptop) creates the storage and an instance role, uploads
the committed source and starts the instance. `node.sh` (on the instance) does
the rest and terminates the instance at the end. Nothing listens on the
instance: everything goes through the store.

## What the instance does (`node.sh`)

1. It records the machine (type, CPU, physical cores, AMX, memory) and
   installs the toolchain of Phase 0 on Amazon Linux 2023 (gcc, Python 3.11,
   numpy, scipy-openblas32), plus boto3.
2. It builds the solver and runs the gates: `make bruteforce resume parallel
   regression dnncheck coord` (5 minutes on this container's 4 cores, 2.5 of
   them the `coord` scenarios). A failed gate stops the instance before any
   unit runs (`--gates` changes the list).
3. If the run is new, `qapctl init` with the recipe: the frontier pass runs
   here, so the whole tree comes from one build (the root bound alone differs
   between Accelerate and OpenBLAS rounding, `experiments/nug30_solve`). The
   recipe as used, with `--no-ub` applied, goes to `runs/RUN/recipe.sh`.
4. The coordinator and one worker with a slot per physical core, until the
   run is finished. Either is restarted (up to 5 times) if it dies first; the
   coordinator replays what was merged. If the coordinator fails a fifth
   time, the worker is drained and the instance stops, rather than idling
   until `--limit`.
5. The check: `bench/parallel.py` with the same recipe on the same cores, then
   `compare.py`, which writes `compare.txt`: with an incumbent, totals equal and
   the optimum; without one, the optimum only (the totals then depend on when
   incumbents arrive). `--no-reference` skips the `parallel.py` run.
6. `out/` (steps, logs, results) goes to `runs/RUN/hosts/HOST/out.tgz` after
   every step and every 5 minutes. At the end the instance shuts down, which
   terminates it. It also shuts down at `--limit` minutes (120) whatever happens.

A Spot interruption notice (the metadata service's `spot/instance-action`,
polled every 5 seconds) drains the worker: each running unit returns its open
nodes, or its last checkpoint if a node does not end within 90 seconds, before
the instance is reclaimed. `launch.sh RUN` again resumes the run on a new
instance.

## Recipes

A recipe (`recipes/NAME.sh`, committed: the instance builds from `HEAD`; or a run's
`runs/RUN/recipe.sh`, which `start.sh` writes and the user data carries) says
what to solve: `INSTANCE`, `PARAMS`, `OPT` (the known optimum the check
expects), `UB` (the incumbent given to the solver), `EXTRA` (solver flags),
`ENV` (the solver's environment) and `COORD` (unit sizing for `qapctl init`).

| recipe | what | CPU [e] |
|---|---|---|
| `nug17` | the smoke test: QPB only, one-second units, a few dozen of them through S3/SQS/DynamoDB | 20 s |
| `nug30` | the acceptance run: the Nov 8 record recipe (rule 9, DNN at depths 1-7) with `ship.bits = 64` (`recipes/nug30.par`), one-minute units | 1,300 s on the M2 Pro; more on Graviton |
| `nug30trace` | `nug30` with `--trace`: every node's line (`QAP_NODE_TRACE`) goes to the store with its unit, the frontier's under `runs/RUN/frontier/`, for the tree viewer; launch with `--no-reference` | as `nug30` plus a few % |
| `tho40sub` | the tho40 smoke test: the depth-3 subtree 26:3,31:7,37:0 of `tho40t.dat` under `experiments/estimates/tho40/tho40.par` (110 nodes, 107 fathomed in the Nov 9 A/B), one-minute units; prices a tho40 node on the instance | 197 s on the M2 Pro |
| `tho40` | the whole tho40 tree from an incumbent of 240,518 (the best known 240,516 + 2, so the run must find 240,516 itself) under `recipes/tho40.par` (the estimate's settings with ship.bits = 32 and the DNN bound at every depth), 20-minute units, with `--trace` (every node's line in the store, a few hundred MB); `--no-reference` and a `--limit` of days (the recipe's comment has the launch line) | 0.9e6 to 4e6 s on the M2 Pro (`experiments/tho40_estimate/`) |

`ship.bits = 64` instead of the record's 32: with 64 a node restored from a
state file starts its DNN from its parent's exact doubles, so a tree cut into
units at other places than `parallel.py` cuts it is the same tree. On nug20
with this recipe both 64 and 32 gave the coordinator and `parallel.py`
identical totals (174 nodes, 163 fathomed, 29,550 FW iterations, the serial
run's), but only 64 is meant to guarantee it. The DNN iterations can still
differ by a few: with the per-level counts of state version 7, nug24 with this
recipe through the coordinator (5-second units) had the serial run's tree at
every depth but 3,165 DNN iterations against 3,175 at depth 3 and 1,406
against 1,404 at depth 4 (Nov 9). Some setting of the recipe keeps state
inside a process that state files do not carry (`make shipcheck`'s smaller
set of settings is exact); `compare.py` checks the tree, not iterations.

## Starting a run (`start.sh`, Nov 16)

**From a spec file (Nov 19):** `start.sh --spec infra/phase1/specs/NAME.spec` takes every answer
from a committed file and asks only "Launch?" (`--yes` answers it); `ARMS` in it makes an A/B, one
run per arm, `RUN-ARM`, identical but for the params file, the difference shown before the launch.
Every new run's record (`runs/RUN/`: `recipe.sh`, `launches.tsv`, `spec.txt`) is committed and pushed
once it is launched (`--no-commit`: not). `infra/phase1/specs/README.md` has the keys and the specs
(`nug30-smoke`, `tai50b-est`, `sko42-zone`). Tested in a throwaway clone with a bare origin, against
the AWS stand-in and moto: the sko42 A/B launched two runs whose recipes differ only in `PARAMS`,
each with its group, and its record was committed and pushed; a spec with an unknown key, or not
committed as it is, is refused.

**Readouts (Nov 19):** each run's `readout.md` (`coord/readout.py`: the verdict, the warnings, the
record) is written by its coordinator at the end. `collect.sh RUN` makes a fresh one into the run's
record and ends with a verdict on itself, in `runs/RUN/collect.txt`: `collect complete`, or
`collect INCOMPLETE: no traces | node trace | readout`. It refuses to start with under 30 minutes
left on the SSO sign-in, and takes the node trace from the trace export instead of reading every
unit a second time (`docs/review_nov18_runs.md`, item 4). With no one collecting,
`.github/workflows/phase1-readouts.yml` writes every run's readout and a board of the runs to the
branch `phase1-readouts` hourly (`infra/phase1/readouts.py`; `git show
origin/phase1-readouts:board.md`), with a role that only reads (`gha_setup.sh`: `qap-gha-readouts`,
environment `aws-phase1-read`, no reviewer). A scheduled workflow runs from the default branch, so
the file must be on `master` too.

**Since Nov 19 a run is a coordinator and a worker group** (`docs/review_nov18_runs.md`, item 3):

- **The coordinator is on demand,** a small instance of the workers' CPU (`c7a.2xlarge`). It
  bounds the root, merges and runs a few units. It is never reclaimed, so a run is never lost with
  its root (tai50b-est lost 13 minutes of a 72-minute root to Spot on Nov 18; tai60b's takes hours).
- **The workers are an EC2 Auto Scaling group,** `qap-RUN` (`launch.sh --workers N`): Spot from
  the cheapest pools with capacity over the three types and every zone of the region.
  - It is created at 0 instances; the coordinator raises it to N once the run exists, so no worker
    waits through the root.
  - A reclaimed worker drains on its 2-minute notice and the group starts another in its place;
    `collect.sh` shows the group, its last activities and the interruptions it answered.
  - A worker whose run is over, or whose build or gates fail, leaves the group (terminating itself
    and decrementing it), so the group does not start another.
  - The coordinator sets the group to 0 when the run ends, and a scheduled action does so at the
    time limit whatever happens.
  - No capacity rebalancing: it would terminate a worker without the notice the worker drains on.
- **The root cache** (`coord/coordinator.py`, `frontier_key`): a frontier pass is stored under
  `runs/.cache/frontier/KEY/` and read back by any later init with the same solver binary, CPU,
  data, flag files, flags, incumbent, depth, traces and resolved settings (the run class aside). A
  relaunch after an interruption during init, or another estimate of the same problem, skips the
  root. `QAP_ROOT_CACHE=off` runs it regardless.
- `cleanup.sh RUN` deletes the group and its launch template before terminating the instances (or
  the group would replace them). A run started before Nov 19 has no group: adding instances to it
  launches `--worker-only` instances as before.
- The instance role may change and leave `qap-*` groups (`autoscaling:UpdateAutoScalingGroup`,
  `SetDesiredCapacity`, `TerminateInstanceInAutoScalingGroup`) and create and delete `qap-units-*`
  queues. The first launch after Nov 19 comes from the laptop, which rewrites the role's policy. The
  GitHub workflow's role can stop and collect such runs once `gha_setup.sh` is run again; it cannot
  create a group, so its launch stays one instance at a time.
- **Tested here** (Nov 19), as the earlier launches were, against a stand-in for the AWS CLI that
  logs every call and keeps the groups' state, and moto: `start.sh --dry-run --yes` (the launch
  template and group shown in full), `launch.sh` creating and then updating a group, and `node.sh`
  in its local mode with a coordinating node and a group worker on nug17 (the group raised to 2
  after init, set to 0 at the end; the worker left it after 17 units; check PASS), a worker that
  booted after its run was over (it left before building), and a worker drained by a Spot notice
  (it did not leave, so the group would replace it; the run finished on the coordinator). Not
  tested: the real Auto Scaling API, its IAM, and a real Spot replacement. The first launch tests
  them; `collect.sh` shows the group's activities.

One command on the laptop starts a run of a QAPLIB problem:

```
cd ~/dev/qap && git pull origin cognac
infra/phase1/start.sh tho40 --dry-run
infra/phase1/start.sh tho40
```

It checks the AWS session first, then asks for each setting with a default (Enter keeps it),
shows the whole launch with the vCPU quota left and today's Spot price, and asks before
launching:

| setting | default |
|---|---|
| problem | the command line's (another name restarts with its defaults) |
| run name | `PROBLEM-YYMMDD`, the UTC date: `tho40-261117` |
| params file | the one `recipes/defaults.txt` names for the problem, its next run's settings (tho40: `recipes/tho40c.par`; tai50b: `recipes/tai50b.par`), else `recipes/PROBLEM.par`, else the estimate's, `experiments/estimates/PROBLEM/PROBLEM.par` (wil50: `experiments/wil50_feasibility/fullrun/fullrun.par`); it must be committed |
| incumbent | the best known (`QAPSOLN/PROBLEM.sln`) + the step of the objective's grid: + 2 where every value is even, else + 1 (`objstep.py`); `none` for no incumbent |
| worker type, market | `c7a.16xlarge`, or the problem's `type=` in `recipes/defaults.txt` (tai60b: `m7a.16xlarge`); Spot |
| worker group's types | the worker type and its siblings on the same CPU (AMD EPYC 9R14): `c7a.16xlarge,m7a.16xlarge,r7a.16xlarge` |
| worker instances | 1 (it shows how many the quota has room for) |
| slots per worker | the type's cores, or as many as its memory holds at the problem's `mem=` GiB a solver process, less 4 GiB (tai60b: 3 GiB, so 64 on an m7a.16xlarge's 256 GiB, 41 on a c7a.16xlarge's 128); `slots=` fixes it. Sized by the first type, which has the least memory of the three |
| coordinator | `c7a.2xlarge` on demand (8 cores, 16 GiB): its cores but one run units, by the same memory rule |
| time limit | 72 hours: the coordinator shuts down, and the group goes to 0 |
| on request | the data file and solver flags (`qaplib/PROBLEMt.dat` and `-s SYM/PROBLEMt.sym` where the grid's symmetry file exists, else `qaplib/PROBLEM.dat`), the disk (100 GB), the unit sizing (tho40-b's 20-minute units, `--trace`), public or private status |

Then it writes the run's recipe to `runs/RUN/recipe.sh` and calls `launch.sh` once per
instance: the first coordinates (`--recipe-file`, `--run-slots` the sum of every instance's
cores, so the queue reserve and steals cover them all), the others are `--worker-only`. All
get `--no-reference`. Commit `runs/RUN/` afterwards (`recipe.sh`, `launches.tsv`; the files
that name the bucket are ignored). For a run that exists, it adds instances with the run's
own settings: the first one coordinates if no instance does (after a Spot interruption or
the time limit), so the same command resumes a run.

- **The incumbent: the next value above the best known**, so the run must find the best
  known solution itself (`status.json`: found by search) and prove that nothing is better.
  The solver fathoms a node whose bound exceeds the incumbent less its tolerance: 1.95 when
  every objective value is an even integer, 0.95 when they are integers (`latest/qap.C`,
  `set_fathom_tolerance`). `objstep.py` makes the same test on the data file (it agrees with
  the solver on all 192 QAPLIB files: 121 even, 70 integers, 1 not integers) and gives the
  step:
  - **+ 2 where every value is even:** tho40, sko42, sko49, wil50 (symmetric matrices, even
    diagonal products), as the 2000 nug30 solve was given the optimum + 2 and tho40-b
    240,518. + 1 there is off the grid, and `bnbdrv` rounds such an incumbent up to the next
    attainable value (`latest/driver.C`), so it would act as + 2: without that, the 1.95
    tolerance would put the cutoff at the best known - 0.95 and could fathom the node holding it.
  - **+ 1 where values are integers but not all even:** the tai-b instances, whose matrices are
    not symmetric (tai50b 458,821,518).
  - Data that are not integers have no grid (tolerance 0): + 1 is a margin there.
  An incumbent typed off the grid (rounded up), at or below the best known, or past the next value is
  pointed out. The estimates
  (`experiments/estimates/*/settings.sh`) keep the best known itself, the tree a run must
  search.
- **Run names:** a run the website shows is named after the QAPLIB problem and the UTC date
  of its launch, `tho40-261117`; a second the same day gets a letter, `tho40-261117b`. The
  website takes the instance's name from the data file, so the name is for people. Earlier
  runs (`tho40-b`, `nug30-a`) keep theirs.
- **AWS credentials:** none are stored. The script runs on the AWS CLI's IAM Identity Center
  (SSO) session: `aws sts get-caller-identity`, and `aws sso login` when there is no session
  or it has expired (once per laptop, `aws configure sso` sets a profile up; `--profile P`
  picks one). The CLI keeps its own short-lived token in `~/.aws/sso/cache`; the instances use
  their role. A sign-in that is not SSO, or keys in the shell's environment, is pointed out.
  The output shows the account as `<account>` (`launch.sh` too), so it can be pasted.
- **A run on several instances:** the coordinator sizes its queue reserve (2 units a slot)
  and its steals by the run's slot count, fixed at `qapctl init`. If an instance is lost, its
  slots stay counted, and the coordinator splits units for them: relaunch with the same
  command. `slots busy` in the results is relative to the run's slots.
- **Tested here** (Nov 16) against a stand-in for the AWS CLI and moto: a dry run with every
  default (`--yes`, signed out first, so through `aws sso login`), and a launch of two
  instances with typed answers: the recipe went in the coordinator's user data unchanged,
  `QAP_RUN_SLOTS=128`, the second instance a worker. `node.sh` in its local mode from a clean
  archive with a recipe file (nug17 from 1,734, 2 slots, the run sized for 8): init at 8
  slots, 1,732 found, the check passed. Dry runs of tho40, sko42, sko49, wil50, tai50b, nug17 and nug30
  give the defaults above (+ 2 on the even ones, + 1 on tai50b), and typed incumbents below,
  at and past the next value give their notes. 119 of the `.sln` files end their lines with `\r`, which first gave
  tai50b no incumbent: the value is now read as digits, and the launch stops if the incumbent
  was never asked.
- **Estimate runs** (Nov 18, `docs/estimates_plan.md`, `coord/dives.py`). To launch one, answer
  "Change the other settings?" with y and give "Estimate run" as `SOLVE_HOURS:DIVE_SECONDS[:B]`,
  e.g. `2:1800`. Every instance then gets `QAP_ESTIMATE` (seconds), and `node.sh`:
  - runs `qapctl estimate-flow` in the coordinator's place: the solve for SOLVE_HOURS of its
    wall time, then the pause, its units returned, and run `RUN-dives` with DIVE_SECONDS of
    survivor dives per open unit;
  - runs `worker --follow` on every instance: the workers move to `RUN-dives` on their own;
  - skips the reference run.

  At the end `out/estimate.txt` (and the store's `runs/RUN-dives/results.txt`) holds the
  estimate: the solve's nodes and CPU so far, the rest estimated, by depth, with a 5-95%
  interval. The solve stays paused, so `qapctl resume RUN` (or a relaunch without the estimate)
  goes on from where it stopped.

  Unless the estimate is not worth its CPU (Nov 19, the decision point: `coord/README.md`): when
  the dives run's pilot (a random tenth of its units) says the dives would cost over 15% of the
  rest of the tree, the dives stop and the solve goes on to the end on the same instances, within
  the time limit; the readout says SOLVING, then COMPLETE. A spec's `ESTIMATE_DECIDE` sets the
  share (`off`: the dives always run to the end), as `QAP_ESTIMATE_DECIDE` on the coordinator. A Spot relaunch resumes the flow at its step. Make the time
  limit cover SOLVE_HOURS plus the dives: the dives take about DIVE_SECONDS per queued unit
  over the slots.

- **Runs from a paused run** (Nov 19, a spec's `FROM=SRC`; `infra/phase1/specs/README.md`). SRC is a
  paused run whose units are back (an estimate run's solve: tho40-est2-261119). Its unit states and
  parent states are read, never changed, so it can still be resumed. Two kinds:
  - `DIVES=DIVE_SECONDS[:B[:SAMPLE[:W]]]`: a dives run of SRC's open nodes, a second estimate with no
    solve phase (`qapctl dives SRC --as RUN` on the coordinating instance, `ESTIMATE_ENV` for the
    estimator's settings). Its readout is its own ("a second estimate of SRC"); `collect.sh RUN`
    fetches its estimate and dives.
  - `SAMPLE=DEPTH:COUNT,...` with `PARAMS` or `ARMS` (and `ARMS_ENV`): a sampled run, one per arm,
    whose first units are a sample of SRC's open nodes (`qapctl init --sample-from SRC`;
    `coord/sample.py`), drawn in proportion to the square root of their subtrees' estimates in a dives
    run of SRC (`SAMPLE_BY`; SRC's own by default) and solved to the end under the arm's settings, the
    same nodes in every arm. `collect.sh RUN-ARM` fetches each arm's `sample.json` and node trace;
    `experiments/estimates_nov19/paired.py runs/RUN-ARM ...` compares them node by node (the CPU of
    arms is comparable only side by side: tai50b's same nodes cost 13-20% less in one run than in
    another). No acceptance check (its results are the sampled subtrees', not a solve's).

  Tested here (Nov 19) with the AWS CLI stand-in and moto: `start.sh --dry-run --yes` on
  `tho40-est3.spec` and `tho40-tune.spec` (the coordinator's environment, the four arms' recipes,
  `cross5`'s own environment), and `node.sh` in its local mode on nug17 in both modes: the dives run
  made from a paused run and coordinated to the end, its readout its own and the source's unchanged;
  the sampled run made from 9 open nodes and solved, with no acceptance check. `tests/sample_test.py`:
  two arms with the same params on a nug24 sample (the DNN at every depth: the parents' states come
  along) keep the same sample and the same subtrees, node for node.

## Running it (laptop)

Once: the AWS CLI and a login (as for Phase 0), and boto3 for `python3`
(`pip3 install boto3`), which `qapctl aws-setup` and the run lookups use.
The instance needs a vCPU quota of at least its own vCPUs in the region (16
for a c7g.4xlarge; Spot quota `L-34B43A08`, on-demand `L-1216C47A`). The
launcher checks it first.

```
git pull origin cognac
infra/phase1/launch.sh nug17-a --recipe nug17 --dry-run    # the AWS calls it would make
infra/phase1/launch.sh nug17-a --recipe nug17               # the smoke test
infra/phase1/collect.sh nug17-a                             # any time: instances, status, steps, results
infra/phase1/launch.sh nug30-a --recipe nug30               # the acceptance run
infra/phase1/launch.sh nug30-b --recipe nug30 --no-ub       # and without the incumbent
git add infra/phase1/runs && git commit -m 'phase 1: ...' && git push
```

- **Options:** `--recipe-file F` (a recipe outside `recipes/`, as `start.sh` writes; it
  goes in the user data, and the files it names must be committed), `--run-slots S` (the
  run's slots when it has more instances than this one), `--type c7g.4xlarge` (16 Graviton3 cores; any type whose
  family ends in `g`, `gd` or `gn` gets the arm64 image), `--region eu-north-1`,
  `--on-demand`, `--worker-only` (more capacity for a run another instance
  coordinates; at most one coordinating instance per run), `--private` (no
  public status), `--no-reference`,
  `--gates "..."` or `--gates none`, `--limit 120` (minutes).
- **Time [e]:** 5-10 minutes of setup and build, about 5 of gates, then the
  run: under a minute for nug17; for nug30 on 16 Graviton3 cores perhaps
  5-10 minutes for the coordinator's run and as long again for `parallel.py`.
- **Cost [e]:** a c7g.4xlarge is a quarter of the plan's c7g.16xlarge: about
  $0.06-0.15 an hour on Spot in eu-north-1 and $0.58 on demand (section 10.1
  snapshot), so cents per run. Storage, queues and tables cost cents a month.
- **What it creates,** once per account and region, kept for later runs:
  - the coordinator's storage, by `qapctl aws-setup`: a private bucket
    `qap-coord-<account>-<region>` (also the source tarballs under `src/`) and
    the DynamoDB tables `qap-runs` and `qap-units` (on-demand capacity). Each run
    has its own SQS queues, `qap-units-RUN-large` and `-small`: created with the
    run, deleted when it finishes or is stopped (kept while it is paused);
  - an IAM role and instance profile `qap-phase1-instance`, which may only make
    the calls `coord/backend.py` makes on those (read and write `runs/`, read
    `src/`, write `public/runs/`, create, use and delete the `qap-units-*` queues,
    items on the two tables) and nothing else. Since Nov 19 the role may create
    and delete queues: the first launch after that comes from the laptop, which
    rewrites the role's policy (the GitHub workflow's role cannot);
  - a security group `qap-phase1` in the default VPC, with no inbound rules.
- **Disk and failures** (after tho40-a, Nov 14, whose 32 GB disk filled in 30
  minutes): `--disk GB` sizes the root disk (default 32; tho40 at 64 workers:
  100). Unit attempts and parent states are deleted as soon as nothing needs
  them (`coord/README.md`, Parent states and disk); a worker takes no unit
  while the work disk has under `QAP_MIN_FREE_GB` (default 2) free; and when
  the coordinator has failed 5 times the instance uploads `out/` and shuts
  down instead of idling. A run given up for good: `qapctl stop RUN`, then
  `qapctl purge-states RUN`.
- **Without an incumbent** (`--no-ub`): a better solution found inside a unit
  reaches the other units when that unit returns and is merged (the
  coordinator verifies and publishes it then), so at most one unit budget
  later: 60 CPU seconds with the nug30 recipe, 30 minutes with the default.
- **Several runs at once:** each run has its own queues (since Nov 19), so runs
  in a region are independent: a worker receives only its run's units. Before,
  every run shared two queues, and another run's large units kept a run's small
  ones from ever being read (the tail of `tho40-est-dives`, Nov 18:
  `docs/review_nov18_runs.md`, item 1). A run created before Nov 19 gets its
  queues when its coordinator next starts, which sends its queued units to them;
  the old shared queues are then unused (`cleanup.sh --storage` deletes them with
  the rest).
- **Interrupted or stopped runs:** `launch.sh RUN` resumes a run that exists;
  `--recipe` is then not needed. Relaunch on the same `--type`: the run's slot
  count is fixed when it is created, and the coordinator sizes its queue
  reserve and its steals by it. If a worker died without draining (the
  instance killed outright), its units come back when their leases run out (30
  minutes), or at once with `python3 -m coord.qapctl requeue-leased RUN
  --backend "aws:qap-coord-<account>-<region>?region=<region>"`.
- **Stopping and removing:** `infra/phase1/cleanup.sh RUN` terminates the run's
  instances. `--all` also deletes the security group, the profile and the
  role; `--all --storage` deletes the queues, the tables and the bucket too,
  after asking (every coordinator run in the region, including any made by
  hand with `qapctl`). The CloudFront distribution of `publish_setup.sh` stays:
  disable and delete it in the CloudFront console (a deletion waits for the
  disabled distribution to deploy).

## From GitHub (Nov 14)

`.github/workflows/phase1.yml` launches, collects or stops a run from GitHub
Actions, so a run can be started and its records committed without the laptop:
from the Actions tab, with `gh workflow run`, or from a Claude session that has
the repository. AWS gives each job short-lived credentials through OpenID
Connect, so no AWS key is stored anywhere.

```
gh workflow run phase1 --ref cognac -f action=launch -f run=tho40-b -f args="--recipe tho40 --on-demand --type c7a.16xlarge --no-reference --limit 4320 --disk 100"
gh workflow run phase1 --ref cognac -f action=collect -f run=tho40-b
gh workflow run phase1 --ref cognac -f action=stop -f run=tho40-b
```

- **Once, on the laptop:** `infra/phase1/gha_setup.sh` (`--dry-run` first),
  signed in to AWS as for `launch.sh` and to GitHub with `gh`. It creates the
  OIDC identity provider, the role `qap-gha-phase1`, the GitHub environment
  `aws-phase1` with you as required reviewer and only `cognac` allowed, the
  environment secret `AWS_ROLE_ARN`, and 14-day log retention. Run it again
  after a change to it.
- **Every job waits for your approval** (the environment's required
  reviewer), from the Actions tab or GitHub's phone app; nothing reaches AWS
  before that. Only a job of this repository in that environment can assume the
  role, so a fork or a pull request cannot.
- **What the role may do:** the calls `launch.sh`, `collect.sh` and
  `cleanup.sh RUN` make, and nothing else. It reads the qap storage and writes
  only `src/`; it launches instances only of the families `c6g c6i c7a c7g c7gd
  c7i c8g`, only tagged `qap-phase1`, and terminates only those; it passes the
  instance role but cannot change any IAM policy. So the laptop still does what
  needs more: the first launch in an account or region (the storage, the
  instance role, the security group), and a launch after a change to the
  instance role's policy in `launch.sh` (which now writes that policy only
  when it differs).
- **What is public** (the repository is): the job's log, with the account
  number masked, the records `collect` and `launch` commit under
  `runs/RUN/`, and a collect's trace artifact (solver output only; its manifest
  names no bucket). Any file there that names the account gets a placeholder
  before the commit, and the commit stops if the number is still in it.
- **Tested here** against stand-ins for `aws` and `gh` only (the calls in
  order, both policies as valid JSON, the input checks, the placeholder). The
  first real run tests the policy: a missing permission shows in the log as
  `AccessDenied` naming the action; add it to `gha_setup.sh` and run that again.

## The public status, for brixius-web

Every 5 minutes, and once at the end, the coordinator publishes the run under
`public/runs/RUN/` in the bucket (`coord/publish.py`): `status.json` (status,
units, incumbent, CPU, open nodes by depth), `totals.state` (the merged
statistics, state version 7: every per-level table, the DNN ones too),
`report.txt` (the solver's report for them), `history.json` (one sample per
publish) and `params.par` (the params file the run was given, for the run
page's branching rules). A `--trace` run (the tho40 recipes) also publishes
its tree under `tree/` (`coord/tree.py`): every node merged so far, for the
website's tree viewer, which reads it through the site, never from the
bucket. `--private` on a new run publishes nothing; `--publish-every` in a
recipe's `COORD` changes the interval (0: only at the end).

The bucket stays private. Once per account and region:

```
infra/phase1/publish_setup.sh --dry-run    # the calls it would make
infra/phase1/publish_setup.sh              # prints the URL base
```

**Nov 13, on real AWS:** `CreateDistribution` was refused with "Your account must be
verified before you can add new CloudFront resources": a new account needs an AWS
Support case (the error message says so) before any distribution exists. The
script now stops there; the origin access control it had created stays, unused
and free, and a rerun after verification completes the setup.

**Meanwhile, straight from S3** (Nov 13, while the Support case is open):
`infra/phase1/public_s3.sh on` grants public read on `public/*` and nothing else
(a bucket policy with `Principal: *` on that prefix; it turns off the two
policy-related Block Public Access flags and leaves the two ACL ones on, which
is all that is needed since object ownership is BucketOwnerEnforced). The base
URL is then
`https://qap-coord-<account>-eu-north-1.s3.eu-north-1.amazonaws.com/public`
(`public_s3.sh url`), so a run's record gets
`"live": {"url": "<base>/runs/RUN/"}`. Once the distribution exists, run
`publish_setup.sh` (it rewrites the policy to the CloudFront-only statement),
then `public_s3.sh off` to put all four blocks back, and point new records at
the CloudFront domain. S3 serves the objects with their `Content-Type` and
`Cache-Control: max-age=60` as the coordinator uploaded them; no CORS is
needed, since the run page fetches the status server-side.

It creates a CloudFront distribution that serves `public/` and nothing else
(origin path `/public`, an origin access control, and a bucket policy that lets
only that distribution read `public/*`), so
`https://<domain>/runs/RUN/status.json` is the run's status. The objects carry
`Cache-Control: max-age=60`. A new distribution takes 5 to 15 minutes to
deploy; it stays within CloudFront's free tier at a dashboard's traffic.
`python3 -m coord.qapctl publish RUN --backend ...` publishes at once (a run
started before the distribution existed has its bundle already: the
distribution only serves it).

## Output

`collect.sh RUN` prints the instances, the coordinator's `status.json` (units by
status, nodes, CPU, incumbent), each host's last steps, and at the end
`results.txt` and `compare.txt` (or `results_partial.txt`: the RESULTS block of a solve paused for an
estimate or stopped, as far as it went; with it, as with `results.txt`, the node trace goes in the
record), and keeps each host's `out/` under
`runs/RUN/hosts/`: `steps.txt`, `node.log`, `coordinator.log`, `worker.log`,
`init.log`, `reference.log`, `gate_*.log`, `sysinfo.txt`. `runs/RUN/launches.tsv`
records each launch (role, type, market, Spot price, instance, zone, commit,
recipe). `runs/RUN/bucket.txt` and `userdata_*.sh` stay on the laptop: they name
the bucket, whose name holds the account number, so `.gitignore` keeps them out
of the repository (removed from it Nov 14), and `collect.sh` finds the bucket
from the signed-in account when `bucket.txt` is missing. A `--trace` run's
traces (`qapctl export`: `nodes.txt.gz`, `dnn.txt.gz`, `manifest.json`) go to
`infra/phase1/exports/RUN/`, out of git; from GitHub Actions they are the job's
artifact `traces-RUN`, kept 90 days. The store keeps them too. A collect exports
them finished or not (since Nov 15): part-way, still running or ended by its
`--limit` (which writes no `results.txt`), they hold the units merged so far, and
`manifest.json` counts the units not yet merged.

## On real EC2 (Nov 13, eu-north-1, on demand; the quota arrived that day)

The first launches of `launch.sh` on EC2 proper (until then only moto and a
logged CLI stand-in). Runs `nug17-a` to `nug17-s` in `runs/`, c7g.xlarge and
c7i.xlarge, cents each.

- **The path works.** User data, the instance role, the Amazon Linux 2023
  toolchain, the build (one minute), uploads to the store, self-termination;
  `nug17-d` on x86 (c7i.xlarge, 2 cores) ran the six gates and the recipe to
  `PASS: totals equal to parallel.py and the optimum` (26 units). Two bugs on
  the way, fixed: `--on-demand` tripped bash 3.2's `set -u` on an empty array
  in `launch.sh`; and a gate that hangs idled the instance until `--limit`
  with an empty, block-buffered log (`nug17-a`), so gates now run `python -u`
  under a 30-minute `timeout` (`nug17-b`: `parallel TIMEOUT 1800s`, stopped).
- **Graviton3's OpenBLAS kernels break the DNN certificate.** On c7g.xlarge the
  DNN gate failed (`nug17-c`): the certified bound is `-inf` on nug12's root
  (prototype 567.98) and on every node with m >= 9, and the parallel gate's
  DNN check never ends because nothing fathoms. The QPB path, the state files
  and the parallel driver passed every check on the same instance. With
  `OPENBLAS_CORETYPE=ARMV8` (`nug17-e`) or `NEOVERSEN1` (`nug17-f`) the gate
  passes with x86's values to the digit and nug17 runs to PASS; `node.sh` now
  sets `NEOVERSEN1` on aarch64 (`--node-env OPENBLAS_CORETYPE=...` overrides it,
  `auto` leaves the choice to OpenBLAS). The probes (`QAP_PROBE`, below;
  `nug17-r`, `nug17-s`, scipy-openblas32 0.3.34, core `neoversev1`) call the
  library directly on random SPD matrices:

  | routine, auto kernels | wrong at n (residual 20 to 104, `info` 0) | right at n |
  |---|---|---|
  | `spotrf` lower | 48, 64, 65, 96, 128, 129, 145, 200 | 16, 32, 256, 300, 400 |
  | `dpotrf` lower | 48, 64, 65, 96 | 16, 32, 128 to 400 |
  | `spotrf` upper | none | all of them |
  | `dtrsm` right side, transposed | 400 (2.79) | 64, 145 |
  | `strsm` (both sides), `ssyrk`, `dsyrk`, `sgemm`, `sgemv`, `strtri`, `ssyevd`, `dsyevd` | none | 64 to 400 |

  The solver's certificate (`dpotrf` of the shifted matrix, lower) and the
  subspace eigensolver's CholQR (`spotrf` then `strtri`, lower) are exactly
  the broken calls; the `-P` prototype certificate prints `chol FAILED` at
  every check. Under the N1 kernels `dnntest` on nug12 costs the same per
  iteration (0.61 against 0.59 ms), so nothing is lost by the override.
  Worth an upstream report (OpenBLAS, target NEOVERSEV1, `potrf_L` at
  mid sizes).
- **nug30 acceptance on Graviton (`nug30-a`, c7g.4xlarge, 16 cores, the override
  on):** six gates passed; the coordinator's tree equals `bench/parallel.py`
  on the same instance (3,688 nodes, 3,547 fathomed, 380,512 FW iterations,
  3,687 DNN nodes), optimum 6124, 412 units, slots busy 72%: PASS. 4,961 CPU s for 3,688 nodes
  against the laptop's 1,423 for 3,716 (`nug30-local-a`, Accelerate): 1.35 s a
  node against 0.38, so a Graviton3 core (N1 kernels, no SVE) is about 3.5x
  slower than a laptop core on this recipe.
- **The tho40 smoke subtree on both families (`tho40sub-g`, `tho40sub-i`):** the
  coordinator gives the Nov 9 A/B's tree on each (110 nodes, 107 fathomed,
  14,278 FW iterations), so the tho40 recipe, its symmetry file and its
  settings file work through the store. Their first launch failed in the
  frontier pass: the source archive lacked `experiments/estimates`, which now
  travels with the settings files. The reference run failed on both because
  `bench/parallel.py` passed the recipe's `-F` to its workers, which the solver
  rejects with `-r`; it now strips it as the coordinator does (checked on a
  nug12 subtree: 72 nodes, 28,439 FW iterations, equal to the serial run).
  CPU per tho40 node, the price of the full run:

  | machine | slots | CPU s | per node |
  |---|---|---|---|
  | M2 Pro, local coordinator | 4 | 127 | 1.15 s |
  | c7i.4xlarge (Sapphire Rapids, 8 cores, 16 vCPUs) | 8 | 351 | 3.19 s |
  | c7g.4xlarge (Graviton3, 16 cores, N1 kernels) | 16 | 797 | 7.25 s |

  At the Nov 8 estimate of 0.6 to 3 million nodes: a c7i.16xlarge (32 cores,
  $2.86/h) 17 to 83 hours and $50 to $240 [e]; a c7g.16xlarge (64 cores,
  $2.32/h) 19 to 94 hours and $44 to $220 [e]. Either is the whole 64-vCPU
  on-demand quota.
- **`QAP_PROBE`** (`--node-env "QAP_PROBE='...'"`): a shell command `node.sh`
  runs in `latest/` after the build, its output uploaded as `out/probe.log`,
  for a diagnosis on the instance's CPU without logging in (there is no SSH
  and the role has no SSM). `make solver` does not build `dnntest`; a probe
  that needs it runs `make dnntest` first (`nug17-p` found out).
- **CloudFront** refused the distribution (the account must be verified by
  AWS Support, above); the live page waits on that, the runs do not.

## Tested here

There is no AWS access from the cloud session, so:

- `node.sh` ran in its local mode (`QAP_LOCAL=1`) from a clean `git archive`
  of the tree, against a local moto server for S3, SQS and DynamoDB:
  - nug17 end to end: toolchain record, build (with the pip OpenBLAS, as on
    the instance), the six gates, init, coordinator and worker (56 units),
    `parallel.py`, `compare.py`: PASS. `collect.sh` then fetched and showed it
    from the same store (the real AWS CLI against moto for S3);
  - Spot: nug24 with the nug30 DNN recipe (249 nodes, 62 CPU s) and a
    stand-in for the metadata service, switched to an interruption notice 23
    seconds in: node A (1 slot) drained in a second after 7 units and stopped;
    node B (a new tree and host, 4 slots) resumed the run from the store and
    finished it, and the totals equalled `parallel.py`'s (249 nodes, 237
    fathomed, 26,580 FW iterations, 248 DNN nodes);
  - a worker-only node started 10 seconds after the coordinating one, 2 slots
    each: the worker ran 34 of the 78 units, and the check passed;
  - a coordinator that fails at its first merge every time
    (`QAP_COORD_CRASH=marked:1`): five restarts, then the worker was stopped,
    `out/` uploaded and the node exited 1.
- `launch.sh` ran against a stand-in for the AWS CLI that logs every call
  (with moto behind `qapctl aws-setup` and the run lookups): a dry run, a new
  run whose first zone had no capacity (the next zone was used), a quota below
  the instance, a second coordinating instance for a run (refused), a new run
  without a recipe (refused), an existing run (resumed without one, x86 image
  for c7i, on-demand), worker-only with and without a coordinator, and another
  run coordinating in the region (noted).
- `tests/coord_test.py` has a new scenario, `tworuns` (two runs on one
  backend at once), which failed before the worker fix (both hung) and passes
  after; all eleven scenarios pass, `aws` included.
- Not tested: the Amazon Linux packages, the IAM policy (moto does not enforce
  IAM), the Spot launch and a real interruption notice. The first launch tests
  them; each failure shows in `steps.txt`, uploaded after every step.

## Run on the laptop (Nov 9)

The account's EC2 vCPU quota is 0 in every region tried (eu-north-1, eu-west-1,
eu-central-1, us-east-1, us-east-2, us-west-2; Spot `L-34B43A08` and
on-demand `L-1216C47A` both, against AWS's default of 5): a new account not yet
cleared for EC2. The Nov 5 Spot request for 256 vCPUs was closed on Nov 8 with
the quota unchanged, the on-demand request for 64 is still open, and a Spot
request for 64 in eu-north-1 was filed Nov 9 (no support plan, so the case
notes cannot be read from the CLI). `launch.sh` stops at the quota check until
one of them is granted. So the acceptance ran on the M2 Pro laptop with the
`local:` backend (SQLite and files in place of S3, SQS and DynamoDB), with the
same init, coordinator and worker, `parallel.py` and `compare.py` steps as
`node.sh` (`runs/nug30-local-a/run.sh`), 8 slots, binaries at ae992f4 with
Accelerate:

- **`nug30-local-a`, with the incumbent: PASS.** 118 one-minute units, wall
  286 s, CPU 1,423 s, slots busy 81.6%; nodes 3,716, fathomed 3,574, FW
  iterations 382,572, DNN nodes 3,715 and fathomed 3,574, all equal to
  `parallel.py`'s on the same settings (CPU 1,493 s, wall 252 s), objective
  6124 (`runs/nug30-local-a/compare.txt`).
- **`nug30-local-b`, without the incumbent: stopped after an hour, not a
  pass.** The solver has no root heuristic inside a branch-and-bound run (`-h`
  runs the heuristic instead of the search) and `qapctl init` does not call
  it, so the first 8 units started with no upper bound at all and each spent
  its 60 CPU s enumerating leaves (5.95 million nodes before the first
  incumbent, 7242, was merged). From then on every improvement came from the
  depth-first dive and was verified and passed through the store to the
  following units as designed, 56 times in the hour, 7242 down to 6336, still
  3.5% above 6124, with 393 units merged, 95 queued, 19.8 million nodes and
  22,649 CPU s. Reaching 6124 this way takes hours, and the `parallel.py`
  reference would repeat the descent. The incumbent pass-through works; the
  test as specified needs a starting incumbent. Proposed: `qapctl init` runs
  `bnbdrv -h` once and records its solution as the run's initial incumbent,
  then the `--no-ub` acceptance tests the pass-through from a realistic start
  (`runs/nug30-local-b/coordinator.log`).
- **`nug20-local-noub`, the no-incumbent check moved to nug20: PASS.** The
  nug30 recipe's settings on nug20 with no upper bound: 40 units, wall 95 s,
  CPU 415 s; five improving solutions found inside units (2602, 2600, 2594,
  2574, 2570), each verified by the coordinator and passed to the following
  units, and the run finished on 2570, the optimum. `parallel.py` without an
  incumbent also finished on 2570 with a different tree (25,101 nodes against
  23,310), as expected when incumbents arrive in a timing-dependent order
  (`runs/nug20-local-noub/compare.txt`). The nug30 variant without an
  incumbent is dropped: a real run always starts from the QAPLIB best known
  solution, and the pass-through it was meant to exercise is shown here.
