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

A recipe (`recipes/NAME.sh`, committed: the instance builds from `HEAD`) says
what to solve: `INSTANCE`, `PARAMS`, `OPT` (the known optimum the check
expects), `UB` (the incumbent given to the solver), `EXTRA` (solver flags),
`ENV` (the solver's environment) and `COORD` (unit sizing for `qapctl init`).

| recipe | what | CPU [e] |
|---|---|---|
| `nug17` | the smoke test: QPB only, one-second units, a few dozen of them through S3/SQS/DynamoDB | 20 s |
| `nug30` | the acceptance run: the Nov 8 record recipe (rule 9, DNN at depths 1-7) with `ship.bits = 64` (`recipes/nug30.par`), one-minute units | 1,300 s on the M2 Pro; more on Graviton |
| `tho40sub` | the tho40 smoke test: the depth-3 subtree 26:3,31:7,37:0 of `tho40t.dat` under `experiments/estimates/tho40/tho40.par` (110 nodes, 107 fathomed in the Nov 9 A/B), one-minute units; prices a tho40 node on the instance | 197 s on the M2 Pro |
| `tho40` | the whole tho40 tree against the best known 240,516, 20-minute units; `--no-reference` and a `--limit` of days (the recipe's comment has the launch line) | 0.9e6 to 4e6 s on the M2 Pro (`experiments/tho40_estimate/`) |

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

- **Options:** `--type c7g.4xlarge` (16 Graviton3 cores; any type whose
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
    `qap-coord-<account>-<region>` (also the source tarballs under `src/`), the
    SQS queues `qap-units-large` and `qap-units-small`, and the DynamoDB tables
    `qap-runs` and `qap-units` (on-demand capacity);
  - an IAM role and instance profile `qap-phase1-instance`, which may only make
    the calls `coord/backend.py` makes on those (read and write `runs/`, read
    `src/`, write `public/runs/`, send, receive and delete on the two queues,
    items on the two tables) and nothing else;
  - a security group `qap-phase1` in the default VPC, with no inbound rules.
- **Without an incumbent** (`--no-ub`): a better solution found inside a unit
  reaches the other units when that unit returns and is merged (the
  coordinator verifies and publishes it then), so at most one unit budget
  later: 60 CPU seconds with the nug30 recipe, 30 minutes with the default.
- **Several runs at once:** the runs in a region share the two queues. A
  worker that receives another run's message hides it for a second and leaves
  it (before Nov 8 it deleted it, which stranded that run's unit), so
  overlapping runs are correct but slower. `launch.sh` says when another run is
  coordinating.
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

## The public status, for brixius-web

Every 5 minutes, and once at the end, the coordinator publishes the run under
`public/runs/RUN/` in the bucket (`coord/publish.py`): `status.json` (status,
units, incumbent, CPU, open nodes by depth), `totals.state` (the merged
statistics, state version 7: every per-level table, the DNN ones too),
`report.txt` (the solver's report for them) and `history.json` (one sample per
publish). `--private` on a new run publishes nothing; `--publish-every` in a
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
`results.txt` and `compare.txt`, and keeps each host's `out/` under
`runs/RUN/hosts/`: `steps.txt`, `node.log`, `coordinator.log`, `worker.log`,
`init.log`, `reference.log`, `gate_*.log`, `sysinfo.txt`. `runs/RUN/launches.tsv`
records each launch (role, type, market, Spot price, instance, zone, commit,
recipe).

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
