# Phase 0: which CPU to rent (Nov 5)

Phase 0 of `docs/distributed/plan.md`: the DNN bound's speed and dollars
per node on the CPU families a fleet could use, before anything bigger is
built on one of them. One Spot instance per family. Each builds the solver
from the committed tree and checks it. It times the 32 wil50 nodes of
`tuning/giveup/all32.txt` with the production DNN flags, uploads its results
to a private S3 bucket and terminates itself.

| family | CPU | why it is in the list |
|---|---|---|
| c7i | Xeon Sapphire Rapids | AMX: the `DNN_AMX=1` kernel (the cloud box is this CPU at 2.1 GHz) |
| c7a | EPYC Genoa | x86 without AMX; one vCPU is a whole core |
| c7g | Graviton3 | about a third of x86's Spot price per core in the plan's Oct 31 snapshot; the first ARM Linux build |
| c8g | Graviton4 | the same, newer |

## What each instance does (`bench.sh`)

1. It records the machine: the instance type and the CPU model, its
   physical cores, whether it has AMX, and memory.
2. It builds on Amazon Linux 2023: gcc, Python 3.11, numpy and the
   scipy-openblas32 wheel (the BLAS the solver links on Linux), then
   `make solver dnntest`.
3. It runs the gates: `make bruteforce resume parallel regression dnncheck`
   (about 4 minutes). A failed gate is recorded and the benchmark goes on.
4. **alone**: nodes 19, 21, 1 and 11 one after another, pinned to one core,
   on an otherwise idle machine (the single-core speed).
5. **loaded**: all 32 nodes, one process per physical core, from a shared
   queue. Once the queue is empty, the idle cores run short filler solves,
   so the machine stays fully loaded until the last real node ends.
6. `summarize.py instance` writes `out/summary.json`: ms per iteration
   alone and loaded, CPU s per node, and the fathom set against the cloud
   box's run of the same 32 nodes (18 fathomed). Iteration counts can
   differ at the edges on another CPU; they are listed, not hidden.

It uploads `out/` after every step, so a failed build or an interrupted
instance still reports how far it got. It shuts itself down when it ends.
Separately, it shuts itself down at `--limit` minutes (75) whatever
happens, and the instance terminates when it shuts down.

Flags: `-T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 2500 -u 48816`,
`DNN_CHOLQR_PASSES=1 DNN_AMX=1` (bundle D; the AMX switch falls back by
itself on CPUs without AMX), one BLAS thread.

## Running it (laptop)

Once: the AWS CLI (`brew install awscli`) and a login (`aws configure sso`;
`aws sso login` when it expires).

```
git pull origin cognac
infra/phase0/launch.sh --dry-run                  # the AWS calls it would make, none made
infra/phase0/launch.sh                            # eu-north-1; c7i c7a c7g c8g, 4xlarge Spot
infra/phase0/collect.sh                           # any time: instance states, results so far, the table
git add infra/phase0/runs && git commit -m 'phase 0: ...' && git push
```

- **Options:** `--region us-west-2`, `--families "c7i c8g"`,
  `--size 2xlarge`, `--on-demand`, `--serial` (one instance at a time),
  `--limit 75`.
- **Quota:** the launcher reads the Spot vCPU quota first. Four 4xlarge
  instances need 64 vCPUs. If the quota covers only one, it runs them one
  at a time; if it covers none, it stops and says so (use `--size 2xlarge`,
  or `--on-demand` if that quota is higher). A new account's default can be
  as low as a few vCPUs.
- **Time:** about 30-45 minutes per instance [e]: 5-10 minutes of setup and
  build, 4 of gates, 5 alone, 5-10 loaded. All four run at the same time
  unless the quota says otherwise.
- **Cost:** four 4xlarge Spot instances for under an hour: about $1-3 in
  eu-north-1 at the plan's snapshot prices [e], plus a few cents of S3.
- **What it creates,** once per account and region, kept for later runs:
  - a private bucket `qap-phase0-<account>-<region>`, with public access blocked;
  - an IAM role and instance profile `qap-phase0-instance`, which may only
    read `src/` and write `results/` in `qap-phase0-<account>-*` buckets;
  - a security group `qap-phase0` in the default VPC, with no inbound rules
    (nobody can connect to the instances).
- **Stopping and removing:** `infra/phase0/cleanup.sh` terminates any
  instance still running. `--all` also deletes the security group, the
  profile and the role; `--all --bucket` deletes the bucket too, after
  asking.

## Output

`collect.sh` writes `runs/<run>/summary.md`, one row per family:

- the CPU, its cores and whether it has AMX;
- ms per iteration alone and loaded, and their ratio (the memory-bandwidth
  loss at full load, risk 8 of the plan);
- CPU s per node under load;
- the Spot price at launch, the price per core-hour, and dollars per 1,000
  DNN nodes on this sample;
- the gates and the fathom set.

The 32 nodes are a fixed sample, not the production mix of nodes. Its
dollars per node rank the families; they are not the attempt's cost per
node.

## Tested here

There is no AWS access from the cloud session, so:

- `launch.sh` ran against a stand-in for the AWS CLI that logs every call,
  including a dry run, a full launch, and a zone with no capacity (the next
  zone is tried). It also ran with quotas below the total (switches to one
  at a time) and below one instance (stops). The arm64 families get the
  arm64 image.
- `bench.sh` ran in its local mode (`QAP_LOCAL=1`) on the cloud box's Ubuntu,
  from a clean `git archive` of exactly what the instances download: the
  build, the gates, one node alone, four nodes on two cores with the queue
  and fillers, and the summary (4 of 4 fathomed, as in the reference).
  `summarize.py table` was run on that output plus a family whose build
  failed (reported, not dropped).
- **Correction (Nov 5, late):** in that run and in the first gcc 11 check,
  the regression gate passed without checking anything. The tarball had no
  `QAPSOLN/`, so `bench.py check` skipped every instance, and it reported
  success when nothing ran. Both are fixed: the tarball carries `QAPSOLN/`,
  and `bench.py check` now fails when no instance runs. The gcc 11 build
  (stage 0 makefile) then passes the regression check on all 7 instances,
  as do brute force (4,536 of 4,536) and `dnncheck`.
- Meschach now builds with the system compiler and `-ffp-contract=off`
  (stage 0 of `docs/meschach_removal.md`), so the instances need no clang,
  and their x86 and ARM builds round alike. Since Nov 6 (stages 2 and 3)
  there is no Meschach to build: the solver needs only `make`, `g++`, `cc`
  for `amx_gemm.c` and LAPACK, and the tarball has no `src/`.
- Not tested: the Amazon Linux packages, the ARM Linux build, the IAM
  role's S3 access from an instance, and the Spot launch itself. Those are
  what the first run tests. Each failure shows in `steps.txt`, which every
  instance uploads after each step.
